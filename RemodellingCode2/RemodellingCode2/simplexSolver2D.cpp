/* simplexSolver.cpp : Solves the Stress Control Solution for Displacement (stretch) Periodic Boundary Condition Problem
with a Periodic Network.

SUMMARY: This function solves the principal stretches necessary to yield the specified stresses using the Nelder-Mead Simplex
algorithm. The convergence state is 0.5% of the maximum stress as the magnitude of the error or when the simplex edge length sum
is below 1e-6.

CALLED FROM: Main.cpp

CALLS ON: networkSolver.cpp

INPUTS: nodes - Nx3 matrix of nodal coordinates in an RVE
		fibers - Mx5 matrix of fiber connectivity, fibers (m,1) and fibers(m,2) represent the number of the two nodes that are connected,
			fibers(m,3:5) are the x,y,z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
			the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
		fib_type - Mx1 vector of fiber type (1 for collagen, 2 for elastin, 3 for failed fibers)
		fib_rads - Mx1 vector of fiber radii (in meters)
		init_lens - Mx1 vector of fiber initial lengths
		fiber_vol_fract - double value or actual tissue fiber volume fraction (typically 10-20% fiber, 80% water/goop)
		F - 3x1 vector of principal network stretches
		num_nodes - int N  number of nodes
		num_fibers - int M number of fibers
		net_stress - 3x3 tensor of volume averaged stress
		guess - bool, 1 if nodes_n contains the solution from a previous iteration, 0 if nodes_n is zeros

OUTPUTS: None

NOTE: This code modifies 'nodes_n', bnd_nodes','fibers_n', 'net_stress', 'fib_forces', 'fib_stress', and 'F'

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-27-2018

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"
#include <vector>

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;


void simplexSolver2D(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess, double fib_vol_m, double phi_m, double P)
{
	double c = 0.25; // parameter for generating initial simplex

	double alpha = 1; double gamma = 2; double beta = 0.5;  double rho = 0.5; // parameters manipulating simplex

	double b = c / (3 * sqrt(2)); double a = b + c / sqrt(2); // values for making equal side length tetrahedron 

	MatrixXd x(3, 3); x = MatrixXd::Zero(3, 3); // matrix containing tetrahedron coordinates (lam1, lam2, lam3, err)

	double err_thresh = 1e-6*appliedStress.maxCoeff(); // error threshold (0.5% of maximum stress)
	double simp_size_thresh = 1e-6; // converged size of simplex
	double simp_size = 1000; // size of simplex

	if (err_thresh > 1000)
	{
		err_thresh = 1000;
	}

	// generating tetrahedron initial guess values
	x(0, 0) = x_guess(0);
	x(0, 1) = x_guess(1);
	x(0, 2) = fib_vol_m / (phi_m*x(0, 0)*x(0, 1)*x_scale*x_scale);


	x(1, 0) = x(0, 0) + a;
	x(1, 1) = x(0, 1) + b;
	x(1, 2) = fib_vol_m / (phi_m*x(1, 0)*x(1, 1)*x_scale*x_scale);


	x(2, 0) = x(0, 0) + b;
	x(2, 1) = x(0, 1) + a;
	x(2, 2) = fib_vol_m / (phi_m*x(2, 0)*x(2, 1)*x_scale*x_scale);

	double err = 1e20; // initializing error
	int flag = 0; // flag to signify first simplex within error tolerance

	//#pragma omp parallel for
	for (int j = 0; j < 3; j++) // generating first error values
	{
		F(0,0) = x(j, 0);
		F(1,1) = x(j, 1);
		F(2,2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

		solve(num_nodes, num_fibers, x_scale, F, network0,
			network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

		x(j, 2) = err;

		if (err < err_thresh)
		{
			flag = 1;
			break;
		}
	}

	if (flag == 0)
	{
		VectorXd xa(3); xa = VectorXd::Zero(3); // vector containing average coordinates (lam1, lam2, lam3, err)
		VectorXd xr(3); xr = VectorXd::Zero(3); // vector containing reflected coordinates (lam1, lam2, lam3, err)
		VectorXd xe(3); xe = VectorXd::Zero(3); // vector containing expanded coordinates (lam1, lam2, lam3, err)
		VectorXd xic(3); xic = VectorXd::Zero(3); // vector containing inside contraction coordinates (lam1, lam2, lam3, err)
		VectorXd xoc(3); xoc = VectorXd::Zero(3); // vector containing outside contraction coordinates (lam1, lam2, lam3, err)

		int num_iter = 200; // number of iterations

		for (int n = 0; n < num_iter; n++) // looping through simplex algorithm
		{

			//cout << "Collagen Fraction: " << phi_C << endl << endl;

			bubbleSort(x);

			simp_size = 0;
			for (int w = 0; w < 2; w++)
			{
				simp_size = simp_size + sqrt((x(w, 0) - x(2,  0))*(x(w, 0) - x(2, 0)) + (x(w, 1) - x(2, 1))*(x(w, 1) - x(2, 1)));
			}

			//cout << "Best Stretches: "<< x(0,0)<<", "<<x(0,1) << ", "<<x(0,2) << endl;
			//cout << "Error: " << x(0, 3) << endl << endl;
			//cin.get();
			
			if (x(0, 2) < err_thresh || simp_size < simp_size_thresh)
			{
				break;
			}

			if (n == num_iter - 1)
			{
				cout << "Simplex Failed to Converge!" << endl;
			}

			// average position excluding worst point
			xa(0) = (x(0, 0) + x(1, 0)) / 2;
			xa(1) = (x(0, 1) + x(1, 1)) / 2;

			xr(0) = xa(0) + alpha * (xa(0) - x(2, 0));
			xr(1) = xa(1) + alpha * (xa(1) - x(2, 1));

			F(0,0) = xr(0);
			F(1,1) = xr(1);
			F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

			solve(num_nodes, num_fibers, x_scale, F, network0,
				network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

			xr(2) = err;

			if (xr(2) < x(0, 2)) // if reflected is better than best point
			{
				// expand point beyond reflected coordinates
				xe(0) = xa(0) + gamma * (xr(0) - xa(0));
				xe(1) = xa(1) + gamma * (xr(1) - xa(1));

				F(0,0) = xe(0);
				F(1,1) = xe(1);
				F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

				solve(num_nodes, num_fibers, x_scale, F, network0,
					network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

				xe(2) = err;

				if (xe(2) < xr(2)) // keep expanded
				{
					x(2, 0) = xe(0);
					x(2, 1) = xe(1);
					x(2, 2) = xe(2);
				}
				else // keep reflected
				{
					x(2, 0) = xr(0);
					x(2, 1) = xr(1);
					x(2, 2) = xr(2);
				}
			}
			else // if reflected is worse than best point
			{
				if (xr(2) > x(2, 2)) // if reflected is worse than worst point
				{
					// inside contraction on tetrahedron
					xic(0) = xa(0) - beta * (xa(0) - x(2, 0));
					xic(1) = xa(1) - beta * (xa(1) - x(2, 1));

					F(0,0) = xic(0);
					F(1,1) = xic(1);
					F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

					solve(num_nodes, num_fibers, x_scale, F, network0,
						network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

					xic(2) = err;

					if (xic(2) > x(2, 2)) // if inside contraction is worse than worst point
						for (int k = 0; k < 2; k++) // performing shrink operation on tetrahedron
						{
							x(k + 1, 0) = x(0, 0) + rho * (x(k + 1, 0) - x(0, 0)); // shrink
							x(k + 1, 1) = x(0, 1) + rho * (x(k + 1, 1) - x(0, 1)); // shrink

							F(0,0) = x(k + 1, 0);
							F(1,1) = x(k + 1, 1);
							F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

							solve(num_nodes, num_fibers, x_scale, F, network0,
								network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

							x(k + 1, 2) = err;
						}
					else // if inside contraction is better than worst point
					{
						x(2, 0) = xic(0);
						x(2, 1) = xic(1);
						x(2, 2) = xic(2);
					}
				}
				else // if reflected point is better than worst point
				{
					xoc(0) = xa(0) + beta * (xr(0) - xa(0)); // outside contraction
					xoc(1) = xa(1) + beta * (xr(1) - xa(1));

					F(0,0) = xoc(0);
					F(1,1) = xoc(1);
					F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

					solve(num_nodes, num_fibers, x_scale, F, network0,
						network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

					xoc(2) = err;

					if (xoc(2) > x(2, 2)) // outside contraction worse than worst point
					{
						for (int k = 0; k < 2; k++) // shrinking tetrahedron
						{
							x(k + 1, 0) = x(0, 0) + rho * (x(k + 1, 0) - x(0, 0)); // shrink
							x(k + 1, 1) = x(0, 1) + rho * (x(k + 1, 1) - x(0, 1)); // shrink

							F(0,0) = x(k + 1, 0);
							F(1,1) = x(k + 1, 1);
							F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

							solve(num_nodes, num_fibers, x_scale, F, network0,
								network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

							x(k + 1, 2) = err;
						}
					}
					else // outside contracted better than worst point
					{
						if (xoc(2) < xr(2)) // if outside contraction is better than reflected
						{
							x(2, 0) = xoc(0); // keep outside contracted point
							x(2, 1) = xoc(1);
							x(2, 2) = xoc(2);
						}
						else // if outside contraction is worse than reflected
						{
							x(2, 0) = xr(0); // keep reflected point
							x(2, 1) = xr(1);
							x(2, 2) = xr(2);
						}
					}
				}

			//cout << "Simplex Error: " << x(0,3) << endl;
			//cout << "Pressure: " << pressure << endl;
			//cout << "Fiber Volume Fraction: " << phi << endl;
			//cout << "Fixed Charge Density: " << c_fcd << endl << endl;
			}
		}

		bubbleSort(x);

		F(0,0) = x(0, 0);
		F(1,1) = x(0, 1);
		F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

		solve(num_nodes, num_fibers, x_scale, F, network0,
			network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);
	}

	cout << "Simplex Stress Error: " << err << endl << endl;
}