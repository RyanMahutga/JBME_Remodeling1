/*
SUMMARY: This function solveSimps the principal stretches necessary to yield the specified stresses using the Nelder - Mead Simplex
algorithm.The convergence state is 0.5% of the maximum stress as the magnitude of the error or when the simplex edge length sum
is below 1e-6.

CALLED FROM : Main.cpp

CALLS ON : networksolveSimpr.cpp

INPUTS : nodes - Nx3 matrix of nodal coordinates in an RVE
fibers - Mx5 matrix of fiber connectivity, fibers(m, 1) and fibers(m, 2) represent the number of the two nodes that are connected,
fibers(m, 3:5) are the x, y, z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
fib_type - Mx1 vector of fiber type(1 for collagen, 2 for elastin, 3 for failed fibers)
fib_rads - Mx1 vector of fiber radii(in meters)
init_lens - Mx1 vector of fiber initial lengths
fiber_vol_fract - double value or actual tissue fiber volume fraction(typically 10 - 20 % fiber, 80 % water / goop)
F - 3x1 vector of principal network stretches
num_nodes - int N  number of nodes
num_fibers - int M number of fibers
net_stress - 3x3 tensor of volume averaged stress
guess - bool, 1 if nodes_n contains the solution from a previous iteration, 0 if nodes_n is zeros

OUTPUTS : None

NOTE : This code modifies 'nodes_n', bnd_nodes','fibers_n', 'net_stress', 'fib_forces', 'fib_stress', and 'F'

Author : Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update : 05 - 27 - 2018

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

void jacobianSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P, int type_flag, double& WSS)
{
	int check1 = 0;
	double dlam = 1e-4, mag_err = 1e6, err_thresh = 100, newflag = -1.0, cond = 0.0;
	Matrix3d dF, J, F_o;
	Vector3d epsilon, init_err, err;
	
	if (type_flag == 0)
	{
		dlam = 1e-5;
	}

	F_o = F;

	for (int n = 0; n < 50; n++)
	{
		// clear vars
		J = Matrix3d::Zero();
		check1 = 0;

		// initial solve
		solveJac(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
			guess, stab_node, appliedStress, net_stress, init_err, P, fib_vol_m, WSS);

		if (net_stress(2, 2) < 0)
		{
			check1 = 1;
		}

		mag_err = init_err.norm();

		//cout << "Error: " << mag_err << endl;

		if (mag_err < err_thresh)
		{
			netResults_n.F = F;
			cout << "Converged! Error: " << mag_err << endl << endl;
			break;
		}
		if (n == 49)
		{
			cout << " Failed to Converge - Error: " << mag_err << endl << endl;
		}

		for (int i = 0; i < 3; i++)
		{
				dF = F;
				dF(i, i) = F(i, i) + (dlam);

				solveJac(num_nodes, num_fibers, x_scale, dF, network0,
					network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, 
					net_stress, err, P, fib_vol_m, WSS);

				J(0, i) = (err(0) - init_err(0)) / (dlam);
				J(1, i) = (err(1) - init_err(1)) / (dlam);
				J(2, i) = (err(2) - init_err(2)) / (dlam);

				if (net_stress(2, 2) < 0)
				{
					check1 = 1;
				}
		}

		LDLT<Matrix3d> LDLT = LDLT.compute(J);

		cond = 1 / LDLT.rcond();
			/*<< endl << J(0,0) << endl;
		if (LDLT.rcond() < 1e-12*abs(J(0, 0)) && n > 0)
		{
			epsilon = 0.5*newflag*epsilon;
			newflag = 1.0;
		}
		else
		{
		*/
			epsilon = LDLT.solve(-1.0*init_err);
		//	newflag = -1.0;
			//epsilon = J.jacobiSvd().solve(-1.0*init_err);
//		}
		/*
		if (type_flag == 0)
		{
			cout << "Deformation Guess: " << endl << F << endl << "Epsilon: " << epsilon << endl << endl;
		}
		*/
		//cout << "External Condition #: " << cond << endl << endl;


				if (epsilon(0) < -0.005)
				{
					epsilon(0) = -0.005;
				}
				else if (epsilon(0) > 0.005)
				{
					epsilon(0) = 0.005;
				}
				if (epsilon(1) < -0.005)
				{
					epsilon(1) = -0.005;
				}
				else if (epsilon(1) > 0.005)
				{
					epsilon(1) = 0.005;
				}
				if (epsilon(2) < -0.005)
				{
					epsilon(2) = -0.005;
				}		
				else if (epsilon(2) > 0.005)
				{
					epsilon(2) = 0.005;
				}

		if ( n > 5 && n < 35)
		{
			F(0, 0) += epsilon(0);
			F(1, 1) += epsilon(1);
			F(2, 2) += epsilon(2);
		}
		else
		{
			F(0, 0) += 0.25*epsilon(0);
			F(1, 1) += 0.25*epsilon(1);
			F(2, 2) += 0.25*epsilon(2);
		}

		cout << "Deformation: " << endl << F << endl << "Wall Shear Stress: " << WSS << endl<< endl;

		netResults_n.F = F;
	}
}