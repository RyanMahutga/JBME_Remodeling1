/* networkSolver.cpp : Solves the Displacement (stretch) Periodic Boundary Condition Problem for a Periodic Network.

SUMMARY: This function solves the force balance on nodes in a periodic fiber RVE subjected to a user defined stretch.
The solution can use any of three constitutive equations for fibers (linear, exponential, of helical spring), the
network equilibrium is found through Netwon's Method relying on a Conjugate Gradient solver in the Eigen Library.

CALLED FROM: simplexSolver.cpp

CALLS ON: calcForcesPeriodic.cpp, calcJacobian.cpp, calcBndNodes.cpp, calcNetStress.cpp

INPUTS: nodes - Nx3 matrix of nodal coordinates in an RVE
		fibers - Mx5 matrix of fiber connectivity, fibers (m,1) and fibers(m,2) represent the number of the two nodes that are connected,
			fibers(m,3:5) are the x,y,z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
			the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
		fib_type - Mx1 vector of fiber type (1 for collagen, 2 for elastin, 3 for ?)
		fib_rads - Mx1 vector of fiber radii (in meters)
		init_lens - Mx1 vector of fiber initial lengths
		fiber_vol_fract - double value or actual tissue fiber volume fraction (typically 10-20% fiber, 80% water/goop)
		RVE_stretch - 3x1 vector of principal network stretches
		num_nodes - int N  number of nodes
		num_fibers - int M number of fibers
		net_stress - 3x3 tensor of volume averaged stress
		guess - bool, 1 if nodes_n contains the solution from a previous iteration, 0 if nodes_n is zeros

OUTPUTS: None

NOTE: This code modifies 'nodes_n', bnd_nodes','fibers_n', 'net_stress', 'fib_forces', 'fib_stress', and 'rve_stretch'

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-07-2018

*/
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>
#include <ctime>
#include <Eigen/PardisoSupport>
//#include "boost/numeric/odeint.hpp"
#include "dataStructureDefinitions.h"

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;
//using namespace boost::numeric::odeint;
/*
template<size_t N>
using vector = Matrix<double, N, 1>;

typedef vector<3> state;

// state_type = double
typedef runge_kutta_dopri5<state, double, state, double, vector_space_algebra> stepper_type;

void rhs( MatrixXd& node_forces, const int& num_fibers, const int& num_nodes)
{
	double C = 1e6;
	VectorXd dxdt(3 * num_nodes);
	//calcForcesPeriodic(nodes_n, nodes_mapped, fibers_n, init_lens, nodes0, fib_type, fib_areas, F,
		//fib_forces, fib_stretch, node_forces, num_fibers, num_nodes);
	dxdt = node_forces / C;
}
*/
void networkSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0, 
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess, VectorXi stab_node, double WSS)
{
	// check for valid stretches
	if (F(0,0) < 0.1 || F(1,1) < 0.1 || F(2,2) < 0.1)
	{
		cout << "RVE stretch too low!" << endl;
		cout << "Possible instability!" << endl;
		cin.get();
		exit(2); // terminate with error
	}

	int iter_max = 200; // maximum Newton Iterations 

	VectorXd residual(3 * num_nodes); residual = VectorXd::Zero(3 * num_nodes); // residual nodal forces
	VectorXd node_forces(3*num_nodes); node_forces = VectorXd::Zero(3 * num_nodes); // nodal forces 
	
	MatrixXd nodes_n = network_n.nodes;
	MatrixXd nodes = network0.nodes;

	MatrixXd fibers = network0.fibers;
	MatrixXd fibers_n = network_n.fibers;

	int flag = 0;

	double err = 1e-6,
		err_tol = 1e-15,
		rel_err = 10,
		max_err = 1,
		min_err = -1,
		max_err_thresh = 1e-15,
		err_init = 0,
		err_fract = 1,
		err_fract_tol = 1e-8,
		damp_fact = 1.0,
		cond = 0.0,
		move_thresh = 1e-4;

	MatrixXd vect(num_fibers, 3); vect = MatrixXd::Zero(num_fibers, 3); // fiber vector from node1 to real position of node2
	VectorXd lambda(num_fibers); lambda = VectorXd::Zero(num_fibers); // fiber stretch
	VectorXd dFdlam(num_fibers); dFdlam = VectorXd::Zero(num_fibers); // derivative of force with respect to lambda (analytically calculated)
	VectorXd epsilon(3 * num_nodes); epsilon = VectorXd::Zero(3 * num_nodes);
	
	MatrixXd nodes_mapped(num_nodes, 3); nodes_mapped = nodes; // nodal positions mapped back to undeformed unit cube (update every iteration)
	
	MatrixXd revert_nodes = nodes_n;
	MatrixXd revert_fibs = fibers_n;

	MatrixXd F_inv(3, 3), I(3, 3), T(3,3), T2(3,3); 
	
	F_inv = F.inverse(); // inverse deformation mapping
	I << 1, 0, 0, 0, 1, 0, 0, 0, 1; // creating identity matrix
	T = 0.5*(F-I); // rigid translation to maintain lower left corner at -0.5*(1 1 1)
	T2 = 0.5*(F_inv-I); // inverse of the rigid translation

	SparseMatrix<double> J(3 * num_nodes, 3 * num_nodes); //J = MatrixXd::Zero(3 * num_nodes, 3 * num_nodes); // initialize Jacobian


		for (int k = 0; k < num_nodes; k++) // initial guesses for deformed nodal coordinates
		{
			if (guess == 0)
			{
				nodes_n(k, 0) = F(0, 0)*(nodes(k, 0)) + F(0, 1)*(nodes(k, 1)) + F(0, 2)*(nodes(k, 2));
					//+ (T(0, 0) + T(0, 1) + T(0, 2)); // Make initial guess for deformed nodal coordinate x
					nodes_n(k, 1) = F(1, 0)*(nodes(k, 0)) + F(1, 1)*(nodes(k, 1)) + F(1, 2)*(nodes(k, 2));
					//+ (T(1, 0) + T(1, 1) + T(1, 2)); // Make initial guess for deformed nodal coordinate y
					nodes_n(k, 2) = F(2, 0)*(nodes(k, 0)) + F(2, 1)*(nodes(k, 1)) + F(2, 2)*(nodes(k, 2));
					//+ (T(2, 0) + T(2, 1) + T(2, 2)); // Make initial guess for deformed nodal coordinate z	
			}
			//nodes_mapped(k, 0) = F_inv(0, 0)*nodes_n(k, 0) + F_inv(0, 1)*nodes_n(k, 1) + F_inv(0, 2)*nodes_n(k, 2);
				//+ (T2(0, 0) + T2(0, 1) + T2(0, 2));;
				//nodes_mapped(k, 1) = F_inv(1, 0)*nodes_n(k, 0) + F_inv(1, 1)*nodes_n(k, 1) + F_inv(1, 2)*nodes_n(k, 2);
				//+ (T2(1, 0) + T2(1, 1) + T2(1, 2));
				//nodes_mapped(k, 2) = F_inv(2, 0)*nodes_n(k, 0) + F_inv(2, 1)*nodes_n(k, 1) + F_inv(2, 2)*nodes_n(k, 2);
					//+ (T2(2, 0) + T2(2, 1) + T2(2, 2));			
		}
	//fibers_n = fibers; // setting deformed fibers as intial fibers

	//cout << "NODES: " << endl << nodes << endl << endl;
	//cout << "NODES MAPPED: " << endl << nodes_mapped << endl << endl;
	//cin.get();

	//clock_t start_loop = clock();
	//double dur_loop;

	//int its = 0;
	//double t = 0.0;

	updateCrossing(nodes_mapped, nodes_n, fibers_n, num_nodes, num_fibers, F); // update crossing conditions

	network_n.nodes = nodes_n;
	network_n.fibers = fibers_n;
	network0.nodes = nodes_mapped;
	network0.fibers = fibers_n;

	int m = 0;

	int idx1i, idx2i,idx1j, idx2j;
	idx1i = idx1j = idx2i = idx2j = 56;

	if (0 == 1) // solve with Dormand-Prince (ODE45) using ODEint
	{
		//integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()),
			//rhs, nodes_n, 1.0, 10.0, 0.01, nodes_n);

	}
	else // solve via newton
	{
		for (int iter = 0; iter < iter_max; iter++) // Newton's Method
		{
			J.setZero(); // set all values of Jacobian to zero
			node_forces.setZero(); // set node forces to zero

			// Jacobain and Residual Calculation
			calcJacobian2(nodes_n, fibers_n, network_n, F, node_forces, num_fibers, num_nodes, J, stab_node, WSS); // Jacobain calculation muust be done at least once to store Jacobian

			//cout << J << endl;
			//cin.get();

			residual = -1.0*node_forces;

			if (iter == 0)
			{
				err_init = residual.norm();
			}
			else if (iter > 0)
			{
				rel_err = abs(residual.norm() - err);
				err_fract = abs(residual.norm() / err_init);
			}

			err = residual.norm(); // error in the sum of nodal forces
			max_err = residual.maxCoeff(&idx1i, &idx1j); // maximum error
			min_err = residual.minCoeff(&idx2i, &idx2j);

			//cout << "Max Min Residual: " << max_err << "  " << min_err << endl << "Nodal Positions: " << int(idx1i / 3) << "  " << int(idx2i / 3) << endl << nodes_n.row(int(idx1i/3)) <<"  " << nodes_n.row(int(idx2i / 3)) << endl;
			/*
			if (max_err > 1e-2)
			{
				network_n.nodes = nodes_n;
				network_n.fibers = fibers_n;
				netResults_n.F = F;

				calcForcesPeriodic(nodes_n, nodes_mapped, fibers_n, network_n, F,
					fiberResults_n, netResults_n, num_fibers, num_nodes);

				componentStress(nodes_n, fibers_n, network_n, netResults_n, fiberResults_n, F, num_fibers, x_scale);

				writeData(num_fibers, network_n, fiberResults_n, netResults_n, m);

				//writePeriodicNet2File(network_n, num_fibers, num_nodes, stab_node, m);

				ofstream myfile1;
				string str1 = "nodes_name" + to_string(m) +".txt";
				myfile1.open(str1);
				myfile1 << int(idx1i / 3) << int(idx2i / 3) ;
				myfile1.close();
				
				m++;
			}
			*/
			if (err < err_init)
			{
				revert_nodes = nodes_n;
				revert_fibs = fibers_n;
			}

			//cout << "Newton Iteration: " << iter + 1 << endl << "Error = " << err << endl << endl;
			//cin.get();

			if (rel_err < err_tol || err_fract < err_fract_tol || (err < err_tol && (max_err < max_err_thresh && abs(min_err) < max_err_thresh)))
			{
				//cout << "Converged!" <<endl<< "Newton Iterations: " << iter << endl
					//<< "Absolute Error: " << err << endl << "Relative Error: "<<rel_err<<endl
					//<<"Fractional Error: "<< err_fract << endl<<endl;
				//its = iter
				netResults_n.net_flag = 0;
				break; // break out of for loop if residual error drops below tolerance
			}

			 // Uncomment this to measure matrix condition (it's very slow!!!)
			//MatrixXd Jac = MatrixXd(J); // convert sparse to dense for calculation of jacobian condition number

			//BDCSVD<MatrixXd> svd(Jac); // calculating Jacobian condition number
			//cond = svd.singularValues()(0)
			//	/ svd.singularValues()(svd.singularValues().size() - 1);
			//cout << "Condition Number = " << cond << endl << endl;
			//if (cond > 1e19) // terminates solutuion if the condition number is too high
							 //Potential solutions to this are increaing k in calcJacobian or changing the network to have a higher connectivity.
			//{
				//cout << "Condition Number = " << cond << endl << endl;
				//cout << "Potential Singularity! Terminating Solution..." << endl << endl;
			//	break;
			//}
			
			// clearing epsilon
			epsilon = VectorXd::Zero(3 * num_nodes); // zero out epsilon (shift vector)

			//clock_t start = clock();
			//double duration;

			if (iter > 225)
			{
				// solving using Conjugate Gradient Solver (this should be in the ballpark of the matlab solution)
				//DGMRES<SparseMatrix<double>, Lower | Upper, DiagonalPreconditioner IncompleteCholesky<double>>  solver;
				ConjugateGradient<SparseMatrix<double>, Lower | Upper, DiagonalPreconditioner<double>>  solver;
				//BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>>  solver;
				solver.setTolerance(0.001*err);
				solver.setMaxIterations(3*num_nodes);
				solver.compute(J);
				epsilon = solver.solve(residual);
				//cout << "CG Solver Iterations: " << solver.iterations() << endl;
				//cout << "CG Solver Error: " << solver.error() << endl << endl;
			}
			else
			{
				PardisoLDLT<SparseMatrix<double>> solver;
				solver.pardisoParameterArray()[59] = 1;
				solver.analyzePattern(J);
				solver.factorize(J);
				epsilon = solver.solve(residual); 

			}

			//duration = (clock() - start) / (double)CLOCKS_PER_SEC;
			//t = t + duration;

			//cout << "ET Solver: " << duration << endl << "Error: " << err<< endl;

			if (iter > 4 && iter < 90)
			{
				damp_fact = 1.0;
			}
			else if (iter < 150)
			{
				damp_fact = 0.25;
			}
			else
			{
				damp_fact = 0.01;
			}

			for (int k = 0; k < num_nodes; k++)
			{
				
				/*
				// updating nodal positions from GMRES solution

				if (epsilon(3 * k + 0) > move_thresh)
				{
					epsilon(3 * k + 0) = move_thresh;
				}
				if (epsilon(3 * k + 1) > move_thresh)
				{
					epsilon(3 * k + 1) = move_thresh;
				}
				if (epsilon(3 * k + 2) > move_thresh)
				{
					epsilon(3 * k + 2) = move_thresh;
				}

				if (epsilon(3 * k + 0) < -move_thresh)
				{
					epsilon(3 * k + 0) = -move_thresh;
				}
				if (epsilon(3 * k + 1) < -move_thresh)
				{
					epsilon(3 * k + 1) = -move_thresh;
				}
				if (epsilon(3 * k + 2) < -move_thresh)
				{
					epsilon(3 * k + 2) = -move_thresh;
				}
				*/

				nodes_n(k, 0) = nodes_n(k, 0) + damp_fact*epsilon(3 * k + 0);
				nodes_n(k, 1) = nodes_n(k, 1) + damp_fact*epsilon(3 * k + 1);	
				nodes_n(k, 2) = nodes_n(k, 2) + damp_fact*epsilon(3 * k + 2);

				// writing in the nodes mapped back to the undeformed unit cube
				nodes_mapped(k, 0) = F_inv(0, 0)*nodes_n(k, 0) + F_inv(0, 1)*nodes_n(k, 1) + F_inv(0, 2)*nodes_n(k, 2);
					//+ (T2(0, 0) + T2(0, 1) + T2(0, 2));;
					nodes_mapped(k, 1) = F_inv(1, 0)*nodes_n(k, 0) + F_inv(1, 1)*nodes_n(k, 1) + F_inv(1, 2)*nodes_n(k, 2);
					//+ (T2(1, 0) + T2(1, 1) + T2(1, 2));
					nodes_mapped(k, 2) = F_inv(2, 0)*nodes_n(k, 0) + F_inv(2, 1)*nodes_n(k, 1) + F_inv(2, 2)*nodes_n(k, 2);
					//+ (T2(2, 0) + T2(2, 1) + T2(2, 2));
			}

			updateCrossing(nodes_mapped, nodes_n, fibers_n, num_nodes, num_fibers, F); // update crossing conditions

			network_n.nodes = nodes_n;
			network_n.fibers = fibers_n;
			network0.nodes = nodes_mapped;
			network0.fibers = fibers_n;
			/*
			if (iter == iter_max - 1)
			{
					nodes_n = revert_nodes;
					fibers_n = revert_fibs;

					for (int k = 0; k < num_nodes; k++)
					{
						// writing in the nodes mapped back to the undeformed unit cube
						nodes_mapped(k, 0) = F_inv(0, 0)*nodes_n(k, 0) + F_inv(0, 1)*nodes_n(k, 1) + F_inv(0, 2)*nodes_n(k, 2);
						//+ (T2(0, 0) + T2(0, 1) + T2(0, 2));;
						nodes_mapped(k, 1) = F_inv(1, 0)*nodes_n(k, 0) + F_inv(1, 1)*nodes_n(k, 1) + F_inv(1, 2)*nodes_n(k, 2);
						//+ (T2(1, 0) + T2(1, 1) + T2(1, 2));
						nodes_mapped(k, 2) = F_inv(2, 0)*nodes_n(k, 0) + F_inv(2, 1)*nodes_n(k, 1) + F_inv(2, 2)*nodes_n(k, 2);
						//+ (T2(2, 0) + T2(2, 1) + T2(2, 2));
					}

					updateCrossing(nodes_mapped, nodes_n, fibers_n, num_nodes, num_fibers, F); // update crossing conditions

					network_n.nodes = nodes_n;
					network_n.fibers = fibers_n;
					network0.nodes = nodes_mapped;
					network0.fibers = fibers_n;
	
				flag = 1;
			}*/

		}
	}

	//dur_loop = (clock() - start_loop) / (double)CLOCKS_PER_SEC;

	//cout << "Total Loop Time: " << dur_loop << endl << "Iterations:" << its << endl << "Time per Loop:" << t / its << endl <<  endl;
	/*
	MatrixXd Jac = MatrixXd(J); // convert sparse to dense for calculation of jacobian condition number

	BDCSVD<MatrixXd> svd(Jac); // calculating Jacobian condition number
	cond = svd.singularValues()(0)
		/ svd.singularValues()(svd.singularValues().size() - 1);

	cout << "Internal Condition #: " << cond/1e7 << " x 1E7" << endl << endl;
	*/

	// Final Force Calculation
	calcForcesPeriodic(nodes_n, nodes_mapped, fibers_n, network_n, F,
		fiberResults_n, netResults_n, num_fibers, num_nodes, WSS);
	
	componentStress(nodes_n, fibers_n, network_n, netResults_n, fiberResults_n, F, num_fibers, x_scale);

	netResults_n.Jacobian = MatrixXd(J);
	
	if (flag == 1)
	{
		netResults_n.net_flag = 1;
		cout << "Failed to Converge!" << endl << "Norm of Residuals: " << err << endl << "Maximum Force: " << fiberResults_n.fib_forces.maxCoeff()<<"  " << fiberResults_n.fib_forces.minCoeff() << endl << endl;
		cout << "lambdas: " << fiberResults_n.fib_stretch.maxCoeff() << "  " << fiberResults_n.fib_stretch.minCoeff() << endl << endl;
		cout << "Crossings: " << fibers_n.col(2).maxCoeff() << "   " << fibers_n.col(3).maxCoeff() << "  "<< fibers_n.col(4).maxCoeff()  << endl << endl;
		cout << "Crossings: " << fibers_n.col(2).minCoeff() << "   " << fibers_n.col(3).minCoeff() << "  " << fibers_n.col(4).minCoeff() << endl << endl;
	}
	/*if (flag == 0)
	{
		cout << "Converged!" << endl << "Norm of Residuals: " << err << endl << "Maximum Force: " << fiberResults_n.fib_forces.maxCoeff() << "  " << fiberResults_n.fib_forces.minCoeff() << endl << endl;
		cout << "Lambdas: " << fiberResults_n.fib_stretch.maxCoeff() << "  " << fiberResults_n.fib_stretch.minCoeff() << endl << endl;
	}*/
}
