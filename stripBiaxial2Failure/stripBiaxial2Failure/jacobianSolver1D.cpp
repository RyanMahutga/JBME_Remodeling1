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
#include "dataStructureDefinitions.h"
#include <vector>

using namespace std;
using namespace Eigen;

// this has been modified to do 2D with constant phi 
void jacobianSolver1D(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P, double& pressure, double& WSS)
{
	double dlam = 1e-6, mag_err = 1000, err_thresh = 50;
	Matrix3d dF;
	Vector3d init_err, err;
	double epsilon, J;

	for (int n = 0; n < 100; n++)
	{
		J = 0.0;

		// initial solve
		solveJac(num_nodes, num_fibers, x_scale, F, network0,
			network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, init_err, P, fib_vol_m, WSS);

		//cout << init_err << endl << pressure << endl<<endl;

		mag_err = abs(init_err(2));

		//cout << "Error: " << mag_err << endl;

		if (mag_err < err_thresh)
		{
			cout << "Converged! Error: " << mag_err << endl << endl;
			break;
		}

		dF = F;
		dF(2,2) = F(2, 2) + dlam;

		solveJac(num_nodes, num_fibers, x_scale, dF, network0,
			network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P, fib_vol_m, WSS);

		J = (err(2) - init_err(2)) / dlam;

		epsilon = -1.0*init_err(2)/J;

		F(2, 2) += epsilon;

	}
	cout << "Deformation: " << endl << F << endl << "Wall Shear Stress: " << WSS << endl << endl;

	netResults_n.F = F;
}