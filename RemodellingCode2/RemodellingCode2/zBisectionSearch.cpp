/*
zBisectionSearch.cpp : Solves the stretch of the free surface to get a 'zero' stress condition.

SUMMARY: This function solves the z-direction stretch that gives <1% of the maximum stress on the free surface. It then updates
F, net_stress, fib_forces, fib_stress, nodes_n, fibers, and bnd_nodes.

CALLED FROM: main.cpp

CALLS ON: calcForcesPeriodic.cpp, calcJacobian.cpp, calcBndNodes.cpp, calcNetStress.cpp, networkSolver.cpp

INPUTS: nodes - Nx3 matrix of nodal coordinates in an RVE
		fibers - Mx5 matrix of fiber connectivity, fibers (m,1) and fibers(m,2) represent the number of the two nodes that are connected,
			fibers(m,3:5) are the x,y,z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
			the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
		fib_type - Mx1 vector of fiber type (1 for collagen, 2 for elastin, 3 for ?)
		fib_rads - Mx1 vector of fiber radii (in meters)
		init_lens - Mx1 vector of fiber initial lengths
		fiber_vol_fract - double value or actual tissue fiber volume fraction (typically 10-20% fiber, 80% water/goop)
		F - 3x1 vector of principal network stretches
		num_nodes - int N  number of nodes
		num_fibers - int M number of fibers
		net_stress - 3x3 tensor of volume averaged stress

OUTPUTS: None

NOTE: This code modifies 'nodes_n', bnd_nodes','fibers_n', 'net_stress', 'fib_forces', 'fib_stress', and 'F'

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-07-2018

*/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

void zBisectionSearch(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess, double fib_vol_m, double phi_m, double P)
{

	double err = 1e6;

	F(0, 0) = x_guess(0);
	F(1, 1) = 1;
	F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

	solveZ(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

	bool flag_high = 0;
	bool flag_low = 0;
	bool flag_bisect = 0;

	double high_val = 0;
	double low_val = 0;

	double err_thresh = 1e-2;

	for (int j = 0; j < 500; j++)
	{
		err_thresh = 1; //0.0001*(net_stress.maxCoeff());

		if (abs(err) < err_thresh && flag_bisect==1)
		{
			break;
		}

		if (flag_high == 1 && flag_low == 1)
		{
			if (flag_bisect == 1)
			{
				if (err >= 0)
				{
					high_val = F(0,0);
				}
				else
				{
					low_val = F(0,0);
				}
			}
			F(0,0) = (high_val + low_val) / 2;
			flag_bisect = 1;
		}
		else if (err >= 0)
		{
			flag_high = 1;
			high_val = F(0,0);
			F(0,0) = F(0,0) - 0.1;
		}
		else if (err < 0)
		{
			flag_low = 1;
			low_val = F(0,0);
			F(0,0) = F(0,0) + 0.1;
		}

		F(2, 2) = fib_vol_m / (phi_m*F(0, 0)*F(1, 1)*x_scale*x_scale);

		solveZ(num_nodes, num_fibers, x_scale, F, network0,
			network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress, net_stress, err, P);

		//cout << "Bisect Iteration: " << j + 1 << endl;
		cout<<"Bisection Error:" <<err<<"   Z Stretch: " << F(2,2)<<endl<< endl;

	}

}