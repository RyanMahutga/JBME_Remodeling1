
/* calc_net_stress_periodic.m calculates the volume averaged network stress for a periodic network

SUMMARY: This function calculates the volume averaged stress tensor from the boundary nodes of a
periodic network. It is dependent of the boundary node coordinates 'bnd_nodes', the actual fiber
volume fraction 'fiber_vol_fract' (typically 10-20% fiber, 80% water/goop), and the fiber volume
in the RVE 'real_fib_vol'.

CALLED FROM: networkSolver.cpp

CALLS ON: None

INPUTS: nodes - Nx3 array of interior node coordinates(xyz)
bnd_node_nums - Bx3 array of boundary node coordinates(xyz)
fibers - Mx5 matrix of nodal coordinates(columns 1 and 2) and
crossing data(x, y, z in columns 3, 4, 5)
fib_forces - Mx1 array of fiber forces between nodes

OUTPUTS: None

NOTE: This function modifies net_stress the 3x3 array for the volume averaged stress tensor

Created by : Ryan Mahutga - Barocas Research Group - University of Minnesota
Date Modified : 04-23-18

*/

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"

using namespace std;
using namespace Eigen;

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, netResults& netResults_n, fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale )
{

	MatrixXd comp_stress(12, 3);

	VectorXi fib_type = network_n.fib_type;
	MatrixXd bnd_nodes = network_n.bnd_nodes;
	VectorXd fib_forces = fiberResults_n.fib_forces;
	VectorXd init_lens = network_n.init_lens;
	Vector3d unit, vect, realNode2; 

	double V = F.determinant();

	Matrix3d sum0, sum1, sum2, sum3, sum_tot, cStress = Matrix3d::Zero(3, 3);
	sum0 = sum1 = sum2 = sum3 = sum_tot = Matrix3d::Zero();

	int node1, node2;
	double vect_len = 0.0;

	for (int f = 0; f < num_fibs; f++)
	{
		node1 = (int)(fibers_n(f, 0) + 0.3);
		node2 = (int)(fibers_n(f, 1) + 0.3);

		//Generating fiber vectors
		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes_n(node2, 0) + fibers_n(f, 2)*F(0, 0) + fibers_n(f, 3)*F(0, 1) + fibers_n(f, 4)*F(0, 2);
		realNode2(1) = nodes_n(node2, 1) + fibers_n(f, 2)*F(1, 0) + fibers_n(f, 3)*F(1, 1) + fibers_n(f, 4)*F(1, 2);
		realNode2(2) = nodes_n(node2, 2) + fibers_n(f, 2)*F(2, 0) + fibers_n(f, 3)*F(2, 1) + fibers_n(f, 4)*F(2, 2);

		vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (fib_type(f) == 1)
				{
					sum1(i, j) = sum1(i,j) + fib_forces(f)*init_lens(f)*unit(i)*unit(j);
				}
				else if (fib_type(f) == 2)
				{
					sum2(i, j) = sum2(i, j) + fib_forces(f)*init_lens(f)*unit(i)*unit(j);
				}
				else if (fib_type(f) == 3)
				{
					sum3(i, j) = sum3(i, j) + fib_forces(f)*init_lens(f)*unit(i)*unit(j);
				}
				else if (fib_type(f) == 0)
				{
					sum0(i, j) = sum0(i, j) + fib_forces(f)*init_lens(f)*unit(i)*unit(j);
				}
			}
		}
	}

		sum_tot = sum0 + sum1 + sum2 + sum3;

		netResults_n.net_stress = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum_tot;

		comp_stress.block(0, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum1;
		comp_stress.block(3, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum2;
		comp_stress.block(6, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum3;
		comp_stress.block(9, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum0;

		netResults_n.component_stress = comp_stress;
		
}
