
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

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale, MatrixXd& comp_stress )
{
	VectorXi fib_type = network_n.fib_type;
	MatrixXd bnd_nodes = network_n.bnd_nodes;
	VectorXd fib_forces = fiberResults_n.fib_forces;

	VectorXi num(3); num = VectorXi::Zero(3);
	VectorXi bnd(3); bnd = VectorXi::Zero(3);
	// calculate number of each fiber type
	for (int i = 0; i < num_fibs; i++)
	{
		num(fib_type(i) - 1)= num(fib_type(i) - 1)+1;

		for (int j = 2; j < 5; j++)
		{
			for (int bt = 0; bt < int(abs(fibers_n(i, j))+0.3); bt++)
			{
				bnd(fib_type(i) - 1) = bnd(fib_type(i) - 1) + 2;
			}
		}
	}

	MatrixXd fibers1(num(0), 5);
	MatrixXd fibers2(num(1), 5);
	MatrixXd fibers3(num(2), 5);

	VectorXd ffib1(num(0)); ffib1 = VectorXd::Zero(num(0));
	VectorXd ffib2(num(1)); ffib2 = VectorXd::Zero(num(1));
	VectorXd ffib3(num(2)); ffib3 = VectorXd::Zero(num(2));

	MatrixXd bnodes1(bnd(0), 3); bnodes1 = MatrixXd::Zero(bnd(0), 3);
	MatrixXd bnodes2(bnd(1), 3); bnodes2 = MatrixXd::Zero(bnd(1), 3);
	MatrixXd bnodes3(bnd(2), 3); bnodes3 = MatrixXd::Zero(bnd(2), 3);

	VectorXi count(3); count = VectorXi::Zero(3);

	Matrix3d cStress = Matrix3d::Zero(3, 3);

	int c1 = 0, c2 = 0, c3 = 0, b=0;
	for (int f = 0; f < num_fibs; f++)
	{
		if (fib_type(f) == 1)
		{
			fibers1(c1,0) = fibers_n(f,0);
			fibers1(c1, 1) = fibers_n(f, 1);
			fibers1(c1, 2) = fibers_n(f, 2);
			fibers1(c1, 3) = fibers_n(f, 3);
			fibers1(c1, 4) = fibers_n(f, 4);
			ffib1(c1) = fib_forces(f);
			for (int v = 2; v < 5; v++)
			{
				for (int p = 0; p < int(abs(fibers_n(f, v)) + 0.3); p++)
				{
					bnodes1(2 * count(0), 0) = bnd_nodes(2 * b, 0);
					bnodes1(2 * count(0), 1) = bnd_nodes(2 * b, 1);
					bnodes1(2 * count(0), 2) = bnd_nodes(2 * b, 2);

					bnodes1(2 * count(0) +1, 0) = bnd_nodes(2 * b+1, 0);
					bnodes1(2 * count(0) +1, 1) = bnd_nodes(2 * b+1, 1);
					bnodes1(2 * count(0) +1, 2) = bnd_nodes(2 * b+1, 2);
					b++;
					count(0)++;
				}
			}
			c1++;
		}
		else if (fib_type(f) == 2)
		{
			fibers2(c2,0) = fibers_n(f,0);
			fibers2(c2, 1) = fibers_n(f, 1);
			fibers2(c2, 2) = fibers_n(f, 2);
			fibers2(c2, 3) = fibers_n(f, 3);
			fibers2(c2, 4) = fibers_n(f, 4);
			ffib2(c2) = fib_forces(f);
			for (int v = 2; v < 5; v++)
			{
				for (int p = 0; p < int(abs(fibers_n(f, v))+0.3); p++)
				{
					bnodes2(2 * count(1), 0) = bnd_nodes(2 * b, 0);
					bnodes2(2 * count(1), 1) = bnd_nodes(2 * b, 1);
					bnodes2(2 * count(1), 2) = bnd_nodes(2 * b, 2);

					bnodes2(2 * count(1) + 1, 0) = bnd_nodes(2 * b + 1, 0);
					bnodes2(2 * count(1) + 1, 1) = bnd_nodes(2 * b + 1, 1);
					bnodes2(2 * count(1) + 1, 2) = bnd_nodes(2 * b + 1, 2);
					b++;
					count(1)++;
				}
			}
			c2++;
		}
		else if (fib_type(f) == 3)
		{
			fibers3(c3,0) = fibers_n(f,0);
			fibers3(c3, 1) = fibers_n(f, 1);
			fibers3(c3, 2) = fibers_n(f, 2);
			fibers3(c3, 3) = fibers_n(f, 3);
			fibers3(c3, 4) = fibers_n(f, 4);
			ffib3(c3) = fib_forces(f);
			for (int v = 2; v < 5; v++)
			{
				for (int p = 0; p < int(abs(fibers_n(f, v)) + 0.3); p++)
				{
					bnodes3(2 * count(2), 0) = bnd_nodes(2 * b, 0);
					bnodes3(2 * count(2), 1) = bnd_nodes(2 * b, 1);
					bnodes3(2 * count(2), 2) = bnd_nodes(2 * b, 2);

					bnodes3(2 * count(2) + 1, 0) = bnd_nodes(2 * b + 1, 0);
					bnodes3(2 * count(2) + 1, 1) = bnd_nodes(2 * b + 1, 1);
					bnodes3(2 * count(2) + 1, 2) = bnd_nodes(2 * b + 1, 2);
					b++;
					count(2)++;
				}
			}
			c3++;
		}
	}

	for (int s = 0; s < 3; s++)
	{
		if (s == 0)
		{
			calcNetStress2(nodes_n, fibers1, ffib1, 
				cStress, bnodes1, F, num(0), x_scale);
		}
		else if (s == 1)
		{
			calcNetStress2(nodes_n, fibers2, ffib2, 
				cStress, bnodes2, F, num(1), x_scale);
		}
		else if (s == 2)
		{
			calcNetStress2(nodes_n, fibers3, ffib3, 
				cStress, bnodes3, F, num(2), x_scale);
		}
		else
		{
			cout << "Error in Code!";
			cin.get();
			exit(4);
		}

		comp_stress(s * 3, 0) = cStress(0, 0);
		comp_stress(s * 3+1, 0) = cStress(1, 0);
		comp_stress(s * 3+2, 0) = cStress(2, 0);
		comp_stress(s * 3, 1) = cStress(0, 1);
		comp_stress(s * 3+1, 1) = cStress(1, 1);
		comp_stress(s * 3+2, 1) = cStress(2, 1);
		comp_stress(s * 3, 2) = cStress(0, 2);
		comp_stress(s * 3+1, 2) = cStress(1, 2);
		comp_stress(s * 3+2, 2) = cStress(2, 2);
	}

}