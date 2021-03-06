
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

void calcNetStress2(MatrixXd& nodes_n, MatrixXd& fibers_n, VectorXd& fib_forces, Matrix3d& net_stress, MatrixXd& bnd_nodes, Matrix3d F, int num_fibs, double x_scale)
{
	// Finding Fiber Crossing Forces
	int node1 = 0;
	int node2 = 0;

	double fiber_current_length = 0;
	Vector3d realnode2 = Vector3d::Zero(3); // real node2 coordinate
	Vector3d vect = Vector3d::Zero(3); // vector from node1 to node2
	Vector3d unit = Vector3d::Zero(3); // unit vector from node1 to node2

	int num_bnd = (int)(bnd_nodes.rows()+0.1); // number of boundary nodes

	MatrixXd bnd_coords(num_bnd, 3);  bnd_coords = MatrixXd::Zero(num_bnd, 3); // boundary coordinates (x,y,z)
	MatrixXd bnd_forces(num_bnd, 3);  bnd_forces = MatrixXd::Zero(num_bnd, 3); // boundary node forces (Fx,Fy,Fz)
	
	double new_volume = 0; // rve volume
	Matrix3d stress = Matrix3d::Zero(3, 3); // unscaled stress
	
	int m = 0; // counter
	int cross = 0;
	for (int n = 0; n < num_fibs; n++)
	{
		node1 = (int)(fibers_n(n, 0) + 0.3);
		node2 = (int)(fibers_n(n, 1) + 0.3);

		//Generating fiber vectors
		realnode2(0) = nodes_n(node2, 0) + fibers_n(n, 2)*F(0, 0) + fibers_n(n, 3)*F(0, 1) + fibers_n(n, 4)*F(0, 2);
		realnode2(1) = nodes_n(node2, 1) + fibers_n(n, 2)*F(1, 0) + fibers_n(n, 3)*F(1, 1) + fibers_n(n, 4)*F(1, 2);
		realnode2(2) = nodes_n(node2, 2) + fibers_n(n, 2)*F(2, 0) + fibers_n(n, 3)*F(2, 1) + fibers_n(n, 4)*F(2, 2);

		vect(0) = realnode2(0) - nodes_n(node1, 0);
		vect(1) = realnode2(1) - nodes_n(node1, 1);
		vect(2) = realnode2(2) - nodes_n(node1, 2);

		fiber_current_length = vect.norm();

		unit = vect / fiber_current_length;

		for (int i = 0; i < 3; i++)
		{
			cross = int(abs(fibers_n(n, 2 + i)) + 0.3);
			for (int k = 0; k < cross; k++)
			{
				bnd_coords(m, 0) = bnd_nodes(m, 0);
				bnd_coords(m, 1) = bnd_nodes(m, 1);
				bnd_coords(m, 2) = bnd_nodes(m, 2);
				bnd_forces(m, 0) = unit(0)*fib_forces(n);
				bnd_forces(m, 1) = unit(1)*fib_forces(n);
				bnd_forces(m, 2) = unit(2)*fib_forces(n);
				m++;

				bnd_coords(m, 0) = bnd_nodes(m, 0);
				bnd_coords(m, 1) = bnd_nodes(m, 1);
				bnd_coords(m, 2) = bnd_nodes(m, 2);
				bnd_forces(m, 0) = -unit(0)*fib_forces(n);
				bnd_forces(m, 1) = -unit(1)*fib_forces(n);
				bnd_forces(m, 2) = -unit(2)*fib_forces(n);
				m++;
			}
		}
	}

	double xfx = 0; double yfx = 0; double zfx = 0;
	double xfy = 0; double yfy = 0; double zfy = 0;
	double xfz = 0; double yfz = 0; double zfz = 0;

	// Sum xi * fi through matrix multiplication
	for (int k = 0; k < num_bnd; k++)
	{
		xfx = xfx + (bnd_coords(k, 0) * bnd_forces(k,0));
			yfx = yfx + (bnd_coords(k, 1) * bnd_forces(k,0));
				zfx = zfx + (bnd_coords(k, 2) * bnd_forces(k,0));

					xfy = xfy +(bnd_coords(k, 0) * bnd_forces(k,1));
						yfy = yfy + (bnd_coords(k, 1) * bnd_forces(k,1));
							zfy = zfy + (bnd_coords(k, 2) * bnd_forces(k,1));

								xfz = xfz + (bnd_coords(k, 0) * bnd_forces(k,2));
									yfz = yfz + (bnd_coords(k, 1) * bnd_forces(k,2));
										zfz = zfz + (bnd_coords(k, 2) * bnd_forces(k,2));
	}

	stress << xfx, yfx, zfx,
		xfy, yfy, zfy,
		xfz, yfz, zfz;

	/*
	 Formula for this is Sij = 1 / V * Sum(xiFj); The off diagnonal terms are
	 symmetric, e.g.xfy = yfx.However, there is a very small difference in
	 these numbers due to numerical issues, so we just average them.The
	 difference seems to be at the 10 - 14 decimal place, so probably could
	 ignore.

	 compact notation[S11 S12 S13 S22 S23 S33]
	 stress = [xfx; 0.5*(xfy + yfx); 0.5*(zfx + xfz); yfy; 0.5*(yfz + zfy); zfz];

	 full notation[S11 S12 S13; S21 S22 S23; S31 S32 S33]
	 stress = [xfx 0.5*(xfy + yfx) 0.5*(zfx + xfz);
			   0.5*(xfy + yfx) yfy 0.5*(yfz + zfy);
			   0.5*(zfx + xfz) 0.5*(yfz + zfy) zfz];
	*/

	// calculate individual fiber stresses

	// calculate stress scaling factor

	new_volume = F.determinant(); // assuming the rve is initially 1x1x1 

	net_stress =  1/(x_scale*x_scale)*(1/new_volume * stress); // rescaled network stress

}