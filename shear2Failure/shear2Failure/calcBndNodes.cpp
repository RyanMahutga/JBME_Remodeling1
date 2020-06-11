/* calcBndNodes.cpp updates the fiber vector according to how nodes enter and exit the RVE and finds the boundary nodes

SUMMARY: This function finds boundary nodes based on fibers that cross the RVE boundary. What is does is find the vector
from node1 to node2 along a fiber, then determine what boundaries it crosses. It then orders the boundary crossings using
a bubble sort algorithm to determine where the boundary nodes are in successive order (i.e. if it crosses the x and z 
boundaries, but crosses z first then x this algorithm accounts for that).

CALLED FROM: networkSolver.cpp

CALLS ON: updateCrossing.cpp

INPUTS: nodes_n - Nx3 matrix of nodal coordinates
fiber - Mx5 matrix of fiber connectivity
num_fibers - M number of fibers
rve_stretch - 3x1 vector of primary RVE stretches

OUTPUTS: None

NOTE: This function does modify the values for nodes_n and fibers in main.

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota
Contact: mahut005@umn.edu

Last Update: 04-02-2018

*/

#include <iostream>
#include "Eigen/Eigen"
#include <complex>
#include <cmath>
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void sortFunction(VectorXd& dL, VectorXi& dir, int& num_cross, int& leng) // this is a bubblesort algorithm
{
	double val_holder; // holder for dl value
	int idx_holder; // holder for crossing index 

	for (int i = 0; i < num_cross - 1; i++)
	{
		leng++;

		if (dL(i + 1) == 0) // break if no values are trailing
		{
			break;
		}

		for (int j = 0; j < num_cross - 1; j++)
		{
			if (dL(j + 1) == 0) // break if no values are trailing
			{
				break;
			}

			if (dL(j) > dL(j + 1))
			{
				val_holder = dL(j);
				idx_holder = dir(j);
				dL(j) = dL(j + 1);
				dL(j + 1) = val_holder;
				dir(j) = dir(j + 1);
				dir(j + 1) = idx_holder;
			}
		}
	}
}

void calcBndNodes(MatrixXd nodes_mapped, MatrixXd nodes_n, MatrixXd fibers_n, int num_nodes, int num_fibers, Matrix3d F, MatrixXd& bnd_nodes)
{
	// defining local values
	double vect_len = 0; // fiber vector length
	int num_cross = 0; // counter for number of fiber crossings
	int node1 = 0; // static node in RVE on fiber from node1 to node2
	int node2 = 0; // dynamic fiber node (i.e. node that moves to be the nearest neighbor to node1 along the fiber
	int leng = 0; // length of boundary crossings

	Matrix3d T(3, 3); // Translation matrix
	Matrix3d I(3, 3); I << 1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0;

	T = 0.5*(F-I); // rigid translation to move corner to -0.5*(1 1 1)

	 // defining local data structures
	MatrixXd nodes2 = nodes_mapped; // unchanging node positions
	Vector3d realNode2(3); realNode2 = Vector3d::Zero(3); // real node2 corrdinates to nearest connected neighbor to node1
	Vector3d unit(3); unit = Vector3d::Zero(3); // fiber unit vector
	Vector3d vect; vect = Vector3d::Zero(3); // fiber vector

	// finding boundaries
	double up_bnd = 0.5 ; // upper boundaries of the RVE
	double lo_bnd = -0.5 ; // lower bound in all cases

	for (int n = 0; n < num_fibers; n++)
	{
		num_cross = num_cross + int((abs(fibers_n(n, 2)) + abs(fibers_n(n, 3)) + abs(fibers_n(n, 4))+0.3));
	}

	// initializing values and data structures for boundary node calculation
	int num_bnd_nodes = 2 * num_cross; // boundary nodes are two per boundary crossing
	VectorXd dL(num_cross); dL = VectorXd::Zero(num_cross); // vector of fiber parts from node to boundary (undeformed coordinates)
	VectorXi dir(num_cross); dir = VectorXi::Zero(num_cross); // direction to boundary
	
	int m = 0;

	// boundary node coordinates
	Vector3d bnd_node1(3); bnd_node1 = Vector3d::Zero(3);
	Vector3d bnd_node2(3); bnd_node2 = Vector3d::Zero(3);

	// resize bounday node matrix
	bnd_nodes.resize(num_bnd_nodes, 3); bnd_nodes = MatrixXd::Zero(num_bnd_nodes, 3);

	int bnd = 0;
	// looping fibers to find boundary nodes
	for (int n = 0; n < num_fibers; n++)
	{
		node1 = (int)(fibers_n(n, 0) + 0.3); // first (static) node of fiber 
		node2 = (int)(fibers_n(n, 1) + 0.3); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,3:5))

		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes_n(node2, 0) + fibers_n(n, 2)*F(0, 0) + fibers_n(n, 3)*F(0, 1) + fibers_n(n, 4)*F(0, 2);
		realNode2(1) = nodes_n(node2, 1) + fibers_n(n, 2)*F(1, 0) + fibers_n(n, 3)*F(1, 1) + fibers_n(n, 4)*F(1, 2);
		realNode2(2) = nodes_n(node2, 2) + fibers_n(n, 2)*F(2, 0) + fibers_n(n, 3)*F(2, 1) + fibers_n(n, 4)*F(2, 2);

		vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		// clearing tracking values
		m = 0;
		dL = VectorXd::Zero(num_cross);
		dir = VectorXi::Zero(num_cross);

		for (int i = 0; i < 3; i++) // i==0 is x, i==1 is y i==2 is z direction
		{
			for (int c = 0; c < (int)(abs(fibers_n(n, 2 + i)) + 0.3); c++)
			{
				if (fibers_n(n, 2 + i) > 0) // decides direction of fiber
				{
					dL(m) = abs(((double(c) + 0.5) - nodes_mapped(node1, i)) / unit(i)); // length of fiber part to intersection
					dir(m) = i; // direction of fiber intersection
					m++;
				}
				else if (fibers_n(n, 2 + i) < 0)
				{
					dL(m) = abs(((-0.5 - double(c)) - nodes_mapped(node1, i)) / unit(i)); // length of fiber part to intersection
					dir(m) = i; // direction of fiber intersection
					m++;
				}
			}
		}

		if (dL(0) != 0) // if the matrix isn't empty, sort it, then find the boundary nodes
		{
			leng = 0; // reset vector length

			sortFunction(dL, dir, num_cross, leng); // sort fiber RVE intersections by length

			for (int j = 0; j < leng; j++)
			{
				// calculating boundary node
				bnd_node1(0) = dL(j)*unit(0) + nodes_mapped(node1,0);
				bnd_node1(1) = dL(j)*unit(1) + nodes_mapped(node1, 1);
				bnd_node1(2) = dL(j)*unit(2) + nodes_mapped(node1, 2);

				for (int i = 0; i < 3; i++) // used if there is more than one boundary crossing for a fiber
				{
					if (bnd_node1(i) > up_bnd)
					{
						bnd_node1(i) = bnd_node1(i) - 1.0;
					}
					else if (bnd_node1(i) < lo_bnd )
					{
						bnd_node1(i) = bnd_node1(i) + 1.0;
					}
				}

				// finding the boundary node on the opposite face from bnd_node1
				if (fibers_n(n, 2 + dir(j)) > 0)
				{
					bnd_node2 = bnd_node1;
					bnd_node2(dir(j)) = bnd_node2(dir(j)) - 1.0;
				}
				else if (fibers_n(n, 2 + dir(j)) < 0)
				{
					bnd_node2 = bnd_node1;
					bnd_node2(dir(j)) = bnd_node2(dir(j)) + 1.0;
				}

				// setting + boundary node (mapping from undeformed to deformed boundary)
				bnd_nodes(2 * bnd, 0) = F(0, 0)*bnd_node1(0) + F(0, 1)*bnd_node1(1) + F(0, 2)*bnd_node1(2) + T(0, 0) + T(0, 1) + T(0, 2);
				bnd_nodes(2 * bnd, 1) = F(1, 0)*bnd_node1(0) + F(1, 1)*bnd_node1(1) + F(1, 2)*bnd_node1(2) + T(1, 0) + T(1, 1) + T(1, 2);
				bnd_nodes(2 * bnd, 2) = F(2, 0)*bnd_node1(0) + F(2, 1)*bnd_node1(1) + F(2, 2)*bnd_node1(2) + T(2, 0) + T(2, 1) + T(2, 2);

				// setting - boundary node (mapping from undeformed to deformed boundary)
				bnd_nodes(2 * bnd + 1, 0) = F(0, 0)*bnd_node2(0) + F(0, 1)*bnd_node2(1) + F(0, 2)*bnd_node2(2) + T(0, 0) + T(0, 1) + T(0, 2);
				bnd_nodes(2 * bnd + 1, 1) = F(1, 0)*bnd_node2(0) + F(1, 1)*bnd_node2(1) + F(1, 2)*bnd_node2(2) + T(1, 0) + T(1, 1) + T(1, 2);
				bnd_nodes(2 * bnd + 1, 2) = F(2, 0)*bnd_node2(0) + F(2, 1)*bnd_node2(1) + F(2, 2)*bnd_node2(2) + T(2, 0) + T(2, 1) + T(2, 2);
				bnd++;
			}
		}
	}
}