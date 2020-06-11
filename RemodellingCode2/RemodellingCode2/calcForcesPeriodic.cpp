/* calcForcesPeriodic.cpp calculate nodal and fiber forces for the periodic RVE

SUMMARY: This function calculates the fiber and nodal forces in the network given 
the current nodal positions. This function is used within the Newton loop to 
equilibrate the network, and to calculate the stress.

CALLED FROM: networkSolver.cpp & calcJacobian.cpp
CALLS ON: updatCrossing.cpp & fiberConstEqn.cpp

INPUTS: nodes_n - Nx3 matrix of nodal coordinates
		nodes0 - Nx3 matrix of nodal coordinates for tethering
		fibers - Mx5 fiber matrix relating fiber end nodes to one another
		init_lens - Mx1 vector of fiber initial lengths
		fib_type - Mx1 vector of fiber types
		fib_areas - Mx1 vector of fiber cross-sectional areas
		F - 3x3 deformation matrix of RVE 
		fib_forces - Mx1 vector of fiber forces (initially empty)
		node_forces - Nx3 matrix of x,y,z components of forces on nodes (initially empty)
		num_fibers - M integer number of fibers
		num_nodes - N integer number of nodes
		bnd_nodes - Px3 coordinates of boundary nodes

OUTPUTS: None

NOTE: nodes_n, fibers, fib_forces, node_forces, and bnd_nodes are modified by this code.

DEPENDS ON: updateCrossing.cpp, fiberConstEqn.cpp

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-04-2018

*/

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"

using namespace std;
using namespace Eigen;

void calcForcesPeriodic(MatrixXd nodes_n, MatrixXd nodes_mapped, MatrixXd fibers_n, network& network_n,
	Matrix3d F, fiberResults& fiberResults_n, netResults& netResults_n, int num_fibers, int num_nodes, double WSS)
{

	VectorXd init_lens = network_n.init_lens;
	VectorXi fib_type = network_n.fib_type;
	VectorXd fib_areas = network_n.fib_areas;
	VectorXd fib_stretch(num_fibers), fib_stress(num_fibers), fib_forces(num_fibers);

	VectorXd node_forces(3*num_nodes); node_forces = VectorXd::Zero(3*num_nodes);

	int node1 = 0;
	int node2 = 0;
	double lambda = 0;
	double fib_force = 0;
	double dFdlam = 0;

	Vector3d node_force1(3); node_force1 = Vector3d::Zero(3);
	Vector3d node_force2(3); node_force2 = Vector3d::Zero(3);
	Vector3d vect(3); vect = Vector3d::Zero(3);
	Vector3d realNode2(3); realNode2 = VectorXd::Zero(3); // real node 2 position

	double vect_len = 0;
	VectorXd unit(3); unit = VectorXd::Zero(3); // unit vector from node1 to node2

	for (int n = 0; n<num_fibers; n++) // loop through fibers to calculate forces
	{
		node1 = (int)(fibers_n(n, 0) + 0.1); // first (static) node of fiber 
		node2 = (int)(fibers_n(n, 1) + 0.1); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,2:4))

		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes_n(node2, 0) + fibers_n(n, 2)*F(0, 0) + fibers_n(n, 3)*F(0, 1) + fibers_n(n, 4)*F(0, 2);
		realNode2(1) = nodes_n(node2, 1) + fibers_n(n, 2)*F(1, 0) + fibers_n(n, 3)*F(1, 1) + fibers_n(n, 4)*F(1, 2);
		realNode2(2) = nodes_n(node2, 2) + fibers_n(n, 2)*F(2, 0) + fibers_n(n, 3)*F(2, 1) + fibers_n(n, 4)*F(2, 2);

		vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		//cout << vect_len << "  " << init_lens(n) << endl;

		fiberConstEqn(vect, node_force1, node_force2, lambda, init_lens(n), fib_type(n), fib_areas(n), fib_force, dFdlam, WSS);

		fib_stretch(n) = lambda; // fiber stretch
		fib_forces(n) = fib_force; // store fiber force
		fib_stress(n) = fib_force / fib_areas(n); // fiber PK1 stress

		// storing nodal forces
		node_forces(3 * node1 + 0) = node_forces(3 * node1 + 0) + node_force1(0);
		node_forces(3 * node1 + 1) = node_forces(3 * node1 + 1) + node_force1(1);
		node_forces(3 * node1 + 2) = node_forces(3 * node1 + 2) + node_force1(2);

		node_forces(3 * node2 + 0) = node_forces(3 * node2 + 0) + node_force2(0);
		node_forces(3 * node2 + 1) = node_forces(3 * node2 + 1) + node_force2(1);
		node_forces(3 * node2 + 2) = node_forces(3 * node2 + 2) + node_force2(2);

	}
	// store fiber results
	fiberResults_n.fib_forces = fib_forces; 
	fiberResults_n.fib_stretch = fib_stretch;
	fiberResults_n.fib_stress = fib_stress;
	netResults_n.node_forces = node_forces;
}