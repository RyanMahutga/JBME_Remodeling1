/*
writeData.cpp is a dumb code that writes things to files

SUMMARY: This function writes data to files with the names of the variables and the integer t appended.

CALLED FROM: Main.cpp
CALLS ON: None

INPUTS:fibers - Mx5 matrix of fiber connectivity, fibers (m,1) and fibers(m,2) represent the number of the two nodes that are connected,
		fibers(m,3:5) are the x,y,z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
		the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
	fib_areas - Mx1 vector of fiber cross sectional areas (in m)
	init_lens - Mx1 vector of fiber initial lengths
	RVE_stretch - 3x1 vector of principal network stretches
	fib_forces - Mx1 vector of fiber forces 
	fib_stress - Mx1 vector of fiber stresses
	bnd_nodes - Px3 coordinates of boundary nodes
	net_stress - 3x3 tensor of volume averaged network stresses
	t - int used to delineate written files from one another (i.e. time steps/remodelling steps)

OUTPUTS: None

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-07-2018

fiberSave = cols 
			1 & 2 node numbers
			3,4,5 connectivity (through what bounds (x,y,z +or-)
			6 - fib type
			7,8 - fib radius and area
			9 - fib length
			10,11,12 - fib stretch, force, and stress
*/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"
#include <vector>

using namespace std;
using namespace Eigen;

void writeData(int num_fibers, network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, int t)
{
		// node save results
		int num_bnd = network_n.num_bnd;
		int num_nodes = network_n.num_nodes;

		MatrixXd bnd_nodes = network_n.bnd_nodes;
		MatrixXd nodes = network_n.nodes;

		MatrixXd netSave(4 + num_nodes + num_bnd, 3);

		Matrix3d F = netResults_n.F;

		netSave.block(0,0,3,3) = F;

		netSave(3, 0) = num_nodes; netSave(3, 1) = num_bnd; // numer of interior and bnd nodes
		netSave(3, 2) = 0; // just a place holder

		netSave.block(4, 0, num_nodes, 3) = nodes; // store interior nodes
		//cout << netSave << endl << endl;
		netSave.block(num_nodes+4, 0, num_bnd, 3) = bnd_nodes; // store boundary nodes

		// fiber save results
		MatrixXd fiberSave(num_fibers, 12);

		MatrixXd fibers = network_n.fibers;

		VectorXi fib_type = network_n.fib_type;

		VectorXd fibType = fib_type.cast<double>();

		fiberSave.topLeftCorner(num_fibers, 5) = fibers;

		fiberSave.col(5) = fibType;
		fiberSave.col(6) = network_n.fib_rads;
		fiberSave.col(7) = network_n.fib_areas;
		fiberSave.col(8) = network_n.init_lens;
		fiberSave.col(9) = fiberResults_n.fib_stretch;
		fiberSave.col(10) = fiberResults_n.fib_forces;
		fiberSave.col(11) = fiberResults_n.fib_stress;


			// saving fiber results
			ofstream myfile1;
			string str1 = "fiberData" + to_string(t) + ".txt";
			myfile1.open(str1);
			myfile1 << fiberSave;
			myfile1.close();

			// saving network node results
			ofstream myfile2;
			string str2 = "nodeData" + to_string(t) + ".txt";
			myfile2.open(str2);
			myfile2 << netSave;
			myfile2.close();

}