/*


*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"


const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void writePeriodicNet2File(network &network_n, int num_fibers, int num_nodes, Vector3i stab_node, int t)
{
	// creating output file
	ofstream outFile;
	string str = "PeriodicNetwork_" + to_string(t) + ".txt";
	outFile.open(str); // open output file for editing

	// writing out first line (number of fibers, number of nodes, fiber volume fraction, stabilizing nodes)
	outFile << num_fibers << "  "<< num_nodes<<"  " << network_n.fiber_vol_fract <<"  " << stab_node(0) <<"  "<< stab_node(1) <<"  "<< stab_node(2)<<endl;

	// writing fiber data then node data
	MatrixXd fibers_n = network_n.fibers;
	VectorXi fib_type = network_n.fib_type;
	VectorXd init_lens = network_n.init_lens;
	VectorXd fib_rads = network_n.fib_rads;
	for (int p = 0; p < num_fibers; p++)
	{
		outFile << fibers_n(p,0)<<"  "<< fibers_n(p,1) << "  "<< fibers_n(p, 2)<<"  "<< fibers_n(p, 3) << "  " << fibers_n(p, 4) << "  " << fib_type(p) << "  " << init_lens(p) << "   " << fib_rads(p) << endl;
	}

	outFile << network_n.nodes;
}