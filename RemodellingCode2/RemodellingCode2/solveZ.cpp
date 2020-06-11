#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;


void solveZ(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, double& err, double P)
{
	double R = 500e-6, H = 50e-6; // aorta radius and wall thickness unloaded

	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, guess, stab_node);

	net_stress = netResults_n.net_stress;
	
	double pressure = net_stress(2, 2);

	netResults_n.pressure = pressure;

	net_stress(0, 0) = net_stress(0, 0) - pressure;
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	// use law of laplace to find bcs
	appliedStress(0) = P * (R*F(0, 0)) / (H*F(2, 2));

	err = (net_stress(0, 0)- appliedStress(0));

}