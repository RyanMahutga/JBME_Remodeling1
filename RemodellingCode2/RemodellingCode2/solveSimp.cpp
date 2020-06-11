#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

const double PI = 3.141592653589793238462643383279502884;
const double R = 8.314; // universal gas constant

using namespace std;
using namespace Eigen;


void solveSimp(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& err, double P, double fib_vol_m)
{
	double R0 = 600e-6, H0 = 100e-6; // aorta radius and wall thickness unloaded

	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, guess, stab_node);

	net_stress = netResults_n.net_stress;

	double c_star = 150; // M - bathing solution concentration
	double T = 310; // temperature K
					//double psi = 0.5; // Ateshian et al. (2009) J. Biomech. Eng. 131
					//double c1 = 381; // Pa*m^3/Eq
					//double c2 = 0.241; // Pa*m^6/Eq^2
					//double c_fcd_0 = 100; // meq/L (roughly gives 0.5M at 0.5 fiber volume fraction see Porterfield paper)
					//double c_pos = 0.5 * (c_fcd + sqrt(c_fcd * c_fcd + 4 * c_star*c_star));
	double k = 280, k2=381, k3 = 0.241; //k gives cfcd give 60mEq/L at 30% actin, k2 is from Ateshian (2009) JBME 131

	double J = F.determinant();

	double phim = fib_vol_m / (x_scale*x_scale*J);

	double c_fcd = k * phim; //k*phi*sum_fib_len / fib_vol; //c_fcd_0 * phi; //c_fcd_0 * (1 - phi); // calculation of fixed charge density from fiber volume fractions
												   //c_neg = 0.5 * (-c_fcd + sqrt(c_fcd * c_fcd + 4 * c_star*c_star)); // negative ion charge density
	//cout << "Phi m: " << phi1<< "   " << "FCD: " << c_fcd << endl;

	double pressure = R * T * (sqrt(c_fcd*c_fcd + 4 * c_star*c_star) - 2 * c_star) + k2*phim*c_fcd + k3*phim*phim*c_fcd*c_fcd; 
	//0.5e6*(1 / ((c_fcd_0 - c_fcd)*(c_fcd_0 - c_fcd)) - 1 / (c_fcd_0 * c_fcd_0))+3000; // +c1 * psi*c_fcd + c2 * psi*psi*c_fcd*c_fcd;

	netResults_n.pressure = pressure;

	net_stress(0, 0) = net_stress(0, 0) - pressure;
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	// use law of laplace to find bcs
	appliedStress(0) = P * (R0*F(0, 0)) / (H0*F(2, 2));
	appliedStress(1) = 0.5* appliedStress(0);
	appliedStress(2) = 0;

	err = ((appliedStress(0) - net_stress(0, 0))*(appliedStress(0) - net_stress(0, 0)) +
		(appliedStress(1) - net_stress(1, 1))*(appliedStress(1) - net_stress(1, 1)) + (appliedStress(2) - net_stress(2, 2))*(appliedStress(2) - net_stress(2, 2)));

}