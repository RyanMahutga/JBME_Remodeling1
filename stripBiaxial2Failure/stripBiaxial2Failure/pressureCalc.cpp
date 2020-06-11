#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

const double PI = 3.141592653589793238462643383279502884;
const double R = 8.3144598; // universal gas constant [J/mol*K]

void pressureCalc(double& phi, double& c_fcd, double& c_neg, 
	double& pressure, MatrixXd& F, double& fib_vol, 
	double& sum_fib_len, double& x_scale)
{
	double c_star = 170;// 154; // M - bathing solution concentration
	double T = 295; // temperature K
	//double psi = 0.5; // Ateshian et al. (2009) J. Biomech. Eng. 131
	//double c1 = 381; // Pa*m^3/Eq
	//double c2 = 0.241; // Pa*m^6/Eq^2
	//double c_fcd_0 = 100; // meq/L (roughly gives 0.5M at 0.5 fiber volume fraction see Porterfield paper)
	//double c_pos = 0.5 * (c_fcd + sqrt(c_fcd * c_fcd + 4 * c_star*c_star));
	double k = 12e-13; //1e-13; //0.9e-6; // mEq/L*m^2 1e-13 for full, 2e-13 for sparse

	double J = F.determinant();

	c_fcd = k * sum_fib_len / (x_scale*x_scale*J); //k*phi*sum_fib_len / fib_vol; //c_fcd_0 * phi; //c_fcd_0 * (1 - phi); // calculation of fixed charge density from fiber volume fractions
		//c_neg = 0.5 * (-c_fcd + sqrt(c_fcd * c_fcd + 4 * c_star*c_star)); // negative ion charge density
/*	
if (phi > 0.999)
	{
		phi = 0.999;
	}
	*/
	pressure = R * T * (sqrt(c_fcd*c_fcd + 4 * c_star*c_star) - 2 * c_star);// +50 / ((1 - phi)*(1 - phi));
		//0.5e6*(1 / ((c_fcd_0 - c_fcd)*(c_fcd_0 - c_fcd)) - 1 / (c_fcd_0 * c_fcd_0))+3000; // +c1 * psi*c_fcd + c2 * psi*psi*c_fcd*c_fcd;
																																						  
// pressure calculation based on Donnan Osmotic Pressure and entropic pressure (virial expansion)
}