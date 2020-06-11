/* fiberConstEqn.cpp calculates nodal anf fiber forces for a number of constitutive equations

SUMMARY: fiberConstEqn.cpp is a function with several fiber constitutive equations implemented.
To change the fiber constitutive equation simply change the value of 'const_eqn' in this function.
This function also calculates dF/dlambda for the constitutive equations for use in the analytical 
Jacobian calculation. 

CALLED FROM: calcForcesPeriodic.cpp
CALLS ON: None

INPUTS: vect1 - 3x1 vector along fibers
		node_force1 - 3x1 vector of forces acting on node1 of fiber
		node_force2 - 3x1 vector of forces acting on node2 of fiber
		init_len - double value of fiber initial length
		fibtype - int value of fiber type
		fib_area - double value of fiber cross-sectional area
		fib_force - double value of fiber force
		dFdlam - double value of dF/dlambda for the constitutive model

OUTPUTS: None

NOTE: node_force1, node_force2, fib_force, and dFdlam are modified by this code.

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 05-04-2018

*/

#include <iostream>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"

using namespace Eigen;
using namespace std;

const double PI = 3.141592653589793238462643383279502884;

void fiberConstEqn(Vector3d vect1, Vector3d& node_force1, Vector3d& node_force2, double& lambda, const double init_len, const int fibtype, const double fib_area, double& fib_force, double& dFdlam)
{
	double GS = 0; // Green Strain
	double fib_mod = 0; // fiber modulus
	lambda = 0; // fiber stretch
	double fib_length = 0; // fiber length

	fib_length = vect1.norm(); // fiber length

	lambda = fib_length / init_len; // fiber stretch

	if (fibtype == 3 || fibtype == 4) // fiber force from helical model in Freed, A.D. & Doehring, T.C. (2005) J. Biomech. Eng. 127  
	{ // (1 for collagen, 4 for active contractile fibers)
		fib_mod = 700e6;

		double R0 = 5.8; // [nm] roughly radius of collagen microfibril
		double r0 = 1.6; // [nm] roughly a collagen triple helix
		double H0 = 67.4; //[nm] d-pattern banding in collagen molecules

		double L0 = sqrt((2 * PI*R0)*(2 * PI*R0) + H0 * H0);

		double lambda_bar = L0 / H0;

		double E_bar = fib_mod * H0*H0 / (L0*(H0 + (1 + 37 / (6 * PI*PI) + 2 * (L0*L0) / ((PI*r0)*(PI*r0)))*(L0 - H0)));

		double sigma = 0;

		if (lambda <= lambda_bar) // see above paper for more information
		{
			double H = lambda * H0;
			double R = sqrt(L0*L0 - H * H) / (2 * PI);
			double eta = (R*R + H * H) / (L0*H*(1 + 4 * R*R / (r0*r0) + 6 * (20 / 9 + R * R / (r0*r0))*R*R / (H*H)));

			double dHdlam = H0;
			double dRdH = -H / (2 * PI*sqrt((L0*L0 - H * H)));
			double dEtadH = -(H*H + R * R) / ((H*H) * L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +
				2 / (L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +
				12 * (R*R) * ((H*H) + (R*R))*((R*R) / (r0*r0) + 20 / 9) / ((H*H*H*H) * L0*((6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
			double dEtadR = 2 * R / (H*L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) -
				((H*H) + (R*R))*(12 * (R*R) / ((r0*r0) * (H*H)) + 12 * R*((R*R) / (r0*r0) + 20 / 9) / (H*H) + 8 * R / (r0*r0)) / (H*L0*((6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
			double dEtadlam = dEtadH * dHdlam + dEtadR * dRdH*dHdlam;

			if (lambda >= 1)
			{
				sigma = eta * E_bar*(lambda - 1);
				dFdlam = fib_area * (eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
			}
			else // make compression very low
			{
				sigma = 1e-12*eta * E_bar*(lambda - 1);
				dFdlam = fib_area * 1e-12*(eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
			}
		}
		else
		{
			sigma = E_bar * (lambda_bar - 1) + fib_mod * (lambda / lambda_bar - 1);
			dFdlam = fib_mod * fib_area / lambda_bar;
		}

			fib_force = fib_area * sigma;
	}
	else if ( fibtype == 2 ) // fiber force from F = k(GS) where GS = Green strain = 0.5(lambda^2-1) 
	{//(use 2 for elastin, 5 for failed fibers)
		GS = lambda - 1.0;

		if (fibtype == 2)
		{
			if (lambda >= 1.0)
			{
				fib_mod = 1e6; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
				dFdlam = fib_mod * fib_area*lambda;
			}
			else
			{
				fib_mod = 1e6/1000; // 1% of elastin stiffness for failed fibers
				dFdlam = fib_mod * fib_area*lambda*5.0;
			}
		}
		else if (fibtype == 11) // actin
		{
			if (lambda >= 1.0) //  give 12pN force at lam = 1.05
			{
				fib_mod = 13e6; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
				dFdlam = fib_mod * fib_area*lambda;
			}
			else
			{
				fib_mod = 13e6/1000; // 1% of elastin stiffness for failed fibers
				dFdlam = fib_mod * fib_area*lambda*5.0;
			}
		}
		else // collagen
		{
			if (lambda >= 1.0)
			{
				fib_mod = 700e6; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
				dFdlam = fib_mod * fib_area*lambda;
			}
			else
			{
				fib_mod = 700e6/1000.0; // 1% of elastin stiffness for failed fibers
				dFdlam = fib_mod * fib_area*lambda;
			}
		}
		fib_force = fib_mod * fib_area*GS;
	}

	else if (fibtype == 1) // actin
	{
		GS = lambda - 1.0;

		// passive contribution
		if (lambda >= 1.0) //  give 12pN force at lam = 1.05
		{
			fib_mod = 13e6; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
			dFdlam = fib_mod * fib_area*lambda;
		}
		else
		{
			fib_mod = 13e6 / 1000; // 1% of elastin stiffness for failed fibers
			dFdlam = fib_mod * fib_area*lambda*5.0;
		}
		//active contribution
			double fib_S = 100e3, l_min = 0.65, l_max = 1.4; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
			
			double f_lam = 1 - (l_max-lambda)*(l_max-lambda) / ((l_max - l_min)*(l_max-l_min));

			double fib_force_active = fib_S*f_lam*fib_area;

			double dFdlamS = fib_area * 2.0 * fib_S * (l_max-lambda) / ((l_max-l_min)*(l_max-l_min)); 

			fib_force = fib_mod * fib_area*GS + fib_force_active;

			dFdlam = dFdlam + dFdlamS;
		
	}
	else if (fibtype==0)
	{
		GS = lambda - 1.0;

		fib_mod = 100; // failed fiber modulus
		dFdlam = fib_mod * fib_area*lambda;
		fib_force = fib_mod * fib_area*GS;
	}
	else if (fibtype == 4) // fiber force from F = EA*(exp(GS)-1) where GS = Green strain = 0.5(lambda^2-1)
	{

		double lambda_limit = 1.15;
		double B;

		fib_mod = 2569; //  fit to Ruberti single fiber experiment
		B = 77.2;

		GS = 0.5*(lambda*lambda - 1);  // Green Strain

		if (lambda >= lambda_limit)
		{
			GS = 0.5*(lambda_limit*lambda_limit - 1);  // Green Strain at limit
			double force_exp = fib_mod * fib_area*(exp(GS*B) - 1) / B;
			double slope_at_limit = fib_mod * fib_area*lambda_limit*exp(B*GS);
			fib_force = force_exp + slope_at_limit * (lambda - lambda_limit);
			dFdlam = slope_at_limit;
		}
		else if (lambda >= 1)
		{
			fib_force = fib_mod * fib_area*(exp(B*GS) - 1) / B;
			dFdlam = fib_mod * fib_area*lambda*exp(B*GS);
		}
		else if (lambda < 1)
		{
			fib_force = 1e-12*fib_mod * fib_area*(exp(B*GS) - 1) / B;
			dFdlam = 1e-12*fib_mod * fib_area*lambda*exp(B*GS);
		}
	}

	// angles for forces
	Vector3d cosine = vect1 / fib_length;

	// nodal forces due to stretched fibers
	node_force1 = fib_force*cosine;
	node_force2 = -1 * fib_force*cosine;

}