/* Main.cpp : Solves the Displacement (stretch) Periodic Boundary Condition Problem for a Periodic Network.

NOTE: To compile this code you will need to have the Eigen library. (https://eigen.tuxfamily.org/)

DETAILED DESCRIPTION: This code relies on reading in two text files that contains the nodal coordinates and the fiber data. These files
can be generated using the appropriate matlab network generator for periodic fiber networks. The solution can use any of three constitutive
equations for fibers (linear, exponential, of helical spring), the network equilibrium is found through Netwon's Method relying on a GMRES
solver in the Eigen Library.

INPUTS: nodes - Nx3 matrix of nodal coordinates in an RVE
		fibers - Mx5 matrix of fiber connectivity, fibers (m,1) and fibers(m,2) represent the number of the two nodes that are connected,
			fibers(m,3:5) are the x,y,z crossing, i.e. {5, 6, +1, 0, -1) means that nodes 5 and 6 are connected, and, from node 5,
			the fiber crosses the posivie x boundary, does not cross the y boundary, and crosses the negative z boundary.
		fib_type - Mx1 vector of fiber type (1 for collagen, 2 for elastin, 3 for ?)
		fib_rads - Mx1 vector of fiber radii (in meters)
		init_lens - Mx1 vector of fiber initial lengths
		F - 3x3 deformation matrix

OUTPUTS:

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Contact: mahut005@umn.edu

Last Modifed: 04/29/2018

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "math.h"
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>
#include <cstdio>
#include <ctime>
#include "dataStructureDefinitions.h"

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

// Main function
int main()
{
	// setting up parallelization
	//Eigen::initParallel();
	//Eigen::setNbThreads(4);
	//int nthreads = Eigen::nbThreads();
	//cout << "THREADS = " << nthreads << endl << endl;

	//omp_set_num_threads(1);

	network network0; // underformed network structure definition
	network network_n; // deformed/modified network structure definition
	network network_n0;
	fiberResults fiberResults_n; // fiber esults
	netResults netResults_n; // network results
	fiberResults fiberResults_n0; // fiber esults
	netResults netResults_n0; // network results

	VectorXi stab_node(3); stab_node = VectorXi::Zero(3); // stabalization node

	Matrix3d F(3, 3); F = Matrix3d::Zero(3, 3); // network deformation gradient

	// reading in data from PeriodicNetwork.txt
	ifstream inFile;

	inFile.open("PeriodicNetworke.txt");
	if (inFile.fail()) { // if can't open file
		cout << "Unable to Open File!" << endl;
		cin.get();
		exit(1); // terminate with error
	}
	int num_fibers = 0;
	int num_nodes = 0;
	// reading in first line (number of fibers, number of nodes)
	inFile >> num_fibers >> num_nodes >> network0.fiber_vol_fract >> stab_node(0) >> stab_node(1) >> stab_node(2);

	network0.num_fibers = num_fibers;
	network0.num_nodes = num_nodes;

	// initialize values
	MatrixXd fibers = MatrixXd::Zero(num_fibers, 5);
	MatrixXd nodes = MatrixXd::Zero(num_nodes, 3);
	VectorXd fib_rads = VectorXd(num_fibers);
	VectorXd init_lens = VectorXd(num_fibers);
	VectorXi fib_type = VectorXi::Zero(num_fibers);

	int count = -1; // counter for fibers
	int count2 = 0; // counter for nodes
	while (!inFile.eof()) // read file until end
	{
		if (count < num_fibers && count > -1)
		{
			inFile >> fibers(count, 0) >> fibers(count, 1) >> fibers(count, 2) >> fibers(count, 3)
				>> fibers(count, 4) >> fib_type(count) >> init_lens(count) >> fib_rads(count);
		}
		else if (count >= num_fibers && count < num_fibers + num_nodes)
		{
			inFile >> nodes(count2, 0) >> nodes(count2, 1) >> nodes(count2, 2);
			count2++;
		}
		else if (count > -1)
		{
			cout << "Problem Reading in Data!" << endl << endl;
			cout << "Check to ensure there are no blank lines in input file." << endl;
			cin.get();
			exit(1);
		}
		count++;
	}

	inFile.close();

	network0.nodes = nodes;
	network0.fibers = fibers;
	network0.fib_type = fib_type;
	network0.init_lens = init_lens;
	network0.fib_rads = fib_rads;

	//-----------------------------------------------------------------------------------------------------------

	// defining passed variables

	cout.precision(9); // set output precision

	double fib_vol_m = 0, fib_vol_e = 0, fib_vol_c = 0, phi_e = 0, phi_c, phi_m = 0;

	VectorXd fib_areas = PI * fib_rads.cwiseProduct(fib_rads);

	for (int f = 0; f < num_fibers; f++)
	{
		if (fib_type(f) == 1)
		{
			fib_rads(f) = 8.0*fib_rads(f);
			fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
			fib_vol_m = fib_vol_m + fib_areas(f)*init_lens(f);
		}
		else if (fib_type(f) == 2)
		{
			fib_rads(f) = 24.0*fib_rads(f);
			fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
			fib_vol_e = fib_vol_e + fib_areas(f)*init_lens(f);
		}
		else
		{
			fib_rads(f) = 8.0*fib_rads(f);
			fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
			fib_vol_c = fib_vol_c + fib_areas(f)*init_lens(f);
		}
	}

	network0.fib_type = fib_type;
	network0.init_lens = init_lens;
	network0.fib_rads = fib_rads;
	network0.fib_areas = fib_areas;

	F << 1.20, 0, 0, 0, 1.10, 0, 0, 0, 0.9; 

	double V = F.determinant();

	double x_scale = 10.0e-6;//fib_areas.dot(init_lens) / (fiber_vol_fract*V);

	double fiber_vol_fract = fib_areas.dot(init_lens) / (x_scale*x_scale*V);

	network0.fiber_vol_fract = fiber_vol_fract;

	VectorXd fib_rads0 = network0.fib_rads;
	VectorXd init_lens0 = network0.init_lens;

	network_n = network0;
	network_n0 = network0;

	double pressure = 0;
	Vector3d x_guess;
	Matrix3d net_stress;

	MatrixXd nodes_n = network_n.nodes;
	MatrixXd fibers_n = network_n.fibers;

	MatrixXd fibers_n0 = fibers_n;
	MatrixXd nodes_n0 = nodes_n;

	/*
	MatrixXd F_inv(3, 3); F_inv = F.inverse(); // inverse deformation mapping

	MatrixXd I(3, 3); I << 1, 0, 0, 0, 1, 0, 0, 0, 1; // creating identity matrix

	MatrixXd T2(3, 3); T2 = 0.5*(F_inv - I); // inverse of the rigid translation

	
	for (int k = 0; k < num_nodes; k++) // map back to undeformed
	{
		nodes(k, 0) = F_inv(0, 0)*nodes_n(k, 0) + F_inv(0, 1)*nodes_n(k, 1) + F_inv(0, 2)*nodes_n(k, 2)
			+ (T2(0, 0) + T2(0, 1) + T2(0, 2));;
		nodes(k, 1) = F_inv(1, 0)*nodes_n(k, 0) + F_inv(1, 1)*nodes_n(k, 1) + F_inv(1, 2)*nodes_n(k, 2)
			+ (T2(1, 0) + T2(1, 1) + T2(1, 2));
		nodes(k, 2) = F_inv(2, 0)*nodes_n(k, 0) + F_inv(2, 1)*nodes_n(k, 1) + F_inv(2, 2)*nodes_n(k, 2)
			+ (T2(2, 0) + T2(2, 1) + T2(2, 2));
	}
	*/

	network0.nodes = nodes;

	network network00 = network0;

	phi_m = fib_vol_m / (x_scale*x_scale*V);
	phi_e = fib_vol_e / (x_scale*x_scale*V);
	phi_c = fib_vol_c / (x_scale*x_scale*V);

	clock_t start;
	double duration;

	start = clock();

	netResults_n.F = F;
	MatrixXd component_stress(12, 3); component_stress = MatrixXd::Zero(12, 3);

	int time_steps = 450;
	double dt = 0, dL = 0, dR = 0;
	Vector3d appliedStress, zeroStress; appliedStress << 0, 0, 0; zeroStress = appliedStress;

	VectorXd time(time_steps), volume(time_steps),press0(time_steps), c2ksi(6), c2theta(6);
	MatrixXd stretch(time_steps ,3), phi(time_steps,3), Orient(time_steps,18),stretch0(time_steps, 3), c2t(time_steps, 6), c2k(time_steps, 6);
	MatrixXd xstress(time_steps,4), ystress(time_steps,4), zstress(time_steps,4), xzstress(time_steps,4);
	time = VectorXd::Zero(time_steps + 1);

	VectorXd tau_wall = time;

	VectorXd fib_stress(num_fibers), fib_stretch(num_fibers);

	Matrix3d F0, net_stress0;

	VectorXd Or(18); 

	F0 = Matrix3d::Zero();
	net_stress0 = F0;
	
	double tau_m = 4.5, tau_c = 90.0, EigRatio = 0, tau_mL = 120.0, tau_cL = 240.0, tau_L = 360.0;
	double sig_m_targ = 750e3, sig_c_targ = 200e3; // changed from 750/30

	string str;
	double P = 0.0, pressure0 = 0.0, WSS=0.0, WSS0=0.0;

	time(0) = 0;

	double M = 1.0, Mc=1.0;

	double ctime = 0.0;
	Vector3d old_F0, del_F0;

	cout << "Initial Remodelling: " << endl;

	for (int t = 0; t < time_steps; t++)
	{
		// setting dt
		if (t < 200)
		{
			dt = 0.05;
		}
		else if (t < 350)
		{
			dt = 0.1;
		}
		else if (t < 450)
		{
			dt = 0.5;
		}
	/*	else if (t < 400)
		{
			dt = 1.0;
		}
		else if (t < 500)
		{
			dt = 1.0;
		}
		*/
		time(t + 1) = time(t) + dt;
		ctime = time(t);

		x_guess(0) = F(0, 0);
		x_guess(1) = F(1, 1);
		x_guess(2) = F(2, 2);

		if (t < 100)
		{
			P = 7000.0*(double(t + 1.0)) / 100.0 + 6500.0;
		}
		else
		{
			P = 13500.0; // MAP in Pa for 100mHg
		}

		jacobianSolver(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
			stab_node, appliedStress, net_stress, pressure, x_guess, fib_vol_m, phi_m, P, 0, WSS);

		tau_wall(t) = WSS;

		fibers_n = network_n.fibers;
		nodes_n = network_n.nodes;
		F = netResults_n.F;

		if (t>-1)
		{
			if (t==0)
			{
				F0(0, 0) = 0.997;//0.931;
				F0(1, 1) = 0.998;//0.917;
				F0(2, 2) = 1.026;//0.916;

				// zero-stress state
				jacobianSolver(num_nodes, num_fibers, x_scale, F0, network00, network_n0, fiberResults_n0, netResults_n0, 0,
					stab_node, zeroStress, net_stress0, pressure0, x_guess, fib_vol_m, phi_m, 500.0, 1, WSS0);
				cout << net_stress0 << endl << netResults_n0.pressure << endl;

			}
			else
			{
				F0(0, 0) += 0.001; //del_F0(0)*dt;
				F0(1, 1) += 0.001;//del_F0(1)*dt;
				F0(2, 2) += 0.001;//del_F0(1)*dt;;

				jacobianSolver(num_nodes, num_fibers, x_scale, F0, network00, network_n0, fiberResults_n0, netResults_n0, 0,
					stab_node, zeroStress, net_stress0, pressure0, x_guess, fib_vol_m, phi_m, 500.0, 1, WSS0);
			}
			F0 = netResults_n0.F;
		}

		press0(t) = netResults_n0.pressure;

		fib_type = network_n.fib_type;

		calcPeriodicOrientation2(Or, F0, nodes_n0, fibers_n0, num_fibers, fib_areas, init_lens, fib_type, c2ksi, c2theta);

		Orient.row(t) = Or;
		c2k.row(t) = c2ksi;
		c2t.row(t) = c2theta;

		if (t == 0 || t == time_steps - 1 || t % 100 == 0)
		{
			writeData(num_fibers, network_n, fiberResults_n, netResults_n, t);

			writePeriodicNet2File(network_n, num_fibers, num_nodes, stab_node, t);

			if (t == time_steps - 1)
			{
				writePeriodicNet2File(network00, num_fibers, num_nodes, stab_node, t*10);
			}
		}

		stretch(t, 0) = F(0, 0);
		stretch(t, 1) = F(1, 1);
		stretch(t, 2) = F(2, 2);

		stretch0(t, 0) = F0(0, 0);
		stretch0(t, 1) = F0(1, 1);
		stretch0(t, 2) = F0(2, 2);

		pressure = netResults_n.pressure;

		xstress(t, 0) = net_stress(0, 0);
		ystress(t, 0) = net_stress(1, 1);
		zstress(t, 0) = net_stress(2, 2);
		xzstress(t, 0) = net_stress(0, 2);

		component_stress = netResults_n.component_stress;

		xstress(t, 1) = component_stress(0, 0);
		ystress(t, 1) = component_stress(1, 1);
		zstress(t, 1) = component_stress(2, 2);
		xzstress(t, 1) = component_stress(0, 2);

		xstress(t, 2) = component_stress(3, 0);
		ystress(t, 2) = component_stress(4, 1);
		zstress(t, 2) = component_stress(5, 2);
		xzstress(t, 2) = component_stress(3, 2);

		xstress(t, 3) = component_stress(6, 0);
		ystress(t, 3) = component_stress(7, 1);
		zstress(t, 3) = component_stress(8, 2);
		xzstress(t, 3) = component_stress(6, 2);

		V = F0.determinant();

		phi_e = fib_vol_e / (x_scale*x_scale*V);
		phi_c = fib_vol_c / (x_scale*x_scale*V);
		phi_m = fib_vol_m / (x_scale*x_scale*V);

		phi(t, 0) = phi_m;
		phi(t, 1) = phi_e;
		phi(t, 2) = phi_c;

		appendResults(stretch.row(t), stretch0.row(t), xstress.row(t), ystress.row(t), zstress.row(t), xzstress.row(t),
			phi.row(t), time(t), Orient.row(t), c2k.row(t), c2t.row(t), press0(t), tau_wall(t));

		cout << "Time: " << t << endl << "Fractions: " << phi.row(t) << endl << "Stress: " << endl << net_stress << endl
			<< "Pressure: " << pressure << endl << "Component Stresses: " << endl << component_stress << endl << endl
			<< "Growth: " << stretch0.row(t) << endl << "Rest Pressure: " << press0(t) << endl << endl;;

		fib_stress = fiberResults_n.fib_stress;
		fib_stretch = fiberResults_n.fib_stretch;

		if (t < time_steps - 1)
		{
			fib_vol_m = 0;
			fib_vol_c = 0;
			for (int f = 0; f < num_fibers; f++)
			{
				if (fib_type(f) == 1)
				{
					if (fib_stress(f) > 0)
					{
							dR = 1.0 / tau_m * (fib_stress(f) / sig_m_targ - 1.0)*fib_rads(f)*dt;
					}
					else
					{
							dR = -1.0 / tau_m * fib_rads(f)* dt;
					}

						if (dR > 1e-18)
						{
							//dL = 1.0 / tau_mL * (fib_stress(f) / sig_m_targ - 1.0)*init_lens(f)*dt;
							init_lens(f) = (init_lens(f) * (((fib_rads(f)*fib_rads(f) / M + 2.0*fib_rads(f)*dR + dR * dR)) / (fib_rads(f)*fib_rads(f) / M + (2.0 * fib_rads(f)*dR + dR * dR) / fib_stretch(f))));
							//init_lens(f) = init_lens(f) + dL;
						}

					fib_rads(f) = fib_rads(f) + dR;

					if (fib_rads(f) < 1e-10)
					{
						fib_rads(f) = 1e-10;
					}

					fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
					fib_vol_m = fib_vol_m + fib_areas(f)*init_lens(f);
				}
				else if (fib_type(f) == 3)
				{
					if (fib_stress(f) > 0)
					{		
						dR = 1.0 / tau_c * (fib_stress(f) / sig_c_targ - 1.0)*fib_rads(f)*dt;
					}
					else
					{
						dR = -1.0 / tau_c * fib_rads(f) * dt;
					}

					if (dR > 1e-18) // there is radial growth
					{
						//dL = 1.0 / tau_cL * (fib_stress(f) / sig_c_targ - 1.0)*init_lens(f)*dt;
						init_lens(f) = (init_lens(f) * (((fib_rads(f)*fib_rads(f) / M + 2.0*fib_rads(f)*dR + dR * dR)) / (fib_rads(f)*fib_rads(f) / M + (2.0 * fib_rads(f)*dR + dR * dR) / fib_stretch(f))));
						//init_lens(f) = init_lens(f) + dL;
					}

					fib_rads(f) = fib_rads(f) + dR;

					if (fib_rads(f) < 1e-10)
					{
						fib_rads(f) = 1e-10;
					}

					fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
					fib_vol_c = fib_vol_c + fib_areas(f)*init_lens(f);
				}
			}

			network_n.init_lens = init_lens;
			network_n.fib_rads = fib_rads;
			network_n.fib_areas = fib_areas;

			network_n0.init_lens = init_lens;
			network_n0.fib_rads = fib_rads;
			network_n0.fib_areas = fib_areas;
		}
	}

	duration = (clock() - start) / (double)CLOCKS_PER_SEC;

	cout << "Deformation Gradient: " << endl << F << endl << endl;
	cout << "Stress: " << endl << netResults_n.net_stress << endl << endl;
	cout << "Component Stress: " << endl<<component_stress << endl << endl;
	cout << "Elapsed time : " << duration << endl << endl; ;

	return 0;
}