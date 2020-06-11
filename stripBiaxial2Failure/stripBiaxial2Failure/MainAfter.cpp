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

#define EIGEN_USE_MKL_ALL

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
#include <direct.h>

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
	fiberResults fiberResults_n; // fiber esults
	netResults netResults_n; // network results

	VectorXi stab_node(3); stab_node = VectorXi::Zero(3); // stabalization node

	Matrix3d F(3, 3); F = Matrix3d::Zero(3, 3); // network deformation gradient

	// reading in data from PeriodicNetwork.txt
	ifstream inFile1, inFile;

	inFile1.open("nodeData399.txt");
	if (inFile1.fail()) { // if can't open file
		cout << "Unable to Open Node Data File!" << endl;
		cin.get();
		exit(1); // terminate with error
	}
	inFile1 >> F(0,0) >> F(0,1) >> F(0,2) >> F(1,0) >> F(1,1) >> F(1,2) >>F(2,0) >> F(2,1) >> F(2,2);
	inFile1.close();

	inFile.open("PeriodicNetwork_399.txt");
	if (inFile.fail()) { // if can't open file
		cout << "Unable to Open Network File!" << endl;
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
			fib_vol_m = fib_vol_m + fib_areas(f)*init_lens(f);
		}
		else if (fib_type(f) == 2)
		{
			fib_vol_e = fib_vol_e + fib_areas(f)*init_lens(f);
		}
		else
		{
			fib_vol_c = fib_vol_c + fib_areas(f)*init_lens(f);
		}
	}

	network0.fib_type = fib_type;
	network0.init_lens = init_lens;
	network0.fib_rads = fib_rads;
	network0.fib_areas = fib_areas;

	double V = F.determinant();

	double x_scale = 10.0e-6;//fib_areas.dot(init_lens) / (fiber_vol_fract*V);

	double fiber_vol_fract = fib_areas.dot(init_lens) / (x_scale*x_scale*V);

	network0.fiber_vol_fract = fiber_vol_fract;

	VectorXd fib_rads0 = network0.fib_rads;
	VectorXd init_lens0 = network0.init_lens;

	network_n = network0;

	double pressure = 0;
	Vector3d x_guess;
	Matrix3d net_stress;

	MatrixXd nodes_n = network_n.nodes;
	MatrixXd fibers_n = network_n.fibers;

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

	network0.nodes = nodes;

	phi_m = fib_vol_m / (x_scale*x_scale*V);
	phi_e = fib_vol_e / (x_scale*x_scale*V);
	phi_c = fib_vol_c / (x_scale*x_scale*V);

	clock_t start;
	double duration;

	start = clock();

	netResults_n.F = F;
	MatrixXd component_stress(9, 3); component_stress = MatrixXd::Zero(9, 3);

	int time_steps = 350;
	double dt = 0, dL = 0, dR = 0, EigRatio = 0;
	Vector3d appliedStress; appliedStress << 0, 0, 0;

	VectorXd time(time_steps), volume(time_steps);
	MatrixXd stretch(time_steps ,3), phi(time_steps,3), Orient(time_steps,4);
	MatrixXd xstress(time_steps,4), ystress(time_steps,4), zstress(time_steps,4), xzstress(time_steps,4);
	time = VectorXd::Zero(time_steps + 1);

	VectorXd fib_stress(num_fibers), fib_stretch(num_fibers);
	VectorXd R(num_fibers), L(num_fibers);
	R = VectorXd::Constant(num_fibers,1);
	L = VectorXd::Constant(num_fibers,1);

	Matrix3d Or;

	double tau_m = 6.0, tau_c = 48.0;
	double sig_m_targ = 750e3, sig_c_targ=30e3; // changed from 500/30

	string str, strfile;
	double P = 0, P_Given, Shear_Given;

	time(0) = 0;

	double M = 1.5, Mc=1.0;

	int simtype = 4;

	if (simtype == 1)
	{
		strfile = "Shear";
		cout << "Shear ... " << endl << endl;
		Shear_Given = tan(0.3491); // 20 degrees from current deformed state
		P_Given = 13000.0;
	}
	else if (simtype == 2)
	{
		strfile = "Underload";
		cout << "Underload ... " << endl << endl;
		Shear_Given = 0.0;
		P_Given = 9750.0; //0.75 of normal
	}
	else if (simtype == 3)
	{
		strfile = "Overload";
		cout << "Overload ... " << endl << endl;
		Shear_Given = 0.0;
		P_Given = 19500.0; // 1.5 of normal
	}
	else if (simtype == 4)
	{
		strfile = "Overload_Shear";
		cout << "Overload + Shear ... " << endl << endl;
		Shear_Given = tan(0.3491); // 20 degrees from current deformed state
		P_Given = 19500.0; // 1.5 of normal
	}
	else
	{
		cout << "Error! Inproper selection fo Simulation Type Choose 1-4." << endl;
		cin.get();
		exit(1);
	}

	_mkdir(strfile.c_str());
	_chdir(strfile.c_str());

	for (int t = 0; t < time_steps; t++)
	{
		// setting dt
		if (t < 50)
		{
			dt = 0.01;
		}
		else if (t < 150) 
		{
			dt = 0.1;
		}
		else if (t < 350)
		{
			dt = 0.5;
		}

		time(t + 1) = time(t) + dt;

		x_guess(0) = F(0, 0);
		x_guess(1) = F(1, 1);
		x_guess(2) = F(2, 2);

		P = P_Given;
		F(0, 2) = F(2,2)*Shear_Given; 

		jacobianSolver(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
			stab_node, appliedStress, net_stress, pressure, x_guess, fib_vol_m, phi_m, P);

		fibers_n = network_n.fibers;
		nodes_n = network_n.nodes;
		F = netResults_n.F;

		calcPeriodicOrientation(Or, EigRatio, F, nodes_n, fibers_n, num_fibers, fib_areas, init_lens);

		Orient(t, 0) = Or(0, 0);
		Orient(t, 1) = Or(1, 1);
		Orient(t, 2) = Or(2, 2);
		Orient(t, 3) = EigRatio;

		if (t == 0 || t == time_steps - 1 )
		{
			writeData(num_fibers, network_n, fiberResults_n, netResults_n, t);

			writePeriodicNet2File(network_n, num_fibers, num_nodes, stab_node, t);
		}

		stretch(t, 0) = F(0, 0);
		stretch(t, 1) = F(1, 1);
		stretch(t, 2) = F(2, 2);

		pressure = netResults_n.pressure;

		xstress(t, 0) = net_stress(0, 0);
		ystress(t, 0) = net_stress(1, 1);
		zstress(t, 0) = net_stress(2, 2);
		xzstress(t, 0) = net_stress(0, 2);

		componentStress(nodes_n, fibers_n, network_n, fiberResults_n, F, num_fibers, x_scale, component_stress);

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

		V = F.determinant();

		phi_e = fib_vol_e / (x_scale*x_scale*V);
		phi_c = fib_vol_c / (x_scale*x_scale*V);
		phi_m = fib_vol_m / (x_scale*x_scale*V);

		phi(t, 0) = phi_m;
		phi(t, 1) = phi_e;
		phi(t, 2) = phi_c;

		cout << "Time: " << t << endl << "Fractions: " << phi.row(t) << endl << "Stress: " << endl << net_stress << endl
			<< "Pressure: " << pressure << endl << "Component Stresses: " << endl << component_stress << endl << endl;

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

						if (dR > 0)
						{
							init_lens(f) = (init_lens(f) * (((fib_rads(f)*fib_rads(f) / M + 2.0*fib_rads(f)*dR + dR * dR)) / (fib_rads(f)*fib_rads(f) / M + (2.0 * fib_rads(f)*dR + dR * dR) / fib_stretch(f))));
						}

					fib_rads(f) = fib_rads(f) + dR;

					if (init_lens(f) < 1e-6)
					{
						init_lens(f) = 1e-6;
					}
					if (fib_rads(f) < 1e-12)
					{
						fib_rads(f) = 1e-12;
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

					if (dR > 0) // there is radial growth
					{
						init_lens(f) = (init_lens(f) * (((fib_rads(f)*fib_rads(f) / Mc + 2.0*fib_rads(f)*dR + dR * dR)) / (fib_rads(f)*fib_rads(f) / Mc + (2.0 * fib_rads(f)*dR + dR * dR) / fib_stretch(f))));
					}

					fib_rads(f) = fib_rads(f) + dR;

					if (fib_rads(f) < 1e-12)
					{
						fib_rads(f) = 1e-12;
					}

					fib_areas(f) = PI * fib_rads(f)*fib_rads(f);
					fib_vol_c = fib_vol_c + fib_areas(f)*init_lens(f);
				}
			}

			network_n.init_lens = init_lens;
			network_n.fib_rads = fib_rads;
			network_n.fib_areas = fib_areas;
		}
	}

	// saving stretch results
	ofstream myfile1;
	str = "stretch.txt";
	myfile1.open(str);
	myfile1 << stretch;
	myfile1.close();

	// saving stress results
	ofstream myfile2;
	str = "xstress.txt";
	myfile2.open(str);
	myfile2 << xstress;
	myfile2.close();

	// saving stress results
	ofstream myfile3;
	str = "ystress.txt";
	myfile3.open(str);
	myfile3 << ystress;
	myfile3.close();

	// saving stress results
	ofstream myfile4;
	str = "zstress.txt";
	myfile4.open(str);
	myfile4 << zstress;
	myfile4.close();

	// saving stress results
	ofstream myfile5;
	str = "xzstress.txt";
	myfile5.open(str);
	myfile5 << xzstress;
	myfile5.close();

	ofstream myfile6;
	str = "phi.txt";
	myfile6.open(str);
	myfile6 << phi;
	myfile6.close();

	ofstream myfile7;
	str = "time.txt";
	myfile7.open(str);
	myfile7 << time;
	myfile7.close();

	ofstream myfile8;
	str = "orientation.txt";
	myfile8.open(str);
	myfile8 << Orient;
	myfile8.close();
	
	duration = (clock() - start) / (double)CLOCKS_PER_SEC;

	MatrixXd comp_stress(15, 3); comp_stress = MatrixXd::Zero(15, 3);

	componentStress(network_n.nodes, network_n.fibers, network_n, fiberResults_n, F, num_fibers, x_scale, comp_stress);

	cout << "Deformation Gradient: " << endl << F << endl << endl;
	cout << "Stress: " << endl << netResults_n.net_stress << endl << endl;
	cout << "Component Stress: " << endl<<comp_stress << endl << endl;
	cout << "Elapsed time : " << duration << endl << endl; ;

	//cout << "Press ENTER to exit." << endl;
	//cin.get();
	return 0;
}