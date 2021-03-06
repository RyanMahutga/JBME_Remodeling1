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
#include "omp.h"
#include <direct.h>

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

#define EIGEN_USE_MKL_ALL

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

	inFile1.open("nodeData349.txt");
	if (inFile1.fail()) { // if can't open file
		cout << "Unable to Open Node Data File!" << endl;
		cin.get();
		exit(1); // terminate with error
	}
	inFile1 >> F(0, 0) >> F(0, 1) >> F(0, 2) >> F(1, 0) >> F(1, 1) >> F(1, 2) >> F(2, 0) >> F(2, 1) >> F(2, 2);
	inFile1.close();

	inFile.open("PeriodicNetwork_349.txt");
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

	int lam_steps = 250;
	double dtau = 0.05, tau=0.0;

	VectorXd volume(lam_steps);
	MatrixXd stretch(lam_steps ,3), phi(lam_steps,3), Orient(lam_steps,4);
	VectorXd shear(lam_steps); shear(0) = 0.0;
	MatrixXd xstress(lam_steps,4), ystress(lam_steps,4), zstress(lam_steps,4), xzstress(lam_steps,4);

	VectorXd fib_stress(num_fibers), fib_stretch(num_fibers);
	Matrix3d Or;

	string str;
	double P = 0, EigRatio=0;

	Vector3d appliedStress; appliedStress << 0.0, 0.0, 0.0;
	
	Vector3d lam_crit; lam_crit << 2.0, 2.35, 1.14;

	cout << "Simple Shear to Failure: " << endl;

	x_guess << F(0),F(1),F(2);

	// find zero stress state
	jacobianSolver(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
		stab_node, appliedStress, net_stress, pressure, x_guess, fib_vol_m, phi_m, P);

	Vector3d lam0; lam0 << F(0,0), F(1,1),F(2,2);
	double tau0 = 0.0;
	Vector3i type; type << 0, 0, 0;
	MatrixXi total_failed_types(lam_steps, 3);
	cout << "Zero-Stress Stretches: " <<endl << lam0 << endl<<endl << "Stress: " << net_stress<<endl<<endl;

	// equibiaxial to failure
	int exit_flag = 0, flag2 = 0, break_flag1=0, break_flag2=0, k = 0, num_fail = 5, tot_fail = 0;

	double drop1=0.0, drop2 = 0.0;

	F(0, 0) = 1.4*lam0(0);
	F(1, 1) = 1.25*lam0(1);

	netResults_n.F = F;

	while (exit_flag == 0)
	{
		num_fail = 200;
		flag2 = 0;
		while (num_fail > 1 && flag2==0)
		{
			if (k > 0)
			{
				F(0, 2) = shear(k) + dtau;
			}
			//jacobianSolver2D(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
				//stab_node, appliedStress, net_stress, pressure, x_guess, fib_vol_m, phi_m, P);

			solveConst(num_nodes, num_fibers, x_scale, F, network0,
				network_n, fiberResults_n, netResults_n, 0, stab_node, appliedStress, net_stress, P, fib_vol_m);

			fib_stretch = fiberResults_n.fib_stretch;

			num_fail = 0;
			for (int m = 0; m < num_fibers; m++)
			{
				if (fib_type(m) == 1 && fib_stretch(m) > lam_crit(0))
				{
					num_fail++;
				}
				else if (fib_type(m) == 2 && fib_stretch(m) > lam_crit(1))
				{
					num_fail++;
				}
				else if (fib_type(m) == 3 && fib_stretch(m) > lam_crit(2))
				{
					num_fail++;
				}
			}

			if (num_fail > 10 && dtau > 0.01)
			{
				dtau = dtau - 0.01;
			}
			else if (num_fail == 0 && dtau < 0.1)
			{
				dtau = dtau + 0.05;
			}
			else
			{
				flag2 = 1;
			}

			if (flag2 == 1)
			{
				for (int m = 0; m < num_fibers; m++)
				{
					if (fib_type(m) == 1 && fib_stretch(m) > lam_crit(0))
					{
						fib_type(m) = 0;
						type(0) = type(0) + 1;
					}
					else if (fib_type(m) == 2 && fib_stretch(m) > lam_crit(1))
					{
						fib_type(m) = 0;
						type(1) = type(1) + 1;
					}
					else if (fib_type(m) == 3 && fib_stretch(m) > lam_crit(2))
					{
						fib_type(m) = 0;
						type(2) = type(2) + 1;
					}
				}

				network_n.fib_type = fib_type;

				//jacobianSolver2D(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
					//stab_node, appliedStress, net_stress, pressure, x_guess, fib_vol_m, phi_m, P);

				solveConst(num_nodes, num_fibers, x_scale, F, network0,
					network_n, fiberResults_n, netResults_n, 0, stab_node, appliedStress, net_stress, P, fib_vol_m);

				fib_stretch = fiberResults_n.fib_stretch;
			}
		}

			shear(k+1) = F(0, 2);

			tot_fail = tot_fail + num_fail;

			total_failed_types.row(k) = type;
			/*
			if (k > 1)
			{
				drop1 = net_stress(2, 2) - zstress(k-1,0);
				drop2 = net_stress(0, 2) - shear(k - 2);
				if (drop1 < -50e3 || drop2 < -30e3)
				{
					break_flag1 = 1;
				}
			}
			*/
			if (tot_fail > 0.5 * num_fibers || shear(k) > 2.5 )//|| break_flag1)
			{
				exit_flag = 1;
			}

		fibers_n = network_n.fibers;
		nodes_n = network_n.nodes;
		F = netResults_n.F;

		cout << "Failed Fibers: "<< tot_fail << endl<< endl << "Deformation: " << endl<< F << endl << endl
			<< "Stress: "<<endl<< net_stress << endl << endl;//<< component_stress << endl << endl;;

		//calcPeriodicOrientation(Or, EigRatio, F, nodes_n, fibers_n, num_fibers, fib_areas, init_lens);

		//Orient(k, 0) = Or(0, 0);
		//Orient(k, 1) = Or(1, 1);
		//Orient(k, 2) = Or(2, 2);
		//Orient(k, 3) = EigRatio;

		if (tot_fail > 0.1*num_fibers)
		{
			writeData(num_fibers, network_n, fiberResults_n, netResults_n, k);

			writePeriodicNet2File(network_n, num_fibers, num_nodes, stab_node, k);
		}

		stretch(k, 0) = F(0, 0);
		stretch(k, 1) = F(1, 1);
		stretch(k, 2) = F(2, 2);

		pressure = netResults_n.pressure;

		xstress(k, 0) = net_stress(0, 0);
		ystress(k, 0) = net_stress(1, 1);
		zstress(k, 0) = net_stress(2, 2);
		xzstress(k, 0) = net_stress(0, 2);

		componentStress(nodes_n, fibers_n, network_n, fiberResults_n, F, num_fibers, x_scale, component_stress);

		xstress(k, 1) = component_stress(0, 0);
		ystress(k, 1) = component_stress(1, 1);
		zstress(k, 1) = component_stress(2, 2);
		xzstress(k, 1) = component_stress(0, 2);

		xstress(k, 2) = component_stress(3, 0);
		ystress(k, 2) = component_stress(4, 1);
		zstress(k, 2) = component_stress(5, 2);
		xzstress(k, 2) = component_stress(3, 2);

		xstress(k, 3) = component_stress(6, 0);
		ystress(k, 3) = component_stress(7, 1);
		zstress(k, 3) = component_stress(8, 2);
		xzstress(k, 3) = component_stress(6, 2);

		V = F.determinant();

		phi_e = fib_vol_e / (x_scale * x_scale * V);
		phi_c = fib_vol_c / (x_scale * x_scale * V);
		phi_m = fib_vol_m / (x_scale * x_scale * V);

		phi(k, 0) = phi_m;
		phi(k, 1) = phi_e;
		phi(k, 2) = phi_c;

		fib_stress = fiberResults_n.fib_stress;
		fib_stretch = fiberResults_n.fib_stretch;

		network_n.init_lens = init_lens;
		network_n.fib_rads = fib_rads;
		network_n.fib_areas = fib_areas;

		//append results
		appendResults(stretch.row(k), xstress.row(k), ystress.row(k), zstress.row(k), xzstress.row(k),phi.row(k), shear(k), total_failed_types.row(k))

		k++; // increment k
	}

	
	duration = (clock() - start) / (double)CLOCKS_PER_SEC;

	MatrixXd comp_stress(15, 3); comp_stress = MatrixXd::Zero(15, 3);

	componentStress(network_n.nodes, network_n.fibers, network_n, fiberResults_n, F, num_fibers, x_scale, comp_stress);

	cout << "Deformation Gradient: " << endl << F << endl << endl;
	cout << "Stress: " << endl << netResults_n.net_stress << endl << endl;
	cout << "Component Stress: " << endl<<comp_stress << endl << endl;
	cout << "Elapsed time : " << duration << endl << endl; ;

	cout << "Press ENTER to exit." << endl;
	cin.get();
	return 0;
}