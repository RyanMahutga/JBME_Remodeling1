/* calcForcePeriodic.h is the header file for solving the stretch control problem for periodic RVEs

Created by : Ryan Mahutga - Barocas Research Group - University of Minnesota
Date Modified : 04-23-18
*/

#define EIGEN_USE_MKL_ALL

#pragma once
#include "Eigen/Eigen"
#include "dataStructureDefinitions.h"

using namespace Eigen;

void networkSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess, VectorXi stab_node); // Cholesky LL^T direct solve

void calcForcesPeriodic(MatrixXd nodes_n, MatrixXd nodes_mapped, MatrixXd fibers_n, network& network_n, Matrix3d F,
	fiberResults& fiberResults_n, netResults& netResults_n, int num_fibers, int num_nodes); // calculating all forces on fibers

void updateCrossing(MatrixXd& nodes_mapped, MatrixXd& nodes_n, MatrixXd& fibers_n, int num_nodes, int num_fibers, Matrix3d F);
// updating crossing conditions if fibers leave RVE boundary

void fiberConstEqn(Vector3d vect, Vector3d& node_force1, Vector3d& node_force2, double& lambda, double init_len, int fibtype, 
	double fib_area, double& fib_force, double& dFdlam); // fiber Constitutive models

void calcJacobian2(MatrixXd& nodes_n, MatrixXd& fibers_n, network& network_n, Matrix3d F, VectorXd& node_forces,
	int num_fibers, int num_nodes, SparseMatrix<double>& J, VectorXi stab_node); // Jacobain and Residual calculation

//void calcJacobian(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, Matrix3d F, int num_fibers,
	// int num_nodes, SparseMatrix<double>& J, VectorXi stab_node); // calculate the Jacobain dF_i/dx_j

void calcNetStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, fiberResults& fiberResults_n, netResults& netResults_n,
	Matrix3d F, int num_fibs, double x_scale); //calculate network stresses

void calcBndNodes(MatrixXd nodes_mapped, MatrixXd nodes_n, MatrixXd fibers_n, int num_nodes, int num_fibers, Matrix3d F, MatrixXd& bnd_nodes);
// calculates the nodal positions of the boundary nodes

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale, MatrixXd& comp_stress);
// component stress calculation

void calcNetStress2(MatrixXd& nodes_n, MatrixXd& fibers_n, VectorXd& fib_forces, Matrix3d& net_stress, MatrixXd& bnd_nodes, Matrix3d F, int num_fibs, double x_scale);

//3D simplex 
//void simplexSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	//network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	//VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess, double fib_vol_m, double phi_m, double P);

//void solveSimp(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	//network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	//VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& err, double P, double fib_vol_m);

// 3D Newton Solver
void jacobianSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d& F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess, double fib_vol_m, double phi_m, double P);

void jacobianSolver2D(int num_nodes, int num_fibers, double x_scale, Matrix3d& F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P);

void solveJac(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, Vector3d& err, double P, double fib_vol_m);

void solveConst(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double P, double fib_vol_m);
//void bubbleSort(MatrixXd& x);

void calcPeriodicOrientation(Matrix3d& R, double& EigRatio, Matrix3d& F, MatrixXd& nodes_n, MatrixXd& fibers_n, int num
	, VectorXd& fib_areas, VectorXd& init_lens);

void writeData(int num_fibers, network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, int t);

void writePeriodicNet2File(network &network_n, int num_fibers, int num_nodes, Vector3i stab_node, int t);