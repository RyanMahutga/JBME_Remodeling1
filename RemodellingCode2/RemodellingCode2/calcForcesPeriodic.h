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
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess, 
	VectorXi stab_node, double WSS); // Cholesky LL^T direct solve

void calcForcesPeriodic(MatrixXd nodes_n, MatrixXd nodes_mapped, MatrixXd fibers_n, network& network_n, Matrix3d F,
	fiberResults& fiberResults_n, netResults& netResults_n, int num_fibers, int num_nodes, double WSS); // calculating all forces on fibers

void updateCrossing(MatrixXd& nodes_mapped, MatrixXd& nodes_n, MatrixXd& fibers_n, int num_nodes, int num_fibers, Matrix3d F);
// updating crossing conditions if fibers leave RVE boundary

void fiberConstEqn(Vector3d vect, Vector3d& node_force1, Vector3d& node_force2, double& lambda, double init_len, int fibtype, 
	double fib_area, double& fib_force, double& dFdlam, double WSS); // fiber Constitutive models

void calcJacobian2(MatrixXd& nodes_n, MatrixXd& fibers_n, network& network_n, Matrix3d F, VectorXd& node_forces,
	int num_fibers, int num_nodes, SparseMatrix<double>& J, VectorXi stab_node, double WSS); // Jacobain and Residual calculation

void calcBndNodes(MatrixXd nodes_mapped, MatrixXd nodes_n, MatrixXd fibers_n, int num_nodes, int num_fibers, 
	Matrix3d F, MatrixXd& bnd_nodes);
// calculates the nodal positions of the boundary nodes

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, netResults& netResults_n, 
	fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale);
// component stress calculation

// 3D Newton Solver
void jacobianSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P, int type_flag, double& WSS, double Q);

void solveJac(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, Vector3d appliedStress, Matrix3d& net_stress, Vector3d& err, double P, 
	double fib_vol_m, double& WSS, double Q);

//void bubbleSort(MatrixXd& x);
void calcPeriodicOrientation2(VectorXd& R, Matrix3d& F, MatrixXd& nodes_n, MatrixXd& fibers_n, int num
	, VectorXd& fib_areas, VectorXd& init_lens, VectorXi fib_type, VectorXd& c2ksi, VectorXd& c2theta);

void writeData(int num_fibers, network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, int t);

void writePeriodicNet2File(network &network_n, int num_fibers, int num_nodes, Vector3i stab_node, int t);

void appendResults(Vector3d stretch, Vector3d stretch0, VectorXd xstress, VectorXd ystress, VectorXd zstress, VectorXd xzstress,
	Vector3d phi, double time, VectorXd Orient, VectorXd c2k, VectorXd c2t, double press0, double WSS);