
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

void calcPeriodicOrientation(Matrix3d& R, double& EigRatio, Matrix3d& F, MatrixXd& nodes_n, MatrixXd& fibers_n, int num
,VectorXd& fib_areas, VectorXd& init_lens)

{
	double om11 = 0.0;
	double om12 = 0.0;
	double om13 = 0.0;
	double om22 = 0.0;
	double om23 = 0.0;
	double om33 = 0.0;

	double fiber_length = 0;

	double del_x = 0;
	double del_y = 0;
	double del_z = 0;

	double cosa = 0;
	double cosb = 0;
	double cosg = 0;
	double fib_vol = 0, total_vol = 0;

	int a = 0;
	int b = 0;

	Vector3d real_node2, vals;

	for (int i = 0; i < num; i++)
	{

		a = fibers_n(i, 0); // node 1 num
		b = fibers_n(i, 1); // node 2 num

		real_node2(0) = nodes_n(b, 0) + fibers_n(i, 2)*F(0, 0) + fibers_n(i, 3)*F(0, 1) + fibers_n(i, 4)*F(0, 2);
		real_node2(1) = nodes_n(b, 1) + fibers_n(i, 2)*F(1, 0) + fibers_n(i, 3)*F(1, 1) + fibers_n(i, 4)*F(1, 2);
		real_node2(2) = nodes_n(b, 2) + fibers_n(i, 2)*F(2, 0) + fibers_n(i, 3)*F(2, 1) + fibers_n(i, 4)*F(2, 2);

		del_x = nodes_n(a, 0) - real_node2(0); // x1 - x2
		del_y = nodes_n(a, 1) - real_node2(1); // y1 - y2
		del_z = nodes_n(a, 2) - real_node2(2); // z1 - z2

		fiber_length = sqrt(del_x*del_x + del_y*del_y + del_z*del_z);

		fib_vol = fib_areas(i)*init_lens(i);

		total_vol = total_vol + fib_vol;

		// direction cosines for current fiber

		cosa = del_x / fiber_length;
		cosb = del_y / fiber_length;
		cosg = del_z / fiber_length;

		om11 = om11 + fib_vol * cosa*cosa;
		om12 = om12 + fib_vol * cosa*cosb;
		om13 = om13 + fib_vol * cosa*cosg;
		om22 = om22 + fib_vol * cosb*cosb;
		om23 = om23 + fib_vol * cosb*cosg;
		om33 = om33 + fib_vol * cosg*cosg;
	}

	om11 = om11 / total_vol;
	om12 = om12 / total_vol;
	om13 = om13 / total_vol;
	om22 = om22 / total_vol;
	om23 = om23 / total_vol;
	om33 = om33 / total_vol;

	R << om11, om12, om13,
		om12, om22, om23,
		om13, om23, om33;

	SelfAdjointEigenSolver<Matrix3d> EigenSolver(R);
	vals = EigenSolver.eigenvalues();
	EigRatio = vals.maxCoeff() / vals.minCoeff();
}