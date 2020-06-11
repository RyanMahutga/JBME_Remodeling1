
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

void calcPeriodicOrientation2(VectorXd& R, Matrix3d& F, MatrixXd& nodes_n, MatrixXd& fibers_n, int num
		,VectorXd& fib_areas, VectorXd& init_lens, VectorXi fib_type, VectorXd& c2ksi, VectorXd& c2theta)

{
	VectorXd om1(6), om2(6), om3(6); om1 = VectorXd::Zero(6);
	om2 = om1; om3 = om1;

	double fiber_length = 0;

	double del_x = 0;
	double del_y = 0;
	double del_z = 0;

	double cosa = 0;
	double cosb = 0;
	double cosg = 0;
	double sin2ksi = 0;
	double fib_vol = 0;
	
	Vector3d total_vol; total_vol << 0.0, 0.0, 0.0;

	int a = 0, b = 0, m1 = 0, m2 = 0, m3 = 0;

	Vector3d real_node2, vals1, vals2, vals3;

	for (int i = 0; i < num; i++)
	{
		if (fib_type(i) == 1)
		{
			m1++;
		}
		else if (fib_type(i) == 2)
		{
			m2++;
		}
		else if (fib_type(i) == 3)
		{
			m3++;
		}
	}

	VectorXd cos2theta1(m1), cos2ksi1(m1), cos2theta2(m2), cos2ksi2(m2), cos2theta3(m3), cos2ksi3(m3);

	m1 = m2 = m3 = 0;

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

		// direction cosines for current fiber

		cosa = del_x / fiber_length;
		cosb = del_y / fiber_length;
		cosg = del_z / fiber_length;

		if (fib_type(i)==1)
		{ 
			om1(0) += fib_vol * cosa*cosa * fib_vol;
			om1(1) += fib_vol * cosa*cosb* fib_vol;
			om1(2) += fib_vol * cosa*cosg* fib_vol;
			om1(3) += fib_vol * cosb*cosb* fib_vol;
			om1(4) += fib_vol * cosb*cosg* fib_vol;
			om1(5) += fib_vol * cosg*cosg* fib_vol;

			total_vol(0) += fib_vol*fib_vol;

			cos2ksi1(m1) = cosg * cosg;
			sin2ksi = 1 - cos2ksi1(m1);
			cos2theta1(m1) = (cosa*cosa) / sin2ksi;

			cos2theta1(m1) = cos2theta1(m1)*fib_vol*fib_vol;
			cos2ksi1(m1) = cos2ksi1(m1)*fib_vol*fib_vol;
			
			m1++;
		}
		else if (fib_type(i)==2)
		{
			om2(0) += fib_vol * cosa*cosa* fib_vol;
			om2(1) += fib_vol * cosa*cosb* fib_vol;
			om2(2) += fib_vol * cosa*cosg* fib_vol;
			om2(3) += fib_vol * cosb*cosb* fib_vol;
			om2(4) += fib_vol * cosb*cosg* fib_vol;
			om2(5) += fib_vol * cosg*cosg* fib_vol;

			total_vol(1) += fib_vol * fib_vol;

			cos2ksi2(m2) = cosg * cosg;
			sin2ksi = 1 - cos2ksi2(m2);
			cos2theta2(m2) = (cosa*cosa) / sin2ksi;

			cos2theta2(m2) = cos2theta2(m2)*fib_vol*fib_vol;
			cos2ksi2(m2) = cos2ksi2(m2)*fib_vol*fib_vol;

			m2++;
		}
		else if (fib_type(i) == 3)
		{
			om3(0) += fib_vol * cosa*cosa* fib_vol;
			om3(1) += fib_vol * cosa*cosb* fib_vol;
			om3(2) += fib_vol * cosa*cosg* fib_vol;
			om3(3) += fib_vol * cosb*cosb* fib_vol;
			om3(4) += fib_vol * cosb*cosg* fib_vol;
			om3(5) += fib_vol * cosg*cosg* fib_vol;

			total_vol(2) += fib_vol * fib_vol;

			cos2ksi3(m3) = cosg * cosg;
			sin2ksi = 1 - cos2ksi3(m3);
			cos2theta3(m3) = (cosa*cosa) / sin2ksi;

			cos2theta3(m3) = cos2theta3(m3)*fib_vol*fib_vol;
			cos2ksi3(m3) = cos2ksi3(m3)*fib_vol*fib_vol;

			m3++;
		}
	}

	om1 = om1 / total_vol(0);
	om2 = om2 / total_vol(1);
	om3 = om3 / total_vol(2);

	R.segment(0, 6) = om1;
	R.segment(6, 6) = om2;
	R.segment(12, 6) = om3;

	cos2theta1 = cos2theta1.array() / (total_vol(0));
	cos2theta2 = cos2theta2.array() / (total_vol(1));
	cos2theta3 = cos2theta3.array() / (total_vol(2));

	cos2ksi1 = cos2ksi1.array() / (total_vol(0));
	cos2ksi2 = cos2ksi2.array() / (total_vol(1));
	cos2ksi3 = cos2ksi3.array() / (total_vol(2));

	c2ksi(0) = cos2ksi1.mean();
	c2ksi(1) = sqrt((cos2ksi1.array() - cos2ksi1.mean()).square().sum() / (cos2ksi1.size() - 1));
	c2ksi(2) = cos2ksi2.mean();
	c2ksi(3) = sqrt((cos2ksi2.array() - cos2ksi2.mean()).square().sum() / (cos2ksi2.size() - 1));
	c2ksi(4) = cos2ksi3.mean();
	c2ksi(5) = sqrt((cos2ksi3.array() - cos2ksi3.mean()).square().sum() / (cos2ksi3.size() - 1));

	c2theta(0) = cos2theta1.mean();
	c2theta(1) = sqrt((cos2theta1.array() - cos2theta1.mean()).square().sum() / (cos2theta1.size() - 1));
	c2theta(2) = cos2theta2.mean();
	c2theta(3) = sqrt((cos2theta2.array() - cos2theta2.mean()).square().sum() / (cos2theta2.size() - 1));
	c2theta(4) = cos2theta3.mean();
	c2theta(5) = sqrt((cos2theta3.array() - cos2theta3.mean()).square().sum() / (cos2theta3.size() - 1));

}