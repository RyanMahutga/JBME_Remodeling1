
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

using namespace std;
using namespace Eigen;

void appendResults(Vector3d stretch, VectorXd xstress, VectorXd ystress, VectorXd zstress, VectorXd xzstress,
	Vector3d phi, double shear, VectorXi total_failed_types)
{
	string str;
	// saving stretch results
	ofstream myfile1;
	str = "stretch_f.txt";
	myfile1.open(str, ofstream::app);
	myfile1 << stretch.transpose() << endl;
	myfile1.close();

	// saving stress results
	ofstream myfile2;
	str = "xstress_f.txt";
	myfile2.open(str, ofstream::app);
	myfile2 << xstress.transpose() << endl;
	myfile2.close();

	// saving stress results
	ofstream myfile3;
	str = "ystress_f.txt";
	myfile3.open(str, ofstream::app);
	myfile3 << ystress.transpose() << endl;
	myfile3.close();

	// saving stress results
	ofstream myfile4;
	str = "zstress_f.txt";
	myfile4.open(str, ofstream::app);
	myfile4 << zstress.transpose() << endl;
	myfile4.close();

	// saving stress results
	ofstream myfile5;
	str = "xzstress_f.txt";
	myfile5.open(str, ofstream::app);
	myfile5 << xzstress.transpose() << endl;
	myfile5.close();

	ofstream myfile6;
	str = "phi_f.txt";
	myfile6.open(str, ofstream::app);
	myfile6 << phi.transpose() << endl;
	myfile6.close();

	ofstream myfile7;
	str = "shear_f.txt";
	myfile7.open(str, ofstream::app);
	myfile7 << shear << endl;;
	myfile7.close();

	ofstream myfile8;
	str = "failed_f.txt";
	myfile8.open(str, ofstream::app);
	myfile8 << total_failed_types.transpose()<<endl;
	myfile8.close();
}