
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

void appendResults(Vector3d stretch, Vector3d stretch0, VectorXd xstress, VectorXd ystress, VectorXd zstress, VectorXd xzstress,
	Vector3d phi,double time, VectorXd Orient, VectorXd c2k, VectorXd c2t,double press0, double tau_wall)
{
	string str;

	// saving stretch results
	ofstream myfile1;
	str = "stretch.txt";
	myfile1.open(str, ofstream::app);
	myfile1 << stretch.transpose() << endl;
	myfile1.close();

	ofstream myfile9;
	str = "stretch0.txt";
	myfile1.open(str, ofstream::app);
	myfile1 << stretch0.transpose() << endl;
	myfile1.close();

	// saving stress results
	ofstream myfile2;
	str = "xstress.txt";
	myfile2.open(str, ofstream::app);
	myfile2 << xstress.transpose() << endl;
	myfile2.close();

	// saving stress results
	ofstream myfile3;
	str = "ystress.txt";
	myfile3.open(str, ofstream::app);
	myfile3 << ystress.transpose() << endl;
	myfile3.close();

	// saving stress results
	ofstream myfile4;
	str = "zstress.txt";
	myfile4.open(str, ofstream::app);
	myfile4 << zstress.transpose() << endl;
	myfile4.close();

	// saving stress results
	ofstream myfile5;
	str = "xzstress.txt";
	myfile5.open(str, ofstream::app);
	myfile5 << xzstress.transpose() << endl;
	myfile5.close();

	ofstream myfile6;
	str = "phi.txt";
	myfile6.open(str, ofstream::app);
	myfile6 << phi.transpose() << endl;
	myfile6.close();

	ofstream myfile7;
	str = "time.txt";
	myfile7.open(str, ofstream::app);
	myfile7 << time << endl;
	myfile7.close();

	ofstream myfile8;
	str = "orientation.txt";
	myfile8.open(str, ofstream::app);
	myfile8 << Orient.transpose() << endl;
	myfile8.close();

	ofstream myfile11;
	str = "cos2ksi.txt";
	myfile11.open(str, ofstream::app);
	myfile11 << c2k.transpose() << endl;
	myfile11.close();

	ofstream myfile12;
	str = "cos2theta.txt";
	myfile12.open(str, ofstream::app);
	myfile12 << c2t.transpose() << endl;
	myfile12.close();

	ofstream myfile10;
	str = "pressure0.txt";
	myfile8.open(str, ofstream::app);
	myfile8 << press0<<endl;
	myfile8.close();

	ofstream myfile13;
	str = "WSS.txt";
	myfile13.open(str, ofstream::app);
	myfile13 << tau_wall << endl;
	myfile13.close();
}