#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include "networkDefinitions.h"
#include "dataStructureDefinition.h"
#include <vector>

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void bubbleSort(MatrixXd& x) // this is a bubblesort algorithm
{
	int r = x.rows();

	VectorXd val_holder(4); val_holder = VectorXd::Zero(4); // holder for x values

	for (int i = 0; i < r-1; i++)
	{
		for (int j = 0; j < r-1; j++)
		{
			if (x(j, r-1) > x(j + 1, r-1))
			{
				val_holder(0) = x(j, 0);
				val_holder(1) = x(j, 1);
				val_holder(2) = x(j, 2);
				if (r == 4)
				{
					val_holder(3) = x(j, 3);
				}

				x(j, 0) = x(j + 1, 0);
				x(j, 1) = x(j + 1, 1);
				x(j, 2) = x(j + 1, 2);
				if (r == 4)
				{
					x(j, 3) = x(j + 1, 3);
				}
			
				x(j + 1, 0) = val_holder(0);
				x(j + 1, 1) = val_holder(1);
				x(j + 1, 2) = val_holder(2);
				if (r == 4)
				{
					x(j + 1, 3) = val_holder(3);
				}
			}
		}
	}
}