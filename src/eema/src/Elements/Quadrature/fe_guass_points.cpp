#include"functions.h"

using namespace Eigen;

VectorXd guass_points(int n){

	VectorXd points(n);
	double tmp;
	if (n==2){
		tmp = -1/sqrt(3);
		points(0) = tmp;
		tmp = -tmp;
		points(1) = tmp;
	}

	return points;
}
