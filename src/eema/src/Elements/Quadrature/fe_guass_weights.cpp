#include"functions.h"

using namespace Eigen;

VectorXd guass_weights(int n){

	VectorXd weights(n);
	double tmp;

	if(n==2){
		tmp = 1.0;
		weights(0) = tmp;
		tmp = 1.0;
		weights(1) = tmp;
	}
	return weights;
}
