#include"functions.h"

using namespace Eigen;

MatrixXd guass_weights_3d(int ndof, int nx,int ny,int nz){
	int n; // n should be the highest value of nx,ny and nz.
	if(nx>=ny){
		if(nx>=nz){
			n = nx;
		}
		else{
			n = nz;
		}
	}
	else{
		if(ny>=nz){
			n = ny;
		}
		else{
		n = nz;
		}
	}

	MatrixXd weights_3d(ndof,n);
	VectorXd tmp(n);

	if(n==2){
		for(int i=0;i<ndof;i++){
			weights_3d(i,0) = guass_weights(2)(0);
			weights_3d(i,1) = guass_weights(2)(1);
		}
	}
	return weights_3d;
}
