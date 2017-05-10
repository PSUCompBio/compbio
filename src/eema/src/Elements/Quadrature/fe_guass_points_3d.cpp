#include"functions.h"
using namespace Eigen;

MatrixXd guass_points_3d(int nx,int ny,int nz){
	int n; // n is the highest value of nx,ny and nz.
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
	MatrixXd points_3d(3,2);
	double tmp;
	if(n==2){
		for(int i=0;i<3;i++){
			tmp = -1/sqrt(3);
			points_3d(i,0) = tmp;
			tmp = -tmp;
			points_3d(i,1) = tmp;
		}
	}
	return points_3d;
}
