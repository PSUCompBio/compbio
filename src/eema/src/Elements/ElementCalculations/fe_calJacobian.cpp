

#include"functions.h"

using namespace Eigen;

MatrixXd fe_calJacobian(int dim, int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord){
	
	MatrixXd jacobian = MatrixXd::Zero(dim,dim);

	if(dim==3){
		for(int k=0;k<nnel;k++){
			jacobian(0,0) = jacobian(0,0) + dndr(k)*xcoord(k);
			jacobian(0,1) = jacobian(0,1) + dndr(k)*ycoord(k);
			jacobian(0,2) = jacobian(0,2) + dndr(k)*zcoord(k);
			jacobian(1,0) = jacobian(1,0) + dnds(k)*xcoord(k);
			jacobian(1,1) = jacobian(1,1) + dnds(k)*ycoord(k);
			jacobian(1,2) = jacobian(1,2) + dnds(k)*zcoord(k);
			jacobian(2,0) = jacobian(2,0) + dndt(k)*xcoord(k);
			jacobian(2,1) = jacobian(2,1) + dndt(k)*ycoord(k);
			jacobian(2,2) = jacobian(2,2) + dndt(k)*zcoord(k);
		}			
	}

	return jacobian;
}
