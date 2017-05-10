
#include"functions.h"

using namespace Eigen;

VectorXd fe_dndx_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian){

	VectorXd dndx(nnel);
	for(int i=0;i<nnel;i++){
		dndx(i) = (invJacobian(0,0)*dndr(i)) + (invJacobian(0,1)*dnds(i)) + (invJacobian(0,2)*dndt(i));
		//dndy(i) = (invJacobian(1,0)*dndr(i)) + (invJacobian(1,1)*dnds(i)) + (invJacobian(1,2)*dndt(i));
		//dndz(i) = (invJacobian(2,0)*dndr(i)) + (invJacobian(2,1)*dnds(i)) + (invJacobian(2,2)*dndt(i));
	}
	return dndx;
}

VectorXd fe_dndy_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian){

	VectorXd dndy(nnel);
	for(int i=0;i<nnel;i++){
		//dndx(i) = (invJacobian(0,0)*dndr(i)) + (invJacobian(0,1)*dnds(i)) + (invJacobian(0,2)*dndt(i));
		dndy(i) = (invJacobian(1,0)*dndr(i)) + (invJacobian(1,1)*dnds(i)) + (invJacobian(1,2)*dndt(i));
		//dndz(i) = (invJacobian(2,0)*dndr(i)) + (invJacobian(2,1)*dnds(i)) + (invJacobian(2,2)*dndt(i));
	}
	return dndy;
}

VectorXd fe_dndz_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian){

	VectorXd dndz(nnel);
	for(int i=0;i<nnel;i++){
		//dndx(i) = (invJacobian(0,0)*dndr(i)) + (invJacobian(0,1)*dnds(i)) + (invJacobian(0,2)*dndt(i));
		//dndy(i) = (invJacobian(1,0)*dndr(i)) + (invJacobian(1,1)*dnds(i)) + (invJacobian(1,2)*dndt(i));
		dndz(i) = (invJacobian(2,0)*dndr(i)) + (invJacobian(2,1)*dnds(i)) + (invJacobian(2,2)*dndt(i));
	}
	return dndz;
}
