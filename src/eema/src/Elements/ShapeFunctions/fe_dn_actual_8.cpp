
#include"functions.h"

using namespace Eigen;

void fe_dndx_8_pbr(VectorXd& dndx, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian) {
	for (int i = 0; i < nnel; i++) {
		dndx(i) = (invJacobian(0, 0) * dndr(i)) + (invJacobian(0, 1) * dnds(i)) + (invJacobian(0, 2) * dndt(i));
		//dndy(i) = (invJacobian(1,0)*dndr(i)) + (invJacobian(1,1)*dnds(i)) + (invJacobian(1,2)*dndt(i));
		//dndz(i) = (invJacobian(2,0)*dndr(i)) + (invJacobian(2,1)*dnds(i)) + (invJacobian(2,2)*dndt(i));
	}
}

void fe_dndy_8_pbr(VectorXd& dndy, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian) {
	for (int i = 0; i < nnel; i++) {
		//dndx(i) = (invJacobian(0,0)*dndr(i)) + (invJacobian(0,1)*dnds(i)) + (invJacobian(0,2)*dndt(i));
		dndy(i) = (invJacobian(1, 0) * dndr(i)) + (invJacobian(1, 1) * dnds(i)) + (invJacobian(1, 2) * dndt(i));
		//dndz(i) = (invJacobian(2,0)*dndr(i)) + (invJacobian(2,1)*dnds(i)) + (invJacobian(2,2)*dndt(i));
	}
}

void fe_dndz_8_pbr(VectorXd& dndz, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian) {
	for (int i = 0; i < nnel; i++) {
		//dndx(i) = (invJacobian(0,0)*dndr(i)) + (invJacobian(0,1)*dnds(i)) + (invJacobian(0,2)*dndt(i));
		//dndy(i) = (invJacobian(1,0)*dndr(i)) + (invJacobian(1,1)*dnds(i)) + (invJacobian(1,2)*dndt(i));
		dndz(i) = (invJacobian(2, 0) * dndr(i)) + (invJacobian(2, 1) * dnds(i)) + (invJacobian(2, 2) * dndt(i));
	}
}



void fe_dndxyz_8_pbr_array(double* dndx, double* dndy, double* dndz, int nnel, double* dndr, double* dnds, double* dndt, double** invJacobian) {
	for (int i = 0; i < nnel; i++) {
		dndx[i] = (invJacobian[0][0] * dndr[i]) + (invJacobian[0][1] * dnds[i]) + (invJacobian[0][2] * dndt[i]);
		dndy[i] = (invJacobian[1][0] * dndr[i]) + (invJacobian[1][1] * dnds[i]) + (invJacobian[1][2] * dndt[i]);
		dndz[i] = (invJacobian[2][0] * dndr[i]) + (invJacobian[2][1] * dnds[i]) + (invJacobian[2][2] * dndt[i]);
	}
}
