#include "functions.h"

using namespace Eigen;

VectorXd fe_saintvenant_elastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt){
	
    VectorXd sigma_local = VectorXd::Zero(6);

	double E = fe_get_mats(opt,1);
	double nu = fe_get_mats(opt,2);
    double mu = fe_get_mats(opt,3);
    double lambda = ((E*nu)/((1+nu)*(1-(2*nu))));

    MatrixXd F = MatrixXd::Zero(3,3); // deformation gradient
	F = fe_calDefGrad(dndx,dndy,dndz,u);
    double defJacobian = F.determinant();
    MatrixXd I = MatrixXd::Identity(3,3);

    MatrixXd GE = ((F.transpose()*F)-I)*(0.5);
    MatrixXd cauchy_sigma = MatrixXd::Zero(3,3);
    cauchy_sigma = ((2*mu*(GE)) + (lambda*(GE.trace())*I));

    MatrixXd pk_S = MatrixXd::Zero(3,3);
	MatrixXd F_inv = F.inverse();
	pk_S = defJacobian*F_inv*cauchy_sigma*F_inv.transpose();

	sigma_local = fe_tensor2voigt(pk_S); /** outputs 2nd cauchy stress tensor in vector form */

	if(return_opt==1){
		sigma_local = fe_tensor2voigt(cauchy_sigma); /** outputs cauchy stress tensor in vector form */
	}

	return sigma_local;
}