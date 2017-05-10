#include "functions.h"

using namespace Eigen;

VectorXd fe_ogden_hyperelastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt){

    VectorXd sigma_local = VectorXd::Zero(ndof*2);

    	MatrixXd F = MatrixXd::Zero(ndof,ndof); // deformation gradient
		F = fe_calDefGrad(dndx,dndy,dndz,u);
		double defJacobian = F.determinant();

		// fe_display_matrix(F);

		MatrixXd I = MatrixXd::Identity(ndof,ndof);

		MatrixXd B = MatrixXd::Zero(ndof,ndof); // Left Cauchy-Green Deformation Tensor
		B = F*F.transpose();
		EigenSolver<MatrixXd> es(B);

		MatrixXd D = es.pseudoEigenvalueMatrix();
		MatrixXd V = es.pseudoEigenvectors();

		double l_1 = D(0,0);
		double l_2 = D(1,1);
		double l_3 = D(2,2);

		VectorXd n_1 = VectorXd::Zero(3);
		n_1 << V(0,0), V(1,0), V(2,0);
		VectorXd n_2 = VectorXd::Zero(3);
		n_2 << V(0,1), V(1,1), V(2,1);
		VectorXd n_3 = VectorXd::Zero(3);
		n_3 << V(0,2), V(1,2), V(2,2);


		double D_1 = fe_get_mats(3,1)/2;
		double p = -2*D_1*(defJacobian-1);

		double n = fe_get_mats(opt,3);
		int counter = 4;
		double W_1 = 0;
		double W_2 = 0;
		double W_3 = 0;
		for(double i=0.0;i<n;){
			double mu = fe_get_mats(opt,counter);
			double alpha = fe_get_mats(3,counter+1);
			W_1 = W_1 + (mu*pow(l_1,(alpha-1)));
			W_2 = W_2 + (mu*pow(l_2,(alpha-1)));
			W_3 = W_3 + (mu*pow(l_3,(alpha-1)));
			i = i+1.0;
			counter = counter+2;
		}
		MatrixXd cauchy_sigma = MatrixXd::Zero(ndof,ndof);

		cauchy_sigma = -p*I + (l_1/(l_1*l_2*l_3))*W_1*(n_1*(n_1.transpose())) + (l_2/(l_1*l_2*l_3))*W_2*(n_2*(n_2.transpose())) + (l_3/(l_1*l_2*l_3))*W_3*(n_3*(n_3.transpose()));

		MatrixXd pk_S = MatrixXd::Zero(3,3);
		MatrixXd F_inv = F.inverse();
		pk_S = defJacobian*F_inv*cauchy_sigma*F_inv.transpose();

		sigma_local = fe_tensor2voigt(pk_S);

		if(return_opt==1){
			sigma_local = fe_tensor2voigt(cauchy_sigma);
		}

    return sigma_local;
}

