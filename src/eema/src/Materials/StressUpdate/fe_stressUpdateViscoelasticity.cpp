#include "functions.h"

using namespace Eigen;

void fe_stressUpdateViscoelasticity_pbr(VectorXd& instantStress, double dT, MatrixXd& defGrad, MatrixXd invDefGrad, double defJacobian, int i_normal, int intx_normal, int inty_normal, int intz_normal, int return_opt)
{

	// Reference: "ABAQUS 6.14 THEORY GUIDE," Section 4.8.2 "FINITE-STRAIN VISCOELASTICITY"

	// define Prony series parameters (Garimella, 2016)
	double g_1 = 0.65425; // i = 1
	double g_2 = 0.0149; // i = 2
	double tau_1 = 0.0066940; // i = 1
	double tau_2 = 0.15642; // i = 2

	MatrixXd instantPk1Stress = MatrixXd::Zero(ndof, ndof);
	MatrixXd instantPk2Stress = MatrixXd::Zero(ndof, ndof);
	MatrixXd instantCauchyStress = MatrixXd::Zero(ndof, ndof);
	MatrixXd invDistGrad_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd devInstantPk1Stress_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd devInternalPk1Stress1_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd devInternalPk1Stress2_prev = MatrixXd::Zero(ndof, ndof);

	// return_opt = 0 for quadrature points used in getforce subroutine, input and output stress are PK2
	if (return_opt == 0) {
		instantPk2Stress = fe_voigt2tensor(instantStress); // convert instantaneous PK2 stress from vector to matrix

		std::cout << "BEFORE: instantPk2Stress = " << '\n' << instantPk2Stress << '\n';

		instantPk1Stress = instantPk2Stress * defGrad.transpose(); // instantaneous PK1 stress
		// retrive data from previous time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				invDistGrad_prev(i, j) = invDistGrad_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // inverse distortion deformation gradient
				devInstantPk1Stress_prev(i, j) = devInstantPk1Stress_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // deviatoric part of instantaneous PK1 stress
				devInternalPk1Stress1_prev(i, j) = devInternalPk1Stress1_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // deviatoric internal PK1 stress for i = 1
				devInternalPk1Stress2_prev(i, j) = devInternalPk1Stress2_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // deviatoric internal PK1 stress for i = 2
			}
		}
	}

	// return_opt = 1 for centroid stress calculation, input and output stress are Cauchy
	if (return_opt == 1) {
		instantCauchyStress = fe_voigt2tensor(instantStress); // convert instantaneous Cauchy stress from vector to matrix
		instantPk1Stress = defJacobian * invDefGrad * instantCauchyStress; // instantaneous PK1 stress
		// retrive data from previous time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				invDistGrad_prev(i, j) = invDistGrad_prev_centroid_store[i_normal][i][j]; // inverse distortion deformation gradient
				devInstantPk1Stress_prev(i, j) = devInstantPk1Stress_prev_centroid_store[i_normal][i][j]; // deviatoric part of instantaneous PK1 stress
				devInternalPk1Stress1_prev(i, j) = devInternalPk1Stress1_prev_centroid_store[i_normal][i][j]; // deviatoric internal PK1 stress for i = 1
				devInternalPk1Stress2_prev(i, j) = devInternalPk1Stress2_prev_centroid_store[i_normal][i][j]; // deviatoric internal PK1 stress for i = 2
			}
		}
	}

	MatrixXd devInstantPk1Stress = instantPk1Stress - (1.0/3.0) * instantPk1Stress.trace() * I; // deviatoric part of instantaneous PK1 stress
	double gamma_1 = exp(-dT/tau_1);
	double gamma_2 = exp(-dT/tau_2);
	double alpha_1 = 1 - (tau_1/dT) * (1 - gamma_1);
	double alpha_2 = 1 - (tau_2/dT) * (1 - gamma_2);
	double beta_1 = (tau_1/dT) * (1 - gamma_1) - gamma_1;
	double beta_2 = (tau_2/dT) * (1 - gamma_2) - gamma_2;
	MatrixXd distGrad = defGrad / pow(defJacobian, 1.0/3.0); // distortion deformation gradient
	MatrixXd relDistGrad = distGrad * invDistGrad_prev; // distorional deformation gradient relative to previous time step
	MatrixXd devInstantPk1Stress_hat = relDistGrad * devInstantPk1Stress_prev * relDistGrad.transpose(); //deviatoric part of instantaneous PK1 stress hat
	MatrixXd devInternalPk1Stress1_hat = relDistGrad * devInternalPk1Stress1_prev * relDistGrad.transpose(); // deviatoric internal PK1 stress hat for i = 1
	MatrixXd devInternalPk1Stress2_hat = relDistGrad * devInternalPk1Stress2_prev * relDistGrad.transpose(); // deviatoric internal PK1 stress hat for i = 2
	MatrixXd devInternalPk1Stress1 = alpha_1 * g_1 * devInstantPk1Stress + beta_1 * g_1 * devInstantPk1Stress_hat + gamma_1 * devInternalPk1Stress1_hat; // deviatoric internal PK1 stress for i = 1
	MatrixXd devInternalPk1Stress2 = alpha_2 * g_2 * devInstantPk1Stress + beta_2 * g_2 * devInstantPk1Stress_hat + gamma_2 * devInternalPk1Stress2_hat; // deviatoric internal PK1 stress for i = 1

	instantPk1Stress = instantPk1Stress - (devInternalPk1Stress1 + devInternalPk1Stress2); // modified PK1 stress

	MatrixXd invDistGrad = MatrixXd::Zero(ndof, ndof);
	fe_invMatrix_pbr(invDistGrad, distGrad);

	if (return_opt == 0) {
		// calculate output stress
		instantPk2Stress = instantPk1Stress * invDefGrad.transpose(); // modified PK2 stress

		std::cout << "AFTER: instantPk2Stress = " << '\n' << instantPk2Stress << '\n';

		instantStress = fe_tensor2voigt(instantPk2Stress); // convert instantaneous PK2 stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				invDistGrad_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = invDistGrad(i, j); // inverse distortion deformation gradient
				devInstantPk1Stress_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = devInstantPk1Stress(i, j); // deviatoric part of instantaneous PK1 stress
				devInternalPk1Stress1_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = devInternalPk1Stress1(i, j); // deviatoric internal PK1 stress for i = 1
				devInternalPk1Stress2_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = devInternalPk1Stress2(i, j); // deviatoric internal PK1 stress for i = 2
			}
		}
	}

	if (return_opt == 1) {
		// calculate output stress
		instantCauchyStress = (1 / defJacobian) * defGrad * instantPk1Stress; // modified Cauchy stress
		instantStress = fe_tensor2voigt(instantCauchyStress); // convert instantaneous Cauchy stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				invDistGrad_prev_centroid_store[i_normal][i][j] = invDistGrad(i, j);  // inverse distortion deformation gradient
				devInstantPk1Stress_prev_centroid_store[i_normal][i][j] = devInstantPk1Stress(i, j); // deviatoric part of instantaneous PK1 stress
				devInternalPk1Stress1_prev_centroid_store[i_normal][i][j] = devInternalPk1Stress1(i, j); // deviatoric internal PK1 stress for i = 1
				devInternalPk1Stress2_prev_centroid_store[i_normal][i][j] = devInternalPk1Stress2(i, j);  // deviatoric internal PK1 stress for i = 2
			}
		}
	}

}

