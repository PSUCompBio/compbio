#include "functions.h"

using namespace Eigen;

void fe_stressUpdateViscoelasticity_pbr(VectorXd& instantStress, double dT, MatrixXd& defGrad, MatrixXd invDefGrad, double defJacobian, int i_normal, int intx_normal, int inty_normal, int intz_normal, int return_opt)
{

	// Primary Reference: "Formulation and implementation of three-dimensional viscoelasticity at small and finite strains," (Kaliske, 1997).
	// Secondary Reference: "ABAQUS 6.14 THEORY GUIDE," Section 4.8.2 "FINITE-STRAIN VISCOELASTICITY"

	// define Prony series parameters based on (Garimella, 2016)
	double g_1 = 0.65425; // i = 1
	double g_2 = 0.0149; // i = 2
	double tau_1 = 0.0066940; // i = 1
	double tau_2 = 0.15642; // i = 2

	// convert Prony series parameters from ABAQUS to notation used in (Kaliske, 1997)
	double g_infinity = 1.0 - g_1 - g_2;
	double gamma_1 = g_1 / g_infinity;
	double gamma_2 = g_2 / g_infinity;

	MatrixXd instantPk2Stress = MatrixXd::Zero(ndof, ndof);
	MatrixXd instantCauchyStress = MatrixXd::Zero(ndof, ndof);

	MatrixXd modifiedPk1Stress = MatrixXd::Zero(ndof, ndof);
	MatrixXd modifiedPk2Stress = MatrixXd::Zero(ndof, ndof);
	MatrixXd modifiedCauchyStress = MatrixXd::Zero(ndof, ndof);

	MatrixXd internalStressVariable1_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd internalStressVariable2_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd devInstantPk2Stress_prev = MatrixXd::Zero(ndof, ndof);

	// return_opt = 0 for quadrature points used in getforce subroutine, input and output stress are PK2
	if (return_opt == 0) {
		instantPk2Stress = fe_voigt2tensor(instantStress); // convert instantaneous PK2 stress from vector to matrix
		// retrive data from previous time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev(i, j) = internalStressVariable1_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // internal stress variable H for i = 1
				internalStressVariable2_prev(i, j) = internalStressVariable2_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // internal stress variable H for i = 2
				devInstantPk2Stress_prev(i, j) = devInstantPK2Stress_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j]; // deviatoric part of instantaneous PK2 stress
			}
		}
	}

	// return_opt = 1 for centroid stress calculation, input and output stress are Cauchy
	if (return_opt == 1) {
		instantCauchyStress = fe_voigt2tensor(instantStress); // convert instantaneous Cauchy stress from vector to matrix
		instantPk2Stress = defJacobian * invDefGrad * instantCauchyStress * invDefGrad.transpose(); // instantaneous PK2 stress
		// retrive data from previous time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev(i, j) = internalStressVariable1_prev_centroid_store[i_normal][i][j]; // internal stress variable H for i = 1
				internalStressVariable2_prev(i, j) = internalStressVariable2_prev_centroid_store[i_normal][i][j]; // internal stress variable H for i = 2
				devInstantPk2Stress_prev(i, j) = devInstantPK2Stress_prev_centroid_store[i_normal][i][j]; // deviatoric part of instantaneous PK2 stress
			}
		}
	}

	MatrixXd rightCauchyGreen = defGrad.transpose() * defGrad; // right Cauchy-Green deformation tensor C
	MatrixXd invRightCauchyGreen = MatrixXd::Zero(ndof, ndof); // inverse of right Cauchy-Green deformation tensor C
	fe_invMatrix_pbr(invRightCauchyGreen, rightCauchyGreen);
	MatrixXd sphInstantPk2Stress = (1.0 / 3.0) * (rightCauchyGreen * instantPk2Stress).trace() * invRightCauchyGreen; // spherical part of instantaneous PK2 stress in reference configuration
	MatrixXd devInstantPk2Stress = instantPk2Stress - sphInstantPk2Stress; // deviatoric part of instantaneous PK2 stress in reference configuration

	// MatrixXd devInstantPk2Stress = instantPk2Stress - (1.0 / 3.0) * (rightCauchyGreen * instantPk2Stress).trace() * invRightCauchyGreen; // deviatoric part of instantaneous PK2 stress in reference configuration - not used...

	MatrixXd internalStressVariable1 = exp(-dT/tau_1) * internalStressVariable1_prev + gamma_1 * ( (1 - exp(-dT/tau_1)) / (dT/tau_1) ) * (devInstantPk2Stress - devInstantPk2Stress_prev); // internal stress variable H for i = 1
	MatrixXd internalStressVariable2 = exp(-dT/tau_2) * internalStressVariable2_prev + gamma_2 * ( (1 - exp(-dT/tau_2)) / (dT/tau_2) ) * (devInstantPk2Stress - devInstantPk2Stress_prev); // internal stress variable H for i = 2

	MatrixXd modifiedDevPk2Stress = devInstantPk2Stress + internalStressVariable1 + internalStressVariable2; // modified deviatoric part of PK2 stress in reference configuration
	modifiedPk2Stress = sphInstantPk2Stress + modifiedDevPk2Stress; // modified PK2 stress

	// MatrixXd modifiedDevPk1Stress = defGrad * modifiedDevPk2Stress * defGrad.transpose(); // modified deviatoric part of PK1 stress in reference configuration - not used...
	// modifiedPk1Stress = defJacobian * (2.0/0.9091e-9) * (defJacobian - 1) * I + modifiedDevPk1Stress; // modified PK1 stress - not used...

	if (return_opt == 0) {
		// calculate output stress
		// modifiedPk2Stress = modifiedPk1Stress * invDefGrad.transpose(); // modified PK2 stress - not used...
		instantStress = fe_tensor2voigt(modifiedPk2Stress); // convert modified PK2 stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = internalStressVariable1(i, j); // internal stress variable H for i = 1
				internalStressVariable2_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = internalStressVariable2(i, j); // internal stress variable H for i = 2
				devInstantPK2Stress_prev_normal_store[i_normal][intx_normal][inty_normal][intz_normal][i][j] = devInstantPk2Stress(i, j); // deviatoric part of instantaneous PK2 stress
			}
		}
	}

	if (return_opt == 1) {
		// calculate output stress
		// modifiedCauchyStress = (1 / defJacobian) * defGrad * modifiedPk1Stress; // modified Cauchy stress - not used...
		modifiedCauchyStress = (1 / defJacobian) * defGrad * modifiedPk2Stress * defGrad.transpose(); // modified Cauchy stress
		instantStress = fe_tensor2voigt(modifiedCauchyStress); // convert modified Cauchy stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev_centroid_store[i_normal][i][j] = internalStressVariable1(i, j); // internal stress variable H for i = 1
				internalStressVariable2_prev_centroid_store[i_normal][i][j] = internalStressVariable2(i, j); // internal stress variable H for i = 2
				devInstantPK2Stress_prev_centroid_store[i_normal][i][j] = devInstantPk2Stress(i, j); // deviatoric part of instantaneous PK2 stress
			}
		}
	}

}

