#include "functions.h"

using namespace Eigen;

void fe_stressUpdateViscoelasticity_pbr(VectorXd& instantStress, double dT, MatrixXd& defGrad, MatrixXd invDefGrad, double defJacobian, int iterator, int ix, int iy, int iz, int opt, int return_opt)
{

	// Primary Reference: "Formulation and implementation of three-dimensional viscoelasticity at small and finite strains," (Kaliske, 1997).
	// Secondary Reference: "ABAQUS 6.14 THEORY GUIDE," Section 4.8.2 "FINITE-STRAIN VISCOELASTICITY"

	// identify the strain energy function being used
	// extract the necessary material properties
	std::string model;
	model = fe_get_model(opt, "mechanical");
	double D1 = 0.0;
	if (model == "mooney-rivlin_hyperelastic") {
		D1 = 2.0 / fe_get_mats(opt, 1, "mechanical");
	} else if (model == "ogden_hyperelastic") {
		D1 = 2.0 / fe_get_mats(opt, 1, "mechanical");
	} else {
		std::cout << "Material must be Mooney-Rivlin or Ogden for viscoelasticity to be included." << '\n'; // Later on, we can add other material models.
		std::cout << "Simulation terminated." << '\n';
		std::exit(1);
	}

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
	MatrixXd modifiedKStress = MatrixXd::Zero(ndof, ndof);

	MatrixXd internalStressVariable1_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd internalStressVariable2_prev = MatrixXd::Zero(ndof, ndof);
	MatrixXd devInstantPk2Stress_prev = MatrixXd::Zero(ndof, ndof);

	// return_opt = 0 for quadrature points used in getforce subroutine, input and output stress are PK2
	if (return_opt == 0) {
		instantPk2Stress = fe_voigt2tensor(instantStress); // convert instantaneous PK2 stress from vector to matrix
		// retrive data from previous time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev(i, j) = internalStressVariable1_prev_normal_store[iterator][ix][iy][iz][i][j]; // internal stress variable H for i = 1
				internalStressVariable2_prev(i, j) = internalStressVariable2_prev_normal_store[iterator][ix][iy][iz][i][j]; // internal stress variable H for i = 2
				devInstantPk2Stress_prev(i, j) = devInstantPK2Stress_prev_normal_store[iterator][ix][iy][iz][i][j]; // deviatoric part of instantaneous PK2 stress
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
				internalStressVariable1_prev(i, j) = internalStressVariable1_prev_centroid_store[iterator][i][j]; // internal stress variable H for i = 1
				internalStressVariable2_prev(i, j) = internalStressVariable2_prev_centroid_store[iterator][i][j]; // internal stress variable H for i = 2
				devInstantPk2Stress_prev(i, j) = devInstantPK2Stress_prev_centroid_store[iterator][i][j]; // deviatoric part of instantaneous PK2 stress
			}
		}
	}

	MatrixXd rightCauchyGreen = defGrad.transpose() * defGrad; // right Cauchy-Green deformation tensor C
	MatrixXd invRightCauchyGreen = MatrixXd::Zero(ndof, ndof); // inverse of right Cauchy-Green deformation tensor C
	fe_invMatrix_pbr(invRightCauchyGreen, rightCauchyGreen);

	// MatrixXd sphInstantPk2Stress = (1.0 / 3.0) * (rightCauchyGreen * instantPk2Stress).trace() * invRightCauchyGreen; // not used...delete later
	// MatrixXd devInstantPk2Stress = instantPk2Stress - sphInstantPk2Stress; // not used...delete later

	MatrixXd devInstantPk2Stress = instantPk2Stress - (1.0 / 3.0) * (rightCauchyGreen * instantPk2Stress).trace() * invRightCauchyGreen; // deviatoric part of instantaneous PK2 stress in reference configuration
	MatrixXd internalStressVariable1 = exp(-dT/tau_1) * internalStressVariable1_prev + gamma_1 * ( (1 - exp(-dT/tau_1)) / (dT/tau_1) ) * (devInstantPk2Stress - devInstantPk2Stress_prev); // internal stress variable H for i = 1
	MatrixXd internalStressVariable2 = exp(-dT/tau_2) * internalStressVariable2_prev + gamma_2 * ( (1 - exp(-dT/tau_2)) / (dT/tau_2) ) * (devInstantPk2Stress - devInstantPk2Stress_prev); // internal stress variable H for i = 2
	MatrixXd modifiedDevPk2Stress = devInstantPk2Stress + internalStressVariable1 + internalStressVariable2; // modified deviatoric part of PK2 stress in reference configuration
	// modifiedPk2Stress = sphInstantPk2Stress + modifiedDevPk2Stress; // modified PK2 stress
	MatrixXd modifiedDevKStress = defGrad * modifiedDevPk2Stress * defGrad.transpose(); // modified deviatoric part of Kirchhoff stress in reference configuration

	// calculate the modified total Kirchoff stress
	// the equation depends on the strain energy function
	if (model == "mooney-rivlin_hyperelastic") {
		modifiedKStress = defJacobian * (2.0 / D1) * (defJacobian - 1) * I + modifiedDevKStress;
	} else if (model == "ogden_hyperelastic") {
		modifiedKStress = defJacobian * (2.0 / D1) * (defJacobian - 1) * I + modifiedDevKStress;
	} else {
		std::cout << "Material must be Mooney-Rivlin or Ogden for viscoelasticity to be included." << '\n'; // Later on, we can add other material models.
		std::cout << "Simulation terminated." << '\n';
		std::exit(1);
	}

	if (return_opt == 0) {
		// calculate output stress
		modifiedPk2Stress = invDefGrad * modifiedKStress * invDefGrad.transpose(); // modified PK2 stress
		instantStress = fe_tensor2voigt(modifiedPk2Stress); // convert modified PK2 stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev_normal_store[iterator][ix][iy][iz][i][j] = internalStressVariable1(i, j); // internal stress variable H for i = 1
				internalStressVariable2_prev_normal_store[iterator][ix][iy][iz][i][j] = internalStressVariable2(i, j); // internal stress variable H for i = 2
				devInstantPK2Stress_prev_normal_store[iterator][ix][iy][iz][i][j] = devInstantPk2Stress(i, j); // deviatoric part of instantaneous PK2 stress
			}
		}
	}

	if (return_opt == 1) {
		// calculate output stress
		modifiedCauchyStress = (1 / defJacobian) * modifiedKStress; // modified Cauchy stress
		// modifiedCauchyStress = (1 / defJacobian) * defGrad * modifiedPk2Stress * defGrad.transpose(); // not used...delete later
		instantStress = fe_tensor2voigt(modifiedCauchyStress); // convert modified Cauchy stress from matrix to vector
		// store current data for next time Step
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				internalStressVariable1_prev_centroid_store[iterator][i][j] = internalStressVariable1(i, j); // internal stress variable H for i = 1
				internalStressVariable2_prev_centroid_store[iterator][i][j] = internalStressVariable2(i, j); // internal stress variable H for i = 2
				devInstantPK2Stress_prev_centroid_store[iterator][i][j] = devInstantPk2Stress(i, j); // deviatoric part of instantaneous PK2 stress
			}
		}
	}
}
