/* Developed by Harsha Teja Garimella */
#include "functions.h"

using namespace Eigen;

/** \brief This function calculates the updated stress for 3d elements - elastic, hyperelastic material models were implemented so far. */

void fe_stressUpdate_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt) {

	/** This function is  for 2d or 3d elements. For 1d elements, different stress update function was included */

	std::string model;
	model = fe_get_model(opt, "mechanical");

	/** This block develops outputs the updated stress for a 3d elastic material */
	if (model == "simple_elastic") {
		fe_simple_elastic_pbr(sigma, dndx, dndy, dndz, disp_mat, u, opt, return_opt);
	}

	if (model == "saint_venant") {
		fe_saintvenant_elastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block outputs the updated stress for a mooney-rivlin hyperelastic material model */
	if (model == "mooney-rivlin_hyperelastic") {
		fe_mooneyrivlin_hyperelastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block develops outputs the updated stress for a 3d ogden hyperelastic material */
	if (model == "ogden_hyperelastic") {
		// n = number of coefficients - used a double type for consistency
		fe_ogden_hyperelastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

}


void fe_stressUpdate_fiber_pbr(VectorXd& fiber_stress, int opt, double lambda, MatrixXd& T_fiber, MatrixXd& T_fiber_inv, int return_opt) {

	// return_opt = 0 is for PK2 stress with respect to the reference coordinate system
	// return_opt = 1 is for Cauchy stress with respect to the fiber coordinate system

	std::string model;
	model = fe_get_model(opt, "mechanical");
	double principal_cauchy_stress = 0;
	MatrixXd cauchy_stress_fib = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the fiber coordinate system.
	MatrixXd PK2_stress_fib = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the fiber coordinate system.
	MatrixXd cauchy_stress_ref = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the reference coordinate system.
	MatrixXd PK2_stress_ref = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the reference coordinate system.

	if (model == "ogden_hyperelastic") {

	  // The following material properties follow the nomenclature convention used in the ABAQUS USER'S GUIDE VOLUME III: MATERIALS, section 22.5.1.
		double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
		double n = fe_get_mats(opt, 3, "mechanical");
		double mu = fe_get_mats(opt, 4, "mechanical");
		double alpha = fe_get_mats(opt, 5, "mechanical");

		principal_cauchy_stress  = ( 2*mu/alpha )*( pow(lambda,alpha) - pow(lambda,-alpha/2) );

  } else {
    std::cout << "Fiber material must be Ogden for damage to be included." << '\n'; // Later on, we can add the material thresholds for other material models.
    std::cout << "Simulation terminated." << '\n';
    std::exit(1);
  }

	cauchy_stress_fib(0,0) = principal_cauchy_stress;
	PK2_stress_fib(0,0) = principal_cauchy_stress/pow(lambda,2);

	PK2_stress_ref = T_fiber_inv.transpose()*PK2_stress_fib*T_fiber_inv;

	if (return_opt == 0) {
		fiber_stress = fe_tensor2voigt(PK2_stress_ref);
	}
	if (return_opt == 1) {
		fiber_stress = fe_tensor2voigt(cauchy_stress_fib);
	}

}
