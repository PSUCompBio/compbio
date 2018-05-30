#include "functions.h"

using namespace Eigen;

void fe_stressUpdateDamage_pbr(VectorXd& fiber_stress_damaged, int opt, int fib, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double lambda, VectorXd& lambda_min, VectorXd& lambda_max, int return_opt)
{

	// return_opt = 0 is for PK2 stress with respect to the reference coordinate system
	// return_opt = 1 is for Cauchy stress with respect to the fiber coordinate system

	std::string model;
	model = fe_get_model(opt, "mechanical");
	double r = 0;  // material constant
	double m = 0;  // material constant
	double principal_cauchy_stress_primary = 0;
	double principal_cauchy_stress = 0;
	double U_density_primary = 0;
	double U_density_max = 0;
	MatrixXd cauchy_stress_fib = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the fiber coordinate system.
	MatrixXd PK2_stress_fib = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the fiber coordinate system.
	MatrixXd cauchy_stress_ref = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the reference coordinate system.
	MatrixXd PK2_stress_ref = MatrixXd::Zero(ndof, ndof); // Stress components with respect to the reference coordinate system.

	if (model == "ogden_hyperelastic") {
    r = 1.25;
    m = 3e3;

	  // The following material properties follow the nomenclature convention used in the ABAQUS USER'S GUIDE VOLUME III: MATERIALS, section 22.5.1.
		double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
		double n = fe_get_mats(opt, 3, "mechanical");
		double mu = fe_get_mats(opt, 4, "mechanical");
		double alpha = fe_get_mats(opt, 5, "mechanical");

		principal_cauchy_stress_primary  = ( 2*mu/alpha )*( pow(lambda,alpha) - pow(lambda,-alpha/2) );

  } else {
    std::cout << "Fiber material must be Ogden for damage to be included." << '\n'; // Later on, we can add the material thresholds for other material models.
    std::cout << "Simulation terminated." << '\n';
    std::exit(1);
  }

	if ( lambda <= lambda_min(fib) || lambda >= lambda_max(fib) ) {
		principal_cauchy_stress = principal_cauchy_stress_primary;
	}

	if ( lambda_min(fib) < lambda < lambda_max(fib) ) {
		U_density_primary = fe_calStrainEnergyDensity(opt, lambda);
		if (lambda > 1) {
			U_density_max = fe_calStrainEnergyDensity(opt, lambda_max(fib));
		}
		if (lambda < 1) {
			U_density_max = fe_calStrainEnergyDensity(opt, lambda_min(fib));
		}
		double etta = 1 - (1/r)*erf( (1/m)*(U_density_max - U_density_primary) );
		principal_cauchy_stress = etta*principal_cauchy_stress_primary;
	}

	cauchy_stress_fib(0,0) = principal_cauchy_stress;
	PK2_stress_fib(0,0) = principal_cauchy_stress/pow(lambda,2);
	MatrixXd T = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3); // transformation matrix from reference coordinate system to fiber coordinate system
	MatrixXd T_inv = MatrixXd::Zero(ndof, ndof);
	fe_invMatrix_pbr(T_inv, T);
	// cauchy_stress_ref = T_inv.transpose()*cauchy_stress_fib*T_inv;
	PK2_stress_ref = T_inv.transpose()*PK2_stress_fib*T_inv;

	if ( lambda > lambda_max(fib) ) {
		lambda_max(fib) = lambda;
	}
	if ( lambda < lambda_min(fib) ) {
		lambda_min(fib) = lambda;
	}

	if (return_opt == 0) {
		fiber_stress_damaged = fe_tensor2voigt(PK2_stress_ref);
	}
	if (return_opt == 1) {
		fiber_stress_damaged = fe_tensor2voigt(cauchy_stress_fib);
	}

}

double fe_calStrainEnergyDensity(int opt, double lambda)
{

  // Here we assume the fiber is incompressible and subjected only to applied axial stretch.
  double U_density = 0;
  std::string model;
  model = fe_get_model(opt, "mechanical");

  if (model == "ogden_hyperelastic") {
	  // The following material properties follow the nomenclature convention used in the ABAQUS USER'S GUIDE VOLUME III: MATERIALS, section 22.5.1.
    double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
    double n = fe_get_mats(opt, 3, "mechanical");
    double mu = fe_get_mats(opt, 4, "mechanical");
    double alpha = fe_get_mats(opt, 5, "mechanical");

    U_density = ( 2*mu/pow(alpha,2) )*( pow(lambda,alpha) + 2*pow(sqrt(1/lambda),alpha) - 3 ); // Strain energy per unit reference volume.

	} else {
    std::cout << "Fiber material must be Ogden for damage to be included." << '\n'; // Later on, we can add the strain energy equation for other material models.
    std::cout << "Simulation terminated." << '\n';
    std::exit(1);
  }

  return U_density;

}
