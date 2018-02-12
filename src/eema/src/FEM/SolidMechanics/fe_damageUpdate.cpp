#include "functions.h"

using namespace Eigen;

void fe_damageUpdate_pbr(VectorXd& d, int fib, double lambda)
{

  double lambda_LL = 0.5;  // lower limit for stretch, d = 1 for lower values
  double lambda_UL = 1.5;  // upper limit for stretch, d = 1 for higher values
  double tmp = 0;  // d for current value of stretch

  if (lambda == 1) {
     tmp = 0;
  }

  if (lambda < 1) {
    if (lambda < lambda_LL) {
      tmp = 1;
    }
    if (lambda > lambda_LL) {
      tmp = (-1/(1 - lambda_LL))*lambda + (1/(1 - lambda_LL));
    }
  }

  if (lambda > 1) {
    if (lambda > lambda_UL) {
      tmp = 1;
    }
    if (lambda < lambda_UL) {
      tmp = (1/(lambda_UL - 1))*lambda - (1/(lambda_UL - 1));
    }
  }

  if (tmp > d(fib)) {
    d(fib) = tmp;
  }

}

void fe_fatigueDamageUpdate_pbr(int opt, VectorXd& d_fatigue, int fib, double lambda, VectorXd& lambda_min, VectorXd& lambda_max, int& n_load_cycle_full, int& n_load_cycle_partial, double t)
{

  // Function actions:
  // 1. It keeps track of the minimum and maximum stretch experienced by each fiber during the current load cycle.
  // 2. It recognizes when the current load cycle is complete (i.e., when lambda = 1).
  // 3. It calculates the additional fatigue damage associated with the completed load cycle.
  // 4. It calculates the additional fatigue damage associated with the final incomplete load cycle, if applicable.

  double tol = pow(10, -4); // This is the tolerance used to determine the end of a load cycle (i.e., when lambda = 1).
  // We want the tolerance to be small, so that lambda will be close to 1; which marks the true end of the load cycle. However, the tolerance cannot be too small because we need lambda to fall within the tolerance at one of the discrete time steps. This is what triggers the calculation of additional fatigue damage. If we make the tolerance too small, we might jump over lambda = 1 and not recognize that the load cycle is complete.

  double d_fatigue_additional = 0;

  std::string model;
  model = fe_get_model(opt, "mechanical");

  double alpha_fatigue = 0; // material constant governing the amount of damage incurred by a single cycle
  double beta_fatigue = 0;  // material constant governing the amount of damage incurred by a single cycle
  double psi_min = 0;       // minimum equivalent strain to elicit the accumulation of fatigue damage
  double psi_max = 0;       // maximum equivalent strain limit (i.e., complete failure after one cycle)
  double n_max = 0;         // maximum number of fatigue cycles

  if (model == "ogden_hyperelastic") {
    alpha_fatigue = 1;
    beta_fatigue = 25;
    psi_min = 5;
    psi_max = 1200;
    n_max = 500;
  } else {
    std::cout << "Fiber material must be Ogden for damage to be included." << '\n'; // Later on, we can add the material thresholds for other material models.
    std::cout << "Simulation terminated." << '\n';
    std::exit(1);
  }

  // This was the original definition of fatigue damage. We are not using this anymore. However, I will leave the code in place for now. The structure might be useful later on. -JG, 2/8/2018
  // double lambda_LL = 0.5;    // lower limit for stretch, d = 1 for lower values
  // double lambda_UL = 1.5;    // upper limit for stretch, d = 1 for higher values
  // double d_fatigue_CUL = 0.1;  // upper limit for d_fatigue during compressive stretch, occurs at lambda_LL
  // double d_fatigue_TUL = 0.1;  // upper limit for d_fatigue during tensile stretch, occurs at lambda_UL

  if (std::abs(lambda - 1) > tol) {
    if (lambda < lambda_min(fib)) {
      lambda_min(fib) = lambda;
    }
    if (lambda > lambda_max(fib)) {
      lambda_max(fib) = lambda;
    }
  }

  if (std::abs(lambda - 1) <= tol) {
    if (lambda_min(fib) != 1 || lambda_max(fib) != 1) {

      double equivalent_strain_lambda_min = fe_calEquivalentStrain(opt, lambda_min(fib));
      double equivalent_strain_lambda_max = fe_calEquivalentStrain(opt, lambda_max(fib));
      double equivalent_strain_peak = 0;

      if (equivalent_strain_lambda_min >= equivalent_strain_lambda_max) {
        equivalent_strain_peak = equivalent_strain_lambda_min;
      }
      if (equivalent_strain_lambda_min < equivalent_strain_lambda_max) {
        equivalent_strain_peak = equivalent_strain_lambda_max;
      }

      if (equivalent_strain_peak < psi_min) {
        d_fatigue_additional = 0;
      }

      if (psi_min <= equivalent_strain_peak <= psi_max) {
        double tmp1 = (beta_fatigue*(n_max - 1)) / (equivalent_strain_peak - psi_min + beta_fatigue);
        double tmp2 = 1 - exp( alpha_fatigue*(1 - equivalent_strain_peak/psi_max) );
        double tmp3 = 1 - exp( alpha_fatigue*(1 - psi_min/psi_max) );
        double n_tot = tmp1*tmp2/tmp3 + 1;
        d_fatigue_additional = 1/n_tot;
      }

      if (equivalent_strain_peak > psi_max) {
        d_fatigue_additional = 1;
      }

      d_fatigue(fib) = d_fatigue(fib) + d_fatigue_additional;
      if (d_fatigue(fib) > 1) {
        d_fatigue(fib) = 1;
      }

      n_load_cycle_full = n_load_cycle_full + 1;

      lambda_min(fib) = 1;
      lambda_max(fib) = 1;
    }
  }

  if (t == t_end) {
    if (std::abs(lambda - 1) > tol) {
      if (lambda_min(fib) != 1 || lambda_max(fib) != 1) {

        double equivalent_strain_lambda_min = fe_calEquivalentStrain(opt, lambda_min(fib));
        double equivalent_strain_lambda_max = fe_calEquivalentStrain(opt, lambda_max(fib));
        double equivalent_strain_peak = 0;

        if (equivalent_strain_lambda_min >= equivalent_strain_lambda_max) {
          equivalent_strain_peak = equivalent_strain_lambda_min;
        }
        if (equivalent_strain_lambda_min < equivalent_strain_lambda_max) {
          equivalent_strain_peak = equivalent_strain_lambda_max;
        }

        if (equivalent_strain_peak < psi_min) {
          d_fatigue_additional = 0;
        }

        if (psi_min <= equivalent_strain_peak <= psi_max) {
          double tmp1 = (beta_fatigue*(n_max - 1)) / (equivalent_strain_peak - psi_min + beta_fatigue);
          double tmp2 = 1 - exp( alpha_fatigue*(1 - equivalent_strain_peak/psi_max) );
          double tmp3 = 1 - exp( alpha_fatigue*(1 - psi_min/psi_max) );
          double n_tot = tmp1*tmp2/tmp3 + 1;
          d_fatigue_additional = 1/n_tot;
        }

        if (equivalent_strain_peak > psi_max) {
          d_fatigue_additional = 1;
        }

        d_fatigue(fib) = d_fatigue(fib) + d_fatigue_additional;
        if (d_fatigue(fib) > 1) {
          d_fatigue(fib) = 1;
        }

        n_load_cycle_partial = n_load_cycle_partial + 1;

        lambda_min(fib) = 1;
        lambda_max(fib) = 1;
      }
    }
  }

  // This was the original definition of fatigue damage. We are not using this anymore. However, I will leave the code in place for now. The structure might be useful later on. -JG, 2/8/2018
  // if (lambda_min(fib) < 1) {
  //   d_fatigue_additional = (-d_fatigue_CUL/(1 - lambda_LL))*lambda_min(fib) + (d_fatigue_CUL/(1 - lambda_LL));
  // }
  // if (lambda_max(fib) > 1) {
  //   d_fatigue_additional = (d_fatigue_TUL/(lambda_UL - 1))*lambda_max(fib) - (d_fatigue_TUL/(lambda_UL - 1));
  // }

}

double fe_calEquivalentStrain(int opt, double lambda)
{

  // Here we assume the fiber is incompressible and subjected only to applied axial stretch.
  double U_density = 0;
  double equivalent_strain = 0;

  std::string model;
  model = fe_get_model(opt, "mechanical");

  // The following material properties follow the nomenclature convention used in the ABAQUS USER'S GUIDE VOLUME III: MATERIALS, section 22.5.1.
  if (model == "ogden_hyperelastic") {
    double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
    double n = fe_get_mats(opt, 3, "mechanical");
    double mu = fe_get_mats(opt, 4, "mechanical");
    double alpha = fe_get_mats(opt, 5, "mechanical");

    U_density = ( 2*mu/pow(alpha,2) )*( pow(lambda,alpha) + 2*pow(sqrt(1/lambda),alpha) - 3 ); // Strain energy per unit reference volume.
    equivalent_strain = sqrt(2*U_density);

	} else {
    std::cout << "Fiber material must be Ogden for damage to be included." << '\n'; // Later on, we can add the strain energy equation for other material models.
    std::cout << "Simulation terminated." << '\n';
    std::exit(1);
  }

  return equivalent_strain;

}




