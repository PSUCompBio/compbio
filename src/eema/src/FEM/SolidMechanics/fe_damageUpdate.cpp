#include "functions.h"

using namespace Eigen;

void fe_staticDamageUpdate_pbr(VectorXd& d_static, int element_id, VectorXd& lambda_max)
{

  // Function actions:
  // 1. It calculates the static damage associated with a specific stretch value. The damage evolution law is linear.
  // 2. For each fiber element, we use the maximum previous stretch for that specific fiber element.
  // 3. For each host element, we use the maximum previous average fiber stretch for the fibers associated with that specific host element.

  double lambda_initiation = 1.91; // stretch level at which damage initiates
  double lambda_separation = 2.66; // stretch level at which final separation occurs

  if (lambda_max(element_id) > lambda_initiation) {
    if (lambda_max(element_id) <= lambda_separation) {
      d_static(element_id) = (1 / (lambda_separation - lambda_initiation) ) * lambda_max(element_id) + ( -lambda_initiation / (lambda_separation - lambda_initiation) );
    }
    if (lambda_max(element_id) > lambda_separation) {
      d_static(element_id) = 1;
    }
  }

}

void fe_fatigueDamageUpdate_pbr(int opt, VectorXd& d_fatigue, int fib, double lambda, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t)
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

  if (std::abs(lambda - 1) > tol) {
    if (lambda < lambda_min_cycle(fib)) {
      lambda_min_cycle(fib) = lambda;
    }
    if (lambda > lambda_max_cycle(fib)) {
      lambda_max_cycle(fib) = lambda;
    }
  }

  if (std::abs(lambda - 1) <= tol) {
    if (lambda_min_cycle(fib) != 1 || lambda_max_cycle(fib) != 1) {
      double equivalent_strain_lambda_min_cycle = fe_calEquivalentStrain(opt, lambda_min_cycle(fib));
      double equivalent_strain_lambda_max_cycle = fe_calEquivalentStrain(opt, lambda_max_cycle(fib));
      double equivalent_strain_peak = 0;
      if (equivalent_strain_lambda_min_cycle >= equivalent_strain_lambda_max_cycle) {
        equivalent_strain_peak = equivalent_strain_lambda_min_cycle;
      }
      if (equivalent_strain_lambda_min_cycle < equivalent_strain_lambda_max_cycle) {
        equivalent_strain_peak = equivalent_strain_lambda_max_cycle;
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
      n_load_cycle_full(fib) = n_load_cycle_full(fib) + 1;
      lambda_min_cycle(fib) = 1;
      lambda_max_cycle(fib) = 1;
    }
  }

  if (t == t_end) {
    if (std::abs(lambda - 1) > tol) {
      if (lambda_min_cycle(fib) != 1 || lambda_max_cycle(fib) != 1) {
        double equivalent_strain_lambda_min_cycle = fe_calEquivalentStrain(opt, lambda_min_cycle(fib));
        double equivalent_strain_lambda_max_cycle = fe_calEquivalentStrain(opt, lambda_max_cycle(fib));
        double equivalent_strain_peak = 0;
        if (equivalent_strain_lambda_min_cycle >= equivalent_strain_lambda_max_cycle) {
          equivalent_strain_peak = equivalent_strain_lambda_min_cycle;
        }
        if (equivalent_strain_lambda_min_cycle < equivalent_strain_lambda_max_cycle) {
          equivalent_strain_peak = equivalent_strain_lambda_max_cycle;
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
        n_load_cycle_partial(fib) = n_load_cycle_partial(fib) + 1;
        lambda_min_cycle(fib) = 1;
        lambda_max_cycle(fib) = 1;
      }
    }
  }

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




