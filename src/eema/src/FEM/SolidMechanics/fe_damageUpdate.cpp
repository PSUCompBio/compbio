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

void fe_deltaDamageUpdate_pbr(VectorXd& delta_d, int fib, double lambda, VectorXd& lambda_min, VectorXd& lambda_max)
{

  double tol = pow(10, -4);  // tolerance used to determine if lambda = 1
  double delta_d_add = 0;

  double lambda_LL = 0.5;    // lower limit for stretch, d = 1 for lower values
  double lambda_UL = 1.5;    // upper limit for stretch, d = 1 for higher values
  double delta_d_CUL = 0.1;  // upper limit for delta_d during compressive stretch, occurs at lambda_LL
  double delta_d_TUL = 0.1;  // upper limit for delta_d during tensile stretch, occurs at lambda_UL

  if (lambda < lambda_min(fib)) {
    lambda_min(fib) = lambda;
  }

  if (lambda > lambda_max(fib)) {
    lambda_max(fib) = lambda;
  }

  if (std::abs(lambda - 1) < tol) {
    if (lambda_min(fib) < 1) {
      delta_d_add = (-delta_d_CUL/(1 - lambda_LL))*lambda_min(fib) + (delta_d_CUL/(1 - lambda_LL));
    }
    if (lambda_max(fib) > 1) {
      delta_d_add = (delta_d_TUL/(lambda_UL - 1))*lambda_max(fib) - (delta_d_TUL/(lambda_UL - 1));
    }
    lambda_min(fib) = 1;
    lambda_max(fib) = 1;
  }

  delta_d(fib) = delta_d(fib) + delta_d_add;

}
