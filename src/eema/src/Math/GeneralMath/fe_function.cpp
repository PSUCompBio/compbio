#include"functions.h"
#include<cmath>
using namespace Eigen;

double fe_function(double a, std::string b, double time){
  double result;

  double T = 0.1;  // period for harmonic functions (sec)
  double omega = 2*M_PI/T;  // frequency for harmonic functions

  if(b=="RAMP"){
    result = a*(time/t_end);
  }

  if(b=="SIN"){
    result = a*sin(omega*time);
  }

  if(b=="COS"){
    result = a*cos(omega*time);
  }

  if(b=="STEP"){
    result = a;
  }

  return result;
}

double fe_function_derivative(double a, std::string b, double time){
  double result;

  double T = 0.1;  // period for harmonic functions (sec)
  double omega = 2*M_PI/T;  // frequency for harmonic functions

  if(b=="RAMP"){
    result = (a/t_end);
  }

  if(b=="SIN"){
    result = a*omega*cos(omega*time);
  }

  if(b=="COS"){
    result = a*omega*(-1.0)*sin(omega*time);
  }

  if(b=="STEP"){
    result = 0;
  }

  return result;
}

double fe_function_d_derivative(double a, std::string b, double time){
  double result;

  double T = 0.1;  // period for harmonic functions (sec)
  double omega = 2*M_PI/T;  // frequency for harmonic functions

  if(b=="RAMP"){
    result = 0;
  }

  if(b=="SIN"){
    result = a*pow(omega,2)*(-1.0)*sin(omega*time);
  }

  if(b=="COS"){
    result = a*pow(omega,2)*-(1.0)*cos(omega*time);
  }

  if(b=="STEP"){
    result = 0;
  }

  return result;
}
