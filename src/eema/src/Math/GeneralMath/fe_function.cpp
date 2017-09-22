#include"functions.h"
#include<cmath>
using namespace Eigen;

double fe_function(double a, std::string b, double time){
  double result;

  if(b=="RAMP"){
    result = a*(time/t_end);
  }

  if(b=="SIN"){
    result = a*sin((2*M_PI/0.1)*time);
  }

  if(b=="COS"){
    result = a*cos(time);
  }

  if(b=="STEP"){
    result = a;
  }

  return result;
}

double fe_function_derivative(double a, std::string b, double time){
  double result;

  if(b=="RAMP"){
    result = (a/t_end);
  }

  if(b=="SIN"){
    result = a*(2*M_PI/0.1)*cos((2*M_PI/0.1)*time);
  }

  if(b=="COS"){
    result = a*(-1.0)*sin(time);
  }

  if(b=="STEP"){
    result = 0;
  }

  return result;
}

double fe_function_d_derivative(double a, std::string b, double time){
  double result;

  if(b=="RAMP"){
    result = 0;
  }

  if(b=="SIN"){
    result = a*-1.0*pow((2*M_PI/0.1),2)*sin((2*M_PI/0.1)*time);
  }

  if(b=="COS"){
    result = a*-1.0*cos(time);
  }

  if(b=="STEP"){
    result = 0;
  }

  return result;
}
