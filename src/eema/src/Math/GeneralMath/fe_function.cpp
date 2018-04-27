#include"functions.h"
#include<cmath>
using namespace Eigen;

double fe_function(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction){

  double result;
  double T = 0.2;  // period for harmonic functions (sec)
  double omega = 2*M_PI/T;  // frequency for harmonic functions
  if(b=="RAMP"){
    result = a*(time/t_end);
//std::cout<<" ";
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

  if(b=="TABULAR"){
	//std::cout<<time<<" ";
    for(int i = 0; i<n_steps; i++){
	if(time_data[i]<=time && time<=time_data[i+1])
	{
	if( direction == 1)
	result = x_data[i] + (x_data[i+1]-x_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	if( direction == 2)
	result = y_data[i] + (y_data[i+1]-y_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	if( direction == 3)
	result = z_data[i] + (z_data[i+1]-z_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	break;}}
std::cout<<result<<" ";
  return result;
}

}

double fe_function_derivative(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction){
  double result;

  double T = 0.2;  // period for harmonic functions (sec)
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
  
  if(b=="TABULAR"){
    for(int i = 0; i<n_steps; i++)
	if(time_data[i]<=time && time<=time_data[i+1])
	{if( direction == 1)
	result = x_data[i] + (x_data[i+1]-x_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	else if( direction == 2)
	result = y_data[i] + (y_data[i+1]-y_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	else if( direction == 3)
	result = z_data[i] + (z_data[i+1]-z_data[i])*(time-time_data[i])/(time_data[i+1]-time_data[i]);
	break;}

  return result;
}
}

double fe_function_d_derivative(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction){
  double result;

  double T = 0.2;  // period for harmonic functions (sec)
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

  if(b=="TABULAR"){
    result = 0;}

  return result;
}
