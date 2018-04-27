#include "functions.h"
#include "BC.h"

using namespace Eigen;

void BC::readBC(std::string type_user, double input_amp_user, int num_dof_user, VectorXi dof_user, std::string time_curve_user, VectorXd time_data_user, VectorXd x_data_user, VectorXd y_data_user, VectorXd z_data_user, int n_steps_user){
  type = type_user;
  amplitude = input_amp_user;
  num_dof = num_dof_user;
  dof = dof_user;
  time_curve = time_curve_user;
time_data = time_data_user;
x_data = x_data_user;
y_data = y_data_user;  
z_data = z_data_user;
n_steps = n_steps_user;
}

std::string BC::getType(){
  return type;
}

double BC::getAmplitude(){
  return amplitude;
}

int BC::getNumDOF(){
  return num_dof;
}

VectorXi BC::getDOF(){
  return dof;
}

std::string BC::getTimeBehavior(){
  return time_curve;
}

VectorXd BC::getTime_data(){
  return time_data;
}

VectorXd BC::getX_data(){
  return x_data;
}

VectorXd BC::getY_data(){
  return y_data;
}

VectorXd BC::getZ_data(){
  return z_data;
}

int BC::getN_steps(){
  return n_steps;
}
void BC::printInfo(){
  std::cout << "Type of Boundary Condition: " << type << "\n";
  std::cout << "Amplitude: " << amplitude << "\n";
  std::cout << "Number of DOF's of interest: " << num_dof << '\n';
  std::cout << "DOF details: \n" << dof << '\n';
  std::cout << "Time dependent behavior: " << time_curve << '\n';
}
