#ifndef HEADERS_BC_H_
#define HEADERS_BC_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class BC{
private:
    std::string type;
    double amplitude;
    int num_dof, n_steps;
    VectorXi dof;
    VectorXd time_data, x_data, y_data, z_data;
    std::string time_curve;
public:
	void readBC(std::string type_user, double input_amp_user, int num_dof_user, VectorXi dof_user, std::string time_curve_user, VectorXd time_data_user, VectorXd x_data_user, VectorXd y_data_user, VectorXd z_data_user, int n_steps_user);

    std::string getType();
    double getAmplitude();
    int getNumDOF();
    VectorXi getDOF();
    std::string getTimeBehavior();
    VectorXd getTime_data();
    VectorXd getX_data();
    VectorXd getY_data();
    VectorXd getZ_data();
    int getN_steps();
    void printInfo();
};

#endif /* HEADERS_BC_H_ */
