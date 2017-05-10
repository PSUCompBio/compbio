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
    int num_dof;
    VectorXi dof;
    std::string time_curve;
public:
	void readBC(std::string type_user, double input_amp_user, int num_dof_user, VectorXi dof_user, std::string time_curve_user);

    std::string getType();
    double getAmplitude();
    int getNumDOF();
    VectorXi getDOF();
    std::string getTimeBehavior();

    void printInfo();
};

#endif /* HEADERS_BC_H_ */
