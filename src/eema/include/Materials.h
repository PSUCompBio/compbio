/*
 * Materials.h
 *
 *  Created on: Dec 19, 2016
 *      Author: vsg111
 */

#ifndef HEADERS_MATERIALS_H_
#define HEADERS_MATERIALS_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Materials{
private:
	int mat_id;
	std::string mat_model;
	VectorXd mat_properties;
	//VectorXd mat_properties = VectorXd::Zero(100);
	double data;
	friend double fe_get_mats(int matl_code, int obj_interest);
	friend std::string fe_get_model(int matl_code);
public:
	void readMats(int a, std::string b, VectorXd c);
	void printMats();
};

#endif /* HEADERS_MATERIALS_H_ */
