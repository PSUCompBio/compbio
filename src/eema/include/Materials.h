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

class Materials {
private:
	int mat_id;
	std::string mechanical_mat_model;
	std::string electrical_mat_model;
	VectorXd mechanical_mat_properties;
	VectorXd electrical_mat_properties;

	double data;
	friend double fe_get_mats(int matl_code, int obj_interest, std::string type);
	friend std::string fe_get_model(int matl_code, std::string type);
public:
	void readMatId(int a);
	void readMats(std::string b, VectorXd c, std::string type);
	void printInfo(std::string type);
	double getMechMatParam(int obj_interest);
	double getElecMatParam(int obj_interest);
	int getMatID();
};

#endif /* HEADERS_MATERIALS_H_ */
