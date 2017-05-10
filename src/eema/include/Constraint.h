#ifndef HEADERS_CONSTRAINT_H_
#define HEADERS_CONSTRAINT_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Constraint {

private:
	std::string constraint_name;
	int constraint_id;

	// Embedded Element Constraint
	std::string master_name;
	std::string slave_name;
	bool address_volume_redundancy;

public:
	void readConstraints(std::string name, int id, std::string master, std::string slave, bool address_VR);
	void printInfo();

	std::string getName();

	// Embedded Element Constraint Get Functions
	std::string get_EmbedMaster();
	std::string get_EmbedSlave();
	bool get_EmbedAddressVR();
};

#endif /* HEADERS_CONSTRAINT_H_ */