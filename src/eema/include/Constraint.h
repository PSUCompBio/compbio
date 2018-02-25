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

	/** Vector Storing the host elements ids for each embedded element */
	VectorXi embed_map;
	VectorXi* embed_map_pointer;

	bool address_volume_redundancy;
	bool include_damage;
	bool import_damage;

	double time_healing;

public:
	void readConstraints(std::string name, int id, std::string master, std::string slave, bool address_VR, bool include_d, bool import_d, double t_healing);
	void preprocess();
	void printInfo();

	std::string getName();

	// Embedded Element Constraint Get Functions
	std::string get_EmbedMaster();
	std::string get_EmbedSlave();
	bool get_EmbedAddressVR();
	bool get_EmbedIncludeDamage();
	bool get_EmbedImportDamage();
	VectorXi* get_EmbedMapPointer();
	double get_HealingTime();
};

#endif /* HEADERS_CONSTRAINT_H_ */