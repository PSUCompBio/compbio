#include<iostream>

#include "functions.h"
using namespace Eigen;

VectorXd fe_scatter(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof) {

	int edof = (node_list.size()) * ndof;
	VectorXi index = fe_find_index(node_list);

	int ii = 0;

	for (int i = 0; i < edof; i++) {
		ii = index(i);
		global_vec(ii) = global_vec(ii) + local_vec(i);
	}

	return global_vec;

}

void fe_scatter_pbr(VectorXd& global_vec, VectorXd& local_vec, VectorXi node_list, int sdof) {

	// pbr - pass by reference

	int edof = (node_list.size()) * ndof;
	VectorXi index = fe_find_index(node_list);

	int ii = 0;

	for (int i = 0; i < edof; i++) {
		ii = index(i);
		global_vec(ii) = global_vec(ii) + local_vec(i);
	}

}

void fe_scatter_electricalForce(VectorXd& global_vec, VectorXd& local_vec, VectorXi node_list) {

	int edof = node_list.size();

	int ii = 0;

	for (int i = 0; i < edof; i++) {
		ii = node_list(i);
		global_vec(ii) = global_vec(ii) + local_vec(i);
	}

}
