#include "functions.h"

using namespace Eigen;

std::string fe_get_model(int matl_code, std::string type) {
	std::string model;

	if (type == "mechanical") {
		int i = matMap[matl_code];
		model = mat[i].mechanical_mat_model;
	}

	if (type == "electrical") {
		int i = matMap[matl_code];
		model = mat[i].electrical_mat_model;
	}

	return model;
}
