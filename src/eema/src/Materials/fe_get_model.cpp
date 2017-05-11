#include "functions.h"

using namespace Eigen;

std::string fe_get_model(int matl_code, std::string type) {
	std::string model;

	if (type == "mechanical") {
		int num_pointer = 0;
		for (int i = 0; i < material_types; i++) {
			if (matl_code == mat[i].mat_id) {
				num_pointer = i;
				break;
			}
		}
		model = mat[num_pointer].mechanical_mat_model;
	}

	if (type == "electrical") {
		int num_pointer = 0;
		for (int i = 0; i < material_types; i++) {
			if (matl_code == mat[i].mat_id) {
				num_pointer = i;
				break;
			}
		}
		model = mat[num_pointer].mechanical_mat_model;
	}

	return model;
}
