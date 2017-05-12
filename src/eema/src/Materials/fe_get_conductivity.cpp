#include "functions.h"
using namespace Eigen;

void fe_get_conductivity(MatrixXd& conductivity, int mat_id) {

	std::string type = fe_get_model(mat_id, "electrical");

	if (type == "isotropic-conductive") {
		conductivity(0, 0) = fe_get_mats(mat_id, 0, "electrical");
		conductivity(1, 1) = fe_get_mats(mat_id, 0, "electrical");
		conductivity(2, 2) = fe_get_mats(mat_id, 0, "electrical");
	}

	if (type == "anisotropic-conductive") {
		int k = 0;
		for (int i = 0; i < ndof; i++) {
			for (int j = 0; j < ndof; j++) {
				conductivity(i, j) = fe_get_mats(mat_id, k, "electrical");
				k = k + 1;
			}
		}
	}

}