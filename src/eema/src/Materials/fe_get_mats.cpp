#include"functions.h"

using namespace Eigen;

double fe_get_mats(int matl_code, int obj_interest, std::string type) {

	double mat_param = 0;

	if (type == "mechanical") {

		int num_pointer = 0;
		for (int i = 0; i < material_types; i++) {
			if (matl_code == mat[i].mat_id) {
				num_pointer = i;
				break;
			}
		}
		VectorXd m_1 = mat[num_pointer].mechanical_mat_properties;
		mat_param = m_1(obj_interest);
	}

	if (type == "electrical") {
		int num_pointer = 0;
		for (int i = 0; i < material_types; i++) {
			if (matl_code == mat[i].mat_id) {
				num_pointer = i;
				break;
			}
		}
		VectorXd m_1 = mat[num_pointer].electrical_mat_properties;
		mat_param = m_1(obj_interest);
	}

	return mat_param;
}

// obj_interest table:
// 0 - density
// 1 - young's or bulk modulus
// 2 - poisson ratio
// 3 - Mooney Rivlin (C01), Ogden(),
// 4 - Mooney Rivlin (C02), Ogden(),
