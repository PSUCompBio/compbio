#include"functions.h"

using namespace Eigen;

double fe_get_mats(int matl_code, int obj_interest, std::string type) {

	double mat_param = 0;
	int i = 0;

	if (type == "mechanical") {
		i = matMap[matl_code];
		mat_param = mat[i].getMechMatParam(obj_interest);
	}

	if (type == "electrical") {
		i = matMap[matl_code];
		mat_param = mat[i].getElecMatParam(obj_interest);
	}

	return mat_param;
}

// obj_interest table:
// 0 - density
// 1 - young's or bulk modulus
// 2 - poisson ratio
// 3 - Mooney Rivlin (C01), Ogden(),
// 4 - Mooney Rivlin (C02), Ogden(),
