#include "functions.h"
using namespace Eigen;

void fe_form_face_matrix(MatrixXi& face_mat, int tmp, MatrixXi& elements) {

	MatrixXi faces(6, 4);
	faces << 0, 1, 2, 3,
	      5, 1, 0, 4,
	      6, 5, 4, 7,
	      3, 2, 6, 7,
	      7, 4, 0, 3,
	      6, 2, 1, 5;

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < faces.cols(); j++) {
			face_mat(i, j) = elements(tmp, (faces(i, j) + 2));
		}
	}

}