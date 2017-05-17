#include "functions.h"
using namespace Eigen;

void fe_find_common_face(VectorXi& common_face, bool com_face_value, int id_1, int id_2, MatrixXd& nodes_host, MatrixXi& elements_host) {

	common_face = VectorXi::Zero(4);
	com_face_value = 0;

	VectorXi element_list = elements_host.col(0);
	MatrixXi face_mat_1 = MatrixXi::Zero(6, 4);
	MatrixXi face_mat_2 = MatrixXi::Zero(6, 4);

	int tmp1 = fe_find(element_list, id_1);
	int tmp2 = fe_find(element_list, id_2);

	fe_form_face_matrix(face_mat_1, tmp1, elements_host);
	fe_form_face_matrix(face_mat_2, tmp2, elements_host);

	MatrixXi tmp_face_mat = face_mat_1;

	for (int k = 0; k < face_mat_1.rows(); k++) {
		face_mat_1.row(k) = fe_sort_ascendingVectorXi(face_mat_1.row(k));
		face_mat_2.row(k) = fe_sort_ascendingVectorXi(face_mat_2.row(k));
	}

	for (int i = 0; i < face_mat_1.rows(); i++) {
		for (int j = 0; j < face_mat_2.rows(); j++) {
			if (face_mat_1.row(i) == face_mat_2.row(j)) {
				common_face = tmp_face_mat.row(i);
				com_face_value = 1;
				break;
			}
		}
		if (com_face_value == 1) {
			break;
		}
	}
}