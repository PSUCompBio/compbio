#include "functions.h"

using namespace Eigen;

void fe_concatenate_vector2matrix(MatrixXd& A, VectorXd& B, int opt) {

	//opt=1: row-wise
	//opt=2: col-wise

	if (opt == 1) {
		MatrixXd tmp = MatrixXd::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXd::Zero((tmp.rows() + 1), tmp.cols());
		A << tmp, B;
	}

	if (opt == 2) {
		MatrixXd tmp = MatrixXd::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXd::Zero(tmp.rows(), (tmp.cols() + 1));
		A << tmp, B;
	}
}

void fe_insert_vector2matrix(MatrixXd& A, VectorXd& B, int num, int opt) {

	//opt-1: row wise
	//opt-2: col wise

	if (opt == 1) {
		MatrixXd tmp = MatrixXd::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXd::Zero((tmp.rows() + 1), tmp.cols());
		int counter = 0;
		for (int i = 0; i < A.rows(); i++) {
			if (i == num) {
				A.row(i) = B;
			}
			else {
				A.row(i) = tmp.row(counter);
				counter = counter + 1;
			}
		}
	}

}