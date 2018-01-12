#include "functions.h"

using namespace Eigen;

void fe_concatenate_vector2matrix(MatrixXd& A, VectorXd& B, int opt) {

	//opt=1: row-wise
	//opt=2: col-wise

	if (opt == 1) {
		if (A.rows() != 0 && A.cols() != 0) {
			MatrixXd tmp = MatrixXd::Zero(A.rows(), A.cols());
			tmp = A;
			A = MatrixXd::Zero((tmp.rows() + 1), tmp.cols());
			for (int i = 0; i < A.rows(); i++) {
				if (i != A.rows() - 1) {
					A.row(i) = tmp.row(i);
				}
				else
				{
					A.row(i) = B;
				}
			}
		}
		else {
			A = MatrixXd::Zero(1, B.size());
			A.row(0) = B;
		}
	}

	if (opt == 2) {
		MatrixXd tmp = MatrixXd::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXd::Zero(tmp.rows(), (tmp.cols() + 1));
		A << tmp, B;
	}
}

void fe_concatenate_vector2matrix(MatrixXi& A, VectorXi& B, int opt) {

	//opt=1: row-wise
	//opt=2: col-wise

	if (opt == 1) {

		if (A.rows() != 0 && A.cols() != 0) {

			MatrixXi tmp = MatrixXi::Zero(A.rows(), A.cols());
			tmp = A;
			A = MatrixXi::Zero((tmp.rows() + 1), tmp.cols());
			for (int i = 0; i < A.rows(); i++) {
				if (i != A.rows() - 1) {
					A.row(i) = tmp.row(i);
				}
				else
				{
					A.row(i) = B;
				}
			}
		}
		else {
			A = MatrixXi::Zero(1, B.size());
			A.row(0) = B;
		}
	}

	if (opt == 2) {
		MatrixXi tmp = MatrixXi::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXi::Zero(tmp.rows(), (tmp.cols() + 1));
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

void fe_insert_vector2matrix(MatrixXi&A, VectorXi& B, int num, int opt) {

	if (opt == 1) {
		MatrixXi tmp = MatrixXi::Zero(A.rows(), A.cols());
		tmp = A;
		A = MatrixXi::Zero((tmp.rows() + 1), tmp.cols());
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

void fe_concatenate_double2vector(VectorXd& A, double b) {

	if (A.size() != 0) {
		VectorXd tmp = VectorXd::Zero(A.size());
		tmp = A;
		A = VectorXd::Zero(tmp.size() + 1);
		for (int i = 0; i < A.size(); i++) {
			if (i != A.size() - 1) {
				A(i) = tmp(i);
			}
			else {
				A(i) = b;
			}
		}
	}
	else {
		A = VectorXd::Zero(1);
		A(0) = b;
	}

}

void fe_insert_int2vector(VectorXi& A, int b) {

	if (A.size() != 0) {
		VectorXi tmp = VectorXi::Zero(A.size());
		tmp = A;
		A = VectorXi::Zero(tmp.size() + 1);
		for (int i = 0; i < A.size(); i++) {
			if (i != A.size() - 1) {
				A(i) = tmp(i);
			}
			else {
				A(i) = b;
			}
		}
	}
	else {
		A = VectorXi::Zero(1);
		A(0) = b;
	}

}

void fe_insert_int2vector(VectorXi& A, int num, int b) {
	if (A.size() != 0) {

		int size_old = A.size();
		int size_new = 0;

		if (num < size_old) {
			size_new = size_old;
		}
		if (num > (size_old - 1) ) {
			size_new = num + 1;
		}

		VectorXi tmp = VectorXi::Zero(size_old);
		tmp = A;

		A = VectorXi::Zero(size_new);

		int counter = 0;

		for (int i = 0; i < size_new; i++) {
			if (i != num) {
				if (counter > (size_old - 1) ) {
					A(i) = 0;
				}
				if (counter < size_old) {
					A(i) = tmp(i);
					counter = counter + 1;
				}
			}
			else {
				A(i) = b;
			}
		}
	}
	else {
		A = VectorXi::Zero(1);
		A(0) = b;
	}

	// if (A.size() != 0) {
	// 	VectorXi tmp = VectorXi::Zero(A.size());
	// 	tmp = A;
	// 	A = VectorXi::Zero(tmp.size() + 1);
	// 	int counter = 0;
	// 	for (int i = 0; i < A.size(); i++) {
	// 		if (i != num) {
	// 			A(i) = tmp(counter);
	// 			counter = counter + 1;
	// 		}
	// 		else {
	// 			A(i) = b;
	// 		}
	// 	}
	// }
	// else {
	// 	A = VectorXi::Zero(1);
	// 	A(0) = b;
	// }

}
