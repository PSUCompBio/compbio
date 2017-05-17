#include "functions.h"
using namespace Eigen;

VectorXd fe_sort_ascendingVectorXd(VectorXd A) {

	VectorXd tmp = VectorXd::Zero(A.size());
	double low = A.minCoeff() - 10;

	for (int i = 0; i < A.size(); i++) {
		tmp(i) = A.maxCoeff();
		double high = A.maxCoeff();
		for (int j = 0; j < A.size(); j++) {
			if (A(j) == high) {
				A(j) = low;
			}
		}
	}

	A = tmp;
	return A;
}

VectorXi fe_sort_ascendingVectorXi(VectorXi A) {

	VectorXi tmp = VectorXi::Zero(A.size());
	int low = A.minCoeff() - 10;

	for (int i = 0; i < A.size(); i++) {
		tmp(i) = A.maxCoeff();
		int high = A.maxCoeff();
		for (int j = 0; j < A.size(); j++) {
			if (A(j) == high) {
				A(j) = low;
			}
		}
	}

	A = tmp;
	return A;
}