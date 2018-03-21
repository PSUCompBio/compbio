
/** \brief This function calculates the deformation gradient based on the displacements */

#include"functions.h"

using namespace Eigen;

MatrixXd fe_calDefGrad(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u) {

	MatrixXd F = MatrixXd::Zero(3, 3); // Deformation Gradient
	MatrixXd H = MatrixXd::Zero(3, 3); // F  = I+H and H = du/dX
	MatrixXd I = MatrixXd::Identity(3, 3); // Identity Matrix
	int i1, i2, i3;

	for (int i = 0; i < dndx.size(); i++) {
		i1 = i * 3;
		i2 = i1 + 1;
		i3 = i2 + 1;

		// First row
		H(0, 0) = H(0, 0) + (dndx(i) * u(i1));
		H(0, 1) = H(0, 1) + (dndy(i) * u(i1));
		H(0, 2) = H(0, 2) + (dndz(i) * u(i1));

		//Second row
		H(1, 0) = H(1, 0) + (dndx(i) * u(i2));
		H(1, 1) = H(1, 1) + (dndy(i) * u(i2));
		H(1, 2) = H(1, 2) + (dndz(i) * u(i2));

		//Third row
		H(2, 0) = H(2, 0) + (dndx(i) * u(i3));
		H(2, 1) = H(2, 1) + (dndy(i) * u(i3));
		H(2, 2) = H(2, 2) + (dndz(i) * u(i3));
	}

	F = I + H;

	return F;
}

/*************************************/

void fe_calDefGrad_pbr(MatrixXd& F, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u) {

	H_DefGrad = MatrixXd::Zero(3, 3);
	for (i_lbv = 0; i_lbv < dndx.size(); i_lbv++) {

		// First row
		H_DefGrad(0, 0) = H_DefGrad(0, 0) + (dndx(i_lbv) * u(i_lbv * 3));
		H_DefGrad(0, 1) = H_DefGrad(0, 1) + (dndy(i_lbv) * u(i_lbv * 3));
		H_DefGrad(0, 2) = H_DefGrad(0, 2) + (dndz(i_lbv) * u(i_lbv * 3));

		//Second row
		H_DefGrad(1, 0) = H_DefGrad(1, 0) + (dndx(i_lbv) * u((i_lbv * 3) + 1));
		H_DefGrad(1, 1) = H_DefGrad(1, 1) + (dndy(i_lbv) * u((i_lbv * 3) + 1));
		H_DefGrad(1, 2) = H_DefGrad(1, 2) + (dndz(i_lbv) * u((i_lbv * 3) + 1));

		//Third row
		H_DefGrad(2, 0) = H_DefGrad(2, 0) + (dndx(i_lbv) * u((i_lbv * 3) + 2));
		H_DefGrad(2, 1) = H_DefGrad(2, 1) + (dndy(i_lbv) * u((i_lbv * 3) + 2));
		H_DefGrad(2, 2) = H_DefGrad(2, 2) + (dndz(i_lbv) * u((i_lbv * 3) + 2));
	}

	F = I + H_DefGrad;
}

void fe_calDefGrad_fiber_pbr(MatrixXd& F_fiber_ref, MatrixXd& T_fiber, MatrixXd& T_fiber_inv, double lambda) {

	MatrixXd F_fiber_fib = MatrixXd::Zero(ndof, ndof);

	F_fiber_fib(0,0) = lambda;
	F_fiber_fib(1,1) = 1/sqrt(lambda);
	F_fiber_fib(2,2) = 1/sqrt(lambda);

	F_fiber_ref = T_fiber_inv.transpose()*F_fiber_fib*T_fiber_inv;

}
