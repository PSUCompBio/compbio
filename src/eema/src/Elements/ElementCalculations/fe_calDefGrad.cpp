
/** \brief This function calculates the deformation gradient based on the displacements */

#include"functions.h"

using namespace Eigen;

void fe_calDefGrad_pbr(MatrixXd& F, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u) {

	i_lbv = 0;

	H_DefGrad(0, 0) = (dndx(i_lbv) * u(i_lbv * 3));
	H_DefGrad(0, 1) = (dndy(i_lbv) * u(i_lbv * 3));
	H_DefGrad(0, 2) = (dndz(i_lbv) * u(i_lbv * 3));

	H_DefGrad(1, 0) = (dndx(i_lbv) * u((i_lbv * 3) + 1));
	H_DefGrad(1, 1) = (dndy(i_lbv) * u((i_lbv * 3) + 1));
	H_DefGrad(1, 2) = (dndz(i_lbv) * u((i_lbv * 3) + 1));

	H_DefGrad(2, 0) = (dndx(i_lbv) * u((i_lbv * 3) + 2));
	H_DefGrad(2, 1) = (dndy(i_lbv) * u((i_lbv * 3) + 2));
	H_DefGrad(2, 2) = (dndz(i_lbv) * u((i_lbv * 3) + 2));

	for (i_lbv = 1; i_lbv < dndx.size(); i_lbv++) {

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

void fe_calDefGrad_pbr_array(MatrixXd& F, int i, int x, int y, int z, VectorXd& u) {

	i_lbv = 0;

	H_DefGrad(0, 0) = (dndx_store[i][x][y][z][i_lbv] * u(i_lbv * 3));
	H_DefGrad(0, 1) = (dndy_store[i][x][y][z][i_lbv] * u(i_lbv * 3));
	H_DefGrad(0, 2) = (dndz_store[i][x][y][z][i_lbv] * u(i_lbv * 3));

	H_DefGrad(1, 0) = (dndx_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));
	H_DefGrad(1, 1) = (dndy_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));
	H_DefGrad(1, 2) = (dndz_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));

	H_DefGrad(2, 0) = (dndx_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));
	H_DefGrad(2, 1) = (dndy_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));
	H_DefGrad(2, 2) = (dndz_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));

	for (i_lbv = 1; i_lbv < nnel_normal; i_lbv++) {

		// First row
		H_DefGrad(0, 0) = H_DefGrad(0, 0) + (dndx_store[i][x][y][z][i_lbv] * u(i_lbv * 3));
		H_DefGrad(0, 1) = H_DefGrad(0, 1) + (dndy_store[i][x][y][z][i_lbv] * u(i_lbv * 3));
		H_DefGrad(0, 2) = H_DefGrad(0, 2) + (dndz_store[i][x][y][z][i_lbv] * u(i_lbv * 3));

		//Second row
		H_DefGrad(1, 0) = H_DefGrad(1, 0) + (dndx_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));
		H_DefGrad(1, 1) = H_DefGrad(1, 1) + (dndy_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));
		H_DefGrad(1, 2) = H_DefGrad(1, 2) + (dndz_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 1));

		//Third row
		H_DefGrad(2, 0) = H_DefGrad(2, 0) + (dndx_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));
		H_DefGrad(2, 1) = H_DefGrad(2, 1) + (dndy_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));
		H_DefGrad(2, 2) = H_DefGrad(2, 2) + (dndz_store[i][x][y][z][i_lbv] * u((i_lbv * 3) + 2));
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
