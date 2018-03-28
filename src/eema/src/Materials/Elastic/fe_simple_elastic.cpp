#include "functions.h"

using namespace Eigen;

void fe_simple_elastic_pbr(VectorXd& sigma_local, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt) {

	sigma_local = VectorXd::Zero(6);

	MatrixXd matl_mat = MatrixXd::Zero(6, 6);
	double E = fe_get_mats(opt, 1, "mechanical");
	double nu = fe_get_mats(opt, 2, "mechanical");
	fe_calculate_matlmat_pbr(matl_mat, 4, E, nu);
	sigma_local = matl_mat * disp_mat * u;
}

void fe_calculate_matlmat_pbr(MatrixXd& matl_mat, int n, double E, double nu) {

	double cons = (E / ((1 + nu) * (1 - (2 * nu))));

	if (n == 4) { // n==4 3 dimensional analysis
		// int dim = (n*(n+1)/2);

		// 1st Row
		matl_mat(0, 0) = (1 - nu);
		matl_mat(0, 1) = nu;
		matl_mat(0, 2) = nu;
		matl_mat(0, 3) = 0;
		matl_mat(0, 4) = 0;
		matl_mat(0, 5) = 0;

		//2nd Row
		matl_mat(1, 0) = (nu);
		matl_mat(1, 1) = (1 - nu);
		matl_mat(1, 2) = nu;
		matl_mat(1, 3) = 0;
		matl_mat(1, 4) = 0;
		matl_mat(1, 5) = 0;

		//3rd Row
		matl_mat(2, 0) = (nu);
		matl_mat(2, 1) = nu;
		matl_mat(2, 2) = (1 - nu);
		matl_mat(2, 3) = 0;
		matl_mat(2, 4) = 0;
		matl_mat(2, 5) = 0;

		//4th Row
		matl_mat(3, 0) = 0;
		matl_mat(3, 1) = 0;
		matl_mat(3, 2) = 0;
		matl_mat(3, 3) = ((1 - (2 * nu)) / 2);
		matl_mat(3, 4) = 0;
		matl_mat(3, 5) = 0;

		//5th Row
		matl_mat(4, 0) = 0;
		matl_mat(4, 1) = 0;
		matl_mat(4, 2) = 0;
		matl_mat(4, 3) = 0;
		matl_mat(4, 4) = ((1 - (2 * nu)) / 2);
		matl_mat(4, 5) = 0;

		//6th Row
		matl_mat(5, 0) = 0;
		matl_mat(5, 1) = 0;
		matl_mat(5, 2) = 0;
		matl_mat(5, 3) = 0;
		matl_mat(5, 4) = 0;
		matl_mat(5, 5) = ((1 - (2 * nu)) / 2);
	}

	matl_mat = matl_mat * cons;

}
