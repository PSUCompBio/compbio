#include"functions.h"

using namespace Eigen;

MatrixXd fe_strDispMatrix(int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz) {

	MatrixXd disp_mat(6, edof);
	disp_mat = MatrixXd::Zero(6, edof);
	int i1, i2, i3;

	for (int i = 0; i < nnel; i++) {

		i1 = 0 + i * 3;
		i2 = i1 + 1;
		i3 = i2 + 1;
		disp_mat(0, i1) = dndx(i);
		disp_mat(1, i2) = dndy(i);
		disp_mat(2, i3) = dndz(i);
		disp_mat(3, i1) = dndy(i);
		disp_mat(3, i2) = dndx(i);
		disp_mat(4, i2) = dndz(i);
		disp_mat(4, i3) = dndy(i);
		disp_mat(5, i1) = dndz(i);
		disp_mat(5, i3) = dndx(i);

	}
	//std::cout<<"Exiting Strain Displacement Generator Function\n";
	return disp_mat;
}

void fe_strDispMatrix_totalLagrangian_pbr(MatrixXd& disp_mat_modified, int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u) {

	// Refer to Belytschko if you have questions about the formation of this matrix

	fe_calDefGrad_pbr(F_sdm, dndx, dndy, dndz, u);
	FT_sdm = F_sdm.transpose();

	for (i_sdm = 0; i_sdm < nnel; i_sdm++) {

		disp_mat_modified(0, (i_sdm * 3)) = dndx(i_sdm) * FT_sdm(0, 0);
		disp_mat_modified(0, ((i_sdm * 3) + 1)) = dndx(i_sdm) * FT_sdm(0, 1);
		disp_mat_modified(0, (((i_sdm * 3) + 1) + 1)) = dndx(i_sdm) * FT_sdm(0, 2);

		disp_mat_modified(1, (i_sdm * 3)) = dndy(i_sdm) * FT_sdm(1, 0);
		disp_mat_modified(1, ((i_sdm * 3) + 1)) = dndy(i_sdm) * FT_sdm(1, 1);
		disp_mat_modified(1, (((i_sdm * 3) + 1) + 1)) = dndy(i_sdm) * FT_sdm(1, 2);

		disp_mat_modified(2, (i_sdm * 3)) = dndz(i_sdm) * FT_sdm(2, 0);
		disp_mat_modified(2, ((i_sdm * 3) + 1)) = dndz(i_sdm) * FT_sdm(2, 1);
		disp_mat_modified(2, (((i_sdm * 3) + 1) + 1)) = dndz(i_sdm) * FT_sdm(2, 2);


		disp_mat_modified(3, (i_sdm * 3)) = (dndy(i_sdm) * FT_sdm(0, 0)) + (dndx(i_sdm) * FT_sdm(1, 0));
		disp_mat_modified(3, ((i_sdm * 3) + 1)) = (dndy(i_sdm) * FT_sdm(0, 1)) + (dndx(i_sdm) * FT_sdm(1, 1));
		disp_mat_modified(3, (((i_sdm * 3) + 1) + 1)) = (dndy(i_sdm) * FT_sdm(0, 2)) + (dndx(i_sdm) * FT_sdm(1, 2));

		disp_mat_modified(4, (i_sdm * 3)) = (dndz(i_sdm) * FT_sdm(1, 0)) + (dndy(i_sdm) * FT_sdm(2, 0));
		disp_mat_modified(4, ((i_sdm * 3) + 1)) = (dndz(i_sdm) * FT_sdm(1, 1)) + (dndy(i_sdm) * FT_sdm(2, 1));
		disp_mat_modified(4, (((i_sdm * 3) + 1) + 1)) = (dndz(i_sdm) * FT_sdm(1, 2)) + (dndy(i_sdm) * FT_sdm(2, 2));

		disp_mat_modified(5, (i_sdm * 3)) = (dndz(i_sdm) * FT_sdm(0, 0)) + (dndx(i_sdm) * FT_sdm(2, 0));
		disp_mat_modified(5, ((i_sdm * 3) + 1)) = (dndz(i_sdm) * FT_sdm(0, 1)) + (dndx(i_sdm) * FT_sdm(2, 1));
		disp_mat_modified(5, (((i_sdm * 3) + 1) + 1)) = (dndz(i_sdm) * FT_sdm(0, 2)) + (dndx(i_sdm) * FT_sdm(2, 2));
	}
}


void fe_strDispMatrix_totalLagrangian_fiber_pbr(MatrixXd& disp_mat_modified, int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, MatrixXd& F_fiber_ref) {

	// Refer to Belytschko if you have questions about the formation of this matrix

	MatrixXd FT = MatrixXd::Zero(3, 3);
	FT = F_fiber_ref.transpose();

	int i1, i2, i3;

	for (int i = 0; i < nnel; i++) {

		i1 = 0 + i * 3;
		i2 = i1 + 1;
		i3 = i2 + 1;

		disp_mat_modified(0, i1) = dndx(i) * FT(0, 0);
		disp_mat_modified(0, i2) = dndx(i) * FT(0, 1);
		disp_mat_modified(0, i3) = dndx(i) * FT(0, 2);

		disp_mat_modified(1, i1) = dndy(i) * FT(1, 0);
		disp_mat_modified(1, i2) = dndy(i) * FT(1, 1);
		disp_mat_modified(1, i3) = dndy(i) * FT(1, 2);

		disp_mat_modified(2, i1) = dndz(i) * FT(2, 0);
		disp_mat_modified(2, i2) = dndz(i) * FT(2, 1);
		disp_mat_modified(2, i3) = dndz(i) * FT(2, 2);


		disp_mat_modified(3, i1) = (dndy(i) * FT(0, 0)) + (dndx(i) * FT(1, 0));
		disp_mat_modified(3, i2) = (dndy(i) * FT(0, 1)) + (dndx(i) * FT(1, 1));
		disp_mat_modified(3, i3) = (dndy(i) * FT(0, 2)) + (dndx(i) * FT(1, 2));

		disp_mat_modified(4, i1) = (dndz(i) * FT(1, 0)) + (dndy(i) * FT(2, 0));
		disp_mat_modified(4, i2) = (dndz(i) * FT(1, 1)) + (dndy(i) * FT(2, 1));
		disp_mat_modified(4, i3) = (dndz(i) * FT(1, 2)) + (dndy(i) * FT(2, 2));

		disp_mat_modified(5, i1) = (dndz(i) * FT(0, 0)) + (dndx(i) * FT(2, 0));
		disp_mat_modified(5, i2) = (dndz(i) * FT(0, 1)) + (dndx(i) * FT(2, 1));
		disp_mat_modified(5, i3) = (dndz(i) * FT(0, 2)) + (dndx(i) * FT(2, 2));

	}

}
