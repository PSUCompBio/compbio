#include "functions.h"

using namespace Eigen;

void fe_ogden_hyperelastic_pbr(VectorXd& sigma_local, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt) {

	// References:

	// I. Doghri, Mechanics of deformable solids: linear, nonlinear, analytical and computational aspects. Berlin [u.a.: Springer, 2010.

	// J. C. Simo and R. L. Taylor, “Quasi-incompressible finite elasticity in principal stretches. continuum basis and numerical algorithms,”
	// Comput. Methods Appl. Mech. Eng., vol. 85, no. 3, pp. 273–310, Feb. 1991.

  double one_third  = 1.0 / 3.0;
  double two_third  = 2.0 / 3.0;
  double four_third = 4.0 / 3.0;
  double one_half   = 1.0 / 2.0;

  MatrixXd F = MatrixXd::Zero(ndof, ndof); // deformation gradient
  fe_calDefGrad_pbr(F, dndx, dndy, dndz, u);

  // double defJacobian = F.determinant(); // Jacobian - determinant of deformation gradient
  double defJacobian = fe_detMatrix_pbr(F); // Jacobian - determinant of deformation gradient

	MatrixXd b = MatrixXd::Zero(ndof, ndof);
	b = F * F.transpose();

	EigenSolver<MatrixXd> es(b);

	MatrixXd D = es.pseudoEigenvalueMatrix();
	MatrixXd V = es.pseudoEigenvectors();

	double l_1 = sqrt(D(0, 0));
	double l_2 = sqrt(D(1, 1));
	double l_3 = sqrt(D(2, 2));

	VectorXd n_1 = VectorXd::Zero(3);
	n_1 << V(0, 0), V(1, 0), V(2, 0);
	VectorXd n_2 = VectorXd::Zero(3);
	n_2 << V(0, 1), V(1, 1), V(2, 1);
	VectorXd n_3 = VectorXd::Zero(3);
	n_3 << V(0, 2), V(1, 2), V(2, 2);

	double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
	double n = fe_get_mats(opt, 3, "mechanical");
	double mu = fe_get_mats(opt, 4, "mechanical");
	double alpha = fe_get_mats(opt, 5, "mechanical");

	double dW_1 = ((one_third / l_1) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_1, alpha) - pow(l_2, alpha) - pow(l_3, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_1));
	double dW_2 = ((one_third / l_2) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_2, alpha) - pow(l_3, alpha) - pow(l_1, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_2));
	double dW_3 = ((one_third / l_3) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_3, alpha) - pow(l_1, alpha) - pow(l_2, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_3));

	MatrixXd tau = MatrixXd::Zero(ndof, ndof); // Kirchhoff stress
	MatrixXd pk_S = MatrixXd::Zero(ndof, ndof); // second Piola-Kirchhoff stress
	MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof); // Cauchy stress

	if (l_1 != l_2 && l_1 != l_3 && l_2 != l_3) {
		tau = (l_1 * dW_1 * n_1 * n_1.transpose()) + (l_2 * dW_2 * n_2 * n_2.transpose()) + (l_3 * dW_3 * n_3 * n_3.transpose());
	}

	if (l_1 == l_2 && l_3 != l_1) {
		tau = (l_1 * dW_1 * I) + ((l_3 * dW_3 - l_1 * dW_1) * n_3 * n_3.transpose());
	}

	if (l_1 == l_3 && l_2 != l_1) {
		tau = (l_1 * dW_1 * I) + ((l_2 * dW_2 - l_1 * dW_1) * n_2 * n_2.transpose());
	}

	if (l_2 == l_3 && l_1 != l_2) {
		tau = (l_2 * dW_2 * I) + ((l_1 * dW_1 - l_2 * dW_2) * n_1 * n_1.transpose());
	}

	if (l_1 == l_2 && l_1 == l_3) {
		tau = (l_1 * dW_1 * I);
	}

	MatrixXd F_inv = MatrixXd::Zero(ndof, ndof); // inverse of deformation gradient
	// MatrixXd F_inv = F.inverse();
	fe_invMatrix_pbr(F_inv, F);

	pk_S = F_inv * tau * F_inv.transpose();

	sigma_local = fe_tensor2voigt(pk_S);

	if (return_opt == 1) {
		cauchy_sigma = (1 / defJacobian) * F * pk_S * F.transpose();
		sigma_local = fe_tensor2voigt(cauchy_sigma);
	}

	// double one_half   = 1.0 / 2.0;
	//
	// sigma_local = VectorXd::Zero(ndof * 2);
	//
	// MatrixXd F = MatrixXd::Zero(ndof, ndof); // deformation gradient
	// fe_calDefGrad_pbr(F, dndx, dndy, dndz, u);
	// double defJacobian = F.determinant();
	//
	// // fe_display_matrix(F);
	//
	// MatrixXd I = MatrixXd::Identity(ndof, ndof);
	//
	// MatrixXd B = MatrixXd::Zero(ndof, ndof); // Left Cauchy-Green Deformation Tensor
	// B = F * F.transpose();
	// EigenSolver<MatrixXd> es(B);
	//
	// MatrixXd D = es.pseudoEigenvalueMatrix();
	// MatrixXd V = es.pseudoEigenvectors();
	//
	// double l_1 = D(0, 0);
	// double l_2 = D(1, 1);
	// double l_3 = D(2, 2);
	//
	// VectorXd n_1 = VectorXd::Zero(3);
	// n_1 << V(0, 0), V(1, 0), V(2, 0);
	// VectorXd n_2 = VectorXd::Zero(3);
	// n_2 << V(0, 1), V(1, 1), V(2, 1);
	// VectorXd n_3 = VectorXd::Zero(3);
	// n_3 << V(0, 2), V(1, 2), V(2, 2);
	//
	//
	// double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
	// double p = -2 * D_1 * (defJacobian - 1);
	//
	// double n = fe_get_mats(opt, 3, "mechanical");
	// int counter = 4;
	// double W_1 = 0;
	// double W_2 = 0;
	// double W_3 = 0;
	// for (double i = 0.0; i < n;) {
	// 	double mu = fe_get_mats(opt, counter, "mechanical");
	// 	double alpha = fe_get_mats(opt, counter + 1, "mechanical");
	// 	W_1 = W_1 + (mu * pow(l_1, (alpha - 1)));
	// 	W_2 = W_2 + (mu * pow(l_2, (alpha - 1)));
	// 	W_3 = W_3 + (mu * pow(l_3, (alpha - 1)));
	// 	i = i + 1.0;
	// 	counter = counter + 2;
	// }
	// MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof);
	//
	// cauchy_sigma = -p * I + (l_1 / (l_1 * l_2 * l_3)) * W_1 * (n_1 * (n_1.transpose())) + (l_2 / (l_1 * l_2 * l_3)) * W_2 * (n_2 * (n_2.transpose())) + (l_3 / (l_1 * l_2 * l_3)) * W_3 * (n_3 * (n_3.transpose()));
	//
	// MatrixXd pk_S = MatrixXd::Zero(3, 3);
	// MatrixXd F_inv = F.inverse();
	// pk_S = defJacobian * F_inv * cauchy_sigma * F_inv.transpose();
	//
	// sigma_local = fe_tensor2voigt(pk_S);
	//
	// if (return_opt == 1) {
	// 	sigma_local = fe_tensor2voigt(cauchy_sigma);
	// }

}


void fe_ogden_hyperelastic_pbr_array(VectorXd& sigma_local, int i, int x, int y, int z, VectorXd& u, int opt, int return_opt) {

	// References:

	// I. Doghri, Mechanics of deformable solids: linear, nonlinear, analytical and computational aspects. Berlin [u.a.: Springer, 2010.

	// J. C. Simo and R. L. Taylor, “Quasi-incompressible finite elasticity in principal stretches. continuum basis and numerical algorithms,”
	// Comput. Methods Appl. Mech. Eng., vol. 85, no. 3, pp. 273–310, Feb. 1991.

  double one_third  = 1.0 / 3.0;
  double two_third  = 2.0 / 3.0;
  double four_third = 4.0 / 3.0;
  double one_half   = 1.0 / 2.0;

  MatrixXd F = MatrixXd::Zero(ndof, ndof); // deformation gradient
  fe_calDefGrad_pbr_array(F, i, x, y, z, u);

  // double defJacobian = F.determinant(); // Jacobian - determinant of deformation gradient
  double defJacobian = fe_detMatrix_pbr(F); // Jacobian - determinant of deformation gradient

	MatrixXd b = MatrixXd::Zero(ndof, ndof);
	b = F * F.transpose();

	EigenSolver<MatrixXd> es(b);

	MatrixXd D = es.pseudoEigenvalueMatrix();
	MatrixXd V = es.pseudoEigenvectors();

	double l_1 = sqrt(D(0, 0));
	double l_2 = sqrt(D(1, 1));
	double l_3 = sqrt(D(2, 2));

	VectorXd n_1 = VectorXd::Zero(3);
	n_1 << V(0, 0), V(1, 0), V(2, 0);
	VectorXd n_2 = VectorXd::Zero(3);
	n_2 << V(0, 1), V(1, 1), V(2, 1);
	VectorXd n_3 = VectorXd::Zero(3);
	n_3 << V(0, 2), V(1, 2), V(2, 2);

	double D_1 = 2 / fe_get_mats(opt, 1, "mechanical");
	double n = fe_get_mats(opt, 3, "mechanical");
	double mu = fe_get_mats(opt, 4, "mechanical");
	double alpha = fe_get_mats(opt, 5, "mechanical");

	double dW_1 = ((one_third / l_1) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_1, alpha) - pow(l_2, alpha) - pow(l_3, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_1));
	double dW_2 = ((one_third / l_2) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_2, alpha) - pow(l_3, alpha) - pow(l_1, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_2));
	double dW_3 = ((one_third / l_3) * (2 * mu / alpha) * pow(defJacobian, (-alpha * one_third)) * (2 * pow(l_3, alpha) - pow(l_1, alpha) - pow(l_2, alpha))) + ((2 / D_1) * (defJacobian - 1) * (defJacobian / l_3));

	MatrixXd tau = MatrixXd::Zero(ndof, ndof); // Kirchhoff stress
	MatrixXd pk_S = MatrixXd::Zero(ndof, ndof); // second Piola-Kirchhoff stress
	MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof); // Cauchy stress

	if (l_1 != l_2 && l_1 != l_3 && l_2 != l_3) {
		tau = (l_1 * dW_1 * n_1 * n_1.transpose()) + (l_2 * dW_2 * n_2 * n_2.transpose()) + (l_3 * dW_3 * n_3 * n_3.transpose());
	}

	if (l_1 == l_2 && l_3 != l_1) {
		tau = (l_1 * dW_1 * I) + ((l_3 * dW_3 - l_1 * dW_1) * n_3 * n_3.transpose());
	}

	if (l_1 == l_3 && l_2 != l_1) {
		tau = (l_1 * dW_1 * I) + ((l_2 * dW_2 - l_1 * dW_1) * n_2 * n_2.transpose());
	}

	if (l_2 == l_3 && l_1 != l_2) {
		tau = (l_2 * dW_2 * I) + ((l_1 * dW_1 - l_2 * dW_2) * n_1 * n_1.transpose());
	}

	if (l_1 == l_2 && l_1 == l_3) {
		tau = (l_1 * dW_1 * I);
	}

	MatrixXd F_inv = MatrixXd::Zero(ndof, ndof); // inverse of deformation gradient

	fe_invMatrix_pbr(F_inv, F);

	pk_S = F_inv * tau * F_inv.transpose();

	sigma_local = fe_tensor2voigt(pk_S);

	if (return_opt == 1) {
		cauchy_sigma = (1 / defJacobian) * F * pk_S * F.transpose();
		sigma_local = fe_tensor2voigt(cauchy_sigma);
	}
}
