#include "functions.h"

using namespace Eigen;

void fe_mooneyrivlin_hyperelastic_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt)
{
    double I1, temp1, I2, defJacobian, defJacobian_frac, I1_bar, I2_bar, c1, c2, D, p;

    MatrixXd F = MatrixXd::Zero(ndof, ndof);
    MatrixXd C = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_square = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_bar = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_bar_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof);
    MatrixXd pk_S = MatrixXd::Zero(ndof, ndof);

    fe_calDefGrad_pbr(F, dndx, dndy, dndz, u);

    C = F.transpose() * F;

    fe_invMatrix_pbr(C_inv, C);

    I1         = C.trace();// First Invariant
    C_square = C * C; // Second Invariant
    temp1       = C_square.trace();
    I2         = (((pow(I1, 2) - (temp1)) * (0.5)));

    defJacobian      = fe_detMatrix_pbr(F); // Jacobian - determinant of deformation gradient
    defJacobian_frac = (1 / defJacobian);

    I1_bar  = I1 * (pow(defJacobian, (-1 * (2.0/3.0))));
    I2_bar  = I2 * (pow(defJacobian, (-1 * (4.0/3.0))));
    C_bar = C * (pow(defJacobian, (-1 * (2.0/3.0))));

    fe_invMatrix_pbr(C_bar_inv, C_bar);

    c1 = fe_get_mats(opt, 3, "mechanical");
    c2 = fe_get_mats(opt, 4, "mechanical");
    D  = fe_get_mats(opt, 1, "mechanical") * 0.5;

    p = -2 * D * (defJacobian - 1);

    pk_S = (-p * defJacobian * C_inv) + ( (2 * (pow(defJacobian, (-1 * (2.0/3.0))))) * ( ((c1 + (c2 * I1_bar)) * I) - (c2 * C_bar) - ((1.0/3.0) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * (C_bar_inv)  )));

    sigma = fe_tensor2voigt(pk_S); /** outputs 2nd cauchy stress tensor in vector form */

    if (return_opt == 1) {
        sigma = fe_tensor2voigt(defJacobian_frac * F * pk_S * F.transpose()); /** outputs cauchy stress tensor in vector form */
    }
}


void fe_mooneyrivlin_hyperelastic_pbr_array(VectorXd& sigma, int i, int x, int y, int z, VectorXd& u, int opt, int return_opt)
{
    double I1, temp1, I2, defJacobian, defJacobian_frac, I1_bar, I2_bar, c1, c2, D, p;

    MatrixXd F = MatrixXd::Zero(ndof, ndof);
    MatrixXd C = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_square = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_bar = MatrixXd::Zero(ndof, ndof);
    MatrixXd C_bar_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof);
    MatrixXd pk_S = MatrixXd::Zero(ndof, ndof);

    fe_calDefGrad_pbr_array(F, i, x, y, z, u);

    C = F.transpose() * F;

    fe_invMatrix_pbr(C_inv, C);

    I1         = C.trace();// First Invariant
    C_square = C * C; // Second Invariant
    temp1       = C_square.trace();
    I2         = (((pow(I1, 2) - (temp1)) * (0.5)));

    defJacobian      = fe_detMatrix_pbr(F); // Jacobian - determinant of deformation gradient
    defJacobian_frac = (1 / defJacobian);

    I1_bar  = I1 * (pow(defJacobian, (-1 * (2.0/3.0))));
    I2_bar  = I2 * (pow(defJacobian, (-1 * (4.0/3.0))));
    C_bar = C * (pow(defJacobian, (-1 * (2.0/3.0))));

    fe_invMatrix_pbr(C_bar_inv, C_bar);

    c1 = fe_get_mats(opt, 3, "mechanical");
    c2 = fe_get_mats(opt, 4, "mechanical");
    D  = fe_get_mats(opt, 1, "mechanical") * 0.5;

    p = -2 * D * (defJacobian - 1);

    pk_S = (-p * defJacobian * C_inv) + ( (2 * (pow(defJacobian, (-1 * (2.0/3.0))))) * ( ((c1 + (c2 * I1_bar)) * I) - (c2 * C_bar) - ((1.0/3.0) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * (C_bar_inv)  )));

    sigma = fe_tensor2voigt(pk_S); /** outputs 2nd cauchy stress tensor in vector form */

    if (return_opt == 1) {
        sigma = fe_tensor2voigt(defJacobian_frac * F * pk_S * F.transpose()); /** outputs cauchy stress tensor in vector form */
    }
}
