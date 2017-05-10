#include "functions.h"

using namespace Eigen;

VectorXd fe_mooneyrivlin_hyperelastic(VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt)
{

    VectorXd sigma = VectorXd::Zero(2 * ndof);

    double one_third  = 1.0 / 3.0;
    double two_third  = 2.0 / 3.0;
    double four_third = 4.0 / 3.0;
    double one_half   = 1.0 / 2.0;

    MatrixXd F = MatrixXd::Zero(ndof, ndof); // deformation gradient
    fe_calDefGrad_pbr(F, dndx, dndy, dndz, u);

    MatrixXd I = MatrixXd::Identity(ndof, ndof);

    MatrixXd C = MatrixXd::Zero(ndof, ndof);
    C = F.transpose() * F;

    double I1         = C.trace();// First Invariant
    MatrixXd C_square = C * C; // Second Invariant
    double tmp1       = C_square.trace();
    double I2         = (((pow(I1, 2) - (tmp1)) * (one_half)));
    // double I3 = B.determinant(); // Third Invariant

    double defJacobian      = F.determinant(); // Jacobian - determinant of deformation gradient
    double defJacobian_frac = (1 / defJacobian);

    double I1_bar  = I1 * (pow(defJacobian, (-1 * two_third)));
    double I2_bar  = I2 * (pow(defJacobian, (-1 * four_third)));
    MatrixXd C_bar = C * (pow(defJacobian, (-1 * two_third)));

    MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof);
    MatrixXd pk_S  = MatrixXd::Zero(ndof, ndof);

    double c1 = fe_get_mats(opt, 3);
    double c2 = fe_get_mats(opt, 4);
    double D  = fe_get_mats(opt, 1) / 2;

    double p = -2 * D * (defJacobian - 1);

    //cauchy_sigma = (((-p) * I) + ((defJacobian_frac) * (c1 + (I1_bar * c2)) * B_bar) - ((defJacobian_frac) * c2 * B_bar * B_bar) - ((defJacobian_frac) * (one_third) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * I));
    //MatrixXd F_inv = F.inverse();
    //pk_S = defJacobian * F_inv * cauchy_sigma * F_inv.transpose();

    pk_S = (-p * defJacobian * C.inverse()) + ( (2 * (pow(defJacobian, (-1 * two_third)))) * ( ((c1 + (c2 * I1_bar)) * I) - (c2 * C_bar) - ((one_third) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * (C_bar.inverse())  )));

    sigma = fe_tensor2voigt(pk_S); /** outputs 2nd cauchy stress tensor in vector form */

    if (return_opt == 1) {
        sigma = fe_tensor2voigt(defJacobian_frac * F * pk_S * F.transpose()); /** outputs cauchy stress tensor in vector form */
    }

    return sigma;

} // fe_mooneyrivlin_hyperelastic

/**********************************/

void fe_mooneyrivlin_hyperelastic_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt)
{

    double one_third  = 1.0 / 3.0;
    double two_third  = 2.0 / 3.0;
    double four_third = 4.0 / 3.0;
    double one_half   = 1.0 / 2.0;

    MatrixXd F = MatrixXd::Zero(ndof, ndof); // deformation gradient
    fe_calDefGrad_pbr(F, dndx, dndy, dndz, u);

    MatrixXd I = MatrixXd::Identity(ndof, ndof);

    MatrixXd C = MatrixXd::Zero(ndof, ndof);
    C = F.transpose() * F;

    double I1         = C.trace();// First Invariant
    MatrixXd C_square = C * C; // Second Invariant
    double tmp1       = C_square.trace();
    double I2         = (((pow(I1, 2) - (tmp1)) * (one_half)));
    // double I3 = B.determinant(); // Third Invariant

    double defJacobian      = F.determinant(); // Jacobian - determinant of deformation gradient
    double defJacobian_frac = (1 / defJacobian);

    double I1_bar  = I1 * (pow(defJacobian, (-1 * two_third)));
    double I2_bar  = I2 * (pow(defJacobian, (-1 * four_third)));
    MatrixXd C_bar = C * (pow(defJacobian, (-1 * two_third)));

    MatrixXd cauchy_sigma = MatrixXd::Zero(ndof, ndof);
    MatrixXd pk_S  = MatrixXd::Zero(ndof, ndof);

    double c1 = fe_get_mats(opt, 3);
    double c2 = fe_get_mats(opt, 4);
    double D  = fe_get_mats(opt, 1) / 2;

    double p = -2 * D * (defJacobian - 1);

    //cauchy_sigma = (((-p) * I) + ((defJacobian_frac) * (c1 + (I1_bar * c2)) * B_bar) - ((defJacobian_frac) * c2 * B_bar * B_bar) - ((defJacobian_frac) * (one_third) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * I));
    //MatrixXd F_inv = F.inverse();
    //pk_S = defJacobian * F_inv * cauchy_sigma * F_inv.transpose();

    pk_S = (-p * defJacobian * C.inverse()) + ( (2 * (pow(defJacobian, (-1 * two_third)))) * ( ((c1 + (c2 * I1_bar)) * I) - (c2 * C_bar) - ((one_third) * ((c1 * I1_bar) + (2 * c2 * I2_bar)) * (C_bar.inverse())  )));

    sigma = fe_tensor2voigt(pk_S); /** outputs 2nd cauchy stress tensor in vector form */

    if (return_opt == 1) {
        sigma = fe_tensor2voigt(defJacobian_frac * F * pk_S * F.transpose()); /** outputs cauchy stress tensor in vector form */
    }

}
