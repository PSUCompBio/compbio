#include "functions.h"

using namespace Eigen;

void fe_mooneyrivlin_hyperelastic_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt)
{
    fe_calDefGrad_pbr(F_mrh, dndx, dndy, dndz, u);

    C_mrh = F_mrh.transpose() * F_mrh;

    fe_invMatrix_pbr(C_inv_mrh, C_mrh);

    I1_mrh         = C_mrh.trace();// First Invariant
    C_square_mrh = C_mrh * C_mrh; // Second Invariant
    tmp1_mrh       = C_square_mrh.trace();
    I2_mrh         = (((pow(I1_mrh, 2) - (tmp1_mrh)) * (0.5)));

    defJacobian_mrh      = fe_detMatrix_pbr(F_mrh); // Jacobian - determinant of deformation gradient
    defJacobian_frac_mrh = (1 / defJacobian_mrh);

    I1_bar_mrh  = I1_mrh * (pow(defJacobian_mrh, (-1 * (2.0/3.0))));
    I2_bar_mrh  = I2_mrh * (pow(defJacobian_mrh, (-1 * (4.0/3.0))));
    C_bar_mrh = C_mrh * (pow(defJacobian_mrh, (-1 * (2.0/3.0))));

    fe_invMatrix_pbr(C_bar_inv_mrh, C_bar_mrh);

    c1_mrh = fe_get_mats(opt, 3, "mechanical");
    c2_mrh = fe_get_mats(opt, 4, "mechanical");
    D_mrh  = fe_get_mats(opt, 1, "mechanical") * 0.5;

    p_mrh = -2 * D_mrh * (defJacobian_mrh - 1);

    pk_S_mrh = (-p_mrh * defJacobian_mrh * C_inv_mrh) + ( (2 * (pow(defJacobian_mrh, (-1 * (2.0/3.0))))) * ( ((c1_mrh + (c2_mrh * I1_bar_mrh)) * I) - (c2_mrh * C_bar_mrh) - ((1.0/3.0) * ((c1_mrh * I1_bar_mrh) + (2 * c2_mrh * I2_bar_mrh)) * (C_bar_inv_mrh)  )));

    sigma = fe_tensor2voigt(pk_S_mrh); /** outputs 2nd cauchy stress tensor in vector form */

    if (return_opt == 1) {
        sigma = fe_tensor2voigt(defJacobian_frac_mrh * F_mrh * pk_S_mrh * F_mrh.transpose()); /** outputs cauchy stress tensor in vector form */
    }
}
