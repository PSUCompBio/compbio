#include "functions.h"

using namespace Eigen;

void fe_calCentroidStress_3d_pbr(VectorXd& element_stress, double dT, int nnel, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, VectorXd& u_e, int material_id, int iterator, int ix, int iy, int iz)
{
    int edof = nnel * ndof;

    VectorXd sigma_centroid = VectorXd::Zero(6);
    element_stress = VectorXd::Zero(ndof * ndof);
    double zero = 0.0;

    /** Strain Displacement Matrix of size 6 X edof (= nnel*ndof) */
    MatrixXd disp_mat(6, edof);
    VectorXd dndr(nnel);
    VectorXd dnds(nnel);
    VectorXd dndt(nnel);
    VectorXd dndx(nnel);
    VectorXd dndy(nnel);
    VectorXd dndz(nnel);
    MatrixXd jacobian(ndof, ndof);
    MatrixXd invJacobian(ndof, ndof);
    MatrixXd defGrad_centroid(ndof, ndof);
    MatrixXd invDefGrad_centroid(ndof, ndof);

    fe_dniso_8(dndr, dnds, dndt, zero, zero, zero);
    jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
    fe_invMatrix_pbr(invJacobian, jacobian);
    fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
    fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
    fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);

    fe_calDefGrad_pbr(defGrad_centroid, dndx, dndy, dndz, u_e);
    double defJacobian_centroid = fe_detMatrix_pbr(defGrad_centroid);
    fe_invMatrix_pbr(invDefGrad_centroid, defGrad_centroid);

    fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);

    // Stress- Calculation

    fe_stressUpdate_pbr(sigma_centroid, dndx, dndy, dndz, disp_mat, u_e, material_id, 1);

    if (include_viscoelasticity == 1) {
      // modifies sigma_e_centroid to include viscoelastic effects
      fe_stressUpdateViscoelasticity_pbr(sigma_centroid, dT, defGrad_centroid, invDefGrad_centroid, defJacobian_centroid, iterator, ix, iy, iz, material_id, 1);
    }

    element_stress(0) = sigma_centroid(0);
    element_stress(1) = sigma_centroid(3);
    element_stress(2) = sigma_centroid(5);
    element_stress(3) = sigma_centroid(3);
    element_stress(4) = sigma_centroid(1);
    element_stress(5) = sigma_centroid(4);
    element_stress(6) = sigma_centroid(5);
    element_stress(7) = sigma_centroid(4);
    element_stress(8) = sigma_centroid(2);

} // fe_calCentroidStress_3d_pbr

void fe_calCentroidStrain_3d_pbr(VectorXd& element_strain, int nnel, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, VectorXd& u_e)
{
    int edof = nnel * ndof;

    element_strain = VectorXd::Zero(ndof * ndof);
    double zero = 0.0;
    MatrixXd disp_mat(6, edof);
    VectorXd dndr(nnel);
    VectorXd dnds(nnel);
    VectorXd dndt(nnel);
    VectorXd dndx(nnel);
    VectorXd dndy(nnel);
    VectorXd dndz(nnel);
    MatrixXd jacobian(ndof, ndof);
    MatrixXd invJacobian(ndof, ndof);

    fe_dniso_8(dndr, dnds, dndt, zero, zero, zero);
    jacobian    = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
    // invJacobian = jacobian.inverse();
    fe_invMatrix_pbr(invJacobian, jacobian);
    fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
    fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
    fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);
    fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);

    // Strain - Calculation
    MatrixXd E = MatrixXd::Zero(3, 3);
    MatrixXd F = MatrixXd::Zero(3, 3);
    fe_calDefGrad_pbr(F, dndx, dndy, dndz, u_e);

    /* Polar Deomposition of F */
    JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV); // F = USV
    MatrixXd Usvd    = svd.matrixU();
    MatrixXd Vsvd    = svd.matrixV();
    VectorXd Ssvd    = svd.singularValues();
    MatrixXd Ssvdmat = MatrixXd::Zero(3, 3);
    Ssvdmat(0, 0) = Ssvd(0);
    Ssvdmat(1, 1) = Ssvd(1);
    Ssvdmat(2, 2) = Ssvd(2);
    MatrixXd Rpd = Usvd * Vsvd.adjoint(); //rotation tensor
    MatrixXd Upd = Vsvd * Ssvdmat * Vsvd.adjoint(); //material stretch
    MatrixXd Vpd = Usvd * Ssvdmat * Usvd.adjoint(); //spatial stretch

    E = Vpd.log();
    element_strain(0) = E(0, 0);
    element_strain(1) = E(0, 1);
    element_strain(2) = E(0, 2);
    element_strain(3) = E(1, 0);
    element_strain(4) = E(1, 1);
    element_strain(5) = E(1, 2);
    element_strain(6) = E(2, 0);
    element_strain(7) = E(2, 1);
    element_strain(8) = E(2, 2);

}

void fe_calCentroidStrain_embed_3d_pbr(VectorXd& element_strain, VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double length_old)
{

    element_strain = VectorXd::Zero(ndof * ndof);

    if (xcoord_embed.size() == 2) {
        VectorXd xcoord_embed_new(xcoord_embed.size()), ycoord_embed_new(ycoord_embed.size()), zcoord_embed_new(zcoord_embed.size());
        int counter = 0;
        for (int i = 0; i < xcoord_embed.size(); i++) {
            xcoord_embed_new(i) = xcoord_embed(i) + u_embed(counter);
            ycoord_embed_new(i) = ycoord_embed(i) + u_embed(counter + 1);
            zcoord_embed_new(i) = zcoord_embed(i) + u_embed(counter + 2);
            counter = counter + 3;
        }

        double length_new = fe_calVolume(xcoord_embed_new, ycoord_embed_new, zcoord_embed_new);

        double lambda = length_new / length_old;

        element_strain(0) = log(lambda);

    }

}

void fe_calCentroidStress_embed_3d_pbr(VectorXd& element_stress, int material_id, VectorXd& u_e, VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double length_old, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord) {

    element_stress = VectorXd::Zero(ndof * ndof);
    MatrixXd temp_stress = MatrixXd::Zero(ndof, ndof);

    if (xcoord_embed.size() == 2) {
        VectorXd local_intg_points = fe_findIntgPoints_1d(xcoord_embed, ycoord_embed, zcoord_embed, 0, length_old);
        VectorXd global_intg_points = fe_newtonRhapson(local_intg_points, xcoord, ycoord, zcoord);

        VectorXd sigma_truss = VectorXd::Zero(6);

        int nnel = xcoord.size();
        int edof = nnel*ndof;


        MatrixXd disp_mat(6, edof);

        VectorXd dndr(nnel);
        VectorXd dnds(nnel);
        VectorXd dndt(nnel);
        VectorXd dndx(nnel);
        VectorXd dndy(nnel);
        VectorXd dndz(nnel);
        MatrixXd jacobian(ndof, ndof);
        MatrixXd invJacobian(ndof, ndof);

        fe_dniso_8(dndr, dnds, dndt, global_intg_points(0), global_intg_points(1), global_intg_points(2));
        jacobian    = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
        // invJacobian = jacobian.inverse();
        fe_invMatrix_pbr(invJacobian, jacobian);
        fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
        fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
        fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);
        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);

        // Procedure - 1: (Same Deformation Gradient - Because No Slip)
        fe_stressUpdate_pbr(sigma_truss, dndx, dndy, dndz, disp_mat, u_e, material_id, 1);

        temp_stress = fe_voigt2tensor(sigma_truss);

        MatrixXd T = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3);
        temp_stress = T.transpose() * temp_stress * T;
        element_stress(0) = temp_stress(0,0);
    }

}
