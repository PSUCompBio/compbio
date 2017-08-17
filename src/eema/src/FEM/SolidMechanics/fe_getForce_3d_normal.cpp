#include "functions.h"

using namespace Eigen;

void fe_getForce_3d_normal(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, VectorXd& u_prev, double dT, VectorXd& fvd)
{

    MatrixXd* nodes_host    = mesh[host_id].getNewNodesPointer();
    MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();

    // Variables - Mesh Details
    int nel   = mesh[host_id].getNumElements();       /*! number of elements */
    int nnel  = mesh[host_id].getNumNodesPerElement(); // number of nodes per element
    int nnode = mesh[host_id].getNumNodes();          // number of nodes
    int sdof  = nnode * ndof;          // system degrees of freedom
    int edof  = nnel * ndof;           // element degrees of freedom

    VectorXd element_stress_host_local = VectorXd::Zero(nel * 9);
    VectorXd element_strain_host_local = VectorXd::Zero(nel * 9);
    VectorXd tmp_storage = VectorXd::Zero(ndof * ndof);

    // Element Data
    VectorXd xcoord      = VectorXd::Zero(nnel);
    VectorXd ycoord      = VectorXd::Zero(nnel);
    VectorXd zcoord      = VectorXd::Zero(nnel);

    for (int i = 0; i < nel; i++) {

        // std::cout << "Here...." << i << "\n";

        for (int j = 0; j < nnel; j++) {
            int g = (*elements_host)(i, j + 2);
            xcoord(j)      = (*nodes_host)(g, 1);
            ycoord(j)      = (*nodes_host)(g, 2);
            zcoord(j)      = (*nodes_host)(g, 3);
        }

        VectorXd u_e = VectorXd::Zero(edof); // element displacements
        fe_gather_pbr(u, u_e, (*elements_host).block<1, 8>(i, 2), sdof);

        VectorXd u_e_prev = VectorXd::Zero(edof); // previous element displacements
        fe_gather_pbr(u_prev, u_e_prev, (*elements_host).block<1, 8>(i, 2), sdof);

        VectorXd f_ext_e = VectorXd::Zero(edof);
        fe_gather_pbr(fext, f_ext_e, (*elements_host).block<1, 8>(i, 2), sdof); // element external nodal forces

        VectorXd f_int_e = VectorXd::Zero(edof);
        VectorXd f_tot_e = VectorXd::Zero(edof);
        VectorXd f_vd_e = VectorXd::Zero(edof);

        int nglx = 2;
        int ngly = 2;
        int nglz = 2;

        MatrixXd disp_mat(6, edof);

        VectorXd dndr(nnel);
        VectorXd dnds(nnel);
        VectorXd dndt(nnel);
        VectorXd dndx(nnel);
        VectorXd dndy(nnel);
        VectorXd dndz(nnel);
        MatrixXd jacobian(ndof, ndof);
        MatrixXd invJacobian(ndof, ndof);
        VectorXd sigma_e = VectorXd::Zero(6);
        VectorXd pressure_e = VectorXd::Zero(6);

        // MatrixXd points_3d = guass_points_3d(nglx,ngly,nglz);
        // MatrixXd weights_3d = guass_weights_3d(ndof,nglx,ngly,nglz);

        VectorXd points  = guass_points(nglx);
        VectorXd weights = guass_weights(nglx);
        int node_counter = 0;

        if (time_step_counter != 0) { // if this is not the first time step the go into the loop
            for (int intx = 0; intx < nglx; intx++) {
                double x   = points(intx);
                double wtx = weights(intx);
                for (int inty = 0; inty < ngly; inty++) {
                    double y   = points(inty);
                    double wty = weights(inty);
                    for (int intz = 0; intz < nglz; intz++) {
                        double z   = points(intz);
                        double wtz = weights(intz);

                        fe_dniso_8(dndr, dnds, dndt, x, y, z);

                        jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                        double detJacobian = jacobian.determinant();
                        invJacobian = jacobian.inverse();

                        fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);
                        // disp_mat = fe_strDispMatrix(edof, nnel, dndx, dndy, dndz);

                        fe_stressUpdate_pbr(sigma_e, dndx, dndy, dndz, disp_mat, u_e, (*elements_host)(i, 1), 0);

                        // f_int_e = f_int_e + ((disp_mat.transpose())*sigma_e*wtx*wty*wtz*detJacobian); (previous correct)
                        // std::cout<<k<<std::endl;

                        f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx * wty * wtz * detJacobian);

                        fe_getPressure_lbv_pbr(pressure_e, dndx, dndy, dndz, u_e, u_e_prev, dT, xcoord, ycoord, zcoord, (*elements_host)(i, 1));

                        f_vd_e = f_vd_e + ((disp_mat.transpose()) * (-1 * pressure_e) * wtx * wty * wtz * detJacobian);
                    }
                }
            }

            fe_calCentroidStress_3d_pbr(tmp_storage, nnel, xcoord, ycoord, zcoord, u_e, (*elements_host)(i, 1));
            element_stress_host_local.segment<9>(i * 9) = tmp_storage;
            fe_calCentroidStrain_3d_pbr(tmp_storage, nnel, xcoord, ycoord, zcoord, u_e);
            element_strain_host_local.segment<9>(i * 9) = tmp_storage;
        }
        f_tot_e = f_ext_e - f_int_e - f_vd_e;

        fe_scatter_pbr(f_tot, f_tot_e, (*elements_host).block<1, 8>(i, 2), sdof);
        fe_scatter_pbr(f_vd, f_vd_e, (*elements_host).block<1, 8>(i, 2), sdof);
    }

    mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);

    nodes_host = NULL;
    elements_host = NULL;
} // fe_getForce_3d_normal
