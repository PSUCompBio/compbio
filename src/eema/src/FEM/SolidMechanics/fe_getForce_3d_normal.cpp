#include "functions.h"

using namespace Eigen;

void fe_getForce_3d_normal(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, VectorXd& u_prev, double dT, VectorXd& f_damp, int t_plot)
{

    MatrixXd* nodes_host    = mesh[host_id].getNewNodesPointer();
    MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();

    // Variables - Mesh Details
    nel_normal   = mesh[host_id].getNumElements();       /*! number of elements */
    nnel_normal  = mesh[host_id].getNumNodesPerElement(); // number of nodes per element
    nnode_normal = mesh[host_id].getNumNodes();          // number of nodes
    sdof_normal  = nnode_normal * ndof;          // system degrees of freedom
    edof_normal  = nnel_normal * ndof;           // element degrees of freedom

    for (i_normal = 0; i_normal < nel_normal; i_normal++) {

        for (j_normal = 0; j_normal < nnel_normal; j_normal++) {
            g_normal = (*elements_host)(i_normal, j_normal + 2);
            xcoord_normal(j_normal)      = (*nodes_host)(g_normal, 1);
            ycoord_normal(j_normal)      = (*nodes_host)(g_normal, 2);
            zcoord_normal(j_normal)      = (*nodes_host)(g_normal, 3);
        }

        fe_gather_pbr(u, u_e_normal, (*elements_host).block<1, 8>(i_normal, 2), sdof_normal);

        fe_gather_pbr(u_prev, u_e_prev_normal, (*elements_host).block<1, 8>(i_normal, 2), sdof_normal);

        fe_gather_pbr(fext, f_ext_e_normal, (*elements_host).block<1, 8>(i_normal, 2), sdof_normal); // element external nodal forces

        f_ext_e_sum_normal = f_ext_e_normal.sum();

        VectorXd f_int_e = VectorXd::Zero(edof_normal);
        VectorXd f_tot_e = VectorXd::Zero(edof_normal);
        VectorXd f_damp_e = VectorXd::Zero(edof_normal);

        if (time_step_counter != 0) { // if this is not the first time step the go into the loop
            for (intx_normal = 0; intx_normal < 2; intx_normal++) {
                x_normal   = points_normal(intx_normal);
                wtx_normal = weights_normal(intx_normal);
                for (inty_normal = 0; inty_normal < 2; inty_normal++) {
                    y_normal   = points_normal(inty_normal);
                    wty_normal = weights_normal(inty_normal);
                    for (intz_normal = 0; intz_normal < 2; intz_normal++) {
                        z_normal   = points_normal(intz_normal);
                        wtz_normal = weights_normal(intz_normal);

                        fe_dniso_8_array(dndr_store, dnds_store, dndt_store, x_normal, y_normal, z_normal, intx_normal, inty_normal, intz_normal);
                        jacobian_normal = fe_calJacobian_array(ndof, nnel_normal, dndr_store[intx_normal][inty_normal][intz_normal], dnds_store[intx_normal][inty_normal][intz_normal], dndt_store[intx_normal][inty_normal][intz_normal], xcoord_normal, ycoord_normal, zcoord_normal);
                        detJacobian_normal = fe_detMatrix_pbr(jacobian_normal);
                        fe_invMatrix_pbr(invJacobian_normal, jacobian_normal);

                        fe_dndx_8_pbr_array(dndx_normal, nnel_normal, dndr_store[intx_normal][inty_normal][intz_normal], dnds_store[intx_normal][inty_normal][intz_normal], dndt_store[intx_normal][inty_normal][intz_normal], invJacobian_normal);
                        fe_dndy_8_pbr_array(dndy_normal, nnel_normal, dndr_store[intx_normal][inty_normal][intz_normal], dnds_store[intx_normal][inty_normal][intz_normal], dndt_store[intx_normal][inty_normal][intz_normal], invJacobian_normal);
                        fe_dndz_8_pbr_array(dndz_normal, nnel_normal, dndr_store[intx_normal][inty_normal][intz_normal], dnds_store[intx_normal][inty_normal][intz_normal], dndt_store[intx_normal][inty_normal][intz_normal], invJacobian_normal);

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat_normal, edof_normal, nnel_normal, dndx_normal, dndy_normal, dndz_normal, u_e_normal);

                        fe_stressUpdate_pbr(sigma_e_normal, dndx_normal, dndy_normal, dndz_normal, disp_mat_normal, u_e_normal, (*elements_host)(i_normal, 1), 0);

                        f_int_e = f_int_e + ((disp_mat_normal.transpose()) * sigma_e_normal * wtx_normal * wty_normal * wtz_normal * detJacobian_normal);

                        if (f_ext_e_sum_normal < 1e-18) { // only include damping when no external forces act on the element

                          // calculate bulk viscosity pressure that is linear in the volumetric strain rate
                          fe_getPressure_lbv_pbr(pressure_e_normal, dndx_normal, dndy_normal, dndz_normal, u_e_normal, u_e_prev_normal, dT, xcoord_normal, ycoord_normal, zcoord_normal, (*elements_host)(i_normal, 1));

                          // calculate internal damping force resulting from bulk viscosity pressure
                          f_damp_e = f_damp_e + ((disp_mat_normal.transpose()) * pressure_e_normal * wtx_normal * wty_normal * wtz_normal * detJacobian_normal);

                        }
                    }
                }
            }

            if (t_plot == 1) {
              fe_calCentroidStress_3d_pbr(tmp_storage_normal, nnel_normal, xcoord_normal, ycoord_normal, zcoord_normal, u_e_normal, (*elements_host)(i_normal, 1));
              element_stress_host_local_normal.segment<9>(i_normal * 9) = tmp_storage_normal;

              fe_calCentroidStrain_3d_pbr(tmp_storage_normal, nnel_normal, xcoord_normal, ycoord_normal, zcoord_normal, u_e_normal);
              element_strain_host_local_normal.segment<9>(i_normal * 9) = tmp_storage_normal;
            }
        }

        f_tot_e = f_ext_e_normal - f_int_e - f_damp_e;

        fe_scatter_pbr(f_tot, f_tot_e, (*elements_host).block<1, 8>(i_normal, 2), sdof_normal);
        fe_scatter_pbr(f_damp, f_damp_e, (*elements_host).block<1, 8>(i_normal, 2), sdof_normal);
    }

    if (t_plot == 1) {
        mesh[host_id].readElementStressStrain(element_stress_host_local_normal, element_strain_host_local_normal);
    }

    nodes_host = NULL;
    elements_host = NULL;
} // fe_getForce_3d_normal
