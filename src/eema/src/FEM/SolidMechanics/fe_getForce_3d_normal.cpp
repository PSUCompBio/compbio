#include "functions.h"
#include <thread>
#include <mutex>

using namespace Eigen;

std::mutex m;
VectorXd global_tot, global_damp;


void task(int start, int total, VectorXd u, VectorXd fext, int time_step_counter, int host_id, VectorXd u_prev, double dT, double t, int t_plot, VectorXd element_stress_host_local, VectorXd element_strain_host_local) {
    int i, j, g, ix, iy, iz;
    double x, y, z;

    VectorXd tmp_storage = VectorXd::Zero(ndof * ndof);
    VectorXd xcoord = VectorXd::Zero(nnel_normal);
    VectorXd ycoord = VectorXd::Zero(nnel_normal);
    VectorXd zcoord = VectorXd::Zero(nnel_normal);

    MatrixXd* nodes_host    = mesh[host_id].getNewNodesPointer();
    MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();

    VectorXd u_e = VectorXd::Zero(edof_normal);
    VectorXd u_e_prev = VectorXd::Zero(edof_normal);
    VectorXd f_ext_e = VectorXd::Zero(edof_normal);
    VectorXd sigma_e = VectorXd::Zero(6);
    VectorXd pressure_e = VectorXd::Zero(6);
    VectorXd local_total = VectorXd::Zero(sdof_normal);
    VectorXd damp_total = VectorXd::Zero(sdof_normal);
    VectorXd f_tot_e = VectorXd::Zero(edof_normal);

    MatrixXd defGrad = MatrixXd::Identity(ndof, ndof);
    MatrixXd invDefGrad = MatrixXd::Identity(ndof, ndof);
    MatrixXd disp_mat = MatrixXd::Zero(6, edof_normal);

    for (i = start; i < start + total; i++) {
        double f_ext_e_sum, defJacobian;

        VectorXd f_int_e = VectorXd::Zero(edof_normal);
        VectorXd f_damp_e = VectorXd::Zero(edof_normal);

        if (t == 0) {
            for (j = 0; j < nnel_normal; j++) {
                g = (*elements_host)(i, j + 2);
                x_store[i][j] = (*nodes_host)(g, 1);
                y_store[i][j] = (*nodes_host)(g, 2);
                z_store[i][j] = (*nodes_host)(g, 3);
            }

            for (ix = 0; ix < 2; ix++) {
                x = points_normal(ix);
                wtx_normal[ix] = weights_normal(ix);

                for (iy = 0; iy < 2; iy++) {
                    y = points_normal(iy);
                    wty_normal[ix][iy] = weights_normal(iy);

                    for (iz = 0; iz < 2; iz++) {
                        z = points_normal(iz);
                        wtz_normal[ix][iy][iz] = weights_normal(iz);

                        fe_dniso_8_array(dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], x, y, z, ix, iy, iz);

                        fe_calJacobian_array(jacobian_store[i][ix][iy][iz], nnel_normal, dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], x_store[i], y_store[i], z_store[i]);

                        det_store[i][ix][iy][iz] = fe_detMatrix_pbr_array(jacobian_store[i][ix][iy][iz]);

                        fe_invMatrix_pbr_array(invJacobian_store[i][ix][iy][iz], jacobian_store[i][ix][iy][iz], det_store[i][ix][iy][iz]);

                        fe_dndxyz_8_pbr_array(dndx_store[i][ix][iy][iz], dndy_store[i][ix][iy][iz], dndz_store[i][ix][iy][iz], nnel_normal, dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], invJacobian_store[i][ix][iy][iz]);
                    }
                }
            }
        }

        for (j = 0; j < nnel_normal; j++) {
            xcoord(j) = x_store[i][j];
            ycoord(j) = y_store[i][j];
            zcoord(j) = z_store[i][j];
        }

        fe_gather_pbr(u, u_e, ((*elements_host).block<1, 8>(i, 2)), sdof_normal);

        fe_gather_pbr(u_prev, u_e_prev, (*elements_host).block<1, 8>(i, 2), sdof_normal);

        fe_gather_pbr(fext, f_ext_e, (*elements_host).block<1, 8>(i, 2), sdof_normal);

        f_ext_e_sum = f_ext_e.sum();

        if (time_step_counter != 0) {
            for (ix = 0; ix < 2; ix++) {
                for (iy = 0; iy < 2; iy++) {
                    for (iz = 0; iz < 2; iz++) {

                        fe_calDefGrad_pbr_array(defGrad, i, ix, iy, iz, u_e);

                        defJacobian = fe_detMatrix_pbr(defGrad);

                        fe_invMatrix_pbr(invDefGrad, defGrad);

                        fe_strDispMatrix_totalLagrangian_pbr_array(disp_mat, edof_normal, nnel_normal, i, ix, iy, iz, u_e, dndx_store[i][ix][iy][iz], dndy_store[i][ix][iy][iz], dndz_store[i][ix][iy][iz]);

                        fe_stressUpdate_pbr_array(sigma_e, i, ix, iy, iz, disp_mat, u_e, (*elements_host)(i, 1), 0);

                        if (include_viscoelasticity == 1)
                            fe_stressUpdateViscoelasticity_pbr(sigma_e, dT, defGrad, invDefGrad, defJacobian, i, ix, iy, iz, (*elements_host)(i, 1), 0);

                        f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx_normal[ix] * wty_normal[ix][iy] * wtz_normal[ix][iy][iz] * det_store[i][ix][iy][iz]);

                        if (f_ext_e_sum < 1e-18) {
                            fe_getPressure_lbv_pbr_array(pressure_e, i, ix, iy, iz, u_e, u_e_prev, dT, xcoord, ycoord, zcoord, (*elements_host)(i, 1));
                            f_damp_e = f_damp_e + ((disp_mat.transpose()) * pressure_e * wtx_normal[ix] * wty_normal[ix][iy] * wtz_normal[ix][iy][iz] * det_store[i][ix][iy][iz]);
                        }
                    }
                }
            }

            if (include_viscoelasticity == 0) {
                if (t_plot == 1) {
                    fe_calCentroidStress_3d_pbr(tmp_storage, dT, nnel_normal, xcoord, ycoord, zcoord, u_e, (*elements_host)(i, 1));
                    element_stress_host_local.segment<9>(i * 9) = tmp_storage;

                    fe_calCentroidStrain_3d_pbr(tmp_storage, nnel_normal, xcoord, ycoord, zcoord, u_e);
                    element_strain_host_local.segment<9>(i * 9) = tmp_storage;
                }
            }

            else {
                fe_calCentroidStress_3d_pbr(tmp_storage, dT, nnel_normal, xcoord, ycoord, zcoord, u_e, (*elements_host)(i, 1));
                element_stress_host_local.segment<9>(i * 9) = tmp_storage;

                fe_calCentroidStrain_3d_pbr(tmp_storage, nnel_normal, xcoord, ycoord, zcoord, u_e);
                element_strain_host_local.segment<9>(i * 9) = tmp_storage;
            }
        }

        f_tot_e = f_ext_e - f_int_e - f_damp_e;

        fe_scatter_pbr(local_total, f_tot_e, (*elements_host).block<1, 8>(i, 2), sdof_normal);
        fe_scatter_pbr(damp_total, f_damp_e, (*elements_host).block<1, 8>(i, 2), sdof_normal);
    }

    m.lock();
    global_tot += local_total;
    global_damp += damp_total;
    m.unlock();
}

void fe_getForce_3d_normal(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, VectorXd& u_prev, double dT, VectorXd& f_damp, double t, int t_plot)
{
    int i, number_of_threads = 4;

    VectorXd element_stress_host_local = VectorXd::Zero(nel_normal * 9);
    VectorXd element_strain_host_local = VectorXd::Zero(nel_normal * 9);

    std::thread thread_runner[number_of_threads];

    global_tot = VectorXd::Zero(sdof_normal);
    global_damp = VectorXd::Zero(sdof_normal);

    global_tot += f_tot;
    global_damp += f_damp;

    for (i = 0; i < number_of_threads; i++)
        thread_runner[i] = std::thread(task, i*2, 2, u, fext, time_step_counter, host_id, u_prev, dT, t, t_plot, element_stress_host_local, element_strain_host_local);

    for (i = 0; i < number_of_threads; i++)
		thread_runner[i].join();

    f_tot = global_tot;
    f_damp = global_damp;

    if (t_plot == 1) {
        mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);
    }
}
