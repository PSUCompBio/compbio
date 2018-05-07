#include "functions.h"

using namespace Eigen;

void fe_getForce_3d_embed(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, int embed_id, bool address_vr, bool include_d, VectorXi& embed_map, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXd& d_avg, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t, int t_plot)
{
    int iterator, j, g;

    MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();
    MatrixXi* elements_embed = mesh[embed_id].getNewElementsPointer();

    MatrixXd* nodes_host = mesh[host_id].getNewNodesPointer();
    MatrixXd* nodes_embed = mesh[embed_id].getNewNodesPointer();

    VectorXd* element_characteristic_embed = mesh[embed_id].getElementCharacteristicPointer();

    VectorXd element_stress_embed_local = VectorXd::Zero((*elements_embed).rows() * 9);
    VectorXd element_strain_embed_local = VectorXd::Zero((*elements_embed).rows() * 9);

    MatrixXd defGrad = MatrixXd::Identity(ndof, ndof);
    MatrixXd invDefGrad = MatrixXd::Identity(ndof, ndof);
    MatrixXd disp_mat = MatrixXd::Zero(6, edof_normal);

    VectorXd element_stress_host_local = VectorXd::Zero(nel_normal * 9);
    VectorXd element_strain_host_local = VectorXd::Zero(nel_normal * 9);

    VectorXd xcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd ycoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd zcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd u_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);
    VectorXd u_embed_local = VectorXd::Zero(((*elements_embed).cols() - 2) * ndof);
    VectorXd v_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);
    VectorXd a_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);

    VectorXd tmp_storage = VectorXd::Zero(ndof * ndof);

    VectorXd xcoord = VectorXd::Zero(nnel_normal);
    VectorXd ycoord = VectorXd::Zero(nnel_normal);
    VectorXd zcoord = VectorXd::Zero(nnel_normal);

    VectorXd u_e = VectorXd::Zero(edof_normal);
    VectorXd u_e_prev = VectorXd::Zero(edof_normal);
    VectorXd f_ext_e = VectorXd::Zero(edof_normal);

    double nfib_el      = 0; // number of fibers associated with an individual host element
    double d_tot_sum_el = 0; // sum of d_tot for all fibers associated with an individual host element

    for (iterator = 0; iterator < nel_normal; iterator++) {

        double f_ext_e_sum, defJacobian;

        for (j = 0; j < nnel_normal; j++) {
            g = (*elements_host)(iterator, j + 2);
            xcoord(j) = (*nodes_host)(g, 1);
            ycoord(j) = (*nodes_host)(g, 2);
            zcoord(j) = (*nodes_host)(g, 3);
        }

        fe_gather_pbr(u, u_e, (*elements_host).block<1, 8>(iterator, 2), sdof_normal);

        fe_gather_pbr(u_prev, u_e_prev, (*elements_host).block<1, 8>(iterator, 2), sdof_normal);

        fe_gather_pbr(fext, f_ext_e, (*elements_host).block<1, 8>(iterator, 2), sdof_normal); // element external nodal forces

        f_ext_e_sum = f_ext_e.sum();

        VectorXd f_int_e = VectorXd::Zero(edof_normal);
        VectorXd f_tot_e = VectorXd::Zero(edof_normal);
        VectorXd f_damp_e = VectorXd::Zero(edof_normal);

        MatrixXd disp_mat_fiber(6, edof_normal); // testing a concept...

        VectorXd dndr(nnel_normal);
        VectorXd dnds(nnel_normal);
        VectorXd dndt(nnel_normal);
        VectorXd dndx(nnel_normal);
        VectorXd dndy(nnel_normal);
        VectorXd dndz(nnel_normal);
        MatrixXd jacobian(ndof, ndof);
        MatrixXd invJacobian(ndof, ndof);
        VectorXd sigma_e = VectorXd::Zero(6);
        VectorXd pressure_e = VectorXd::Zero(6);

        VectorXd points  = guass_points(2);
        VectorXd weights = guass_weights(2);
        int node_counter = 0;
        if (time_step_counter != 0) { // if this is not the first time step the go into the loop

            VectorXd points_embed  = guass_points(2);
            VectorXd weights_embed = guass_weights(2);
            VectorXd sigma_embed = VectorXd::Zero(6);
            VectorXd fiber_stress_cauchy = VectorXd::Zero(6);   // Cauchy fiber stress with respect to the fiber coordinate system in voigt notation
            VectorXd fiber_stress_PK2 = VectorXd::Zero(6);  // PK2 fiber stress with respect to the reference coordinate system in voigt notation
            VectorXd fiber_stress_cauchy_damaged = VectorXd::Zero(6);   // stress softened Cauchy fiber stress with respect to the fiber coordinate system in voigt notation
            VectorXd fiber_stress_PK2_damaged = VectorXd::Zero(6);  // stress softened PK2 fiber stress with respect to the reference coordinate system in voigt notation
            VectorXd sigma_correction = VectorXd::Zero(6);
            d_tot_sum_el = 0;
            nfib_el = 0;

            for (int fib = 0; fib < (*elements_embed).rows(); fib++) {
                if (embed_map(fib) == (*elements_host)(iterator, 0)) {

                    nfib_el = nfib_el + 1; // count number of fibers associated with current host elment

                    for (j = 0; j < ((*elements_embed).cols() - 2); j++) {
                        g = (*elements_embed)(fib, j + 2);
                        xcoord_embed(j)      = (*nodes_embed)(g, 1);
                        ycoord_embed(j)      = (*nodes_embed)(g, 2);
                        zcoord_embed(j)      = (*nodes_embed)(g, 3);

                        VectorXd nat_cooord_nodes = VectorXd::Zero(3);
                        nat_cooord_nodes(0) = xcoord_embed(j);
                        nat_cooord_nodes(1) = ycoord_embed(j);
                        nat_cooord_nodes(2) = zcoord_embed(j);
                        VectorXd iso_coord_nodes = fe_newtonRhapson(nat_cooord_nodes, xcoord, ycoord, zcoord);
                        VectorXd shapes = fe_shapes_8(iso_coord_nodes(0), iso_coord_nodes(1), iso_coord_nodes(2));
                        MatrixXd shape_mat_embed = fe_shapeMatrix(edof_normal, nnel_normal, shapes);
                        u_embed.segment<3>(g * ndof)   = (shape_mat_embed * u_e);
                        u_embed_local.segment<3>(j * ndof) = u_embed.segment<3>(g * ndof);
                    }

                    double length_embed = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed);
                    double length_embed_orig = (*element_characteristic_embed)(fib);
                    double length_embed_curr = fe_calCurrLength_pbr(u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed);
                    double lambda = length_embed_curr/length_embed;
                    MatrixXd F_fiber = MatrixXd::Zero(ndof, ndof);
                    MatrixXd T_fiber = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3); // transformation matrix from reference coordinate system to fiber coordinate system
                    MatrixXd T_fiber_inv = MatrixXd::Zero(ndof, ndof);
                    fe_invMatrix_pbr(T_fiber_inv, T_fiber);

                    if (include_d == 1) {

                      fe_stressUpdateDamage_pbr(fiber_stress_PK2_damaged, (*elements_embed)(fib, 1), fib, xcoord_embed, ycoord_embed, zcoord_embed, lambda, lambda_min, lambda_max, 0); // calculate current fiber stress based on Mullins effect

                      fe_staticDamageUpdate_pbr(d_static, fib, lambda_max); // calculate d_static for current fiber

                      fe_fatigueDamageUpdate_pbr((*elements_embed)(fib, 1), d_fatigue, fib, lambda, lambda_min_cycle, lambda_max_cycle, n_load_cycle_full, n_load_cycle_partial, t); // calculate d_fatigue for current fiber

                      d_tot(fib) = d_static(fib) + d_fatigue(fib);
                      if (d_tot(fib) > 1) {
                        d_tot(fib) = 1;
                      }

                      d_tot_sum_el = d_tot_sum_el + d_tot(fib); // calculate sum of d_tot for fibers associated with current host element
                    }


                    for (int embed_intg = 0; embed_intg < 2; embed_intg++) {

                        VectorXd local_intg_points = fe_findIntgPoints_1d(xcoord_embed, ycoord_embed, zcoord_embed, points_embed(embed_intg), length_embed);
                        VectorXd global_intg_points = fe_newtonRhapson(local_intg_points, xcoord, ycoord, zcoord);

                        double wtt = weights_embed(embed_intg);

                        fe_dniso_8(dndr, dnds, dndt, global_intg_points(0), global_intg_points(1), global_intg_points(2));

                        jacobian    = fe_calJacobian(ndof, nnel_normal, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                        fe_invMatrix_pbr(invJacobian, jacobian);
                        fe_dndx_8_pbr(dndx, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel_normal, dndr, dnds, dndt, invJacobian);

                        fe_calDefGrad_fiber_pbr(F_fiber, T_fiber, T_fiber_inv, lambda); // testing a concept...
                        fe_strDispMatrix_totalLagrangian_fiber_pbr(disp_mat_fiber, edof_normal, nnel_normal, dndx, dndy, dndz, u_e, F_fiber); // testing a concept...
                        fe_stressUpdate_fiber_pbr(fiber_stress_PK2, (*elements_embed)(fib, 1), lambda, T_fiber, T_fiber_inv, 0); // testing a concept...

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof_normal, nnel_normal, dndx, dndy, dndz, u_e);

                        fe_stressUpdate_pbr(sigma_embed, dndx, dndy, dndz, disp_mat, u_e, (*elements_embed)(fib, 1), 0);

                        if (include_d == 0) {
                          // VectorXd f_int_truss = (disp_mat.transpose() * sigma_embed * wtt * (length_embed / 2) * area_truss); // Official version. 2/18/2018
                          VectorXd f_int_truss = (disp_mat_fiber.transpose() * fiber_stress_PK2 * wtt * (length_embed / 2) * area_truss); // testing a concept...
                          f_int_e = f_int_e + f_int_truss;
                        }
                        if (include_d == 1) {
                          // VectorXd f_int_truss = (disp_mat.transpose() * (1 - d_tot(fib)) * fiber_stress_PK2_damaged * wtt * (length_embed / 2) * area_truss); // Official version. 2/18/2018
                          VectorXd f_int_truss = (disp_mat_fiber.transpose() * (1 - d_tot(fib)) * fiber_stress_PK2_damaged * wtt * (length_embed / 2) * area_truss); // testing a concept...
                          f_int_e = f_int_e + f_int_truss;
                        }

                        if (address_vr == true) {
                            fe_stressUpdate_pbr(sigma_correction, dndx, dndy, dndz, disp_mat, u_e, (*elements_host)(iterator, 1), 0);
                            VectorXd f_int_correction = (disp_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss);

                            f_int_e = f_int_e - f_int_correction;
                        }

                    } // end of fiber element quad point loop

                    if (t_plot == 1) {
                        fe_calCentroidStrain_embed_3d_pbr(tmp_storage, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed);
                        element_strain_embed_local.segment<9>(fib * 9) = tmp_storage;

                        fe_calCentroidStress_embed_3d_pbr(tmp_storage, (*elements_embed)(fib, 1), u_e, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, xcoord, ycoord, zcoord);

                        if (include_d == 0) {
                          // element_stress_embed_local.segment<9>(fib * 9) = tmp_storage;  // Official version. 2/18/2018
                          fe_stressUpdate_fiber_pbr(fiber_stress_cauchy, (*elements_embed)(fib, 1), lambda, T_fiber, T_fiber_inv, 1); // testing a concept...
                          tmp_storage = VectorXd::Zero(ndof * ndof); // testing a concept...
                          tmp_storage(0) = fiber_stress_cauchy(0,0); // testing a concept...
                          element_stress_embed_local.segment<9>(fib * 9) = tmp_storage; // testing a concept...
                        }

                        if (include_d == 1) {
                          fe_stressUpdateDamage_pbr(fiber_stress_cauchy_damaged, (*elements_embed)(fib, 1), fib, xcoord_embed, ycoord_embed, zcoord_embed, lambda, lambda_min, lambda_max, 1);
                          tmp_storage = VectorXd::Zero(ndof * ndof);
                          tmp_storage(0) = fiber_stress_cauchy_damaged(0,0);
                          element_stress_embed_local.segment<9>(fib * 9) = (1 - d_tot(fib)) * tmp_storage;
                        }
                    }

                }

            } // end of fiber element loop

            if (nfib_el == 0) {
              d_avg(iterator) = 0;
            }
            if (nfib_el != 0) {
              d_avg(iterator) = d_tot_sum_el/nfib_el; // calculate average damage in current host element
            }

            // *** beginning of moved code section ***

            for (int intx = 0; intx < 2; intx++) {
                double x   = points(intx);
                double wtx = weights(intx);
                for (int inty = 0; inty < 2; inty++) {
                    double y   = points(inty);
                    double wty = weights(inty);
                    for (int intz = 0; intz < 2; intz++) {
                        double z   = points(intz);
                        double wtz = weights(intz);

                        fe_dniso_8(dndr, dnds, dndt, x, y, z);

                        jacobian = fe_calJacobian(ndof, nnel_normal, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                        double detJacobian = fe_detMatrix_pbr(jacobian);
                        fe_invMatrix_pbr(invJacobian, jacobian);

                        fe_dndx_8_pbr(dndx, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel_normal, dndr, dnds, dndt, invJacobian);

                        fe_calDefGrad_pbr(defGrad, dndx, dndy, dndz, u_e); // In the future, reference defGrad in other functions too. We repeat this calculation many times.
                        defJacobian = fe_detMatrix_pbr(defGrad);
                        fe_invMatrix_pbr(invDefGrad, defGrad);

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof_normal, nnel_normal, dndx, dndy, dndz, u_e);

                        fe_stressUpdate_pbr(sigma_e, dndx, dndy, dndz, disp_mat, u_e, (*elements_host)(iterator, 1), 0);

                        if (include_viscoelasticity == 1) {
                          // modifies sigma_e to include viscoelastic effects
                          fe_stressUpdateViscoelasticity_pbr(sigma_e, dT, defGrad, invDefGrad, defJacobian, iterator, intx, inty, intz, (*elements_host)(iterator, 1), 0);
                        }

                        if (include_d == 0) {
                          f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx * wty * wtz * detJacobian);
                        }

                        if (include_d == 1) {
                          f_int_e = f_int_e + ((disp_mat.transpose()) * (1 - d_avg(iterator)) * sigma_e * wtx * wty * wtz * detJacobian);
                        }

                        if (f_ext_e_sum < 1e-18) { // only include damping when no external forces act on the element

                          // calculate bulk viscosity pressure that is linear in the volumetric strain rate
                          fe_getPressure_lbv_pbr(pressure_e, dndx, dndy, dndz, u_e, u_e_prev, dT, xcoord, ycoord, zcoord, (*elements_host)(iterator, 1));

                          // calculate internal damping force resulting from bulk viscosity pressure
                          f_damp_e = f_damp_e + ((disp_mat.transpose()) * pressure_e * wtx * wty * wtz * detJacobian);

                        }
                    }
                }

            } // end of host element quad point loop

            if (t_plot == 1) {
              fe_calCentroidStrain_3d_pbr(tmp_storage, nnel_normal, xcoord, ycoord, zcoord, u_e);
              element_strain_host_local.segment<9>(iterator * 9) = tmp_storage;
            }

            if (include_viscoelasticity == 0) {
              if (t_plot == 1) {
                fe_calCentroidStress_3d_pbr(tmp_storage, dT, nnel_normal, xcoord, ycoord, zcoord, u_e, (*elements_host)(iterator, 1), iterator, 0, 0, 0);
                if (include_d == 0) {
                  element_stress_host_local.segment<9>(iterator * 9) = tmp_storage;
                }
                if (include_d == 1) {
                  element_stress_host_local.segment<9>(iterator * 9) = (1 - d_avg(iterator)) * tmp_storage;
                }
              }
            }

            if (include_viscoelasticity == 1) {
              fe_calCentroidStress_3d_pbr(tmp_storage, dT, nnel_normal, xcoord, ycoord, zcoord, u_e, (*elements_host)(iterator, 1), iterator, 0, 0, 0);
              if (include_d == 0) {
                element_stress_host_local.segment<9>(iterator * 9) = tmp_storage;
              }
              if (include_d == 1) {
                element_stress_host_local.segment<9>(iterator * 9) = (1 - d_avg(iterator)) * tmp_storage;
              }
            }

            // *** end of moved code section ***

        } // end of if statement

        f_tot_e = f_ext_e - f_int_e - f_damp_e;

        fe_scatter_pbr(f_tot, f_tot_e, (*elements_host).block<1, 8>(iterator, 2), sdof_normal);
        fe_scatter_pbr(f_damp, f_damp_e, (*elements_host).block<1, 8>(iterator, 2), sdof_normal);

    } // end of host element loop

    if (t_plot == 1) {
        mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);
        mesh[host_id].readDamage(d_avg);
        mesh[embed_id].readElementStressStrain(element_stress_embed_local, element_strain_embed_local);
        mesh[embed_id].readDamage(d_tot);
        mesh[embed_id].readStretchMin(lambda_min);
        mesh[embed_id].readStretchMax(lambda_max);
        mesh[embed_id].readNodalKinematics(u_embed, v_embed, a_embed);
    }
}
