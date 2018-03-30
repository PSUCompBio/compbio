#include "functions.h"

using namespace Eigen;

void fe_getForce_3d_embed(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, int embed_id, bool address_vr, bool include_d, VectorXi& embed_map, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXd& d_avg, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t, int t_plot)
{
    VectorXd xcoord_embed      = VectorXd::Zero((*elements_embed_normal).cols() - 2);
    VectorXd ycoord_embed      = VectorXd::Zero((*elements_embed_normal).cols() - 2);
    VectorXd zcoord_embed      = VectorXd::Zero((*elements_embed_normal).cols() - 2);
    VectorXd u_embed       = VectorXd::Zero(((*nodes_embed_normal).rows()) * ndof);
    VectorXd u_embed_local = VectorXd::Zero(((*elements_embed_normal).cols() - 2) * ndof);
    VectorXd v_embed       = VectorXd::Zero(((*nodes_embed_normal).rows()) * ndof);
    VectorXd a_embed       = VectorXd::Zero(((*nodes_embed_normal).rows()) * ndof);

    double nfib_el      = 0; // number of fibers associated with an individual host element
    double d_tot_sum_el = 0; // sum of d_tot for all fibers associated with an individual host element

    for (i_normal = 0; i_normal < nel_normal; i_normal++) {
        
        for (j_normal = 0; j_normal < nnel_normal; j_normal++) {
            g_normal = (*elements_host_normal)(i_normal, j_normal + 2);
            xcoord_normal(j_normal) = (*nodes_host_normal)(g_normal, 1);
            ycoord_normal(j_normal) = (*nodes_host_normal)(g_normal, 2);
            zcoord_normal(j_normal) = (*nodes_host_normal)(g_normal, 3);
        }

        fe_gather_pbr(u, u_e_normal, (*elements_host_normal).block<1, 8>(i_normal, 2), sdof_normal);

        fe_gather_pbr(u_prev, u_e_prev_normal, (*elements_host_normal).block<1, 8>(i_normal, 2), sdof_normal);

        fe_gather_pbr(fext, f_ext_e_normal, (*elements_host_normal).block<1, 8>(i_normal, 2), sdof_normal); // element external nodal forces

        f_ext_e_sum_normal = f_ext_e_normal.sum();

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

                        jacobian = fe_calJacobian(ndof, nnel_normal, dndr, dnds, dndt, xcoord_normal, ycoord_normal, zcoord_normal);
                        double detJacobian = fe_detMatrix_pbr(jacobian);
                        fe_invMatrix_pbr(invJacobian, jacobian);

                        fe_dndx_8_pbr(dndx, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel_normal, dndr, dnds, dndt, invJacobian);

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat_normal, edof_normal, nnel_normal, dndx, dndy, dndz, u_e_normal);
                        // disp_mat_normal = fe_strDispMatrix(edof_normal, nnel_normal, dndx, dndy, dndz);

                        fe_stressUpdate_pbr(sigma_e, dndx, dndy, dndz, disp_mat_normal, u_e_normal, (*elements_host_normal)(i_normal, 1), 0);

                        // f_int_e = f_int_e + ((disp_mat_normal.transpose())*sigma_e*wtx*wty*wtz*detJacobian); (previous correct)

                        f_int_e = f_int_e + ((disp_mat_normal.transpose()) * sigma_e * wtx * wty * wtz * detJacobian);

                        if (f_ext_e_sum_normal < 1e-18) { // only include damping when no external forces act on the element

                          // calculate bulk viscosity pressure that is linear in the volumetric strain rate
                          fe_getPressure_lbv_pbr(pressure_e, dndx, dndy, dndz, u_e_normal, u_e_prev_normal, dT, xcoord_normal, ycoord_normal, zcoord_normal, (*elements_host_normal)(i_normal, 1));

                          // calculate internal damping force resulting from bulk viscosity pressure
                          f_damp_e = f_damp_e + ((disp_mat_normal.transpose()) * pressure_e * wtx * wty * wtz * detJacobian);

                        }
                    }
                }
            }

            if (t_plot == 1) {
                fe_calCentroidStress_3d_pbr(tmp_storage_normal, nnel_normal, xcoord_normal, ycoord_normal, zcoord_normal, u_e_normal, (*elements_host_normal)(i_normal, 1));
                element_stress_host_local_normal.segment<9>(i_normal * 9) = tmp_storage_normal;

                fe_calCentroidStrain_3d_pbr(tmp_storage_normal, nnel_normal, xcoord_normal, ycoord_normal, zcoord_normal, u_e_normal);
                element_strain_host_local_normal.segment<9>(i_normal * 9) = tmp_storage_normal;
            }

            /*
             *  EMBEDDED FIBER ANALYSIS - STARTS FROM HERE
             *
             *  Step-1: Find the number of truss elements present this single host element using the mapping vector.
             *  Step-2: For that number of truss elements and using the mapping vector, one should perform the analysis and update the local element stress/strain information.
             */

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

            for (int fib = 0; fib < (*elements_embed_normal).rows(); fib++) {
                if (embed_map(fib) == (*elements_host_normal)(i_normal, 0)) {

                    nfib_el = nfib_el + 1; // count number of fibers associated with current host elment

                    for (j_normal = 0; j_normal < ((*elements_embed_normal).cols() - 2); j_normal++) {
                        g_normal = (*elements_embed_normal)(fib, j_normal + 2);
                        xcoord_embed(j_normal)      = (*nodes_embed_normal)(g_normal, 1);
                        ycoord_embed(j_normal)      = (*nodes_embed_normal)(g_normal, 2);
                        zcoord_embed(j_normal)      = (*nodes_embed_normal)(g_normal, 3);

                        VectorXd nat_cooord_nodes = VectorXd::Zero(3);
                        nat_cooord_nodes(0) = xcoord_embed(j_normal);
                        nat_cooord_nodes(1) = ycoord_embed(j_normal);
                        nat_cooord_nodes(2) = zcoord_embed(j_normal);
                        VectorXd iso_coord_nodes = fe_newtonRhapson(nat_cooord_nodes, xcoord_normal, ycoord_normal, zcoord_normal);
                        VectorXd shapes = fe_shapes_8(iso_coord_nodes(0), iso_coord_nodes(1), iso_coord_nodes(2));
                        MatrixXd shape_mat_embed = fe_shapeMatrix(edof_normal, nnel_normal, shapes);
                        u_embed.segment<3>(g_normal * ndof)   = (shape_mat_embed * u_e_normal);
                        u_embed_local.segment<3>(j_normal * ndof) = u_embed.segment<3>(g_normal * ndof);
                    }

                    //std::cout << "Embedded Displacements: \n" << u_embed << "\n";
                    //std::cout << "Local Embed Displacement: \n" << u_embed_local << "\n";

                    double length_embed = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed);
                    double length_embed_orig = (*element_characteristic_embed_normal)(fib);
                    double length_embed_curr = fe_calCurrLength_pbr(u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed);
                    double lambda = length_embed_curr/length_embed;
                    MatrixXd F_fiber = MatrixXd::Zero(ndof, ndof);
                    MatrixXd T_fiber = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3); // transformation matrix from reference coordinate system to fiber coordinate system
                    MatrixXd T_fiber_inv = MatrixXd::Zero(ndof, ndof);
                    fe_invMatrix_pbr(T_fiber_inv, T_fiber);

                    if (include_d == 1) {

                      // length_embed_orig = (*element_characteristic_embed_normal)(fib);
                      // length_embed_curr = fe_calCurrLength_pbr(u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed);
                      // lambda = length_embed_curr/length_embed;

                      // fe_damageUpdate_pbr(d, fib, lambda); // This was the original definition of overload damage (similar to the Mullins effect). We are not using this anymore. However, I will leave the code in place for now. The structure might be useful later on. -JG, 2/8/2018

                      fe_fatigueDamageUpdate_pbr((*elements_embed_normal)(fib, 1), d_fatigue, fib, lambda, lambda_min_cycle, lambda_max_cycle, n_load_cycle_full, n_load_cycle_partial, t);

                      d_tot(fib) = d(fib) + d_fatigue(fib);
                      if (d_tot(fib) > 1) {
                        d_tot(fib) = 1;
                      }

                      d_tot_sum_el = d_tot_sum_el + d_tot(fib); // calculate sum of d_tot for fibers associated with current host element

                      fe_stressUpdateDamage_pbr(fiber_stress_PK2_damaged, (*elements_embed_normal)(fib, 1), fib, xcoord_embed, ycoord_embed, zcoord_embed, lambda, lambda_min, lambda_max, 0);

                    }


                    for (int embed_intg = 0; embed_intg < 2; embed_intg++) {

                        VectorXd local_intg_points = fe_findIntgPoints_1d(xcoord_embed, ycoord_embed, zcoord_embed, points_embed(embed_intg), length_embed);
                        VectorXd global_intg_points = fe_newtonRhapson(local_intg_points, xcoord_normal, ycoord_normal, zcoord_normal);

                        // VectorXd global_intg_points(3);
                        // global_intg_points(0) = 0;
                        // global_intg_points(1) = 0;
                        // global_intg_points(2) = points_embed(embed_intg);

                        double wtt = weights_embed(embed_intg);

                        fe_dniso_8(dndr, dnds, dndt, global_intg_points(0), global_intg_points(1), global_intg_points(2));

                        jacobian    = fe_calJacobian(ndof, nnel_normal, dndr, dnds, dndt, xcoord_normal, ycoord_normal, zcoord_normal);
                        fe_invMatrix_pbr(invJacobian, jacobian);
                        fe_dndx_8_pbr(dndx, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel_normal, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel_normal, dndr, dnds, dndt, invJacobian);

                        fe_calDefGrad_fiber_pbr(F_fiber, T_fiber, T_fiber_inv, lambda); // testing a concept...
                        fe_strDispMatrix_totalLagrangian_fiber_pbr(disp_mat_fiber, edof_normal, nnel_normal, dndx, dndy, dndz, u_e_normal, F_fiber); // testing a concept...
                        fe_stressUpdate_fiber_pbr(fiber_stress_PK2, (*elements_embed_normal)(fib, 1), lambda, T_fiber, T_fiber_inv, 0); // testing a concept...

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat_normal, edof_normal, nnel_normal, dndx, dndy, dndz, u_e_normal);

                        // Procedure - 1: (Same Deformation Gradient - Because No Slip)
                        fe_stressUpdate_pbr(sigma_embed, dndx, dndy, dndz, disp_mat_normal, u_e_normal, (*elements_embed_normal)(fib, 1), 0);

                        // New concept that I am trying out...
                        // fe_stressModify(sigma_embed, xcoord_embed, ycoord_embed, zcoord_embed, 3);

                        if (include_d == 0) {
                          // VectorXd f_int_truss = (disp_mat_normal.transpose() * sigma_embed * wtt * (length_embed / 2) * area_truss); // Official version. 2/18/2018
                          VectorXd f_int_truss = (disp_mat_fiber.transpose() * fiber_stress_PK2 * wtt * (length_embed / 2) * area_truss); // testing a concept...
                          f_int_e = f_int_e + f_int_truss;
                        }
                        if (include_d == 1) {
                          // VectorXd f_int_truss = (disp_mat_normal.transpose() * (1 - d_tot(fib)) * fiber_stress_PK2_damaged * wtt * (length_embed / 2) * area_truss); // Official version. 2/18/2018
                          VectorXd f_int_truss = (disp_mat_fiber.transpose() * (1 - d_tot(fib)) * fiber_stress_PK2_damaged * wtt * (length_embed / 2) * area_truss);
                          f_int_e = f_int_e + f_int_truss; // testing a concept...
                        }

                        // Procedure - 2: (Same Displacements - Same Deformation Gradient - Transformation Matrix Inside)
                        /* sigma_truss = fe_stressUpdate_1d(elements_embed_normal(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, dndx, dndy, dndz, u_e_normal, 0);
                        VectorXd f_int_truss = (disp_mat_normal.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss); */

                        // Procedure - 3: (Same Displacements - Same Deformation Gradient - Transformation Matrix Outside)
                        /* MatrixXd stress_transformation_mat = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 1);
                        sigma_truss = fe_stressUpdate_1d(elements_embed_normal(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                        VectorXd f_int_truss = (disp_mat_normal.transpose() * stress_transformation_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss);*/

                        if (address_vr == true) {
                            fe_stressUpdate_pbr(sigma_correction, dndx, dndy, dndz, disp_mat_normal, u_e_normal, (*elements_host_normal)(i_normal, 1), 0);
                            VectorXd f_int_correction = (disp_mat_normal.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss);

                            /*sigma_correction = fe_stressUpdate_1d(elements_host_normal(i_normal, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                            VectorXd f_int_correction = (disp_mat_normal.transpose() * stress_transformation_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss );*/

                            f_int_e = f_int_e - f_int_correction;
                        }

                        // f_int_e = f_int_e + f_int_truss - f_int_correction; // no Volume Redundancy
                        // f_int_e = (0.5*f_int_e) + f_int_truss; // no Volume Redundancy using volume fractions
                        // With Volume Redundancy
                    }

                    if (t_plot == 1) {
                        fe_calCentroidStrain_embed_3d_pbr(tmp_storage_normal, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed);
                        element_strain_embed_local_normal.segment<9>(fib * 9) = tmp_storage_normal;

                        fe_calCentroidStress_embed_3d_pbr(tmp_storage_normal, (*elements_embed_normal)(fib, 1), u_e_normal, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, xcoord_normal, ycoord_normal, zcoord_normal);

                        if (include_d == 0) {
                          // element_stress_embed_local_normal.segment<9>(fib * 9) = tmp_storage_normal;  // Official version. 2/18/2018
                          fe_stressUpdate_fiber_pbr(fiber_stress_cauchy, (*elements_embed_normal)(fib, 1), lambda, T_fiber, T_fiber_inv, 1); // testing a concept...
                          tmp_storage_normal = VectorXd::Zero(ndof * ndof); // testing a concept...
                          tmp_storage_normal(0) = fiber_stress_cauchy(0,0); // testing a concept...
                          element_stress_embed_local_normal.segment<9>(fib * 9) = (1 - d_tot(fib)) * tmp_storage_normal; // testing a concept...
                        }

                        if (include_d == 1) {
                          fe_stressUpdateDamage_pbr(fiber_stress_cauchy_damaged, (*elements_embed_normal)(fib, 1), fib, xcoord_embed, ycoord_embed, zcoord_embed, lambda, lambda_min, lambda_max, 1);
                          tmp_storage_normal = VectorXd::Zero(ndof * ndof);
                          tmp_storage_normal(0) = fiber_stress_cauchy_damaged(0,0);
                          element_stress_embed_local_normal.segment<9>(fib * 9) = (1 - d_tot(fib)) * tmp_storage_normal;
                        }
                    }

                }
            }

            if (nfib_el == 0) {
              d_avg(i_normal) = 0;
            }
            if (nfib_el != 0) {
              d_avg(i_normal) = d_tot_sum_el/nfib_el; // calculate avg damage in current host element
            }

        }

        //std::cout << "Host Final - nodal force: " << f_int_e.maxCoeff() << "\n";

        f_tot_e = f_ext_e_normal - f_int_e - f_damp_e;

        fe_scatter_pbr(f_tot, f_tot_e, (*elements_host_normal).block<1, 8>(i_normal, 2), sdof_normal);
        fe_scatter_pbr(f_damp, f_damp_e, (*elements_host_normal).block<1, 8>(i_normal, 2), sdof_normal);
    }

    if (t_plot == 1) {
        mesh[host_id].readElementStressStrain(element_stress_host_local_normal, element_strain_host_local_normal);
        mesh[host_id].readDamage(d_avg);
        mesh[embed_id].readElementStressStrain(element_stress_embed_local_normal, element_strain_embed_local_normal);
        mesh[embed_id].readDamage(d_tot);
        mesh[embed_id].readStretchMin(lambda_min);
        mesh[embed_id].readStretchMax(lambda_max);
        mesh[embed_id].readNodalKinematics(u_embed, v_embed, a_embed);
    }
}
