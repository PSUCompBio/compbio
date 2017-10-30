#include "functions.h"

using namespace Eigen;

void fe_getForce_3d_embed(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, int embed_id, bool address_vr, VectorXi& embed_map, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d, VectorXd& delta_d, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& d_avg)
{

    MatrixXd* nodes_host     = mesh[host_id].getNewNodesPointer();
    MatrixXi* elements_host  = mesh[host_id].getNewElementsPointer();
    MatrixXd* nodes_embed    = mesh[embed_id].getNewNodesPointer();
    MatrixXi* elements_embed = mesh[embed_id].getNewElementsPointer();

    // Variables - Mesh Details
    int nel   = mesh[host_id].getNumElements();       /*! number of elements */
    int nnel  = mesh[host_id].getNumNodesPerElement(); // number of nodes per element
    int nnode = mesh[host_id].getNumNodes();          // number of nodes
    int sdof  = nnode * ndof;          // system degrees of freedom
    int edof  = nnel * ndof;           // element degrees of freedom

    VectorXd element_stress_host_local  = VectorXd::Zero(nel * 9);
    VectorXd element_strain_host_local  = VectorXd::Zero(nel * 9);
    VectorXd element_stress_embed_local = VectorXd::Zero((*elements_embed).rows() * 9);
    VectorXd element_strain_embed_local = VectorXd::Zero((*elements_embed).rows() * 9);
    VectorXd tmp_storage = VectorXd::Zero(ndof * ndof);

    // Element Data
    VectorXd xcoord      = VectorXd::Zero(nnel);
    VectorXd ycoord      = VectorXd::Zero(nnel);
    VectorXd zcoord      = VectorXd::Zero(nnel);

    VectorXd xcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd ycoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd zcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
    VectorXd u_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);
    VectorXd u_embed_local = VectorXd::Zero(((*elements_embed).cols() - 2) * ndof);
    VectorXd v_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);
    VectorXd a_embed       = VectorXd::Zero(((*nodes_embed).rows()) * ndof);

    double nfib_el      = 0; // number of fibers associated with an individual host element
    double d_tot_sum_el = 0; // sum of d_tot for all fibers associated with an individual host element

    for (int i = 0; i < nel; i++) {
        for (int j = 0; j < nnel; j++) {
            //  int g = -1;
            //  for(int f=0;f<nnode;f++){
            //      if(elements(i,j+2)==nodes(f,0)){
            //          g = f;
            //          break;
            //      }
            //  }

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

        double f_ext_e_sum = f_ext_e.sum();

        VectorXd f_int_e = VectorXd::Zero(edof);
        VectorXd f_tot_e = VectorXd::Zero(edof);
        VectorXd f_damp_e = VectorXd::Zero(edof);

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
                        double detJacobian = fe_detMatrix_pbr(jacobian);
                        fe_invMatrix_pbr(invJacobian, jacobian);

                        fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);

                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);
                        // disp_mat = fe_strDispMatrix(edof, nnel, dndx, dndy, dndz);

                        fe_stressUpdate_pbr(sigma_e, dndx, dndy, dndz, disp_mat, u_e, (*elements_host)(i, 1), 0);

                        // f_int_e = f_int_e + ((disp_mat.transpose())*sigma_e*wtx*wty*wtz*detJacobian); (previous correct)

                        f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx * wty * wtz * detJacobian);

                        if (f_ext_e_sum < 1e-18) { // only include damping when no external forces act on the element

                          // calculate bulk viscosity pressure that is linear in the volumetric strain rate
                          fe_getPressure_lbv_pbr(pressure_e, dndx, dndy, dndz, u_e, u_e_prev, dT, xcoord, ycoord, zcoord, (*elements_host)(i, 1));

                          // calculate internal damping force resulting from bulk viscosity pressure
                          f_damp_e = f_damp_e + ((disp_mat.transpose()) * pressure_e * wtx * wty * wtz * detJacobian);

                        }
                    }
                }
            }

            fe_calCentroidStress_3d_pbr(tmp_storage, nnel, xcoord, ycoord, zcoord, u_e, (*elements_host)(i, 1));
            element_stress_host_local.segment<9>(i * 9) = tmp_storage;
            fe_calCentroidStrain_3d_pbr(tmp_storage, nnel, xcoord, ycoord, zcoord, u_e);
            element_strain_host_local.segment<9>(i * 9) = tmp_storage;

            /*
             *  EMBEDDED FIBER ANALYSIS - STARTS FROM HERE
             *
             *  Step-1: Find the number of truss elements present this single host element using the mapping vector.
             *  Step-2: For that number of truss elements and using the mapping vector, one should perform the analysis and update the local element stress/strain information.
             */

            int ngl_embed          = 2;
            VectorXd points_embed  = guass_points(ngl_embed);
            VectorXd weights_embed = guass_weights(ngl_embed);
            VectorXd sigma_embed = VectorXd::Zero(6);
            VectorXd sigma_correction = VectorXd::Zero(6);

            d_tot_sum_el = 0;
            nfib_el = 0;

            for (int fib = 0; fib < (*elements_embed).rows(); fib++) {
                if (embed_map(fib) == (*elements_host)(i, 0)) {

                    nfib_el = nfib_el + 1; // count number of fibers associated with current host elment

                    for (int j = 0; j < ((*elements_embed).cols() - 2); j++) {
                        int g = (*elements_embed)(fib, j + 2);
                        xcoord_embed(j)      = (*nodes_embed)(g, 1);
                        ycoord_embed(j)      = (*nodes_embed)(g, 2);
                        zcoord_embed(j)      = (*nodes_embed)(g, 3);

                        VectorXd nat_cooord_nodes = VectorXd::Zero(3);
                        nat_cooord_nodes(0) = xcoord_embed(j);
                        nat_cooord_nodes(1) = ycoord_embed(j);
                        nat_cooord_nodes(2) = zcoord_embed(j);
                        VectorXd iso_coord_nodes = fe_newtonRhapson(nat_cooord_nodes, xcoord, ycoord, zcoord);
                        VectorXd shapes = fe_shapes_8(iso_coord_nodes(0), iso_coord_nodes(1), iso_coord_nodes(2));
                        MatrixXd shape_mat_embed = fe_shapeMatrix(edof, nnel, shapes);
                        u_embed.segment<3>(g * ndof)   = (shape_mat_embed * u_e);
                        u_embed_local.segment<3>(j * ndof) = u_embed.segment<3>(g * ndof);
                    }



                    //std::cout << "Embedded Displacements: \n" << u_embed << "\n";
                    //std::cout << "Local Embed Displacement: \n" << u_embed_local << "\n";


                    double length_embed = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed);

                    double length_embed_curr = fe_calCurrLength_pbr(u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed);

                    double lambda = length_embed_curr/length_embed;

                    fe_damageUpdate_pbr(d, fib, lambda);

                    fe_deltaDamageUpdate_pbr(delta_d, fib, lambda, lambda_min, lambda_max);

                    d_tot(fib) = d(fib) + delta_d(fib);

                    if (d_tot(fib) > 1) {
                      d_tot(fib) = 1;
                    }

                    d_tot_sum_el = d_tot_sum_el + d_tot(fib); // calculate sum of d_tot for fibers associated with current host element

                    for (int embed_intg = 0; embed_intg < ngl_embed; embed_intg++) {

                        VectorXd local_intg_points = fe_findIntgPoints_1d(xcoord_embed, ycoord_embed, zcoord_embed, points_embed(embed_intg), length_embed);
                        VectorXd global_intg_points = fe_newtonRhapson(local_intg_points, xcoord, ycoord, zcoord);

                        // VectorXd global_intg_points(3);
                        // global_intg_points(0) = 0;
                        // global_intg_points(1) = 0;
                        // global_intg_points(2) = points_embed(embed_intg);

                        double wtt = weights_embed(embed_intg);

                        fe_dniso_8(dndr, dnds, dndt, global_intg_points(0), global_intg_points(1), global_intg_points(2));

                        jacobian    = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                        fe_invMatrix_pbr(invJacobian, jacobian);
                        fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
                        fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);
                        fe_strDispMatrix_totalLagrangian_pbr(disp_mat, edof, nnel, dndx, dndy, dndz, u_e);

                        // Procedure - 1: (Same Deformation Gradient - Because No Slip)
                        fe_stressUpdate_pbr(sigma_embed, dndx, dndy, dndz, disp_mat, u_e, (*elements_embed)(fib, 1), 0);

                        // New concept that I am trying out...
                        // fe_stressModify(sigma_embed, xcoord_embed, ycoord_embed, zcoord_embed, 3);

                        VectorXd f_int_truss = (disp_mat.transpose() * (1 - d_tot(fib)) * sigma_embed * wtt * (length_embed / 2) * area_truss);

                        // Procedure - 2: (Same Displacements - Same Deformation Gradient - Transformation Matrix Inside)
                        /* sigma_truss = fe_stressUpdate_1d(elements_embed(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, dndx, dndy, dndz, u_e, 0);
                        VectorXd f_int_truss = (disp_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss); */

                        // Procedure - 3: (Same Displacements - Same Deformation Gradient - Transformation Matrix Outside)
                        /* MatrixXd stress_transformation_mat = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 1);
                        sigma_truss = fe_stressUpdate_1d(elements_embed(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                        VectorXd f_int_truss = (disp_mat.transpose() * stress_transformation_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss);*/

                        f_int_e = f_int_e + f_int_truss;

                        if (address_vr == true) {
                            fe_stressUpdate_pbr(sigma_correction, dndx, dndy, dndz, disp_mat, u_e, (*elements_host)(i, 1), 0);
                            VectorXd f_int_correction = (disp_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss);

                            /*sigma_correction = fe_stressUpdate_1d(elements_host(i, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                            VectorXd f_int_correction = (disp_mat.transpose() * stress_transformation_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss );*/

                            f_int_e = f_int_e - f_int_correction;
                        }

                        // f_int_e = f_int_e + f_int_truss - f_int_correction; // no Volume Redundancy
                        // f_int_e = (0.5*f_int_e) + f_int_truss; // no Volume Redundancy using volume fractions
                        // With Volume Redundancy
                    }

                    fe_calCentroidStrain_embed_3d_pbr(tmp_storage, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed);
                    element_strain_embed_local.segment<9>(fib * 9) = tmp_storage;
                    fe_calCentroidStress_embed_3d_pbr(tmp_storage, (*elements_embed)(fib, 1), u_e, u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, xcoord, ycoord, zcoord);
                    element_stress_embed_local.segment<9>(fib * 9) = (1 - d_tot(fib)) * tmp_storage;

                }
            }

            d_avg(i) = d_tot_sum_el/nfib_el; // calculate avg damage in current host element

        }

        //std::cout << "Host Final - nodal force: " << f_int_e.maxCoeff() << "\n";

        f_tot_e = f_ext_e - f_int_e - f_damp_e;

        fe_scatter_pbr(f_tot, f_tot_e, (*elements_host).block<1, 8>(i, 2), sdof);
        fe_scatter_pbr(f_damp, f_damp_e, (*elements_host).block<1, 8>(i, 2), sdof);
    }

    mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);
    mesh[host_id].readDamage(d_avg);
    mesh[embed_id].readElementStressStrain(element_stress_embed_local, element_strain_embed_local);
    mesh[embed_id].readDamage(d_tot);
    mesh[embed_id].readNodalKinematics(u_embed, v_embed, a_embed);

    nodes_host = NULL;
    elements_host = NULL;
    nodes_embed = NULL;
    elements_embed = NULL;

} // fe_getForce_3d_embed
