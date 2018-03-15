#include "functions.h"
using namespace Eigen;

double eps_energy = 0.01;
double failure_time_step = 1e-8;

double ****dndr_store;
double ****dnds_store;
double ****dndt_store;

int counter_test = 0;
MatrixXd I = MatrixXd::Identity(3, 3);

double detJacobian_normal = 0;

VectorXd weights_normal = guass_weights(2);

MatrixXd jacobian_normal = MatrixXd::Zero(3, 3);
MatrixXd invJacobian_normal = MatrixXd::Zero(3, 3);
MatrixXd disp_mat_normal = MatrixXd::Zero(6, 24);

VectorXd dndx_normal = VectorXd::Zero(8);
VectorXd dndy_normal = VectorXd::Zero(8);
VectorXd dndz_normal = VectorXd::Zero(8);

VectorXd xcoord_normal = VectorXd::Zero(8);
VectorXd ycoord_normal = VectorXd::Zero(8);
VectorXd zcoord_normal = VectorXd::Zero(8);

void experimental() {

    int nnel = mesh[0].getNumNodesPerElement();
    VectorXd points  = guass_points(2);
    MatrixXd* nodes_host    = mesh[0].getNewNodesPointer();
    MatrixXi* elements_host = mesh[0].getNewElementsPointer();

    dndr_store = new double***[2];
    for (int i = 0; i < 2; i++) {
        dndr_store[i] = new double**[2];
        for (int j = 0; j < 2; j++) {
            dndr_store[i][j] = new double*[2];
            for (int k = 0; k < 2; k++)
                dndr_store[i][j][k] = new double[8];
            }
        }

    dnds_store = new double***[2];
    for (int i = 0; i < 2; i++) {
        dnds_store[i] = new double**[2];
        for (int j = 0; j < 2; j++) {
            dnds_store[i][j] = new double*[2];
            for (int k = 0; k < 2; k++)
                dnds_store[i][j][k] = new double[8];
            }
        }

    dndt_store = new double***[2];
    for (int i = 0; i < 2; i++) {
        dndt_store[i] = new double**[2];
        for (int j = 0; j < 2; j++) {
            dndt_store[i][j] = new double*[2];
            for (int k = 0; k < 2; k++)
                dndt_store[i][j][k] = new double[8];
            }
        }

    for (int j = 0; j < nnel; j++) {
        int g = (*elements_host)(0, j + 2);
        xcoord_normal(j) = (*nodes_host)(g, 1);
        ycoord_normal(j) = (*nodes_host)(g, 2);
        zcoord_normal(j) = (*nodes_host)(g, 3);
    }

    for (int intx = 0; intx < 2; intx++) {
        double x   = points(intx);
        for (int inty = 0; inty < 2; inty++) {
            double y   = points(inty);
            for (int intz = 0; intz < 2; intz++) {
                double z   = points(intz);
                fe_dniso_8_array(dndr_store, dnds_store, dndt_store, x, y, z, intx, inty, intz);

                jacobian_normal = fe_calJacobian_array(3, nnel, dndr_store[intx][inty][intz], dnds_store[intx][inty][intz], dndt_store[intx][inty][intz], xcoord_normal, ycoord_normal, zcoord_normal);
                detJacobian_normal = fe_detMatrix_pbr(jacobian_normal);
                fe_invMatrix_pbr(invJacobian_normal, jacobian_normal);

                fe_dndx_8_pbr_array(dndx_normal, nnel, dndr_store[intx][inty][intz], dnds_store[intx][inty][intz], dndt_store[intx][inty][intz], invJacobian_normal);
                fe_dndy_8_pbr_array(dndy_normal, nnel, dndr_store[intx][inty][intz], dnds_store[intx][inty][intz], dndt_store[intx][inty][intz], invJacobian_normal);
                fe_dndz_8_pbr_array(dndz_normal, nnel, dndr_store[intx][inty][intz], dnds_store[intx][inty][intz], dndt_store[intx][inty][intz], invJacobian_normal);
            }
        }
    }
}

/*! \brief
 * This function carries out the explicit dynamic analysis of the FEM problem.
 */

void
fe_mainEXPLICIT()
{
    experimental();
    // Following variables - Only for Hex Element
    int nnode = mesh[0].getNumNodes();            // number of nodes
    int sdof  = nnode * ndof;                     // system degrees of freedom
    int nel   = mesh[0].getNumElements();         // number of elements

    // Initialization
    double dT             = dt_initial;
    double t              = t_start;
    int time_step_counter = 0; // time step count
    int plot_state_counter = 1;
    int t_plot = 1;

    double output_temp_1  = ((double) (t_end / output_frequency));
    VectorXd A            = VectorXd::Zero(sdof); // Acceleration Vector
    VectorXd V            = VectorXd::Zero(sdof); // Velocity Vector
    VectorXd V_half       = VectorXd::Zero(sdof); // Velocity Vector at n+1/2
    VectorXd U            = VectorXd::Zero(sdof); // Displacement Vector
    VectorXd U_prev       = VectorXd::Zero(sdof); // Displacement Vector used to calculate strain rate for damping force calculation
    VectorXd F_net        = VectorXd::Zero(sdof); // Total Nodal force vector
    VectorXd fe           = VectorXd::Zero(sdof); // External Nodal force vector
    VectorXd fe_prev      = VectorXd::Zero(sdof);

    VectorXd fr_prev      = VectorXd::Zero(sdof); // Reaction Nodal force vector at previous timestep
    VectorXd fr_curr      = VectorXd::Zero(sdof); // Reaction Nodal force vector at current timestep
    VectorXd fi_prev      = VectorXd::Zero(sdof); // Internal Nodal force vector at previous timestep
    VectorXd fi_curr      = VectorXd::Zero(sdof); // Internal Nodal force vector at current timestep
    VectorXd f_damp_curr  = VectorXd::Zero(sdof); // Linear Bulk Viscosity Damping Nodal force vector
    VectorXd f_damp_prev  = VectorXd::Zero(sdof); // Linear Bulk Viscosity Damping Nodal force vector at previous timestep

    double d_mesh_avg;                            // Average of d_avg values for entire host mesh
    double d_avg_tot;                             // Sum of d_tot values for entire host mesh
    VectorXd d_avg        = VectorXd::Zero(nel);  // Average d_tot of truss elements associated with each hex element

    // Following variables - Only for Truss Element
    int nel_truss = 1;                                   // number of truss elements

    if (embedded_constraint == 1) {
      nel_truss = mesh[1].getNumElements();              // number of truss elements
      fe_checkFiberVolumeFraction(mesh[0], mesh[1]);     // Check to confirm that maximum fiber volume fraction is not too large.
    }

    VectorXd d                  = VectorXd::Zero(nel_truss); // not currently used (i.e., d = 0)
    VectorXd d_fatigue          = VectorXd::Zero(nel_truss); // damage variable representing accumulated damage due to repeated loading
    VectorXd d_tot              = VectorXd::Zero(nel_truss); // damage variable representing total damage, d_tot = d + d_fatigue
    VectorXd lambda_min         = VectorXd::Ones(nel_truss); // minimum stretch experienced over lifetime of fiber
    VectorXd lambda_max         = VectorXd::Ones(nel_truss); // maximum stretch experienced over lifetime of fiber
    VectorXd lambda_min_cycle   = VectorXd::Ones(nel_truss); // minimum stretch experienced during current load cycle
    VectorXd lambda_max_cycle   = VectorXd::Ones(nel_truss); // maximum stretch experienced during current load cycle

    VectorXi n_load_cycle_full    = VectorXi::Zero(nel_truss);  // number of recognized full load cycles
    VectorXi n_load_cycle_partial = VectorXi::Zero(nel_truss);  // number of recognized partial load cycles

    bool import_damage = 0;                            // toggle indicating whether or not the user would like to import initial damage values
    bool include_d = 0;                                // toggle indicating whether or not the user would like to include damage in simulation

    if (embedded_constraint == 1) {
      import_damage = cons[0].get_EmbedImportDamage();
      include_d = cons[0].get_EmbedIncludeDamage();

      if (import_damage == 1) {
        double t_healing = cons[0].get_HealingTime();
        std::string damage_variables_import = home_path + "/" + "fiber_damage_input.txt";
        fe_damageVariableImport(damage_variables_import, d, d_fatigue, d_tot, lambda_min, lambda_max, t_healing);
      }

    }

    double energy_int_old = 0;
    double energy_int_new = 0;
    double energy_vd_old  = 0;
    double energy_vd_new  = 0;
    double energy_ext_old = 0;
    double energy_ext_new = 0;
    double energy_kin     = 0;
    double energy_total   = 0;
    double energy_max     = 0;

    double total_applied_fr = 0;

    std::string internal_energy = home_path + "/" + "results/internal_energy_system.txt";
    std::string external_energy = home_path + "/" + "results/external_energy_system.txt";
    std::string kinetic_energy = home_path + "/" + "results/kinetic_energy_system.txt";
    std::string total_energy = home_path + "/" + "results/total_energy_system.txt";
    std::string viscous_dissipation_energy = home_path + "/" + "results/viscous_dissipation_energy_system.txt";
    fe_energyWrite_new(internal_energy, viscous_dissipation_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_vd_new, energy_ext_new, energy_kin, energy_total);

    std::string reaction_forces = home_path + "/" + "results/total_reaction_force.txt";
    fe_singleDoubleWrite_new(reaction_forces, plot_state_counter, t, total_applied_fr);

    std::string damage_variables = home_path + "/" + "results/damage_variables_single_fiber.txt";
    fe_damageVariableWrite_new(damage_variables, plot_state_counter, t, d[0], d_fatigue[0], d_tot[0], lambda_min[0], lambda_max[0]);

    std::string average_damage = home_path + "/" + "results/average_damage_full_model.txt";
    fe_singleDoubleWrite_new(average_damage, plot_state_counter, t, d_mesh_avg);

    // Loading Conditions
    fe_apply_bc_load(fe, t_start);

    // ----------------------------------------------------------------------------
    // Step-1: Calculate the mass matrix similar to that of belytschko.
    VectorXd m_system = VectorXd::Zero(sdof);
    fe_calculateMass(m_system, "direct_lumped");

    std::string mass = home_path + "/" + "results/system_mass.txt";
    new_vector2text(mass.c_str(), m_system, m_system.cols());

    // ----------------------------------------------------------------------------
    // Step-2: getforce step from Belytschko
    fe_getforce(F_net, ndof, U, fe, time_step_counter, U_prev, dT, f_damp_curr, d, d_fatigue, d_tot, lambda_min, lambda_max, lambda_min_cycle, lambda_max_cycle, d_avg, n_load_cycle_full, n_load_cycle_partial, t, t_plot);

    mesh[0].readNodalKinematics(U, V, A);

    for (int i = 0; i < num_meshes; i++) {
        fe_vtuWrite(plot_state_counter - 1, t, mesh[i]);
    }

    dT = fe_getTimeStep();

    // ----------------------------------------------------------------------------
    // Step-3: Calculate accelerations
    fe_calculateAccln(A, m_system, F_net);
    U_prev = U;

    // ----------------------------------------------------------------------------
    // Step-4: Time loop starts....
    time_step_counter = time_step_counter + 1;
    clock_t s, s_prev, ds;
    s = clock();

    while (t < t_end) {

        if ((t + dT) >= t_end) {
            dT = t_end - t;
            if (dT <= 0) {
                break;
            }
        }

        /** Steps - 4,5,6 and 7 from Belytschko Box 6.1 - Update time, velocity and displacements */
        fe_timeUpdate(U, V, V_half, A, t, dT, "newmark-beta-central-difference");

        /** Update Loading Conditions - time dependent loading conditions */
        fe_apply_bc_load(fe, t);

        if (t >= (plot_state_counter * (output_temp_1)))
            t_plot = 1;

        else
          t_plot = 0;

        /** Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
        fe_getforce(F_net, ndof, U, fe, time_step_counter, U_prev, dT, f_damp_curr, d, d_fatigue, d_tot, lambda_min, lambda_max, lambda_min_cycle, lambda_max_cycle, d_avg, n_load_cycle_full, n_load_cycle_partial, t, t_plot); // Calculating the force term.

        /** Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
        fe_calculateAccln(A, m_system, F_net); // Calculating the new accelerations from total nodal forces.
        fe_apply_bc_acceleration(A, t);

        /** Step- 10 from Belytschko Box 6.1 - Second Partial Update of Nodal Velocities */
        fe_timeUpdate_velocity(V, V_half, A, t, dT, "newmark-beta-central-difference");

        fi_curr = fe - F_net;
        fe_calculateFR(fr_curr, fi_curr, m_system, A);

        /** Step - 11 from Belytschko Box 6.1 - Calculating energies and Checking Energy Balance */
        fe_checkEnergies(U_prev, U, fi_prev, fi_curr, f_damp_prev, f_damp_curr, fe_prev, fe, fr_prev, fr_curr, m_system, V, energy_int_old, energy_int_new, energy_vd_old, energy_vd_new, energy_ext_old, energy_ext_new, energy_kin, energy_total, energy_max);

        mesh[0].readNodalKinematics(U, V, A);

        if (t >= (plot_state_counter * (output_temp_1))) {

            for (int i = 0; i < num_meshes; i++) {
                fe_vtuWrite(plot_state_counter, t, mesh[i]);
            }

            plot_state_counter = plot_state_counter + 1;

            std::cout << "-----------------------------------" << "\n";
            std::cout << " Timestep Value = " << std::setw(10) << std::scientific << std::setprecision(5) << dT
                      << "\n Current Time = " << std::setw(5) << std::setprecision(1) << t
                      << "\n Timestep Number = " << (time_step_counter)
                      << "\n Plot State = " << (plot_state_counter - 1)
                      << "\n CPU Time = " << std::setw(5) << std::setprecision(1)
                      << ((float) ds / CLOCKS_PER_SEC) << "s" << "\n";
            // std::cout << std::setw(5) << std::scientific << std::setprecision(5) << "Z Displacement: " << U(20) << "\n";

            if (include_d == 1) {
              std::cout << " Number of Full Load Cycles = " << n_load_cycle_full[0] << "\n";
              std::cout << " Number of Partial Load Cycles = " << n_load_cycle_partial[0] << "\n";
            }

            /*print current frame, current time, and energy to individual .txt files*/
            fe_energyWrite_append(internal_energy, viscous_dissipation_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_vd_new, energy_ext_new, energy_kin, energy_total);

            total_applied_fr = fe_calculateTotalAppliedFR(fr_curr);

            fe_singleDoubleWrite_append(reaction_forces, plot_state_counter, t, total_applied_fr);

            d_avg_tot = 0;
            for (int i = 0; i < nel; i++) {
              d_avg_tot = d_avg_tot + d_avg(i);
            }
            d_mesh_avg = d_avg_tot / nel;

            fe_damageVariableWrite_append(damage_variables, plot_state_counter, t, d[0], d_fatigue[0], d_tot[0], lambda_min[0], lambda_max[0]);
            fe_singleDoubleWrite_append(average_damage, plot_state_counter, t, d_mesh_avg);

        }

        s_prev = s;
        s      = clock();
        ds     = s - s_prev;
        time_step_counter = time_step_counter + 1;

        dT = fe_getTimeStep();

    }

    if (embedded_constraint == 1) {
      std::string damage_variables_export = home_path + "/" + "results/fiber_damage_output.txt";
      fe_damageVariableExport(damage_variables_export, d, d_fatigue, d_tot, lambda_min, lambda_max);
    }

} // fe_mainEXPLICIT
