#include "functions.h"
using namespace Eigen;

double eps_energy = 0.01;
double area_truss = 7.85398e-7;
double failure_time_step = 1e-8;

/*! \brief
 * This function carries out the explicit dynamic analysis of the FEM problem.
 */

void
fe_mainEXPLICIT()
{
    // Following variables - Only for Hex Element
    int nnode = mesh[0].getNumNodes();          // number of nodes
    int sdof  = nnode * ndof;          // system degrees of freedom

    // Initialization
    double dT             = dt_initial;
    double t              = t_start;
    int time_step_counter = 0; // time step count
    int plot_state_counter = 1;
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

    // Following variables - Only for Truss Element
    int nel_truss = mesh[1].getNumElements();          // number of truss elements

    VectorXd d            = VectorXd::Zero(nel_truss); // damage variable representing damage due to single most severe stretch experienced
    VectorXd delta_d      = VectorXd::Zero(nel_truss); // damage variable representing accumulated damage due to repeated loading
    VectorXd lambda_min   = VectorXd::Zero(nel_truss); // minimum stretch experienced during current load cycle
    VectorXd lambda_max   = VectorXd::Zero(nel_truss); // maximum stretch experienced during current load cycle

    double energy_int_old = 0;
    double energy_int_new = 0;
    double energy_vd_old  = 0;
    double energy_vd_new  = 0;
    double energy_ext_old = 0;
    double energy_ext_new = 0;
    double energy_kin     = 0;
    double energy_total   = 0;
    double energy_max     = 0;

    std::string internal_energy = home_path + "/" + "results/internal_energy_system.txt";
    std::string external_energy = home_path + "/" + "results/external_energy_system.txt";
    std::string kinetic_energy = home_path + "/" + "results/kinetic_energy_system.txt";
    std::string total_energy = home_path + "/" + "results/total_energy_system.txt";
    std::string viscous_dissipation_energy = home_path + "/" + "results/viscous_dissipation_energy_system.txt";
    fe_energyWrite_new(internal_energy, viscous_dissipation_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_vd_new, energy_ext_new, energy_kin, energy_total);

    std::string reaction_forces = home_path + "/" + "results/reaction_forces.txt";
    fe_reactionForceWrite_new(reaction_forces, plot_state_counter, t, fr_curr[5], fr_curr[8], fr_curr[17], fr_curr[20]);

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
    fe_getforce(F_net, ndof, U, fe, time_step_counter, U_prev, dT, f_damp_curr, d, delta_d, lambda_min, lambda_max);

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

        /** Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
        fe_getforce(F_net, ndof, U, fe, time_step_counter, U_prev, dT, f_damp_curr, d, delta_d, lambda_min, lambda_max); // Calculating the force term.

        /** Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
        fe_calculateAccln(A, m_system, F_net); // Calculating the new accelerations from total nodal forces.
        fe_apply_bc_acceleration(A, t);

        /** Step- 10 from Belytschko Box 6.1 - Second Partial Update of Nodal Velocities */
        fe_timeUpdate_velocity(V, V_half, A, t, dT, "newmark-beta-central-difference");

        fi_curr = fe - F_net;
        fe_calculateFR(fr_curr, sdof, fi_curr, m_system, A);

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
                      << ((float) ds / CLOCKS_PER_SEC) << "s \n";
            std::cout << std::setw(5) << std::scientific << std::setprecision(5) << "Z Displacement: " << U(20) << "\n";

            /*print current frame, current time, and energy to individual .txt files*/
            fe_energyWrite_append(internal_energy, viscous_dissipation_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_vd_new, energy_ext_new, energy_kin, energy_total);

            fe_reactionForceWrite_append(reaction_forces, plot_state_counter, t, fr_curr[5], fr_curr[8], fr_curr[17], fr_curr[20]);

            std::cout << "d = " << '\n' << d << '\n';
            std::cout << "delta_d = " << '\n' << delta_d << '\n';
            // std::exit(1);

        }

        s_prev = s;
        s      = clock();
        ds     = s - s_prev;
        time_step_counter = time_step_counter + 1;

        dT = fe_getTimeStep();

    }
} // fe_mainEXPLICIT
