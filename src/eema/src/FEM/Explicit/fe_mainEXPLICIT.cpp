#include "functions.h"
using namespace Eigen;

double eps_energy = 0.01;
double failure_time_step = 1e-8;

// set include_viscoelasticity = 0 to ignore viscoelasticity in the hex elements
// set include_viscoelasticity = 1 to include viscoelasticity in the hex elements
bool include_viscoelasticity = 0;

// Global

int counter_test = 0, *matMap;
MatrixXd I;


// fe_getForce_3d_normal

double *****dndr_store, *****dnds_store, *****dndt_store, **x_store, **y_store, **z_store, x_normal, y_normal, z_normal, *wtx_normal, **wty_normal, ***wtz_normal, f_ext_e_sum_normal, defJacobian_normal, ******jacobian_store, ******invJacobian_store, ****det_store, *****dndx_store, *****dndy_store, *****dndz_store, ******internalStressVariable1_prev_normal_store, ******internalStressVariable2_prev_normal_store, ******devInstantPK2Stress_prev_normal_store, ***internalStressVariable1_prev_centroid_store, ***internalStressVariable2_prev_centroid_store, ***devInstantPK2Stress_prev_centroid_store;
int i_normal, j_normal, g_normal, nel_normal, nnel_normal, nnode_normal, sdof_normal, edof_normal, intx_normal, inty_normal, intz_normal;
VectorXd points_normal, weights_normal, dndx_normal, dndy_normal, dndz_normal, xcoord_normal, ycoord_normal, zcoord_normal, element_stress_host_local_normal, element_strain_host_local_normal, tmp_storage_normal, u_e_normal, u_e_prev_normal, f_ext_e_normal, pressure_e_normal, sigma_e_normal, *element_characteristic_embed_normal, element_stress_embed_local_normal, element_strain_embed_local_normal;
MatrixXd disp_mat_normal, defGrad_normal, invDefGrad_normal, *nodes_host_normal, *nodes_embed_normal;
MatrixXi *elements_host_normal, *elements_embed_normal;

// fe_getPressure_lbv_pbr

int nnel_lbv, i_lbv;
double vol_strain_rate_lbv, volume_initial_lbv, volume_current_lbv, lc_lbv, c_wave_lbv, rho_initial_lbv, rho_current_lbv, pressure_scalar_lbv;
MatrixXd F_curr_lbv, F_inv_lbv, F_invT_lbv, F_prev_lbv, F_dot_lbv, F_dotT_lbv, D_lbv, pressure_matrix_lbv;

// fe_calDefGrad

MatrixXd H_DefGrad;

// fe_strDispMatrix_totalLagrangian_pbr

MatrixXd F_sdm, FT_sdm;
int i_sdm;

// fe_mooneyrivlin_hyperelastic_pbr

MatrixXd F_mrh, C_mrh, C_inv_mrh, C_square_mrh, C_bar_mrh, C_bar_inv_mrh, cauchy_sigma_mrh, pk_S_mrh;
double I1_mrh, tmp1_mrh, I2_mrh, defJacobian_mrh, defJacobian_frac_mrh, I1_bar_mrh, I2_bar_mrh, c1_mrh, c2_mrh, D_mrh, p_mrh;

// fe_calWaveSpeed

double E_cws, nu_cws, rho_cws, c_wave_cws;
std::string model_cws;

// fe_calArea_4

Vector3d tr1_side1, tr1_side2, tr2_side1, tr2_side2, area_tr1, area_tr2;
double area;

// fe_calVolume

Vector3d a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, e1, e2, e3;
double volume, tet1_vol, tet2_vol, tet3_vol, tet4_vol, tet5_vol;

void experimental() {

    int i, j, k ,l, m, n;
    nodes_host_normal = mesh[0].getNewNodesPointer();
    elements_host_normal = mesh[0].getNewElementsPointer();
    nel_normal   = mesh[0].getNumElements();
    nnel_normal  = mesh[0].getNumNodesPerElement();
    nnode_normal = mesh[0].getNumNodes();
    sdof_normal  = nnode_normal * ndof;
    edof_normal  = nnel_normal * ndof;
    if (embedded_constraint == 1) {
      nodes_embed_normal    = mesh[1].getNewNodesPointer();
      elements_embed_normal = mesh[1].getNewElementsPointer();
      element_characteristic_embed_normal = mesh[1].getElementCharacteristicPointer();
    }

    // Allocating Memory

    I = MatrixXd::Identity(3, 3);
    points_normal  = guass_points(2);
    weights_normal = guass_weights(2);
    disp_mat_normal = MatrixXd::Zero(6, edof_normal);
    defGrad_normal = MatrixXd::Identity(ndof, ndof);
    invDefGrad_normal = MatrixXd::Identity(ndof, ndof);
    dndx_normal = VectorXd::Zero(nnel_normal);
    dndy_normal = VectorXd::Zero(nnel_normal);
    dndz_normal = VectorXd::Zero(nnel_normal);
    xcoord_normal = VectorXd::Zero(nnel_normal);
    ycoord_normal = VectorXd::Zero(nnel_normal);
    zcoord_normal = VectorXd::Zero(nnel_normal);
    element_stress_host_local_normal = VectorXd::Zero(nel_normal * 9);
    element_strain_host_local_normal = VectorXd::Zero(nel_normal * 9);
    tmp_storage_normal = VectorXd::Zero(ndof * ndof);
    u_e_normal = VectorXd::Zero(edof_normal);
    u_e_prev_normal = VectorXd::Zero(edof_normal);
    f_ext_e_normal = VectorXd::Zero(edof_normal);
    pressure_e_normal = VectorXd::Zero(6);
    sigma_e_normal = VectorXd::Zero(6);
    F_curr_lbv = MatrixXd::Zero(ndof, ndof);
    F_inv_lbv = MatrixXd::Zero(ndof, ndof);
    F_invT_lbv = MatrixXd::Zero(ndof, ndof);
    F_prev_lbv = MatrixXd::Zero(ndof, ndof);
    F_dot_lbv = MatrixXd::Zero(ndof, ndof);
    F_dotT_lbv = MatrixXd::Zero(ndof, ndof);
    D_lbv = MatrixXd::Zero(ndof, ndof);
    pressure_matrix_lbv = MatrixXd::Zero(ndof, ndof);
    H_DefGrad = MatrixXd::Zero(ndof, ndof);
    F_sdm = MatrixXd::Zero(3, 3);
    FT_sdm = MatrixXd::Zero(3, 3);
    F_mrh = MatrixXd::Zero(ndof, ndof);
    C_mrh = MatrixXd::Zero(ndof, ndof);
    C_inv_mrh = MatrixXd::Zero(ndof, ndof);
    C_square_mrh = MatrixXd::Zero(ndof, ndof);
    C_bar_mrh = MatrixXd::Zero(ndof, ndof);
    C_bar_inv_mrh = MatrixXd::Zero(ndof, ndof);
    cauchy_sigma_mrh = MatrixXd::Zero(ndof, ndof);
    pk_S_mrh = MatrixXd::Zero(ndof, ndof);
    tr1_side1 = Vector3d(3);
    tr1_side2 = Vector3d(3);
    tr2_side1 = Vector3d(3);
    tr2_side2 = Vector3d(3);
    area_tr1 = Vector3d(3);
    area_tr2 = Vector3d(3);

    if (embedded_constraint == 1) {
      element_stress_embed_local_normal = VectorXd::Zero((*elements_embed_normal).rows() * 9);
      element_strain_embed_local_normal = VectorXd::Zero((*elements_embed_normal).rows() * 9);
    }

    if (material_types_counter >= matTypeHigh)
        matMap = new int[material_types_counter + 1];
    else
        matMap = new int[matTypeHigh + 1];

    for (i_lbv = 0; i_lbv < material_types_counter; i_lbv++) {
        i_normal = mat[i_lbv].getMatID();
        matMap[i_normal] = i_lbv;
    }

    i_lbv = 0;
    i_normal = 0;

    dndr_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndr_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndr_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndr_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndr_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dnds_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dnds_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dnds_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dnds_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dnds_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndt_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndt_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndt_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndt_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndt_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    x_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        x_store[i] = new double[nnel_normal];
    }

    y_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        y_store[i] = new double[nnel_normal];
    }

    z_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        z_store[i] = new double[nnel_normal];
    }

    wtx_normal = new double[2];

    wty_normal = new double*[2];
    for (i = 0; i < 2; i++)
        wty_normal[i] = new double[2];

    wtz_normal = new double**[2];
    for (i = 0; i < 2; i++) {
        wtz_normal[i] = new double*[2];
        for (j = 0; j < 2; j++) {
            wtz_normal[i][j] = new double[2];
        }
    }

    jacobian_store = new double*****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        jacobian_store[i] = new double****[2];
        for (j = 0; j < 2; j++) {
            jacobian_store[i][j] = new double***[2];
            for (k = 0; k < 2; k++) {
                jacobian_store[i][j][k] = new double**[2];
                for (l = 0; l < 2; l++) {
                    jacobian_store[i][j][k][l] = new double*[ndof];
                    for (m = 0; m < ndof; m++) {
                        jacobian_store[i][j][k][l][m] = new double[ndof];
                    }
                }
            }
        }
    }

    invJacobian_store = new double*****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        invJacobian_store[i] = new double****[2];
        for (j = 0; j < 2; j++) {
            invJacobian_store[i][j] = new double***[2];
            for (k = 0; k < 2; k++) {
                invJacobian_store[i][j][k] = new double**[2];
                for (l = 0; l < 2; l++) {
                    invJacobian_store[i][j][k][l] = new double*[ndof];
                    for (m = 0; m < ndof; m++) {
                        invJacobian_store[i][j][k][l][m] = new double[ndof];
                    }
                }
            }
        }
    }

    det_store = new double***[nel_normal];
    for (j = 0; j < nel_normal; j++) {
        det_store[j] = new double**[2];
        for (k = 0; k < 2; k++) {
            det_store[j][k] = new double*[2];
            for (l = 0; l < 2; l++) {
                det_store[j][k][l] = new double[2];
            }
        }
    }

    dndx_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndx_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndx_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndx_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndx_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndy_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndy_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndy_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndy_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndy_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndz_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndz_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndz_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndz_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndz_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

  internalStressVariable1_prev_normal_store = new double*****[nel_normal];
  internalStressVariable2_prev_normal_store = new double*****[nel_normal];
  devInstantPK2Stress_prev_normal_store = new double*****[nel_normal];
  for (i = 0; i < nel_normal; i++) {
      internalStressVariable1_prev_normal_store[i] = new double****[2];
      internalStressVariable2_prev_normal_store[i] = new double****[2];
      devInstantPK2Stress_prev_normal_store[i] = new double****[2];
      for (j = 0; j < 2; j++) {
          internalStressVariable1_prev_normal_store[i][j] = new double***[2];
          internalStressVariable2_prev_normal_store[i][j] = new double***[2];
          devInstantPK2Stress_prev_normal_store[i][j] = new double***[2];
          for (k = 0; k < 2; k++) {
              internalStressVariable1_prev_normal_store[i][j][k] = new double**[2];
              internalStressVariable2_prev_normal_store[i][j][k] = new double**[2];
              devInstantPK2Stress_prev_normal_store[i][j][k] = new double**[2];
              for (l = 0; l < 2; l++) {
                  internalStressVariable1_prev_normal_store[i][j][k][l] = new double*[ndof];
                  internalStressVariable2_prev_normal_store[i][j][k][l] = new double*[ndof];
                  devInstantPK2Stress_prev_normal_store[i][j][k][l] = new double*[ndof];
                  for (m = 0; m < ndof; m++) {
                      internalStressVariable1_prev_normal_store[i][j][k][l][m] = new double[ndof];
                      internalStressVariable2_prev_normal_store[i][j][k][l][m] = new double[ndof];
                      devInstantPK2Stress_prev_normal_store[i][j][k][l][m] = new double[ndof];
                  }
              }
          }
      }
  }
  for (i = 0; i < nel_normal; i++) {
      for (j = 0; j < 2; j++) {
          for (k = 0; k < 2; k++) {
              for (l = 0; l < 2; l++) {
                  for (m = 0; m < ndof; m++) {
                    for (n = 0; n < ndof; n++) {
                        internalStressVariable1_prev_normal_store[i][j][k][l][m][n] = 0.0;
                        internalStressVariable2_prev_normal_store[i][j][k][l][m][n] = 0.0;
                        devInstantPK2Stress_prev_normal_store[i][j][k][l][m][n] = 0.0;
                    }
                  }
              }
          }
      }
  }

  internalStressVariable1_prev_centroid_store = new double**[nel_normal];
  internalStressVariable2_prev_centroid_store = new double**[nel_normal];
  devInstantPK2Stress_prev_centroid_store = new double**[nel_normal];
  for (i = 0; i < nel_normal; i++) {
      internalStressVariable1_prev_centroid_store[i] = new double*[ndof];
      internalStressVariable2_prev_centroid_store[i] = new double*[ndof];
      devInstantPK2Stress_prev_centroid_store[i] = new double*[ndof];
      for (j = 0; j < ndof; j++) {
          internalStressVariable1_prev_centroid_store[i][j] = new double[ndof];
          internalStressVariable2_prev_centroid_store[i][j] = new double[ndof];
          devInstantPK2Stress_prev_centroid_store[i][j] = new double[ndof];
      }
  }
  for (i = 0; i < nel_normal; i++) {
      for (j = 0; j < ndof; j++) {
        for (k = 0; k < ndof; k++) {
            internalStressVariable1_prev_centroid_store[i][j][k] = 0.0;
            internalStressVariable2_prev_centroid_store[i][j][k] = 0.0;
            devInstantPK2Stress_prev_centroid_store[i][j][k] = 0.0;
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

    // Initialization
    double dT             = dt_initial;
    double t              = t_start;
    int time_step_counter = 0; // time step count
    int plot_state_counter = 1;
    int t_plot = 1;

    double output_temp_1  = ((double) (t_end / output_frequency));
    VectorXd A            = VectorXd::Zero(sdof_normal); // Acceleration Vector
    VectorXd V            = VectorXd::Zero(sdof_normal); // Velocity Vector
    VectorXd V_half       = VectorXd::Zero(sdof_normal); // Velocity Vector at n+1/2
    VectorXd U            = VectorXd::Zero(sdof_normal); // Displacement Vector
    VectorXd U_prev       = VectorXd::Zero(sdof_normal); // Displacement Vector used to calculate strain rate for damping force calculation
    VectorXd F_net        = VectorXd::Zero(sdof_normal); // Total Nodal force vector
    VectorXd fe           = VectorXd::Zero(sdof_normal); // External Nodal force vector
    VectorXd fe_prev      = VectorXd::Zero(sdof_normal);

    VectorXd fr_prev      = VectorXd::Zero(sdof_normal); // Reaction Nodal force vector at previous timestep
    VectorXd fr_curr      = VectorXd::Zero(sdof_normal); // Reaction Nodal force vector at current timestep
    VectorXd fi_prev      = VectorXd::Zero(sdof_normal); // Internal Nodal force vector at previous timestep
    VectorXd fi_curr      = VectorXd::Zero(sdof_normal); // Internal Nodal force vector at current timestep
    VectorXd f_damp_curr  = VectorXd::Zero(sdof_normal); // Linear Bulk Viscosity Damping Nodal force vector
    VectorXd f_damp_prev  = VectorXd::Zero(sdof_normal); // Linear Bulk Viscosity Damping Nodal force vector at previous timestep

    double d_mesh_avg;                            // Average of d_avg values for entire host mesh
    double d_avg_tot;                             // Sum of d_tot values for entire host mesh
    VectorXd d_avg        = VectorXd::Zero(nel_normal);  // Average d_tot of truss elements associated with each hex element

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
    VectorXd m_system = VectorXd::Zero(sdof_normal);
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
                      << ((double) ds / CLOCKS_PER_SEC) << "s" << "\n";
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
            for (int i = 0; i < nel_normal; i++) {
              d_avg_tot = d_avg_tot + d_avg(i);
            }
            d_mesh_avg = d_avg_tot / nel_normal;

            fe_damageVariableWrite_append(damage_variables, plot_state_counter, t, d[0], d_fatigue[0], d_tot[0], lambda_min[0], lambda_max[0]);
            fe_singleDoubleWrite_append(average_damage, plot_state_counter, t, d_mesh_avg);

        }

        s_prev = s;
        s      = clock();
        ds     = s - s_prev;
        time_step_counter = time_step_counter + 1;

        dT = fe_getTimeStep();
//dT = 1e-6;
    }

    if (embedded_constraint == 1) {
      std::string damage_variables_export = home_path + "/" + "results/fiber_damage_output.txt";
      fe_damageVariableExport(damage_variables_export, d, d_fatigue, d_tot, lambda_min, lambda_max);
    }
    std::cout << "\n Counter = " << counter_test << "\n";
} // fe_mainEXPLICIT
