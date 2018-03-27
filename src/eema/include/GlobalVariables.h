/*! \mainpage
 *
 * This is the property of The Penn State Computational Biomechanics Group. \n
 * This code is developed by Harsha Teja Garimella under the supervision of Dr. Reuben H Kraft.
 *
 * ## Motivation:
 * Computational Brain Biomechanics
 *
 * ## Acknowledgements:
 * Funding from CFDRC and Funding from ARL
 *
 * ## Contact:
 * # Harsha T Garimella,
 * \n Ph.D. Candidate, Mechanical Engineering,
 * \n The Pennsylvania State University,
 * \n University Park, Pennsylvania, USA.
 * \n Email: harshatejagarimella@gmail.com
 *
 * # Jesse Gerber
 * \n MS Candidate, Mechanical Engineering,
 * \n The Pennsylvania State University,
 * \n University Park, Pennsylvania, USA.
 * \n Email: jig6@psu.edu
 *
 * # Reuben H. Kraft, Ph.D.
 * \n Shuman Asst. Professor,
 * \n Department of Mechanical Engineering,
 * \n Department of Biomedical Engineering,
 * \n The Pennsylvania State University,
 * \n University Park, Pennsylvania, USA.
 * \n Email: reuben.kraft@psu.edu
 */

#ifndef GLOBALVARIABLES_H_
#define GLOBALVARIABLES_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "unsupported/Eigen/MatrixFunctions"

#include "Mesh.h"
#include "Materials.h"
#include "BC.h"
#include "Constraint.h"

using namespace Eigen;

/*
 * Global Variables
 */

/** Job Folder */
extern std::string home_path;

/** Job Input File Name */
extern std::string job_file;

/** Dimension of the Simulation (3 for 3D, 2 for 2D and 1 for 1D) */
extern int ndof;

/** Number of Meshes in the simulation */
extern int num_meshes;

/** Pointer to the Mesh Data */
extern Mesh *mesh;

/** Start Time of the Simulation */
extern double t_start;

/** End Time of the Simulation */
extern double t_end;

/** Output Frequency - Number of output states desired */
extern int output_frequency;

/** Initial Time Step */
extern double dt_initial;

/** Max Time Step */
extern double dt_min;

/** Reduction Factor for the Time Step */
extern double reduction;

extern double failure_time_step;

/** Number of different material types in the simulation */
extern int material_types;

/** Pointer to the Material Data */
extern Materials *mat;

/** Types of Boundary Conditions */
extern int bc_types;

/** Pointer to the Boundary Condition Data */
extern BC *bc;

/** Number of constraints */
extern int num_constraints;
extern Constraint *cons;
extern bool embedded_constraint;
extern double area_truss;
// extern bool correct_volume_redundancy; // I do not think this is being used. Delete later on. -JG

/** Convergence Criteria for Newton Rhapson Method */
extern double eps_nr;

/** Maximum iterations for Newton Rhapson Method */
extern double iterations_nr;

/** Energy Balance Criterion */
extern double eps_energy;


// Variables added/moved

extern int counter_test, material_types_counter, matTypeHigh;
extern MatrixXd I;
extern int *matMap;

// fe_getForce_3d_normal

extern double *dndr_store, *dnds_store, *dndt_store, **x_store, **y_store, **z_store, x_normal, y_normal, z_normal, *wtx_normal, **wty_normal, ***wtz_normal, f_ext_e_sum_normal, ******jacobian_store, ******invJacobian_store, ****det_store, *****dndx_store, *****dndy_store, *****dndz_store;
extern int i_normal, j_normal, g_normal, nel_normal, nnel_normal, nnode_normal, sdof_normal, edof_normal, intx_normal, inty_normal, intz_normal;
extern VectorXd points_normal, weights_normal, dndx_normal, dndy_normal, dndz_normal, xcoord_normal, ycoord_normal, zcoord_normal, element_stress_host_local_normal, element_strain_host_local_normal, tmp_storage_normal, u_e_normal, u_e_prev_normal, f_ext_e_normal, pressure_e_normal, sigma_e_normal;
extern MatrixXd disp_mat_normal;

// fe_getPressure_lbv_pbr

extern int nnel_lbv, i_lbv;
extern double vol_strain_rate_lbv, volume_initial_lbv, volume_current_lbv, lc_lbv, c_wave_lbv, rho_initial_lbv, rho_current_lbv, pressure_scalar_lbv;
extern MatrixXd F_curr_lbv, F_inv_lbv, F_invT_lbv, F_prev_lbv, F_dot_lbv, F_dotT_lbv, D_lbv, pressure_matrix_lbv;

// fe_calDefGrad

extern MatrixXd H_DefGrad;

// fe_strDispMatrix_totalLagrangian_pbr

extern MatrixXd F_sdm, FT_sdm;
extern int i_sdm;

// fe_mooneyrivlin_hyperelastic_pbr

extern MatrixXd F_mrh, C_mrh, C_inv_mrh, C_square_mrh, C_bar_mrh, C_bar_inv_mrh, cauchy_sigma_mrh, pk_S_mrh;
extern double I1_mrh, tmp1_mrh, I2_mrh, defJacobian_mrh, defJacobian_frac_mrh, I1_bar_mrh, I2_bar_mrh, c1_mrh, c2_mrh, D_mrh, p_mrh;

// fe_calWaveSpeed

extern double E_cws, nu_cws, rho_cws, c_wave_cws;
extern std::string model_cws;

// fe_calArea_4

extern Vector3d tr1_side1, tr1_side2, tr2_side1, tr2_side2, area_tr1, area_tr2;
extern double area;

// fe_calVolume

extern Vector3d a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, e1, e2, e3;
extern double volume, tet1_vol, tet2_vol, tet3_vol, tet4_vol, tet5_vol;

#endif
