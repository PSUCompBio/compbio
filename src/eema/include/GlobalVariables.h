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

extern int counter_test; // Testing loops

extern MatrixXd I; // fe_calDefGrad

// fe_getForce_3d_normal

extern double ****dndr_store;
extern double ****dnds_store;
extern double ****dndt_store;

extern double detJacobian_normal;

extern MatrixXd disp_mat_normal;
extern VectorXd dndx_normal;
extern VectorXd dndy_normal;
extern VectorXd dndz_normal;
extern MatrixXd jacobian_normal;
extern MatrixXd invJacobian_normal;

extern VectorXd xcoord_normal;
extern VectorXd ycoord_normal;
extern VectorXd zcoord_normal;

extern VectorXd weights_normal;

extern int i_normal;
extern int j_normal;
extern int g_normal;
extern int nel_normal;
extern int nnel_normal;
extern int nnode_normal;
extern int sdof_normal;
extern int edof_normal;
extern int intx_normal;
extern int inty_normal;
extern int intz_normal;
extern double wtx_normal;
extern double wty_normal;
extern double wtz_normal;

// fe_getPressure_lbv_pbr

extern MatrixXd F_curr_lbv;
extern MatrixXd F_inv_lbv;
extern MatrixXd F_invT_lbv;
extern MatrixXd F_prev_lbv;
extern MatrixXd F_dot_lbv;
extern MatrixXd F_dotT_lbv;
extern MatrixXd D_lbv;
extern MatrixXd pressure_matrix_lbv;
extern double vol_strain_rate_lbv;
extern double volume_initial_lbv;
extern int nnel_lbv;
extern int i_lbv;
extern double volume_current_lbv;
extern double lc_lbv;
extern double c_wave_lbv;
extern double rho_initial_lbv;
extern double rho_current_lbv;
extern double pressure_scalar_lbv;

#endif
