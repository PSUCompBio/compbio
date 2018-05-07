#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <omp.h>
#include <utmpx.h>
#include <thread>
#include <mutex>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "unsupported/Eigen/MatrixFunctions"

#include "GlobalVariables.h"
#include "Mesh.h"
#include "Materials.h"
#include "BC.h"
#include "Constraint.h"

using namespace Eigen;

/* Input-Output */

/** Reads a text file into matrix */
MatrixXd text2matrix(std::string name, int cols);

/** Writes a matrix into a text file */
void matrix2text(std::string name, MatrixXd new_slave_master, int width);

/** Writes a vector into a text file */
void new_vector2text(std::string name, VectorXd vector, int width);

/** Function appends a double value to a text file */
void append_double2text(std::string name, double a);

/** Function appends a double value with solution times to a text file */
void append_double2textWithTime(std::string name, int frame, double time, double value);

/** Function creates a new text file and writes a double value to it */
void new_double2text(std::string name, double a);

/** Function creates a new text file and writes a double value with solution times to it */
void new_double2textWithTime(std::string name, int frame, double time, double value);

/** Function writes the results into the VTU folder */
void fe_vtuWrite(int time_step, double time, Mesh& mesh1);

/** Function creates new PVD file in the VTU folder */
void fe_pvdNew(std::string output, int time_step, double time);

/** Function appends PVD file in the VTU folder */
void fe_pvdAppend(std::string output, int time_step, double time);

/** Read the input text file -- for a particular job */
void fe_mainRead(std::string file);

/* =================================================================== */
/* Elements */
/* =================================================================== */

/** Creates the shape functions for an 8 noded element */
VectorXd fe_shapes_8(double rvalue, double svalue, double tvalue);

/** dn of isoparametric element calculated for particular r, s, and t */
void fe_dniso_8(VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, double& r, double& s, double& t);
void fe_dniso_8_array(double *dndr, double *dnds, double *dndt, double& rvalue, double& svalue, double& tvalue, int intx, int inty, int intz);

/** Create a guass_point vector of n values */
VectorXd guass_points(int n);

/** Creates a guass_weight vector of n values */
VectorXd guass_weights(int n);

/** Creates a guass point matrix in 3D */
MatrixXd guass_points_3d(int nx, int ny, int nz);

/** Creates a guass weight matrix in 3D */
MatrixXd guass_weights_3d(int ndof, int nx, int ny, int nz);

/** Calculates the jacobian -- using the derivates of shape functions */
MatrixXd fe_calJacobian(int dim, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);
void fe_calJacobian_array(double **jacobian_store, int nnel, double* dndr, double* dnds, double* dndt, double* xcoord, double* ycoord, double* zcoord);

/** dndx, dndy, dndz of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
void fe_dndx_8_pbr(VectorXd& dndx, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian);
void fe_dndy_8_pbr(VectorXd& dndy, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian);
void fe_dndz_8_pbr(VectorXd& dndz, int nnel, VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, MatrixXd& invJacobian);
void fe_dndxyz_8_pbr_array(double* dndx, double* dndy, double* dndz, int nnel, double* dndr, double* dnds, double* dndt, double** invJacobian);

/** Strain displacement matrix B */
MatrixXd fe_strDispMatrix(int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz);

/** Outputs the shape function matrix for an element */
MatrixXd fe_shapeMatrix(int edof, int nnel, VectorXd shapes);

/** Calculates the mass matrix for a hex element */
MatrixXd fe_mass_hex(MatrixXd nodes, VectorXi elements_row);

/** Calculates the mass of a truss element */
MatrixXd fe_mass_truss(double rho, double A_truss, int edof, MatrixXd nodes, MatrixXd elements);

/** Converts a normal mass matrix into a lumped mass matrix */
MatrixXd fe_transformMass(MatrixXd m, int opt);

/** Calculates the transformation matrix - transformation from local (truss) coordinate system to global (3d hex) coordinate system */
MatrixXd fe_calTransformation(VectorXd x, VectorXd y, VectorXd z, int choice);
MatrixXd fe_calSimpTransformation(MatrixXd truss_nodes);

/** Function outputs the minimum length for an element */
double fe_minElementLength(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);

/** Function outputs the maximum length for an element */
double fe_maxElementLength(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);

/** Calculates the strain displacement matrix in total lagrangian system */
void fe_strDispMatrix_totalLagrangian_pbr(MatrixXd& disp_mat_modified, int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u);
void fe_strDispMatrix_totalLagrangian_pbr_array(MatrixXd& disp_mat_modified, int edof, int nnel, int i, int x, int y, int z, VectorXd& u, double* dndx, double* dndy, double* dndz);
void fe_strDispMatrix_totalLagrangian_fiber_pbr(MatrixXd& disp_mat_modified, int edof, int nnel, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, MatrixXd& F_fiber_ref);

/** Calculates the deformation gradient */
void fe_calDefGrad_pbr(MatrixXd& F, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u);
void fe_calDefGrad_pbr_array(MatrixXd& F, int i, int x, int y, int z, VectorXd& u);
void fe_calDefGrad_fiber_pbr(MatrixXd& F_fiber, MatrixXd& T_fiber, MatrixXd& T_fiber_inv, double lambda);

/** Writes a symmetric tensor into a Voigt vector form */
VectorXd fe_tensor2voigt(MatrixXd A);

/** Writes a Voigt vector into a symmetric tensor form */
MatrixXd fe_voigt2tensor(VectorXd B);

/* =================================================================== */
/* Materials */
/* =================================================================== */

/** Create material matrix for isotropic elastic case */
void fe_calculate_matlmat_pbr(MatrixXd& matl_mat, int n, double E, double nu);

/** Updates the stress at each time step based on the material model */
void fe_stressUpdate_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt);
void fe_stressUpdate_pbr_array(VectorXd& sigma, int i, int x, int y, int z, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt);
void fe_stressUpdate_fiber_pbr(VectorXd& fiber_stress, int opt, double lambda, MatrixXd& T_fiber, MatrixXd& T_fiber_inv, int return_opt);

/** Updates the stress at each time step based on the material model for a 1d element */
VectorXd fe_stressUpdate_1d(int material_id, VectorXd u_embed, VectorXd xcoord_embed, VectorXd ycoord_embed, VectorXd zcoord_embed, double length_old, int choice);
VectorXd fe_stressUpdate_1d(int material_id, VectorXd u_embed, VectorXd xcoord_embed, VectorXd ycoord_embed, VectorXd zcoord_embed, double length_old, VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u_e, int choice);

/** Extracts the material parameter values based on the material id */
double fe_get_mats(int matl_code, int obj_interest, std::string type);

/** Extracts the material model name */
std::string fe_get_model(int matl_code, std::string type);

/** Function calculates the stress vector for a mooney-rivlin material */
void fe_mooneyrivlin_hyperelastic_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt);
void fe_mooneyrivlin_hyperelastic_pbr_array(VectorXd& sigma, int i, int x, int y, int z, VectorXd& u, int opt, int return_opt);

/** Function calculates the stress vector for a ogden-hyperelastic material */
void fe_ogden_hyperelastic_pbr(VectorXd& sigma_local, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt);
void fe_ogden_hyperelastic_pbr_array(VectorXd& sigma_local, int i, int x, int y, int z, VectorXd& u, int opt, int return_opt);

/** Function calculates the stress vector for a simple elastic (small strain) material */
void fe_simple_elastic_pbr(VectorXd& sigma_local, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt);
void fe_simple_elastic_pbr_array(VectorXd& sigma_local, int i, int x, int y, int z, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt);

/** Function calculates the stress vector for a saint venant elastic material model */
void fe_saintvenant_elastic_pbr(VectorXd& sigma_local, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, int opt, int return_opt);
void fe_saintvenant_elastic_pbr_array(VectorXd& sigma_local, int i, int x, int y, int z, VectorXd& u, int opt, int return_opt);

/* =================================================================== */
/* Math */
/* =================================================================== */

/** Calculates the area of a face with 4 vertices */
double fe_calArea_4(double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double c1, double c2, double c3, double c4);

/** Calculates the volume of a 3d element */
double fe_calVolume(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);

/** Calculates current length of truss element */
double fe_calCurrLength_pbr(VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed);

/** find the poistion index of a double value in a vector -- analogous to 'find' function in MATLAB */
int fe_find(VectorXd& A, double a);

/** find the poistion index of a integer value in a vector -- analogous to 'find' function in MATLAB */
int fe_find(VectorXd& A, int a);

int fe_find(VectorXi& A, int a);

/** Function calculates vector result using newton rhapson model */
VectorXd fe_newtonRhapson(VectorXd& nat_coord, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);

/** Function outputs the standard curve values */

double
fe_function(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction);
/** Function outputs the derivative of a time dependent variable at a time instant */
double
fe_function_derivative(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction);
/** Function outputs the double derivative of a time dependent variable at a time instant */
double
fe_function_d_derivative(double a, std::string b, double time, int n_steps, VectorXd time_data, VectorXd x_data, VectorXd y_data, VectorXd z_data, int direction);

/** Concatenate a vector to a matrix -- rowwise or coloumn wise */
void fe_concatenate_vector2matrix(MatrixXd& A, VectorXd& B, int opt);
void fe_concatenate_vector2matrix(MatrixXi& A, VectorXi& B, int opt);

/** Function concatenates vector to a matrix */
void fe_insert_vector2matrix(MatrixXd& A, VectorXd& B, int num, int opt);
void fe_insert_vector2matrix(MatrixXi& A, VectorXi& B, int num, int opt);
void fe_insert_double2vector(VectorXd& A, double b);
void fe_insert_int2vector(VectorXi& A, int b);
void fe_insert_int2vector(VectorXi& A, int num, int b);

/** Linear Algebra functions */
double fe_detMatrix_pbr(MatrixXd& A);
double fe_detMatrix_pbr_array(double** A);
void fe_invMatrix_pbr(MatrixXd& A_inv, MatrixXd& A);
void fe_invMatrix_pbr_array(double** A_inv, double** A, double det);

/* =================================================================== */
/* FEM */
/* =================================================================== */

/** Calculates the resultant force vector - Box 6.1 of Belytschko */
void fe_getforce(VectorXd& f_tot, int ndof, VectorXd& u, VectorXd& fext, int time_step_counter, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXd& d_avg, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t, int t_plot);

/** Find the index based on the DOF of a particular node */
VectorXi fe_find_index(VectorXi node_list);

/** Run the finite element analysis using an explicit dynamic method */
void fe_mainEXPLICIT();

/** Function enforces displacement boundary condition */
void fe_apply_bc_displacement(VectorXd& U, double& time);

/** Function enforces velocity boundary condition */
void fe_apply_bc_velocity(VectorXd& V, double& time);

/** Function enforces acceleration boundary condition */
void fe_apply_bc_acceleration(VectorXd& A, double& time);

/** Function updates the applied load */
void fe_apply_bc_load(VectorXd& fe, double& time);

/** Assembles the global mass matrix */
MatrixXd fe_assemble_mass(MatrixXd mm, MatrixXd m, VectorXi node_list, int sdof);

/** Creates element level vector (displacement, velocity, acceleration etc.) from a system level vector */
VectorXd fe_gather(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);
void fe_gather_pbr(VectorXd& global_vec, VectorXd& local_vec, VectorXi node_list, int sdof);

/** Updates a system level vector based on the element level vector */
VectorXd fe_scatter(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);
void fe_scatter_pbr(VectorXd& global_vec, VectorXd& local_vec, VectorXi node_list, int sdof);

/** Updates the nodal coordinates based on the displacements */
MatrixXd fe_updateNodes(MatrixXd nodes, VectorXd displacements);

// void fe_embed_preprocessing_mapping(VectorXi& host_elements_embed_nodes, Mesh& host, Mesh& embed);

int fe_embed_node_host_check(int host_row, int embed_node, Mesh& host, Mesh& embed);

VectorXd fe_embed_element_intersection_point(int host_row, int embed_row, Mesh& host, Mesh& embed, int& intersection_check);

VectorXi fe_embed_preprocessing_host_map(VectorXi embed_map, Mesh& host, Mesh& embed);

VectorXi fe_embed_preprocessing(Mesh& host, Mesh& embed);

void fe_checkFiberVolumeFraction(Mesh& host, Mesh& embed);

/* =================================================================== */
/* TimeStep */
/* =================================================================== */

/** Calculates the wavespeed for a particular material model */
double fe_calWaveSpeed(int material_id, double volume_intial, double volume_current);

/** Outputs the critical time step based on all the elements in a FE analysis */
double fe_getTimeStep(void);

/** Calculates the time step for a single element based on its dimensions and material model */
double fe_calTimeStep(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id, VectorXd& u_e, double volume_initial);

/* =================================================================== */
/* To be Removed */
/* =================================================================== */

/** Writes VTK files for host mesh */
void fe_vtkWrite_host(std::string output, int format_choice, int mesh_choice, int time_step, MatrixXd nodes, MatrixXi elements);

/** Writes VTK files for truss mesh */
void fe_vtkWrite_truss(std::string output, int format_choice, int mesh_choice, int time_step, MatrixXd nodes, MatrixXi elements);

/** Prints out a vector */
void fe_display_vector(VectorXd A);

/** Prints out a matrix */
void fe_display_matrix(MatrixXd A);

/** Internal nodal force vector for a hexahedral element */
MatrixXd fe_stiffness_hex(double E, double nu, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

/** Internal nodal force vector for a truss (1D) element */
MatrixXd fe_stiffness_embed_truss(MatrixXd nodes_truss, MatrixXd elements_truss, double E_truss, double A_truss, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

void fe_writeElementStress(MatrixXd sigma_all, double time);

/* =================================================================== */
/* New Changes */
/* =================================================================== */

void initial_alloc();

void fe_calculateMass(VectorXd& m_system, std::string type);

void fe_calculateMassDirectLumped(VectorXd& m_system, int mesh_id);
void fe_calculateMassDirectLumped_embed(VectorXd& m_system, int host_id, int embed_id, bool address_vr, VectorXi& embed_map);

VectorXd fe_massLumped(MatrixXd* nodes, VectorXi elements_row);

void fe_calculateAccln(VectorXd& A, VectorXd& mm, VectorXd& F_net);

void fe_calculateKE(double& KE, VectorXd& mm, VectorXd& V);

double fe_calculateKE(MatrixXd mm, VectorXd V);

VectorXd text2vector(std::string name);

void fe_getForce_3d_normal(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, VectorXd& u_prev, double dT, VectorXd& f_damp, double t, int t_plot);

void fe_getForce_3d_embed(VectorXd& f_tot, VectorXd& u, VectorXd& fext, int time_step_counter, int host_id, int embed_id, bool address_vr, bool include_d, VectorXi& embed_map, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXd& d_avg, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t, int t_plot);

void fe_stressUpdateViscoelasticity_pbr(VectorXd& instantStress, double dT, MatrixXd& defGrad, MatrixXd invDefGrad, double defJacobian, int iterator, int ix, int iy, int iz, int opt, int return_opt);

void fe_getPressure_lbv_pbr(VectorXd& pressure, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id);
void fe_getPressure_lbv_pbr_array(VectorXd& pressure, int i, int x, int y, int z, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id);

void fe_stressModify(VectorXd& sigma_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, int choice);

void fe_staticDamageUpdate_pbr(VectorXd& d_static, int element_id, VectorXd& lambda_max);

void fe_fatigueDamageUpdate_pbr(int opt, VectorXd& d_fatigue, int fib, double lambda, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t);

double fe_calEquivalentStrain(int opt, double lambda);

void fe_stressUpdateDamage_pbr(VectorXd& fiber_stress_damaged, int opt, int fib, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double lambda, VectorXd& lambda_min, VectorXd& lambda_max, int return_opt);

double fe_calStrainEnergyDensity(int opt, double lambda);

void fe_calCentroidStress_3d_pbr(VectorXd& element_stress, double dT, int nnel, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, VectorXd& u_e, int material_id, int iterator, int ix, int iy, int iz);

void fe_calCentroidStrain_3d_pbr(VectorXd& element_strain, int nnel, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, VectorXd& u_e);

VectorXd fe_findIntgPoints_1d(VectorXd& x, VectorXd& y, VectorXd& z, double point, double length);

void fe_calCentroidStress_embed_3d_pbr(VectorXd& element_stress, int material_id, VectorXd& u_e, VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double length_old, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord);

void fe_calCentroidStrain_embed_3d_pbr(VectorXd& element_strain, VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, double length_old);

void fe_calculateFR(VectorXd& fr, VectorXd& fi_curr, VectorXd& mm, VectorXd& A);

double fe_calculateTotalAppliedFR(VectorXd& fr);

void fe_timeUpdate(VectorXd& U, VectorXd& V, VectorXd& V_half, VectorXd& A, double& t, double& dT, std::string time_integration_method);

void fe_timeUpdate_velocity(VectorXd& V, VectorXd& V_half, VectorXd& A, double& t, double& dT, std::string time_integration_method);

void fe_calculateEW(VectorXd& U_prev, VectorXd& U, VectorXd& fe_prev, VectorXd& fe, VectorXd& fr_prev, VectorXd& fr_curr, double& energy_ext_old, double& energy_ext_new);

void fe_calculateIE(VectorXd& U_prev, VectorXd& U, VectorXd& fi_prev, VectorXd& fi_curr, double& energy_int_old, double& energy_int_new);

void fe_calculateVDE(VectorXd& U_prev, VectorXd& U, VectorXd& f_damp_prev, VectorXd& f_damp_curr, double& energy_vd_old, double& energy_vd_new);

void fe_checkEnergies(VectorXd& U_prev, VectorXd& U, VectorXd& fi_prev, VectorXd& fi_curr, VectorXd& fe_damp_prev, VectorXd& fe_damp_curr, VectorXd& fe_prev, VectorXd& fe, VectorXd& fr_prev, VectorXd& fr_curr, VectorXd& m_system, VectorXd& V, double& energy_int_old, double& energy_int_new, double& energy_vd_old, double& energy_vd_new, double& energy_ext_old, double& energy_ext_new, double& energy_kin, double& energy_total, double& energy_max) ;

void fe_energyWrite_new(std::string& internal_energy, std::string& viscous_dissipation_energy, std::string& external_energy, std::string& kinetic_energy, std::string& total_energy, int plot_state_counter, double& t, double& energy_int_new, double& energy_vd_new, double& energy_ext_new, double& energy_kin, double& energy_total);

void fe_energyWrite_append(std::string& internal_energy, std::string& viscous_dissipation_energy, std::string& external_energy, std::string& kinetic_energy, std::string& total_energy, int plot_state_counter, double& t, double& energy_int_new, double& energy_vd_new, double& energy_ext_new, double& energy_kin, double& energy_total);

void fe_singleDoubleWrite_new(std::string& filename, int plot_state_counter, double& t, double& variable);

void fe_singleDoubleWrite_append(std::string& filename, int plot_state_counter, double& t, double& variable);

void fe_damageVariableWrite_new(std::string& damage_variables, int plot_state_counter, double& t, double& d_static_curr, double& d_fatigue_curr, double& d_tot, double& lambda_min, double& lambda_max);

void fe_damageVariableWrite_append(std::string& damage_variables, int plot_state_counter, double& t, double& d_static_curr, double& d_fatigue_curr, double& d_tot, double& lambda_min, double& lambda_max);

void fe_damageVariableExport(std::string& damage_variables_export, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max);

void fe_damageVariableImport(std::string& damage_variables_import, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, double t_healing);

/* =================================================================== */
/* BioElectroPhysics */
/* =================================================================== */

void fe_electrical_shapeMatrix(MatrixXd& electrical_shape_mat, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz);
void fe_apply_bc_current(VectorXd& I, double &time);
void fe_bioelectrophysics();
void fe_assemble_electricStiffness(MatrixXd& global, MatrixXd& local, VectorXi node_list);
void fe_scatter_electricalForce(VectorXd& global_vec, VectorXd& local_vec, VectorXi node_list);
void fe_apply_bc_potential(VectorXd& VP, double &time);
void fe_apply_bc_potential(MatrixXd& kk, VectorXd& ff, double time);
void fe_get_conductivity(MatrixXd& conductivity, int mat_id);
void fe_electroStatics_normal(double time, int mesh_id);
void fe_electroStatics_embed(double time, int host_id, int embed_id, VectorXi& embed_map);
void fe_electroStatics(double time);

void fe_form_face_matrix(MatrixXi& face_mat, int tmp, MatrixXi& elements);
void fe_find_common_face(VectorXi& common_face, bool com_face_value, int id_1, int id_2, MatrixXd& nodes_host, MatrixXi& elements_host);
void fe_find_intersection(VectorXd& intersection, VectorXi& common_face, MatrixXd& nodes_embed, MatrixXd& nodes_host, int i_1, int i_2);
Vector3d fe_find_normals(VectorXd& c1, VectorXd& c2, VectorXd& c3);
void fe_crossProduct(VectorXd& result, VectorXd& A, VectorXd& B);
VectorXd fe_sort_ascendingVectorXd(VectorXd A);
VectorXi fe_sort_ascendingVectorXi(VectorXi A);


#endif // ifndef FUNCTIONS_H_
