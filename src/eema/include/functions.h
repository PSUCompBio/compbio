#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

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

#include "GlobalVariables.h"

using namespace Eigen;

/* Input-Output */

/** Reads a text file into matrix */
MatrixXd text2matrix(std::string name, int cols);
/** Writes a matrix into a text file */
void matrix2text(std::string name,MatrixXd new_slave_master,int width);
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
void fe_vtuWrite(std::string output, int time_step, Mesh mesh1);
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
/** dndr of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dndr_8(double rvalue, double svalue, double tvalue);
/** dnds of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dnds_8(double rvalue, double svalue, double tvalue);
/** dndt of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dndt_8(double rvalue, double svalue, double tvalue);
/** Create a guass_point vector of n values */
VectorXd guass_points(int n);
/** Creates a guass_weight vector of n values */
VectorXd guass_weights(int n);
/** Creates a guass point matrix in 3D */
MatrixXd guass_points_3d(int nx, int ny, int nz);
/** Creates a guass weight matrix in 3D */
MatrixXd guass_weights_3d(int ndof, int nx, int ny, int nz);
/** Calculates the jacobian -- using the derivates of shape functions */
MatrixXd fe_calJacobian(int dim, int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);
/** dndx of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndx_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);
/** dndy of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndy_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);
/** dndz of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndz_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);
/** Strain displacement matrix B */
MatrixXd fe_strDispMatrix(int edof, int nnel, VectorXd dndx, VectorXd dndy, VectorXd dndz);
/** Outputs the shape function matrix for an element */
MatrixXd fe_shapeMatrix(int edof, int nnel, VectorXd shapes);
/** Calculates the mass matrix for a hex element */
MatrixXd fe_mass_hex(MatrixXd nodes, VectorXi elements_row);
/** Calculates the mass of a truss element */
MatrixXd fe_mass_truss(double rho,double A_truss,int edof,MatrixXd nodes,MatrixXd elements);
/** Converts a normal mass matrix into a lumped mass matrix */
MatrixXd fe_transformMass(MatrixXd m, int opt);
/** Calculates the transformation matrix - transformation from local (truss) coordinate system to global (3d hex) coordinate system */
MatrixXd fe_calTransformation(MatrixXd truss_nodes, int choice);
MatrixXd fe_calSimpTransformation(MatrixXd truss_nodes);
/** Function outputs the minimum length for an element */
double fe_minElementLength(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);
/** Function outputs the maximum length for an element */
double fe_maxElementLength(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);
/** Calculates the strain displacement matrix in total lagrangian system */
MatrixXd fe_strDispMatrix_totalLagrangian(int edof, int nnel, VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u);
/** Calculates the deformation gradient */
MatrixXd fe_calDefGrad(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u);
/** Writes a symmetric tensor into a Voigt vector form */
VectorXd fe_tensor2voigt(MatrixXd A);
/** Writes a Voigt vector into a symmetric tensor form */
MatrixXd fe_voigt2tensor(VectorXd B);


/* =================================================================== */
/* Materials */
/* =================================================================== */

/** Create material matrix for isotropic elastic case */
MatrixXd fe_calculate_matlmat(int n, double E, double nu);
/** Updates the stress at each time step based on the material model */
VectorXd fe_stressUpdate(VectorXd dndx, VectorXd dndy, VectorXd dndz, MatrixXd disp_mat, VectorXd u, int opt, int return_opt);
/** Updates the stress at each time step based on the material model for a 1d element */
VectorXd fe_stressUpdate_1d(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u_e, int opt, MatrixXd nodes);
/** Extracts the material parameter values based on the material id */
double fe_get_mats(int matl_code, int obj_interest);
/** Extracts the material model name */
std::string fe_get_model(int matl_code);
/** Function calculates the stress vector for a mooney-rivlin material */
VectorXd fe_mooneyrivlin_hyperelastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);
/** Function calculates the stress vector for a ogden-hyperelastic material */
VectorXd fe_ogden_hyperelastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);
/** Function calculates the stress vector for a simple elastic (small strain) material */
VectorXd fe_simple_elastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, MatrixXd disp_mat, VectorXd u, int opt, int return_opt);
/** Function calculates the stress vector for a saint venant elastic material model */
VectorXd fe_saintvenant_elastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);


/* =================================================================== */
/* Math */
/* =================================================================== */
/** Calculates the area of a face with 4 vertices */
double fe_calArea_4(double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double c1, double c2, double c3, double c4);
/** Calculates the volume of a 3d element */
double fe_calVolume(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);
/** find the poistion index of a double value in a vector -- analogous to 'find' function in MATLAB */
int fe_find(VectorXd A, double a);
/** find the poistion index of a integer value in a vector -- analogous to 'find' function in MATLAB */
int fe_find(VectorXd A, int a);
/** Function calculates vector result using newton rhapson model */
VectorXd fe_newtonRhapson(VectorXd nat_coord, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);
/** Function outputs the standard curve values */
double fe_function(double a, std::string b, double time);
/** Function outputs the derivative of a time dependent variable at a time instant */
double fe_function_derivative(double a, std::string b, double time);
/** Function outputs the double derivative of a time dependent variable at a time instant */
double fe_function_d_derivative(double a, std::string b, double time);
/** Concatenate a vector to a matrix -- rowwise or coloumn wise */
MatrixXd fe_concatenate_vector2matrix(MatrixXd A, VectorXd B,int opt);
/** Function concatenates vector to a matrix */
MatrixXd fe_insert_vector2matrix(MatrixXd A, VectorXd B, int num, int opt);

/* =================================================================== */
/* FEM */
/* =================================================================== */
/** Calculates the resultant force vector - Box 6.1 of Belytschko */
VectorXd fe_getforce(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext, int size_counter, MatrixXd nodes_truss, MatrixXi elements_truss);
/** Find the index based on the DOF of a particular node */
VectorXi fe_find_index(VectorXi node_list);
/** Run the finite element analysis using an explicit dynamic method */
void fe_mainEXPLICIT();
/** Function enforces displacement boundary condition */
VectorXd fe_apply_bc_displacement(VectorXd U,double time);
/** Function enforces velocity boundary condition */
VectorXd fe_apply_bc_velocity(VectorXd V,double time);
/** Function enforces acceleration boundary condition */
VectorXd fe_apply_bc_acceleration(VectorXd A,double time);
/** Function updates the applied load */
VectorXd fe_apply_bc_load(VectorXd fe,double time);
/** Assembles the global mass matrix */
MatrixXd fe_assemble_mass(MatrixXd mm, MatrixXd m, VectorXi node_list, int sdof);
/** Creates element level vector (displacement, velocity, acceleration etc.) from a system level vector */
VectorXd fe_gather(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);
/** Updates a system level vector based on the element level vector */
VectorXd fe_scatter(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);
/** Updates the nodal coordinates based on the displacements */
MatrixXd fe_updateNodes(MatrixXd nodes,VectorXd displacements);

VectorXd fe_embed_preprocessing_mapping(Mesh host, Mesh embed);
VectorXd fe_embed_preprocessing(Mesh host,Mesh embed);
void fe_embed_preprocessing_length(Mesh host, Mesh embed);
int fe_compute_host(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host_tmp);
MatrixXd fe_create_bbox(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host, double length);

/* =================================================================== */
/* TimeStep */
/* =================================================================== */
/** Calculates the wavespeed for a particular material model */
double fe_calWaveSpeed(int material_id);
/** Outputs the critical time step based on all the elements in a FE analysis */
double fe_getTimeStep(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext);
/** Calculates the time step for a single element based on its dimensions and material model */
double fe_calTimeStep(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord, int material_id);


/* =================================================================== */
/* To be Removed */
/* =================================================================== */

/** Writes VTK files for host mesh */
void fe_vtkWrite_host(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements);
/** Writes VTK files for truss mesh */
void fe_vtkWrite_truss(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements);

/** Prints out a vector */
void fe_display_vector(VectorXd A);
/** Prints out a matrix */
void fe_display_matrix(MatrixXd A);

MatrixXd fe_apply_bc_stiffness(MatrixXd kk, VectorXi bcdof, VectorXd bcval);

/** Internal nodal force vector for a hexahedral element */
MatrixXd fe_stiffness_hex(double E, double nu, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

/** Internal nodal force vector for a truss (1D) element */
MatrixXd fe_stiffness_embed_truss(MatrixXd nodes_truss, MatrixXd elements_truss, double E_truss, double A_truss, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

void fe_writeElementStress(MatrixXd sigma_all, double time);

/* =================================================================== */
/* New Changes */
/* =================================================================== */
VectorXd fe_calculateMass(VectorXd mm, std::string type);
VectorXd fe_calculateMassDirectLumped(VectorXd mm);
VectorXd fe_massLumped(MatrixXd nodes, VectorXi elements_row);

VectorXd fe_calculateAccln(VectorXd mm, VectorXd F_net);

double fe_calculateKE(VectorXd mm, VectorXd V);
double fe_calculateKE(MatrixXd mm, VectorXd V);

VectorXd text2vector(std::string name);

#endif
