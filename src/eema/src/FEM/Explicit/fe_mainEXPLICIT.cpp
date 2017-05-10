#include "functions.h"
using namespace Eigen;

double eps_energy = 0.01;

/*! \brief
 * This function carries out the explicit dynamic analysis of the FEM problem.
 */

void fe_mainEXPLICIT(){

	for(int i=0;i<num_meshes;i++){
		mesh[i].preprocessMesh();
	}

	MatrixXd nodes = mesh[0].getNewNodes();
	MatrixXi elements = mesh[0].getNewElements();
	MatrixXd nodes_truss = mesh[1].getNewNodes();
	MatrixXi elements_truss = mesh[1].getNewElements();

	// Following variables - Only for Hex Element
	int nel = elements.rows(); /*! number of elements */
	int nnel = (elements.cols()-2); // number of nodes per element
	int nnode = nodes.rows(); // number of nodes
	int sdof = nnode * ndof; // system degrees of freedom
	int edof = nnel * ndof; // element degrees of freedom

	// Updated Nodes and Elements
	MatrixXd updated_nodes = nodes;

	// Writing the Read Mesh
	//fe_vtkWrite_host("eem_matrix",1,5,0,updated_nodes,elements);
	//return;

	// Initialization
	double dT = dt_initial;
	double t_half = 0; // t(n + 1/2)
	double t=t_start;
	int size_counter = 0; // time step count
	int time_temp_1 = 1;
	double output_temp_1 = ((double)(t_end/output_frequency));
	VectorXd A = VectorXd::Zero(sdof); // Acceleration Vector
	VectorXd V = VectorXd::Zero(sdof); // Velocity Vector
	VectorXd V_half = VectorXd::Zero(sdof); // Velocity Vector at n+1/2
	VectorXd U = VectorXd::Zero(sdof); // Displacement Vector
	VectorXd F_net = VectorXd::Zero(sdof); // Total Nodal force vector
	VectorXd fe = VectorXd::Zero(sdof); // External Nodal force vector
	VectorXd fe_prev = VectorXd::Zero(sdof);
	VectorXd fe_curr = VectorXd::Zero(sdof);
	VectorXd fi_prev = VectorXd::Zero(sdof); // Internal nodal force vector at previous timestep
	VectorXd fi_curr = VectorXd::Zero(sdof); // Internal Nodal force vector at current timestep
	VectorXd U_prev = VectorXd::Zero(sdof); // Nodal displacements at previous time
	VectorXd U_curr = VectorXd::Zero(sdof); // Nodal displacements at current time
	double energy_int_old = 0;
	double energy_int_new = 0;
	double energy_ext_old = 0;
	double energy_ext_new = 0;
	double energy_kin = 0;
	double energy_total = 0;
	double energy_max = 0;

	std::string internal_energy = home_path + "/" + "results/internal_energy_system.txt";
	new_double2textWithTime(internal_energy,time_temp_1 - 1,t,energy_int_new);
	std::string external_energy = home_path + "/" +"results/external_energy_system.txt";
	new_double2textWithTime(external_energy,time_temp_1 - 1,t,energy_ext_new);
	std::string kinetic_energy = home_path + "/" + "results/kinetic_energy_system.txt";
	new_double2textWithTime(kinetic_energy,time_temp_1 - 1,t,energy_kin);
	std::string total_energy = home_path + "/" +"results/total_energy_system.txt";
	new_double2textWithTime(total_energy,time_temp_1 - 1,t,energy_total);


	// Loading Conditions
	fe = fe_apply_bc_load(fe,0);

	// ----------------------------------------------------------------------------
	//Step-1: Calculate the mass matrix similar to that of belytschko.

	VectorXd mm = VectorXd::Zero(sdof);
	mm = fe_calculateMass(mm,"direct_lumped");

	std::string mass = home_path+"/"+"results/system_mass.txt";
	new_vector2text(mass.c_str(),mm,mm.cols());

	// ----------------------------------------------------------------------------
 //Step-2: getforce step from Belytschko
	F_net = fe_getforce(nodes,elements,ndof,U,V,fe,size_counter,nodes_truss,elements_truss);

	mesh[0].readNodalKinematics(U,V,A);
	fe_vtuWrite("eem_matrix",size_counter,mesh[0]);
	fe_pvdNew("eem_matrix",size_counter,t);

	dT = reduction * fe_getTimeStep(nodes,elements,ndof,U,V,fe);
	if(dT>dt_min){
		dT = dt_min;
	}

	// ----------------------------------------------------------------------------

	//Step-3: Calculate accelerations
	A  = fe_calculateAccln(mm,F_net);
	U_prev = U;

	// ----------------------------------------------------------------------------

	//Step-4: Time loop starts....
	size_counter = size_counter+1;
	clock_t s,s_prev,ds;
	s = clock();

	   while(t<=t_end){

				/** Apply Loading Conditions - time dependent loading conditions */
				fe = fe_apply_bc_load(fe,t);

				/** Calculate the time at half time step */
				t_half = 0.5*(t+t+dT);

				/** Update the time by adding full time step */
				t = t+dT;

				/** Partially Update Nodal Velocities */
				V_half = V + (t_half - (t-dT))*A;

				/** Enforce Velocity Boundary Conditions */
				//V_half = fe_apply_bc_velocity(V_half,t_half);
				//std::cout << "V_half: \n" << V_half << '\n';

				/** Update Nodal Displacements */
				U = U + dT*(V_half);
				U = fe_apply_bc_displacement(U,t);

				F_net = fe_getforce(nodes,elements,ndof,U,V,fe,size_counter,nodes_truss,elements_truss); // Calculating the force term.

				dT = reduction * fe_getTimeStep(nodes, elements, ndof, U, V, fe);
				if(dT>dt_min){
					dT = dt_min;
				}
				if( ((t+dT) > t_end) && (t!=t_end)){
					dT = t_end-t;
				}

				/** Calculate Accelerations */
				A = fe_calculateAccln(mm,F_net); // Calculating the new accelerations from total nodal forces.
				A = fe_apply_bc_acceleration(A,t);

				/** Completely Update the nodal velocities */
				V = V_half + (t - t_half)*A; // Calculating the new velocities.
				V = fe_apply_bc_velocity(V,t);

				/* Calculating the internal energy terms */
				U_curr = U;
				VectorXd del_U = U_curr - U_prev;
				fi_curr = fe - F_net;
				energy_int_new = energy_int_old + 0.5*(del_U.dot(fi_prev + fi_curr));
				fi_prev = fi_curr;
				fi_curr = VectorXd::Zero(sdof);
				U_prev = U_curr;
				energy_int_old = energy_int_new;


				/* Calculating the external energy terms */
				fe_curr = fe ;
				energy_ext_new = energy_ext_old + 0.5*(del_U.dot(fe_prev + fe_curr));
				fe_prev = fe_curr;
				energy_ext_old = energy_ext_new;

				/* Calculating the kinetic energy */
				energy_kin = fe_calculateKE(mm,V);

				/* Calculating the total energy of the system */
				energy_total = std::abs(energy_kin + energy_int_new - energy_ext_new);
				energy_max = std::max(std::max(energy_kin,energy_int_new),energy_ext_new);

				if(energy_total > (eps_energy*(energy_max))){
					std::cout << "**********************************************" << std::endl;
					std::cout << "ALERT: INSTABILITIES IN THE SYSTEM DETECTED \n BASED ON THE ENERGY BALANCE CHECK \n";
					std::cout << "**********************************************" << std::endl;
				}

				/** Projection of displacements to the embedded mesh is needed */
				/* */
				/* */

				if(t >= (time_temp_1 * (output_temp_1))){

					mesh[0].readNodalKinematics(U,V,A);

					fe_vtuWrite("eem_matrix",size_counter,mesh[0]);
					fe_pvdAppend("eem_matrix",size_counter,t);

					time_temp_1 = time_temp_1 + 1;
					//fe_vtkWrite_host("eem_matrix",1,5,size_counter,nodes,elements);
					//fe_vtkWrite_truss("eem_truss",1,5,size_counter,nodes_truss,elements_truss);
	        std::cout <<"Timestep Value = "<<std::setw(5)<<std::scientific<<std::setprecision(1)<<dT
					  								<<"  Current Time = "<<std::setw(5)<<std::setprecision(1)<< t
					  								<<"  Timestep Number = "<<(size_counter)
					  								<<"  CPU Time = " <<std::setw(5)<<std::setprecision(1)
														<< ((float)ds/CLOCKS_PER_SEC) << "s \n";
					//std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Current Precise Time: " << t << "\n";
					//std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Z Displacement: " << U(20) << "\n";
					//std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Internal Energy: " << energy_int_new << "\n";
					//std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"External Work: " << energy_ext_new << "\n";
					//std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Kinetic Energy: " << energy_kin << "\n";

					/*print current frame, current time, and energy to individual .txt files*/
					append_double2textWithTime(internal_energy,time_temp_1 -1,t,energy_int_new);
					append_double2textWithTime(external_energy,time_temp_1 - 1,t,energy_ext_new);
					append_double2textWithTime(kinetic_energy,time_temp_1 - 1,t,energy_kin);
					append_double2textWithTime(total_energy,time_temp_1 - 1,t,energy_total);
			}

			s_prev = s;
			s = clock();
			ds = s - s_prev;
			size_counter = size_counter+1;
		}

}
