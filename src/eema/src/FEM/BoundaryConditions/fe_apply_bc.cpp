#include"functions.h"

using namespace Eigen;

VectorXd fe_apply_bc_displacement(VectorXd U,double time){

 for(int i=0;i<bc_types;i++){
		std::string type = bc[i].getType();
		if(type=="displacement"){
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				U(c) = fe_function(input_disp_amp,disp_curve,time);
			}
		}
	}

	return U;
}

VectorXd fe_apply_bc_velocity(VectorXd V,double time){

 for(int i=0;i<bc_types;i++){
		std::string type = bc[i].getType();

		if(type=="displacement"){
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				V(c) = fe_function_derivative(input_disp_amp,disp_curve,time);
			}
		}

		if(type=="velocity"){
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				V(c) = fe_function(input_vel_amp,vel_curve,time);
			}
		}

	}

	return V;
}

VectorXd fe_apply_bc_acceleration(VectorXd A,double time){

 for(int i=0;i<bc_types;i++){
		std::string type = bc[i].getType();

		if(type=="displacement"){
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				A(c) = fe_function_d_derivative(input_disp_amp,disp_curve,time);
			}
		}

		if(type=="velocity"){
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				A(c) = fe_function_derivative(input_vel_amp,vel_curve,time);
			}
		}

		if(type=="acceleration"){
			double input_acc_amp = bc[i].getAmplitude();
			std::string acc_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				A(c) = fe_function(input_acc_amp,acc_curve,time);
			}
		}

	}

	return A;
}

VectorXd fe_apply_bc_load(VectorXd fe,double time){

	for(int i=0;i<bc_types;i++){
		std::string type = bc[i].getType();
		if(type=="load"){
			double input_load_amp = bc[i].getAmplitude();
			std::string load_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for(int j=0;j<number_of_dof;j++){
				int c = local_dof(j);
				fe(c) = fe_function(input_load_amp,load_curve,time);
			}
		}
	}

	return fe;
}

MatrixXd fe_apply_bc_stiffness(MatrixXd kk,VectorXi bcdof,VectorXd bcval){

	int n = bcdof.size();
	int sdof = kk.rows();

	for(int i=0;i<n;i++){
		int c = bcdof(i);
		for(int j=0;j<sdof;j++){
			kk(c,j) = 0;
		}
		kk(c,c) = 1;
	}
	return kk;
}
