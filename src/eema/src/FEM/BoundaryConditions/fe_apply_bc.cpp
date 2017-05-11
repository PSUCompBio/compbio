#include"functions.h"

using namespace Eigen;

void fe_apply_bc_displacement(VectorXd& U, double& time) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "displacement") {
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				U(c) = fe_function(input_disp_amp, disp_curve, time);
			}
		}
	}

}

void fe_apply_bc_velocity(VectorXd& V, double& time) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();

		if (type == "displacement") {
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				V(c) = fe_function_derivative(input_disp_amp, disp_curve, time);
			}
		}

		if (type == "velocity") {
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				V(c) = fe_function(input_vel_amp, vel_curve, time);
			}
		}

	}

}

void fe_apply_bc_acceleration(VectorXd& A, double& time) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();

		if (type == "displacement") {
			double input_disp_amp = bc[i].getAmplitude();
			std::string disp_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				A(c) = fe_function_d_derivative(input_disp_amp, disp_curve, time);
			}
		}

		if (type == "velocity") {
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				A(c) = fe_function_derivative(input_vel_amp, vel_curve, time);
			}
		}

		if (type == "acceleration") {
			double input_acc_amp = bc[i].getAmplitude();
			std::string acc_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				A(c) = fe_function(input_acc_amp, acc_curve, time);
			}
		}

	}

}

void fe_apply_bc_load(VectorXd& fe, double& time) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "load") {
			double input_load_amp = bc[i].getAmplitude();
			std::string load_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				fe(c) = fe_function(input_load_amp, load_curve, time);
			}
		}
	}

}

void fe_apply_bc_current(VectorXd& I, double &time) {
	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "current") {
			double input_load_amp = bc[i].getAmplitude();
			std::string load_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				I(c) = fe_function(input_load_amp, load_curve, time);
			}
		}
	}
}

void fe_apply_bc_potential(VectorXd& VP, double &time) {
	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "potential") {
			double input_load_amp = bc[i].getAmplitude();
			std::string load_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				VP(c) = fe_function(input_load_amp, load_curve, time);
			}
		}
	}
}

void fe_apply_bc_potential(MatrixXd& kk, VectorXd& ff, double time) {
	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "potential") {
			double input_load_amp = bc[i].getAmplitude();
			std::string load_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				for (int k = 0; k < kk.cols(); k++) {
					kk(c, k) = 0;
				}
				kk(c, c) = 1;
				ff(c) = fe_function(input_load_amp, load_curve, time);
			}
		}
	}
}
