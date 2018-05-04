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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();	
			//std::cout<<x_data<<" "<<y_data<<" "<<z_data<<" "<<time_data<<" "<<direction<<" "<<n_steps<<" ";
			//std::exit(1);		
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				U(c) = fe_function(input_disp_amp, disp_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				V(c) = fe_function_derivative(input_disp_amp, disp_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
			}
		}

		if (type == "velocity") {
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				V(c) = fe_function(input_vel_amp, vel_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				A(c) = fe_function_d_derivative(input_disp_amp, disp_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
			}
		}

		if (type == "velocity") {
			double input_vel_amp = bc[i].getAmplitude();
			std::string vel_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				A(c) = fe_function_derivative(input_vel_amp, vel_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
			}
		}

		if (type == "acceleration") {
			double input_acc_amp = bc[i].getAmplitude();
			std::string acc_curve = bc[i].getTimeBehavior();
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				A(c) = fe_function(input_acc_amp, acc_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				fe(c) = fe_function(input_load_amp, load_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				I(c) = fe_function(input_load_amp, load_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				VP(c) = fe_function(input_load_amp, load_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
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
			int n_steps = bc[i].getN_steps();
			int direction = 0;
			VectorXd time_data = bc[i].getTime_data();
			VectorXd x_data = bc[i].getX_data();
			VectorXd y_data = bc[i].getY_data();
			VectorXd z_data = bc[i].getZ_data();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				if(c%3 == 0)
				direction = 1;
				else if (c%3 == 1)
				direction = 2;
				else if (c%3 == 2) direction = 3;
				for (int k = 0; k < kk.cols(); k++) {
					kk(c, k) = 0;
				}
				kk(c, c) = 1;
				ff(c) = fe_function(input_load_amp, load_curve, time, n_steps, time_data, x_data, y_data, z_data, direction);
			}
		}
	}
}
