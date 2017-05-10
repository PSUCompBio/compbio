#include "functions.h"

using namespace Eigen;

void fe_timeUpdate(VectorXd& U, VectorXd& V, VectorXd& V_half, VectorXd& A, double& t, double& dT, std::string time_integration_method) {

	if (time_integration_method == "newmark-beta-central-difference") {

		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		/** Update the time by adding full time step */
		t = t + dT;

		/** Partially Update Nodal Velocities */
		V_half = V + ((1 - gamma) * dT * A);

		/** Update Nodal Displacements */
		U = U + (dT * V) + ((pow(dT, 2) / 2.0) * (1 - (2 * beta)) * A) + (beta * dT * dT * A);
		fe_apply_bc_displacement(U, t);
	}

	if (time_integration_method == "newmark-beta-undamped-trapezoidal") {
		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		/** Update the time by adding full time step */
		t = t + dT;

		/** Partially Update Nodal Velocities */
		V_half = V + ((1 - gamma) * dT * A);

		/** Update Nodal Displacements */
		U = U + (dT * V) + ((pow(dT, 2) / 2.0) * (1 - (2 * beta)) * A) + (beta * dT * dT * A);
		fe_apply_bc_displacement(U, t);
	}

	if (time_integration_method == "newmark-beta-custom") {

		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		/** Update the time by adding full time step */
		t = t + dT;

		/** Partially Update Nodal Velocities */
		V_half = V + ((1 - gamma) * dT * A);

		/** Update Nodal Displacements */
		U = U + (dT * V) + ((pow(dT, 2) / 2.0) * (1 - (2 * beta)) * A) + (beta * dT * dT * A);
		fe_apply_bc_displacement(U, t);

	}

}

void fe_timeUpdate_velocity(VectorXd& V, VectorXd& V_half, VectorXd& A, double& t, double& dT, std::string time_integration_method) {

	if (time_integration_method == "newmark-beta-central-difference") {
		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		V = V_half + (gamma * dT * A);
		fe_apply_bc_velocity(V, t);
	}

	if (time_integration_method == "newmark-beta-undamped-trapezoidal") {

		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		V = V_half + (gamma * dT * A);
		fe_apply_bc_velocity(V, t);
	}

	if (time_integration_method == "newmark-beta-custom") {
		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

		V = V_half + (gamma * dT * A);
		fe_apply_bc_velocity(V, t);
	}
}