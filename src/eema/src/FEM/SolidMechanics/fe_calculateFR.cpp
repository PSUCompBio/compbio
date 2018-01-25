#include "functions.h"

using namespace Eigen;

void fe_calculateFR(VectorXd& fr, VectorXd& fi_curr, VectorXd& mm, VectorXd& A) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "displacement") {
			int n_dof = bc[i].getNumDOF();
			VectorXi dof = bc[i].getDOF();
			for (int j = 0; j < n_dof; j++) {
				int dof_id = dof(j);
				fr(dof_id) = fi_curr(dof_id) + mm(dof_id) * A(dof_id);
			}
		}
	}

}


double fe_calculateTotalAppliedFR(VectorXd& fr) {

	double total_applied_fr = 0;

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "displacement") {
			double amplitude = bc[i].getAmplitude();
			if (amplitude != 0) {
				int n_dof = bc[i].getNumDOF();
				VectorXi dof = bc[i].getDOF();
				for (int j = 0; j < n_dof; j++) {
					int dof_id = dof(j);
					total_applied_fr = total_applied_fr + fr(dof_id);
				}
			}
		}
	}

	return total_applied_fr;

}
