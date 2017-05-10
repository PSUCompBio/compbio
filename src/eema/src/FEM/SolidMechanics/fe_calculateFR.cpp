#include "functions.h"

using namespace Eigen;

void fe_calculateFR(VectorXd& fr, int& sdof, VectorXd& fi_curr, VectorXd& mm, VectorXd& A) {

	for (int i = 0; i < bc_types; i++) {
		std::string type = bc[i].getType();
		if (type == "displacement") {
			int number_of_dof = bc[i].getNumDOF();
			VectorXi local_dof = bc[i].getDOF();
			for (int j = 0; j < number_of_dof; j++) {
				int c = local_dof(j);
				fr(c) = fi_curr(c) + mm(c) * A(c);
			}
		}
	}

}