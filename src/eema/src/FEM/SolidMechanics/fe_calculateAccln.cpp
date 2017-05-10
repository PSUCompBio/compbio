#include "functions.h"

using namespace Eigen;

void fe_calculateAccln(VectorXd& A, VectorXd& mm, VectorXd& F_net) {

	A = VectorXd::Zero(F_net.size());

	for (int i = 0; i < F_net.size(); i++) {

		A(i) = (F_net(i) / mm(i));

	}

}
