#include "functions.h"

using namespace Eigen;

void fe_calculateIE(VectorXd& U_prev, VectorXd& U, VectorXd& fi_prev, VectorXd& fi_curr, double& energy_int_old, double& energy_int_new) {
	VectorXd del_U = U - U_prev;
	energy_int_new = energy_int_old + 0.5 * (del_U.dot(fi_prev + fi_curr));
	fi_prev        = fi_curr;
	energy_int_old = energy_int_new;
}