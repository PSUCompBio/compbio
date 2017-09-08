#include "functions.h"

using namespace Eigen;

void fe_calculateVDE(VectorXd& U_prev, VectorXd& U, VectorXd& f_damp_prev, VectorXd& f_damp_curr, double& energy_vd_old, double& energy_vd_new) {
	VectorXd del_U = U - U_prev;
	energy_vd_new = energy_vd_old + 0.5 * (del_U.dot(f_damp_prev + f_damp_curr));
	f_damp_prev        = f_damp_curr;
	energy_vd_old = energy_vd_new;
}
