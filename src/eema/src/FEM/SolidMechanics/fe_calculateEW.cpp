#include "functions.h"

using namespace Eigen;

void fe_calculateEW(VectorXd& U_prev, VectorXd& U, VectorXd& fe_prev, VectorXd& fe, VectorXd& fr_prev, VectorXd& fr_curr, double& energy_ext_old, double& energy_ext_new) {

	VectorXd del_U = U - U_prev;
	energy_ext_new = energy_ext_old + 0.5 * (del_U.dot((fe_prev + fr_prev) + (fe + fr_curr)));
	fe_prev        = fe;
	fr_prev        = fr_curr;
	U_prev 		   = U;
	energy_ext_old = energy_ext_new;
}