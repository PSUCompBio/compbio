#include "functions.h"

using namespace Eigen;

VectorXd
fe_findIntgPoints_1d(VectorXd& x, VectorXd& y, VectorXd& z, double point, double length)
{
	VectorXd intg_points = VectorXd::Zero(3);

	intg_points(0) = x(0) + (((1.0 + point) / 2.0) ) * (x(1) - x(0));
	intg_points(1) = y(0) + (((1.0 + point) / 2.0) ) * (y(1) - y(0));
	intg_points(2) = z(0) + (((1.0 + point) / 2.0) ) * (z(1) - z(0));

	return intg_points;
}
