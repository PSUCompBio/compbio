
#include"functions.h"

using namespace Eigen;

void fe_dniso_8(VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, double& rvalue, double& svalue, double& tvalue) {

	// Calculating dndr

	int i = 0;
	double tmp = 0.125 * (-1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;

	// Calculating dnds

	i = 0;
	tmp = 0.125 * (1 - rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

	// Calculating dndt

	i = 0;
	tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;

}


