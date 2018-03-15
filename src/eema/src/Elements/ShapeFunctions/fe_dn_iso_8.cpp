#include"functions.h"

using namespace Eigen;

void fe_dniso_8(VectorXd& dndr, VectorXd& dnds, VectorXd& dndt, double& rvalue, double& svalue, double& tvalue) {

	double nnr = (1 - svalue) * (1 - tvalue);
	double pnr = (1 + svalue) * (1 - tvalue);
	double npr = (1 - svalue) * (1 + tvalue);
	double ppr = (1 + svalue) * (1 + tvalue);

	double nnt = (1 - rvalue) * (1 - svalue);
	double pnt = (1 + rvalue) * (1 - svalue);
	double npt = (1 - rvalue) * (1 + svalue);
	double ppt = (1 + rvalue) * (1 + svalue);

	double nns = (1 - rvalue) * (1 - tvalue);
	double pns = (1 + rvalue) * (1 - tvalue);
	double nps = (1 - rvalue) * (1 + tvalue);
	double pps = (1 + rvalue) * (1 + tvalue);

	double use;
	int i = 0, sign = 1;

	for (i = 0; i < 8; i++) {
		if (i == 1 || i == 2 || i == 5 || i == 6)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 1)
			use = nnr;
		else if (i == 2 || i == 3)
			use = pnr;
		else if (i == 4 || i == 5)
			use = npr;
		else
			use = ppr;

		dndr(i) = 0.125 * sign * use;

		if (i == 2 || i == 3 || i == 6 || i == 7)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 3)
			use = nns;
		else if (i == 1 || i == 2)
			use = pns;
		else if (i == 4 || i == 7)
			use = nps;
		else
			use = pps;

		dnds(i) = 0.125 * sign * use;

		if (i > 3)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 4)
			use = nnt;
		else if (i == 1 || i == 5)
			use = pnt;
		else if (i == 3 || i == 7)
			use = npt;
		else
			use = ppt;

		dndt(i) = 0.125 * sign * use;
	}
}

void fe_dniso_8_array(double ****dndr, double ****dnds, double ****dndt, double& rvalue, double& svalue, double& tvalue, int intx, int inty, int intz) {

	double nnr = (1 - svalue) * (1 - tvalue);
	double pnr = (1 + svalue) * (1 - tvalue);
	double npr = (1 - svalue) * (1 + tvalue);
	double ppr = (1 + svalue) * (1 + tvalue);

	double nnt = (1 - rvalue) * (1 - svalue);
	double pnt = (1 + rvalue) * (1 - svalue);
	double npt = (1 - rvalue) * (1 + svalue);
	double ppt = (1 + rvalue) * (1 + svalue);

	double nns = (1 - rvalue) * (1 - tvalue);
	double pns = (1 + rvalue) * (1 - tvalue);
	double nps = (1 - rvalue) * (1 + tvalue);
	double pps = (1 + rvalue) * (1 + tvalue);

	double use;
	int i = 0, sign = 1;

	for (i = 0; i < 8; i++) {
		if (i == 1 || i == 2 || i == 5 || i == 6)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 1)
			use = nnr;
		else if (i == 2 || i == 3)
			use = pnr;
		else if (i == 4 || i == 5)
			use = npr;
		else
			use = ppr;

		dndr[intx][inty][intz][i] = 0.125 * sign * use;

		if (i == 2 || i == 3 || i == 6 || i == 7)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 3)
			use = nns;
		else if (i == 1 || i == 2)
			use = pns;
		else if (i == 4 || i == 7)
			use = nps;
		else
			use = pps;

		dnds[intx][inty][intz][i] = 0.125 * sign * use;

		if (i > 3)
			sign = 1;
		else
			sign = -1;

		if (i == 0 || i == 4)
			use = nnt;
		else if (i == 1 || i == 5)
			use = pnt;
		else if (i == 3 || i == 7)
			use = npt;
		else
			use = ppt;

		dndt[intx][inty][intz][i] = 0.125 * sign * use;
	}
}
