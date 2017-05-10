

#include"functions.h"

using namespace Eigen;

VectorXd fe_shapes_8(double rvalue,double svalue,double tvalue){

	VectorXd shapes(8);

	int i = 0;
	double tmp = 0.125*(1-rvalue)*(1-svalue)*(1-tvalue);
	shapes(i) = tmp;

	i = i+1;
	tmp = 0.125*(1+rvalue)*(1-svalue)*(1-tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1+rvalue)*(1+svalue)*(1-tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1-rvalue)*(1+svalue)*(1-tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1-rvalue)*(1-svalue)*(1+tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1+rvalue)*(1-svalue)*(1+tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1+rvalue)*(1+svalue)*(1+tvalue);
        shapes(i) = tmp;

	i = i+1;
        tmp = 0.125*(1-rvalue)*(1+svalue)*(1+tvalue);
        shapes(i) = tmp;

	return shapes;
}
