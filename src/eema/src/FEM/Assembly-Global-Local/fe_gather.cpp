#include<iostream>

#include "functions.h"
using namespace Eigen;

VectorXd fe_gather(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof){

	int edof = (node_list.size())*ndof;
	VectorXi index = fe_find_index(node_list);

	int ii = 0;

	for (int i=0; i<edof; i++){
		ii = index(i);
		local_vec(i) = global_vec(ii);
	}

	return local_vec;

}
