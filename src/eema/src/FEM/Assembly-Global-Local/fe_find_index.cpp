/*
 * fe_find_index.cpp
 *
 * Description: This function calculates the index vector for each element.
 * Index vector is a vector consisting of system DOF number for each node of the element.
 *
 *  Created on: Mar 15, 2017
 *  Author: Harsha T. Garimella
 */

#include<iostream>

#include "functions.h"
using namespace Eigen;

VectorXi fe_find_index(VectorXi node_list){

	VectorXi index;
	int edof = (node_list.size()) * ndof;
	int start = 0;
	int k = 0;

	index = VectorXi::Zero(edof);

	for(int i=0; i<node_list.size(); i++){
		start = (node_list(i))*(ndof);
		for (int j=0; j<ndof ; j++){
			index(k) = start + j;
			k = k+1;
		}
	}

	return index;
}



