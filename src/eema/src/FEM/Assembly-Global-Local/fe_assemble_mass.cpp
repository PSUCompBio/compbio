/*
 * fe_assemble_mass.cpp
 *
 * Description: This function takes an element level mass matrix and assembles it into a
 * global mass matrix.
 *
 *  Created on: Mar 15, 2017
 *  Author: Harsha T. Garimella
 *
 *  Algorithm:
 *  1. For the particular element, get the node list.
 *  2. Find the indices corresponding to each node. (finding using fe_find_index.cpp)
 *  	- Each index represents the DOF of each node.
 *  3. Based on the index value, add the local mass matrix to global mass matrix.
 *
 *
 */

#include<iostream>

#include "functions.h"
using namespace Eigen;

MatrixXd fe_assemble_mass(MatrixXd mm, MatrixXd m, VectorXi node_list, int sdof){

	int edof = (node_list.size())*ndof;
	VectorXi index = fe_find_index(node_list);

	int ii = 0;
	int jj = 0;

	for (int i=0; i<edof; i++){
		ii = index(i);
		for (int j=0; j<edof; j++){
			jj = index(j);
			mm(ii,jj) = mm(ii,jj) + m(i,j) ;
		}
	}

	return mm;

}
