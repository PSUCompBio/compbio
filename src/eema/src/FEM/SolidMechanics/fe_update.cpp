#include "functions.h"

using namespace Eigen;

/*
 * This function update the nodal coordinates based on the displacements.
 */

MatrixXd fe_updateNodes(MatrixXd nodes,VectorXd displacements){
	MatrixXd updated_nodes = MatrixXd::Zero(nodes.rows(),nodes.cols());
	int size_counter = 0;

	for(int i=0;i<nodes.rows();i++){
			updated_nodes(i,0) = nodes(i,0);
			updated_nodes(i,1) = nodes(i,1) + displacements(size_counter);
			updated_nodes(i,2) = nodes(i,2) + displacements(size_counter+1);
			updated_nodes(i,3) = nodes(i,3) + displacements(size_counter+2);
			size_counter = size_counter+3;
	}

	return updated_nodes;
}



