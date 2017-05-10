/*
 * Mesh.h
 *
 *  Created on: Dec 15, 2016
 *      Author: vsg111
 */

#ifndef HEADERS_MESH_H_
#define HEADERS_MESH_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Mesh{
	/** This contains the nodal data read from the mesh */
	MatrixXd nodes;
	/** This contains the element data read from the mesh */
	MatrixXi elements;
	/** This contains the nodal data after preprocessing or any other operation */
	MatrixXd nodes_new;
	/** This contains the element data after preprocessing or any other operation */
	MatrixXi elements_new;

	/** Mesh - Nodal Displacement Vector. This Vector stores the nodal displacement at one time point only. */
	VectorXd U;

	/** Mesh - Nodal Velocity Vector. This Vector stores the nodal displacement at one time point only. */
	VectorXd V;

	/** Mesh - Nodal Acceleration Vector. This Vector stores the nodal displacement at one time point only. */
	VectorXd A;

	/** Mesh - Element-wise stress vector */
	MatrixXd stress;

	/** Mesh - Element-wise strain vector */
	MatrixXd strain;

public:
	void readMesh(MatrixXd n, MatrixXi e);
	void readNodalKinematics(VectorXd Usystem, VectorXd Vsystem, VectorXd Asystem);
	void readStressStrain(MatrixXd stress_tmp,MatrixXd strain_tmp);

	MatrixXd getNodes();
	MatrixXi getElements();
	MatrixXd getNewNodes();
	MatrixXi getNewElements();
	VectorXd getMaxCharLength(std::string choice);
	VectorXd getMinCharLength(std::string choice);

	/** Outputs nodal displacement vector */
	VectorXd getNodalDisp();

	/** Outputs nodal velocity vector */
	VectorXd getNodalVel();

	/** Outputs nodal acceleration vector */
	VectorXd getNodalAcc();

	/** Outputs element-wise stress of the mesh */
	MatrixXd getCellStress();

	/** Ouputs element-wise strain of the entire mesh */
	MatrixXd getCellStrain();

	void preprocessMesh();
	void replaceNodes(MatrixXd A, std::string B);
	void replaceElements(MatrixXi A, std::string B);
	void checkMesh();
};

#endif /* HEADERS_MESH_H_ */
