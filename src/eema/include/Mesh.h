/*
 * Mesh.h
 *
 *  Created on: Dec 15, 2016
 *      Author: vsg111
 */

#ifndef HEADERS_MESH_H_
#define HEADERS_MESH_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <list>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Mesh
{

  std::string mesh_name;

  /** This contains the nodal data read from the mesh */
  MatrixXd nodes;

  /** This contains the element data read from the mesh */
  MatrixXi elements;

  /** This contains the element and node processor data read from the mesh and metis file */
  VectorXi pe;
  MatrixXi pn;

  /** This contains the nodal data after preprocessing or any other operation */
  MatrixXd nodes_new;
  MatrixXd *nodes_new_pointer;

  /** This contains the element data after preprocessing or any other operation */
  MatrixXi elements_new;
  MatrixXi *elements_new_pointer;

  /** Mesh - Nodal Displacement Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd U;
  VectorXd *U_pointer;

  /** Mesh - Nodal Velocity Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd V;
  VectorXd *V_pointer;

  /** Mesh - Nodal Acceleration Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd A;
  VectorXd *A_pointer;

  /** Mesh - Nodal Stress */
  VectorXd nodal_stress;

  /** Mesh - Nodal Strain */
  VectorXd nodal_strain;

  /** Mesh - Element-wise stress vector */
  VectorXd element_stress;
  VectorXd* element_stress_pointer;

  /** Mesh - Element-wise strain vector */
  VectorXd element_strain;
  VectorXd* element_strain_pointer;

  /** Mesh - Element-wise damage scalar */
  VectorXd damage;
  VectorXd* damage_pointer;

  /** Mesh - Element-wise stretch scalar */
  VectorXd stretch_min;
  VectorXd* stretch_min_pointer;

  /** Mesh - Element-wise stretch scalar */
  VectorXd stretch_max;
  VectorXd* stretch_max_pointer;

  VectorXd element_charateristic;
  VectorXd* element_charateristic_pointer;

  VectorXd nodal_electric_potential;
  VectorXd* nodal_electric_potential_pointer;

public:
  void readMesh(std::string name, MatrixXd& n, MatrixXi& e, VectorXi& pe, MatrixXi& pn);
  void readNodalKinematics(VectorXd& Usystem, VectorXd& Vsystem, VectorXd& Asystem);
  void readNodalStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp);
  void readElementStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp);
  void readDamage(VectorXd& d_tot);
  void readStretchMin(VectorXd& lambda_min);
  void readStretchMax(VectorXd& lambda_max);
  void readNodalElectroPhysics(VectorXd& VP);

  MatrixXd getNodes();
  MatrixXi getElements();
  VectorXi getElementsProcessorID();
  MatrixXi getNodesProcessorID();
  MatrixXd getNewNodes();
  MatrixXi getNewElements();
  MatrixXd* getNewNodesPointer();
  MatrixXi* getNewElementsPointer();

  std::string getName();
  VectorXd getMaxCharLength(std::string choice);
  VectorXd getMinCharLength(std::string choice);

  /** Outputs nodal displacement vector */
  VectorXd getNodalDisp();
  VectorXd* getNodalDispPointer();

  /** Outputs nodal velocity vector */
  VectorXd getNodalVel();
  VectorXd* getNodalVelPointer();

  /** Outputs nodal acceleration vector */
  VectorXd getNodalAcc();
  VectorXd* getNodalAccPointer();

  /** Outputs element-wise stress of the mesh */
  VectorXd getCellStress();
  VectorXd* getCellStressPointer();

  /** Ouputs element-wise strain of the entire mesh */
  VectorXd getCellStrain();
  VectorXd* getCellStrainPointer();

  /** Ouputs element-wise damage of the entire mesh */
  VectorXd getCellDamage();
  VectorXd* getCellDamagePointer();

  /** Ouputs element-wise stretch of the entire mesh */
  VectorXd getCellStretchMin();
  VectorXd* getCellStretchMinPointer();

  /** Ouputs element-wise stretch of the entire mesh */
  VectorXd getCellStretchMax();
  VectorXd* getCellStretchMaxPointer();

  /** Outputs Nodal Stress of the entire mesh */
  VectorXd getNodalStress();
  /** Outputs Nodal Strain of the entire mesh */
  VectorXd getNodalStrain();

  void calculateElementCharateristic();
  VectorXd* getElementCharacteristicPointer();

  VectorXd getNodalEPotential();
  VectorXd* getNodalEPotentialPointer();

  int getNumNodes();
  int getNumElements();
  int getNumNodesPerElement();

  void preprocessMesh(std::string choice);
  void replaceNodes(MatrixXd A, std::string B);
  void replaceElements(MatrixXi A, std::string B);
  void append2nodes(std::string choice, VectorXd& b);
  void append2elements(std::string choice, VectorXi& b);
  void reNumber();
  void checkMesh();
  void printInfo();
};

#endif /* HEADERS_MESH_H_ */
