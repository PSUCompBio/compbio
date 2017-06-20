/*
 * System.h
 *
 *  Created on: Dec 15, 2016
 *      Author: vsg111
 */

#ifndef HEADERS_SYSTEM_H_
#define HEADERS_SYSTEM_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class System
{

  /** System - This contains the nodal data for the entire model */
  VectorXd coordinates;
  VectorXd *coordinates_pointer;

  /** System - This contains the nodal data for the entire model */
  int n_nodes;

  /** System - This contains the nodal data after preprocessing or any other operation */
  VectorXd coordinates_new;
  VectorXd *coordinates_new_pointer;

  /** System - Nodal Displacement Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd displacement;
  VectorXd *displacement_pointer;

  /** System - Nodal Velocity Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd velocity;
  VectorXd *velocity_pointer;

  /** System - Nodal Acceleration Vector. This Vector stores the nodal displacement at one time point only. */
  VectorXd acceleration;
  VectorXd *acceleration_pointer;

  /** System - Nodal Internal Force Vector. This Vector stores the nodal internal forces at one time point only. */
  VectorXd force_internal;
  VectorXd *force_internal_pointer;

  /** System - Nodal External Force Vector. This Vector stores the nodal external forces at one time point only. */
  VectorXd force_external;
  VectorXd *force_external_pointer;

  /** System - Nodal Reaction Force Vector. This Vector stores the nodal reaction forces at one time point only. */
  VectorXd force_reaction;
  VectorXd *force_reaction_pointer;

  /** System - Nodal Net Force Vector. This Vector stores the nodal net forces at one time point only. */
  VectorXd force_net;
  VectorXd *force_net_pointer;

  /** Mesh - Nodal Stress */
  VectorXd stress;

  /** Mesh - Nodal Strain */
  VectorXd strain;

  VectorXd element_charateristic;
  VectorXd* element_charateristic_pointer;

  VectorXd nodal_electric_potential;
  VectorXd* nodal_electric_potential_pointer;

public:
  void readMesh(std::string name, MatrixXd& n, MatrixXi& e);
  void readNodalKinematics(VectorXd& Usystem, VectorXd& Vsystem, VectorXd& Asystem);
  void readNodalStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp);
  void readElementStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp);
  void readNodalElectroPhysics(VectorXd& VP);

  MatrixXd getNodes();
  MatrixXi getElements();

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
