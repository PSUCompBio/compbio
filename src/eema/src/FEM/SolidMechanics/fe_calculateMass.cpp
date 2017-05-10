#include "functions.h"

using namespace Eigen;

VectorXd fe_calculateMass(VectorXd mm, std::string type){

    if(type=="direct_lumped"){ // Embedded Mass
      mm = fe_calculateMassDirectLumped(mm);
    }

  return mm;
}

VectorXd fe_calculateMassDirectLumped(VectorXd mm){

  MatrixXd nodes = mesh[0].getNewNodes();
  MatrixXi elements = mesh[0].getNewElements();

  int nel = elements.rows(); /*! number of elements */
  int nnel = (elements.cols()-2); // number of nodes per element
  int nnode = nodes.rows(); // number of nodes
  int sdof = nnode * ndof; // system degrees of freedom
  int edof = nnel * ndof; // element degrees of freedom

  for(int i=0;i<nel;i++){

    VectorXi nodes_local = VectorXi::Zero(elements.cols()-2);

		for(int j=0;j<(elements.cols()-2);j++){
			nodes_local(j) = elements(i,j+2);
		}

		VectorXd m_element = VectorXd::Zero(edof); // mass of hex elements
		m_element = fe_massLumped(nodes,elements.row(i));

		/* TRUSS ANALYSIS - WORKS FOR ONLY SINGLE ELEMENT
		// Truss elements
		// by default, mass of the truss is added to the hex element in a lumped form
		double rho_truss = fe_get_mats(elements_truss(i,1),0);
		double A_truss = 0.5;
		MatrixXd m_truss = MatrixXd::Zero(edof,edof);
		m_truss = fe_mass_truss(rho_truss,A_truss,edof,nodes_truss,elements_truss);

		MatrixXd m_correction = MatrixXd::Zero(edof,edof);
		m_correction = fe_mass_truss(rho,A_truss,edof,nodes_truss,elements_truss); // This is the redundant volume

		// mm = fe_assemble_mass(m_hex+m_truss); // This is where the mass magic happens (m_hex + m_truss - m_correction)*/

    mm = fe_scatter(mm,m_element,nodes_local,sdof);
  }

  return mm;

}
