#include "functions.h"

using namespace Eigen;

void fe_calculateMass(VectorXd& m_system, std::string type) {

  if (type == "direct_lumped" && embedded_constraint == true) { // Embedded Mass
    for (int i = 0; i < num_constraints; ++i) {
      if (cons[i].getName() == "embedded") {
        std::string host = cons[i].get_EmbedMaster();
        std::string slave = cons[i].get_EmbedSlave();
        bool correct_vr = cons[i].get_EmbedAddressVR();
        int host_id, embed_id;
        for (int i = 0; i < num_meshes; ++i) {
          std::string name = mesh[i].getName();
          if (name == host) {
            host_id = i;
          }
          if (name == slave) {
            embed_id = i;
          }
        }
        fe_calculateMassDirectLumped_embed(m_system, host_id, embed_id, correct_vr);
      }
    }
  }
  else if (type == "direct_lumped" && embedded_constraint != true) {
    fe_calculateMassDirectLumped(m_system, 0);
  }
  else {
    fe_calculateMassDirectLumped(m_system, 0);
  }

}

void fe_calculateMassDirectLumped(VectorXd& m_system, int mesh_id) {

  MatrixXd* nodes = mesh[mesh_id].getNewNodesPointer();
  MatrixXi* elements = mesh[mesh_id].getNewElementsPointer();

  int nel = mesh[mesh_id].getNumElements();
  int nnel = mesh[mesh_id].getNumNodesPerElement();
  int nnode = mesh[mesh_id].getNumNodes();
  int sdof = nnode * ndof;
  int edof = nnel * ndof;

  for (int i = 0; i < nel; i++) {
    VectorXd m_element = VectorXd::Zero(edof); // mass of hex elements
    m_element = fe_massLumped(nodes, (*elements).row(i));
    fe_scatter_pbr(m_system, m_element, (*elements).block<1, 8>(i, 2), sdof);
  }

  nodes = NULL;
  elements = NULL;

}

void fe_calculateMassDirectLumped_embed(VectorXd& m_system, int host_id, int embed_id, bool address_vr) {

  //std::cout << "Entered Mass - Embedded" << "\n";

  MatrixXd* nodes_host = mesh[host_id].getNewNodesPointer();
  MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();
  MatrixXd* nodes_embed = mesh[embed_id].getNewNodesPointer();
  MatrixXi* elements_embed = mesh[embed_id].getNewElementsPointer();

  int nel = mesh[host_id].getNumElements(); /*! number of elements */
  int nnel = mesh[host_id].getNumNodesPerElement(); // number of nodes per element
  int nnode = mesh[host_id].getNumNodes(); // number of nodes
  int sdof = nnode * ndof; // system degrees of freedom
  int edof = nnel * ndof; // element degrees of freedom

  VectorXd xcoord_embed((*elements_embed).cols() - 2);
  VectorXd ycoord_embed((*elements_embed).cols() - 2);
  VectorXd zcoord_embed((*elements_embed).cols() - 2);

  for (int i = 0; i < nel; i++) {

    VectorXd m_element = VectorXd::Zero(edof); // mass of hex elements
    m_element = fe_massLumped(nodes_host, (*elements_host).row(i));

    for (int fib = 0; fib < (*elements_embed).rows(); fib++) {
      for (int j = 0; j < (*elements_embed).cols() - 2; j++) {
        int g = (*elements_embed)(fib, j + 2);
        xcoord_embed(j) = (*nodes_embed)(g, 1);
        ycoord_embed(j) = (*nodes_embed)(g, 2);
        zcoord_embed(j) = (*nodes_embed)(g, 3);
      }

      double rho_truss = fe_get_mats((*elements_embed)(fib, 1), 0);
      double rho_host = fe_get_mats((*elements_host)(i, 1), 0);

      double m_truss = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed) * area_truss * rho_truss;
      double m_correction = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed) * area_truss * rho_host;

      for (int k = 0; k < edof; k++) {
        m_element(k) = m_element(k) + (m_truss / nnel);
        if (address_vr == true) {
          m_element(k) = m_element(k) - (m_correction / nnel);
        }
      }
    }

    fe_scatter_pbr(m_system, m_element, (*elements_host).block<1, 8>(i, 2), sdof);
  }

  nodes_host = NULL;
  elements_host = NULL;
  nodes_embed = NULL;
  elements_embed = NULL;

}
