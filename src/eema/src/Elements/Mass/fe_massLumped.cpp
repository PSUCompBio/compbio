#include "functions.h"
using namespace Eigen;

VectorXd fe_massLumped(MatrixXd* nodes, VectorXi elements_row) {

  double mass = 0;
  int nnel = elements_row.size() - 2;
  int nnode = (*nodes).rows();
  int edof = nnel * ndof;
  VectorXd m_element = VectorXd::Zero(edof);
  VectorXd xcoord = VectorXd::Zero(nnel);
  VectorXd ycoord = VectorXd::Zero(nnel);
  VectorXd zcoord = VectorXd::Zero(nnel);

  for (int j = 0; j < nnel; j++) {
    int g = elements_row(j + 2);
    xcoord(j) = (*nodes)(g, 1);
    ycoord(j) = (*nodes)(g, 2);
    zcoord(j) = (*nodes)(g, 3);
  }

  double rho = fe_get_mats(elements_row(1), 0); // material density
  double volume_element = fe_calVolume(xcoord, ycoord, zcoord);
  // The above function fe_calVolume gives volume for a 3d element, area for a 2d element and length for a 1d element.

  if (xcoord.size() == 2) {
    volume_element = volume_element * area_truss;
    // 0.5 is area of the single 1D element cross-section
  }

  if (xcoord.size() == 4) {
    volume_element = volume_element * 1;
    // volume_element = volume_element * thickness;
  }

  mass = rho * volume_element;

  for (int i_1 = 0; i_1 < edof; i_1++) {
    m_element(i_1) = (mass / (nnel));
  }

  return m_element;
}
