#include "functions.h"
using namespace Eigen;

int fe_compute_host(VectorXd& A, MatrixXd& nodes_host, MatrixXi& elements_host_tmp) {

  int host_id;

  MatrixXi faces = MatrixXi::Zero(12, 3);
  faces << 0, 1, 2,
        2, 3, 0,
        6, 5, 4,
        4, 7, 6,
        5, 1, 0,
        0, 4, 5,
        6, 7, 3,
        3, 2, 6,
        1, 5, 6,
        6, 2, 1,
        7, 4, 0,
        0, 3, 7;

  for (int i = 0; i < elements_host_tmp.rows(); i++) {

    int counter = 0;
    for (int j = 0; j < faces.rows(); j++) {

      int i_1, i_2, i_3;
      VectorXd nodes_list = nodes_host.col(0);

      i_1 = fe_find(nodes_list, elements_host_tmp(i, (faces(j, 0) + 2)));
      i_2 = fe_find(nodes_list, elements_host_tmp(i, (faces(j, 1) + 2)));
      i_3 = fe_find(nodes_list, elements_host_tmp(i, (faces(j, 2) + 2)));

      VectorXd center = VectorXd::Zero(ndof);
      for (int k = 0; k < ndof; k++) {
        center(k) = ((nodes_host(i_1, k + 1) + nodes_host(i_2, k + 1) + nodes_host(i_3, k + 1)) / 3.0);
      }

      VectorXd edge_1 = VectorXd::Zero(ndof);
      VectorXd edge_2 = VectorXd::Zero(ndof);

      for (int k = 0; k < ndof; k++) {
        edge_1(k) = (nodes_host(i_2, k + 1) - nodes_host(i_1, k + 1));
        edge_2(k) = (nodes_host(i_2, k + 1) - nodes_host(i_3, k + 1));
      }

      VectorXd normal = VectorXd::Zero(ndof);
      fe_crossProduct(normal, edge_1, edge_2);

      VectorXd center_normal = VectorXd::Zero(ndof);
      center_normal << (center(0) - A(1)), (center(1) - A(2)), (center(2) - A(3));

      double angle = center_normal.dot(normal) / (normal.norm() * center_normal.norm());

      if (angle < 0) {
        break;
      }

      counter = counter + 1;
    }

    if (counter == 12) {
      host_id = elements_host_tmp(i, 0);
      break;
    }
  }

  return host_id;

}