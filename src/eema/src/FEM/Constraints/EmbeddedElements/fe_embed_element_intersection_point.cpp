#include "functions.h"
using namespace Eigen;

VectorXd fe_embed_element_intersection_point(int host_row, int embed_row, Mesh& host, Mesh& embed) {

  MatrixXd* nodes_host = host.getNewNodesPointer();
  MatrixXi* elements_host = host.getNewElementsPointer();
  MatrixXd* nodes_embed = embed.getNewNodesPointer();
  MatrixXi* elements_embed = embed.getNewElementsPointer();

  int num_elements_host = host.getNumElements();
  int num_nodes_per_element_host = host.getNumNodesPerElement();
  int num_elements_embed = embed.getNumElements();
  int num_nodes_per_element_embed = embed.getNumNodesPerElement();

  VectorXd intersection = VectorXd::Zero(ndof);

  // Extract embed node ids for embedded element.
  int embed_node_1 = (*elements_embed)(embed_row, 2);
  int embed_node_2 = (*elements_embed)(embed_row, 3);

  // Extract coordinates for embedded nodes.
  VectorXd embed_node_coord_1 = VectorXd::Zero(ndof);
  VectorXd embed_node_coord_2 = VectorXd::Zero(ndof);
  for (int j = 0; j < ndof; j++) {
    embed_node_coord_1(j) = (*nodes_embed)(embed_node_1, j + 1);
    embed_node_coord_2(j) = (*nodes_embed)(embed_node_2, j + 1);
  }

  // Define twelve triangles, based on generic host element nodes.
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

    // Loop over the twelve triangles that define the host element surface.
    // We are looking for intersection point of embedded element and host element.
    int num_faces = faces.rows();
    for (int i = 0; i < num_faces; i++) {

      // Extract generic host element node ids for triangle.
      int host_element_node_1 = faces(i, 0);
      int host_element_node_2 = faces(i, 1);
      int host_element_node_3 = faces(i, 2);

      // Extract actual host element node ids for triangle.
      int triangle_node_1 = (*elements_host)(host_row, host_element_node_1 + 2);
      int triangle_node_2 = (*elements_host)(host_row, host_element_node_2 + 2);
      int triangle_node_3 = (*elements_host)(host_row, host_element_node_3 + 2);

      // Extract coordinates for each triangle node.
      VectorXd triangle_node_coord_1 = VectorXd::Zero(ndof);
      VectorXd triangle_node_coord_2 = VectorXd::Zero(ndof);
      VectorXd triangle_node_coord_3 = VectorXd::Zero(ndof);
      for (int j = 0; j < ndof; j++) {
        triangle_node_coord_1(j) = (*nodes_host)(triangle_node_1, j + 1);
        triangle_node_coord_2(j) = (*nodes_host)(triangle_node_2, j + 1);
        triangle_node_coord_3(j) = (*nodes_host)(triangle_node_3, j + 1);
      }

      // Define two edges of triangle.
      VectorXd edge_1 = VectorXd::Zero(ndof);
      VectorXd edge_2 = VectorXd::Zero(ndof);
      for (int j = 0; j < ndof; j++) {
        edge_1(j) = triangle_node_coord_2(j) - triangle_node_coord_1(j);
        edge_2(j) = triangle_node_coord_2(j) - triangle_node_coord_3(j);
      }

      // Define outward normal of triangle.
      VectorXd normal = VectorXd::Zero(ndof);
      fe_crossProduct(normal, edge_1, edge_2);

      // The plane defined by the triangle can be described using the following equation.
      // A(x - xo) + B(y - yo) + C(z -zo) = 0
      // The embedded element can be described using the following parametric equations.
      // x = x1 + t(x2 - x1), y = y1 + t(y2 - y1), z = z1 + t(z2 - z1)
      // Solve for paramenter t, where the embedded element intersects the plane.

      double tmp1 = -normal(0)*( embed_node_coord_1(0) - triangle_node_coord_2(0) ) - normal(1)*( embed_node_coord_1(1) - triangle_node_coord_2(1) ) - normal(2)*( embed_node_coord_1(2) - triangle_node_coord_2(2) );

      double tmp2 = normal(0)*( embed_node_coord_2(0) - embed_node_coord_1(0) ) + normal(1)*( embed_node_coord_2(1) - embed_node_coord_1(1) ) + normal(2)*( embed_node_coord_2(2) - embed_node_coord_1(2) );

      // If line is parallel with plane, do not continue. Move on to next triangle.
      if (tmp2 != 0) {
        double t = tmp1/tmp2;

        // If infinite line intersects plane, but line segment does not, do not continue. Move on to next triangle.
        if (t > 0 && t < 1) {

          // Calculate the intersection point of line segment and plane.
          VectorXd intersection_tmp = VectorXd::Zero(ndof);
          intersection_tmp(0) = embed_node_coord_1(0) + t*( embed_node_coord_2(0) - embed_node_coord_1(0) );
          intersection_tmp(1) = embed_node_coord_1(1) + t*( embed_node_coord_2(1) - embed_node_coord_1(1) );
          intersection_tmp(2) = embed_node_coord_1(2) + t*( embed_node_coord_2(2) - embed_node_coord_1(2) );

          // Check if intersection point is inside triangle.
          // There are three tests to run. If any test fails, intersection point is not inside triangle. Move on to next triangle.
          VectorXd v_edge = triangle_node_coord_2 - triangle_node_coord_1;
          VectorXd v_ref = triangle_node_coord_3 - triangle_node_coord_1;
          VectorXd v_point = intersection_tmp - triangle_node_coord_1;
          VectorXd cp_point = VectorXd::Zero(ndof);
          VectorXd cp_ref = VectorXd::Zero(ndof);
          fe_crossProduct(cp_point, v_edge, v_point);
          fe_crossProduct(cp_ref, v_edge, v_ref);
          double dp = cp_point.dot(cp_ref);

          // Test 1.
          if (dp > 0 || std::abs(dp) < 1e-18) {
            v_edge = triangle_node_coord_3 - triangle_node_coord_2;
            v_ref = triangle_node_coord_1 - triangle_node_coord_2;
            v_point = intersection_tmp - triangle_node_coord_2;
            fe_crossProduct(cp_point, v_edge, v_point);
            fe_crossProduct(cp_ref, v_edge, v_ref);
            double dp = cp_point.dot(cp_ref);

            // Test 2.
            if (dp > 0 || std::abs(dp) < 1e-18) {
              v_edge = triangle_node_coord_1 - triangle_node_coord_3;
              v_ref = triangle_node_coord_2 - triangle_node_coord_3;
              v_point = intersection_tmp - triangle_node_coord_3;
              fe_crossProduct(cp_point, v_edge, v_point);
              fe_crossProduct(cp_ref, v_edge, v_ref);
              double dp = cp_point.dot(cp_ref);
              // std::cout << "dp = " << dp << '\n';

              // Test 3.
              if (dp > 0 || std::abs(dp) < 1e-18) {
                intersection = intersection_tmp;
                // We prefer an intersection with 0 < t < 1, not t = 0 or t = 1.
                // Therefore, we will break the loop here and return the current intersection point.
                break;
              }
            }
          }
        }

        if (t == 0 || t == 1) {

          // Calculate the intersection point of line segment and plane.
          VectorXd intersection_tmp = VectorXd::Zero(ndof);
          intersection_tmp(0) = embed_node_coord_1(0) + t*( embed_node_coord_2(0) - embed_node_coord_1(0) );
          intersection_tmp(1) = embed_node_coord_1(1) + t*( embed_node_coord_2(1) - embed_node_coord_1(1) );
          intersection_tmp(2) = embed_node_coord_1(2) + t*( embed_node_coord_2(2) - embed_node_coord_1(2) );

          // Check if intersection point is inside triangle.
          // There are three tests to run. If any test fails, intersection point is not inside triangle. Move on to next triangle.
          VectorXd v_edge = triangle_node_coord_2 - triangle_node_coord_1;
          VectorXd v_ref = triangle_node_coord_3 - triangle_node_coord_1;
          VectorXd v_point = intersection_tmp - triangle_node_coord_1;
          VectorXd cp_point = VectorXd::Zero(ndof);
          VectorXd cp_ref = VectorXd::Zero(ndof);
          fe_crossProduct(cp_point, v_edge, v_point);
          fe_crossProduct(cp_ref, v_edge, v_ref);
          double dp = cp_point.dot(cp_ref);
          // Test 1.
          if (dp > 0 || std::abs(dp) < 1e-18) {
            v_edge = triangle_node_coord_3 - triangle_node_coord_2;
            v_ref = triangle_node_coord_1 - triangle_node_coord_2;
            v_point = intersection_tmp - triangle_node_coord_2;
            fe_crossProduct(cp_point, v_edge, v_point);
            fe_crossProduct(cp_ref, v_edge, v_ref);
            double dp = cp_point.dot(cp_ref);
            // Test 2.
            if (dp > 0 || std::abs(dp) < 1e-18) {
              v_edge = triangle_node_coord_1 - triangle_node_coord_3;
              v_ref = triangle_node_coord_2 - triangle_node_coord_3;
              v_point = intersection_tmp - triangle_node_coord_3;
              fe_crossProduct(cp_point, v_edge, v_point);
              fe_crossProduct(cp_ref, v_edge, v_ref);
              double dp = cp_point.dot(cp_ref);
              // Test 3.
              if (dp > 0 || std::abs(dp) < 1e-18) {
                intersection = intersection_tmp;
              }
            }
          }

          // Note that we do not break here. It is good that we found an intersection point.
          // However, we prefer 0 < t < 1, not t = 0 or t = 1.
        }
      }
    }

  return intersection;

}
