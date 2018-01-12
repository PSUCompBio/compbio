#include "functions.h"
using namespace Eigen;

int fe_embed_node_host_check(int host_row, int embed_node, Mesh& host, Mesh& embed) {

  MatrixXd* nodes_host = host.getNewNodesPointer();
  MatrixXi* elements_host = host.getNewElementsPointer();
  MatrixXd* nodes_embed = embed.getNewNodesPointer();
  MatrixXi* elements_embed = embed.getNewElementsPointer();

  int num_elements_host = host.getNumElements();
  int num_nodes_per_element_host = host.getNumNodesPerElement();
  int num_elements_embed = embed.getNumElements();
  int num_nodes_per_element_embed = embed.getNumNodesPerElement();

  // embed_node_check = 0 if embedded node IS NOT inside the host element
  // embed_node_check = 1 if embedded node IS inside the host element
  int embed_node_check = 0;

  // Extract coordinates for embedded node.
  VectorXd embed_node_coord = VectorXd::Zero(ndof);
  for (int j = 0; j < ndof; j++) {
    embed_node_coord(j) = (*nodes_embed)(embed_node, j + 1);
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

    // Define counter to count how many triangles embedded node is behind, relative to triangle outward normal.
    // If embedded node is behind all twelve triangles, the embedded node is inside the host element.
    int counter = 0;

    // Loop over the twelve triangles that define the host element surface.
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

      // Define centroid of triangle.
      VectorXd center = VectorXd::Zero(ndof);
      for (int j = 0; j < ndof; j++) {
        center(j) = ( triangle_node_coord_1(j) + triangle_node_coord_2(j) + triangle_node_coord_3(j) ) / 3.0;
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

      // Define vector from embedded node to triangle centroid.
      VectorXd center_normal = VectorXd::Zero(ndof);
      center_normal << ( center(0) - embed_node_coord(0) ), ( center(1) - embed_node_coord(1) ), ( center(2) - embed_node_coord(2) );

      // Calculate the cosine of the angle between center_normal and normal.
      double angle = center_normal.dot(normal) / (normal.norm() * center_normal.norm());

      // This means the embedded node is outside of the triangle.
      // Note that if angle = 0, the embedded node is on the triangle surface.
      // Note that there were numerical issues when we used "angle < 0".
      // Embedded nodes on the host element surface are considered to be inside the host element.
      if (angle < -1e-9) {
        embed_node_check = 0;
        break;
      }

      counter = counter + 1;
    }

    // This means the embedded node is inside the host element.
    if (counter == 12) {
      embed_node_check = 1;
    }

  return embed_node_check;

}