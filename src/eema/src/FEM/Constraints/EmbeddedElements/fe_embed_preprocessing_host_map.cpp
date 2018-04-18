#include "functions.h"

using namespace Eigen;

VectorXi fe_embed_preprocessing_host_map(VectorXi embed_map, Mesh& host, Mesh& embed) {

    MatrixXd* nodes_host = host.getNewNodesPointer();
    MatrixXi* elements_host = host.getNewElementsPointer();
    MatrixXd* nodes_embed = embed.getNewNodesPointer();
    MatrixXi* elements_embed = embed.getNewElementsPointer();

    int num_elements_host = host.getNumElements();
    int num_nodes_per_element_host = host.getNumNodesPerElement();
    int num_elements_embed = embed.getNumElements();
    int num_nodes_embed = embed.getNumNodes();
    int num_nodes_per_element_embed = embed.getNumNodesPerElement();

    // Loop over host elements.
    for (int host_row = 0; host_row < num_elements_host; host_row++) {

      // Loop over embedded elements.
      for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {

        // Check if embedded element already has host id.
        // If so, move on to next embedded element.
        if (embed_map(embed_row) == 0) {

          // Extract node ids for the embedded element.
          int embed_node_1 = (*elements_embed)(embed_row, 2);
          int embed_node_2 = (*elements_embed)(embed_row, 3);

          // Check if embedded element nodes are inside host element.
          // Being on boundary of host element counts as being inside host element.
          int embed_node_check_1 = fe_embed_node_host_check(host_row, embed_node_1, host, embed);
          int embed_node_check_2 = fe_embed_node_host_check(host_row, embed_node_2, host, embed);

          // Check if both embedded element nodes are inside host element.
          // If so, assign host element id to embedded element. Then, move on to next embedded element.
          if (embed_node_check_1 + embed_node_check_2 == 2) {
            embed_map(embed_row) = host_row + 1;
          }

          // Check if only one embedded element node is inside host element.
          if (embed_node_check_1 + embed_node_check_2 == 1) {

            // Extract coordinates for embedded nodes.
            VectorXd embed_node_coord_1 = VectorXd::Zero(ndof);
            VectorXd embed_node_coord_2 = VectorXd::Zero(ndof);
            for (int j = 0; j < ndof; j++) {
              embed_node_coord_1(j) = (*nodes_embed)(embed_node_1, j + 1);
              embed_node_coord_2(j) = (*nodes_embed)(embed_node_2, j + 1);
            }

            // Extract material id for embedded element.
            int mat_id_embed = (*elements_embed)(embed_row,1);

            // Find intersection of embedded element and host element.
            int intersection_check = 0;
            VectorXd intersection = VectorXd::Zero(ndof);
            intersection = fe_embed_element_intersection_point(host_row, embed_row, host, embed, intersection_check);

            // Ensure that intersection point is not coincident with either embedded element nodes.
            // If not coincident, split embedded element at intersection point.
            if (embed_node_coord_1 != intersection && embed_node_coord_2 != intersection && intersection_check == 1) {

              // Create new embedded element node row. Append new row to embedded nodes matrix.
              VectorXd new_embed_node_row = VectorXd::Zero(4);
              new_embed_node_row << num_nodes_embed, intersection(0), intersection(1), intersection(2);
              embed.append2nodes("new", new_embed_node_row);

              // Create new embedded element element row.
              VectorXi new_embed_element_row = VectorXi::Zero(4);

              if (embed_node_check_1 == 1) {
                (*elements_embed)(embed_row, 3) = num_nodes_embed;
                new_embed_element_row << num_elements_embed + 1, mat_id_embed, num_nodes_embed, embed_node_2;
                embed.append2elements("new", new_embed_element_row);
              }

              if (embed_node_check_2 == 1) {
                (*elements_embed)(embed_row, 2) = num_nodes_embed;
                new_embed_element_row << num_elements_embed + 1, mat_id_embed, embed_node_1, num_nodes_embed;
                embed.append2elements("new", new_embed_element_row);
              }

              // Define host id for embedded element. Increase length of host id array by one.
              embed_map(embed_row) = host_row + 1;
              fe_insert_int2vector(embed_map, 0);

              // Update number of embedded nodes and elements.
              num_elements_embed = embed.getNumElements();
              num_nodes_embed = embed.getNumNodes();

            }
          }
        }
      }
    }

  return embed_map;

}
