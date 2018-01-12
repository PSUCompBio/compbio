#include "functions.h"

using namespace Eigen;

/** This function carries out all the required pre-processing steps for embedded element constraint */

VectorXi fe_embed_preprocessing(Mesh& host, Mesh& embed) {

    // Retrieve number of embedded elements.
    int num_elements_embed = embed.getNumElements();

    // Define array to store host id for each embedded element.
    VectorXi embed_map = VectorXi::Zero(num_elements_embed);

    // Run sub-routine to assign host ids to embedded elements.
    // During process, embedded elements are split at host interface and new embedded elements are generated.
    // As a result, some new embedded elements are not assigned host ids on first pass of sub-routine.
    // Therefore, sub-routine may need to run multiple times.
    embed_map = fe_embed_preprocessing_host_map(embed_map, host, embed);

    int embed_map_size = embed_map.size();
    int host_id_check = 0;

    // Check if any embedded elements have not been assigned a host id.
    for (int i = 0; i < embed_map_size; i++) {
      if (embed_map(i) == 0) {
        host_id_check = 1;
        break;
      }
    }

    int counter = 0;

    // Run sub-routine until all embedded elements have been assigned host ids.
    while (host_id_check == 1) {

      std::cout << "counter = " << counter << '\n';

      embed_map = fe_embed_preprocessing_host_map(embed_map, host, embed);
      embed_map_size = embed_map.size();

      host_id_check = 0;
      for (int i = 0; i < embed_map_size; i++) {
        if (embed_map(i) == 0) {
          host_id_check = 1;

          break;
        }
      }

      counter = counter + 1;

    }

    return embed_map;

}



// VectorXi fe_embed_preprocessing(Mesh& host, Mesh& embed) {
//
//     fe_embed_preprocessing_length(host, embed);
//
//     VectorXi host_elements_embed_nodes = VectorXi::Zero(embed.getNumNodes());
//     VectorXi host_elements_embed_elements;
//     fe_embed_preprocessing_mapping(host_elements_embed_nodes, host, embed);
//
//
//     MatrixXd* nodes_host = host.getNewNodesPointer();
//     MatrixXi* elements_host = host.getNewElementsPointer();
//     MatrixXd* nodes_embed = embed.getNewNodesPointer();
//     MatrixXi elements_embed = embed.getNewElements();
//
//     int num_embed_elements = embed.getNumElements();
//     int num_nodes_per_embed_element = embed.getNumNodesPerElement();
//
//     for (int i = 0; i < num_embed_elements; i++) {
//
//         elements_embed = embed.getNewElements();
//
//         VectorXd node_list = (*nodes_embed).col(0);
//
//         int i_1 = fe_find(node_list, elements_embed(i, 2));
//         int i_2 = fe_find(node_list, elements_embed(i, 3));
//
//         if (host_elements_embed_nodes(i_1) == host_elements_embed_nodes(i_2)) {
//             fe_insert_int2vector(host_elements_embed_elements, i, host_elements_embed_nodes(i_1));
//         }
//         else {
//
//             int id_1 = host_elements_embed_nodes(i_1);
//             int id_2 = host_elements_embed_nodes(i_2);
//
//             VectorXi common_face;
//             bool face_value;
//
//             fe_find_common_face(common_face, face_value, id_1, id_2, *nodes_host, *elements_host);
//
//             VectorXd intersection;
//             fe_find_intersection(intersection, common_face, *nodes_embed, *nodes_host, i_1, i_2);
//
//             int new_node_num = (*nodes_embed)(((*nodes_embed).rows() - 1), 0) + 1 ;
//             VectorXd new_node = VectorXd::Zero(ndof + 1);
//             new_node(0) = new_node_num;
//             new_node(1) = intersection(0);
//             new_node(2) = intersection(1);
//             new_node(3) = intersection(2);
//
//             embed.append2nodes("new", new_node);
//
//             elements_embed(i, 3) = new_node_num;
//             embed.replaceElements(elements_embed, "new");
//
//             fe_insert_int2vector(host_elements_embed_elements, i, id_1);
//
//             VectorXi new_element = VectorXi::Zero(elements_embed.cols());
//             new_element(0) = 0;
//             new_element(1) = elements_embed(i, 1);
//             new_element(2) = new_node_num;
//             new_element(3) = i_2;
//
//             embed.append2elements("new", new_element);
//             embed.reNumber();
//
//             int new_element_num = embed.getNumElements() - 1;
//
//             fe_insert_int2vector(host_elements_embed_elements, new_element_num, id_2);
//         }
//     }
//
//     return host_elements_embed_elements;
//
// }


