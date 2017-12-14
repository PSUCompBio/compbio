#include "functions.h"

using namespace Eigen;

/** This function carries out all the required pre-processing steps for embedded element constraint */

VectorXi fe_embed_preprocessing(Mesh& host, Mesh& embed) {

    fe_embed_preprocessing_length(host, embed);

    VectorXi host_elements_embed_nodes = VectorXi::Zero(embed.getNumNodes());
    VectorXi host_elements_embed_elements;
    fe_embed_preprocessing_mapping(host_elements_embed_nodes, host, embed);


    MatrixXd* nodes_host = host.getNewNodesPointer();
    MatrixXi* elements_host = host.getNewElementsPointer();
    MatrixXd* nodes_embed = embed.getNewNodesPointer();
    MatrixXi elements_embed = embed.getNewElements();

    int num_embed_elements = embed.getNumElements();
    int num_nodes_per_embed_element = embed.getNumNodesPerElement();

    for (int i = 0; i < num_embed_elements; i++) {

        elements_embed = embed.getNewElements();

        VectorXd node_list = (*nodes_embed).col(0);

        int i_1 = fe_find(node_list, elements_embed(i, 2));
        int i_2 = fe_find(node_list, elements_embed(i, 3));

        if (host_elements_embed_nodes(i_1) == host_elements_embed_nodes(i_2)) {
            fe_insert_int2vector(host_elements_embed_elements, host_elements_embed_nodes(i_1));
        }
        else {

            int id_1 = host_elements_embed_nodes(i_1);
            int id_2 = host_elements_embed_nodes(i_2);

            VectorXi common_face;
            bool face_value;

            fe_find_common_face(common_face, face_value, id_1, id_2, *nodes_host, *elements_host);

            VectorXd intersection;
            fe_find_intersection(intersection, common_face, *nodes_embed, *nodes_host, i_1, i_2);

            int new_node_num = (*nodes_embed)(((*nodes_embed).rows() - 1), 0) + 1 ;
            VectorXd new_node = VectorXd::Zero(ndof + 1);
            new_node(0) = new_node_num;
            new_node(1) = intersection(0);
            new_node(2) = intersection(1);
            new_node(3) = intersection(2);
            embed.append2nodes("new", new_node);

            elements_embed(i, 3) = new_node_num;
            embed.replaceElements(elements_embed, "new");
            fe_insert_int2vector(host_elements_embed_elements, i, id_1);

            VectorXi new_element = VectorXi::Zero(elements_embed.cols());
            new_element(0) = 0;
            new_element(1) = elements_embed(i, 1);
            new_element(2) = new_node_num;
            new_element(3) = i_2;
            fe_insert_int2vector(host_elements_embed_elements, id_2);
            embed.append2elements("new", new_element);
            embed.reNumber();
        }
    }

    return host_elements_embed_elements;

}



