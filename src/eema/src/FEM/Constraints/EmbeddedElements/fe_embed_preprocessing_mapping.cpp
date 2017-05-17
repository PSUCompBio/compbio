#include "functions.h"
using namespace Eigen;

void fe_embed_preprocessing_mapping(VectorXi& host_elements_embed_nodes, Mesh& host, Mesh& embed) {

    MatrixXd* nodes_host = host.getNewNodesPointer();
    MatrixXi* elements_host = host.getNewElementsPointer();
    MatrixXd* nodes_embed = embed.getNewNodesPointer();
    MatrixXi* elements_embed = embed.getNewElementsPointer();

    VectorXd lc_host = host.getMaxCharLength("new");

    /** Temporary Matrices useful in the mapping scheme */
    MatrixXi elements_host_tmp;
    int host_id;

    for (int i = 0; i < (*nodes_embed).rows(); i++) {
        /** Creates a bounding box and seperates the mesh in the bouding box from the whole mesh */

        VectorXd node_info = (*nodes_embed).row(i);
        elements_host_tmp = fe_create_bbox(node_info, (*nodes_host), (*elements_host), lc_host(0));

        host_id = fe_compute_host(node_info, (*nodes_host), elements_host_tmp);
        host_elements_embed_nodes(i) = host_id;
    }

}