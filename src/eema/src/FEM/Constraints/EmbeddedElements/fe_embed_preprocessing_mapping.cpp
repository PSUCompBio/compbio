#include "functions.h"
using namespace Eigen;

VectorXd fe_embed_preprocessing_mapping(Mesh host, Mesh embed){

    MatrixXd nodes_host = host.getNewNodes();
    MatrixXi elements_host = host.getNewElements();
    VectorXd lc_host = host.getMaxCharLength("new");
    MatrixXd nodes_embed = embed.getNewNodes();
    MatrixXi elements_embed = embed.getNewElements();

    VectorXi embedded_nodes(nodes_embed.rows());
    VectorXi embedded_elements;

    /** Temporary Matrices useful in the mapping scheme */
    MatrixXi elements_host_tmp1;
    int host_id;

    for(int i=0;i<nodes_embed.rows();i++){

        /** Creates a bounding box and seperates the mesh in the bouding box from the whole mesh */
        elements_host_tmp1 = fe_create_bbox(nodes_embed.row(i),nodes_host,elements_host,lc_host(0));
        host_id = fe_compute_host(nodes_embed.row(i),nodes_host,elements_host_tmp1);

        embedded_nodes(i) << host_id; 

    }

    return embedded_elements;

}