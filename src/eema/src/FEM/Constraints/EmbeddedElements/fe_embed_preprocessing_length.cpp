#include "functions.h"
using namespace Eigen;

void fe_embed_preprocessing_length(Mesh& host, Mesh& embed) {

	/** Checks if the host mesh has zero length edges */
	host.checkMesh();
	/** Checks if the embedded mesh has zero length edges */
	embed.checkMesh();

	VectorXd lc_host = host.getMinCharLength("new");
	VectorXd lc_embed = embed.getMaxCharLength("new");

	MatrixXd nodes_local;
	MatrixXi elements_local;

	while (lc_embed(0) > lc_host(0)) {

		nodes_local = embed.getNewNodes();
		elements_local = embed.getNewElements();

		for (int i = 0; i < pow(10, 9) ; i++) {

			if ((i == (int)lc_embed(1)) && (lc_embed(0) > lc_host(0))) {

				VectorXd node_list = nodes_local.col(0);
				int i_1 = fe_find(node_list, elements_local(i, 2));
				int i_2 = fe_find(node_list, elements_local(i, 3));

				VectorXd mid_point = VectorXd::Zero(ndof);
				for (int k = 0; k < ndof; k++) {
					mid_point(k) = ( nodes_local(i_1, k + 1) + nodes_local(i_2, k + 1) ) * 0.5;
				}

				double new_node_num = nodes_local((nodes_local.rows() - 1), 0) + 1 ;
				VectorXd new_node = VectorXd::Zero(ndof + 1);
				new_node << new_node_num, mid_point(0), mid_point(1), mid_point(2);

				fe_concatenate_vector2matrix(nodes_local, new_node, 1);

				elements_local(i, 3) = new_node_num;

				VectorXi new_element = VectorXi::Zero(elements_local.cols());
				new_element << 0, elements_local(i, 1), new_node_num, i_2;
				fe_concatenate_vector2matrix(elements_local, new_element, 1);

				embed.replaceElements(elements_local, "new");
				embed.replaceNodes(nodes_local, "new");
				embed.reNumber();

				lc_embed = embed.getMaxCharLength("new");

			}

			if (i == embed.getNumElements()) {
				break;
			}

		}

	}

}