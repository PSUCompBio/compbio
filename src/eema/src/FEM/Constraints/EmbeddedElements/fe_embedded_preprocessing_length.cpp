#include "functions.h"
using namespace Eigen;

void fe_embed_preprocessing_length(Mesh& host, Mesh& embed) {

	/** Checks if the host mesh has zero length edges */
	host.checkMesh();
	/** Checks if the embedded mesh has zero length edges */
	embed.checkMesh();

	VectorXd lc_host = host.getMinCharLength("new");
	VectorXd lc_embed = embed.getMaxCharLength("new");

	MatrixXd* nodes_local;
	MatrixXi* elements_local;

	while (lc_embed(0) > lc_host(0)) {

		nodes_local = embed.getNewNodesPointer();
		elements_local = embed.getNewElementsPointer();

		for (int i = 0; i < (*elements_local).rows(); i++) {

			if (i == (int)lc_embed(1)) {
				int i_1 = fe_find((*nodes_local).col(0), (*elements_local)(i, 2));
				int i_2 = fe_find((*nodes_local).col(0), (*elements_local)(i, 3));

				VectorXd mid_point;
				mid_point << (((*nodes_local)(i_1, 1) + (*nodes_local)(i_2, 1)) / 2) , (((*nodes_local)(i_1, 2) + (*nodes_local)(i_2, 2)) / 2), (((*nodes_local)(i_1, 3) + (*nodes_local)(i_2, 3)) / 2);

				int new_node_num = (*nodes_local)(((*nodes_local).rows() - 1), 0) + 1 ;
				VectorXd new_node;
				new_node << new_node_num, mid_point(0), mid_point(1), mid_point(2);
				fe_concatenate_vector2matrix((*nodes_local), new_node, 1);
				(*elements_local)(i, 3) = new_node_num;

				VectorXd new_element;
				new_element << 0, new_node_num, i_2;
				fe_insert_vector2matrix((*elements_local), new_element, i, 1);

				embed.replaceElements((*elements_local), "old");
				embed.replaceNodes((*nodes_local), "old");
				embed.preprocessMesh();

				lc_embed = embed.getMaxCharLength("new");

			}

		}

	}

}