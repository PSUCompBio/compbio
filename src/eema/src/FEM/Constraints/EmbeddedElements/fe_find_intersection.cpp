#include "functions.h"
using namespace Eigen;

void fe_find_intersection(VectorXd& intersection, VectorXi& common_face, MatrixXd& nodes_embed, MatrixXd& nodes_host, int i_1, int i_2) {

	VectorXi ids = VectorXi::Zero(common_face.size());

	VectorXd p_0 = VectorXd::Zero(ndof);

	for (int i = 0; i < common_face.size(); i++) {
		VectorXd node_list = nodes_embed.col(0);
		ids(i) = fe_find(node_list, common_face(i));
	}

	MatrixXi traingles(2, 3);
	traingles << 0, 1, 2,
	          2, 3, 0;

	MatrixXd normals(2, 3);
	VectorXd final_normal;

	for (int i = 0; i < traingles.rows(); i++) {
		VectorXd coord1 = VectorXd::Zero(ndof);
		VectorXd coord2 = VectorXd::Zero(ndof);
		VectorXd coord3 = VectorXd::Zero(ndof);
		VectorXd node_list = nodes_host.col(0);

		for (int j = 0; j < ndof; j++) {
			int tmp = fe_find(node_list, common_face(traingles(i, 0)));
			coord1(j) = nodes_host(tmp, j + 1);
			tmp = fe_find(node_list, common_face(traingles(i, 1)));
			coord2(j) = nodes_host(tmp, j + 1);
			tmp = fe_find(node_list, common_face(traingles(i, 2)));
			coord3(j) = nodes_host(tmp, j + 1);
			if (i == 0) {
				p_0(j) = coord1(j);
			}
		}
		normals.row(i) = fe_find_normals(coord1, coord2, coord3);
	}

	if (normals.row(0) == normals.row(1)) {
		final_normal = normals.row(0);
	}
	else {
		final_normal = 0.5 * (normals.row(0) + normals.row(1) );
	}

	final_normal = final_normal / final_normal.norm();

	VectorXd l_0 = VectorXd::Zero(ndof);
	VectorXd l_1 = VectorXd::Zero(ndof);

	for (int k = 0; k < ndof; k++) {
		l_0(k) = nodes_embed(i_1, k + 1);
		l_1(k) = nodes_embed(i_2, k + 1);
	}

	VectorXd l = l_1 - l_0 ;

	double d = ((p_0 - l_0).dot(final_normal)) / (l.dot(final_normal));

	intersection = l_0 + (d * l);

}

Vector3d fe_find_normals(VectorXd& c1, VectorXd& c2, VectorXd& c3) {

	Vector3d edge_1;
	edge_1 << (c1(0) - c2(0)), (c1(1) - c2(1)), (c1(2) - c2(2));
	Vector3d edge_2;
	edge_2 << (c1(0) - c3(0)), (c1(1) - c3(1)), (c1(2) - c3(2));

	Vector3d normal = edge_1.cross(edge_2);

	return normal;

}