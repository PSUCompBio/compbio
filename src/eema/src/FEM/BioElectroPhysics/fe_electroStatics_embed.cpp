#include "functions.h"

using namespace Eigen;

void fe_electroStatics_embed(double time, int host_id, int embed_id, VectorXi& embed_map) {

	MatrixXd* nodes_host = mesh[host_id].getNewNodesPointer();
	MatrixXi* elements_host = mesh[host_id].getNewElementsPointer();
	MatrixXd* nodes_embed = mesh[embed_id].getNewNodesPointer();
	MatrixXi* elements_embed = mesh[embed_id].getNewElementsPointer();

	int nel = mesh[host_id].getNumElements();
	int nnode = mesh[host_id].getNumNodes();
	int nnel = mesh[host_id].getNumNodesPerElement();
	int sdof = nnode;

	MatrixXd electrical_kk = MatrixXd::Zero(nnode, nnode);
	MatrixXd electrical_kk_element = MatrixXd::Zero(nnel, nnel);
	VectorXd electrical_force = VectorXd::Zero(nnode);
	VectorXd electrical_force_element = VectorXd::Zero(nnel);
	MatrixXd conductivity = MatrixXd::Zero(ndof, ndof);
	VectorXd VP = VectorXd::Zero(nnode);
	VectorXd I = VectorXd::Zero(nel);

	// Element Data
	VectorXd xcoord      = VectorXd::Zero(nnel);
	VectorXd ycoord      = VectorXd::Zero(nnel);
	VectorXd zcoord      = VectorXd::Zero(nnel);
	VectorXd xcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
	VectorXd ycoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);
	VectorXd zcoord_embed      = VectorXd::Zero((*elements_embed).cols() - 2);

	//fe_apply_bc_current(I, time);

	for (int i = 0; i < nel; i++) {

		for (int j = 0; j < nnel; j++) {
			int g = (*elements_host)(i, j + 2);
			xcoord(j)      = (*nodes_host)(g, 1);
			ycoord(j)      = (*nodes_host)(g, 2);
			zcoord(j)      = (*nodes_host)(g, 3);
		}

		//VectorXd I_element = VectorXd::Zero(nnel);
		//fe_gather_pbr(I, I_element, (*elements).block<1, 8>(i, 2), sdof);

		int nglx = 2;
		int ngly = 2;
		int nglz = 2;
		VectorXd points  = guass_points(nglx);
		VectorXd weights = guass_weights(nglx);

		VectorXd dndr(nnel);
		VectorXd dnds(nnel);
		VectorXd dndt(nnel);
		VectorXd dndx(nnel);
		VectorXd dndy(nnel);
		VectorXd dndz(nnel);
		MatrixXd jacobian(ndof, ndof);
		MatrixXd invJacobian(ndof, ndof);
		VectorXd shapes(nnel);

		fe_get_conductivity(conductivity, (*elements_host)(i, 1));

		for (int fib = 0; fib < (*elements_embed).rows(); fib++) {
			if (embed_map(fib) == (*elements_host)(i, 0)) {
				for (int j = 0; j < ((*elements_embed).cols() - 2); j++) {
					int g = (*elements_embed)(fib, j + 2);
					xcoord_embed(j)      = (*nodes_embed)(g, 1);
					ycoord_embed(j)      = (*nodes_embed)(g, 2);
					zcoord_embed(j)      = (*nodes_embed)(g, 3);
				}

				MatrixXd T = MatrixXd::Zero(3, 3);
				T = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3);

				MatrixXd fiber_conductivity = MatrixXd::Zero(3, 3);
				fiber_conductivity(0, 0) = 10;

				fiber_conductivity = T * fiber_conductivity * T.transpose();
				conductivity = conductivity + fiber_conductivity;
				I(i) = I(i) + 10;
			}
		}


		MatrixXd electrical_shape_mat = MatrixXd(nnel, ndof);

		for (int intx = 0; intx < nglx; intx++) {
			double x = points(intx);
			double wtx = weights(intx);
			for (int inty = 0; inty < ngly; inty++) {
				double y = points(inty);
				double wty = weights(inty);
				for (int intz = 0; intz < nglz; intz++) {
					double z = points(intz);
					double wtz = weights(intz);

					fe_dniso_8(dndr, dnds, dndt, x, y, z);
					jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
					double detJacobian = jacobian.determinant();
					invJacobian = jacobian.inverse();
					fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
					fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
					fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);
					fe_electrical_shapeMatrix(electrical_shape_mat, dndx, dndy, dndz);

					electrical_kk_element = electrical_kk_element + (electrical_shape_mat * conductivity * electrical_shape_mat.transpose() * wtx * wty * wtz * detJacobian);

					shapes = fe_shapes_8(x, y, z);
					electrical_force_element = electrical_force_element + (I(i) * shapes * wtx * wty * detJacobian);
				}
			}
		}

		fe_assemble_electricStiffness(electrical_kk, electrical_kk_element, (*elements_host).block<1, 8>(i, 2));
		fe_scatter_electricalForce(electrical_force, electrical_force_element, (*elements_host).block<1, 8>(i, 2));
	}

	fe_apply_bc_potential(electrical_kk, electrical_force, time);
	VP = electrical_kk.inverse() * electrical_force ;

	std::cout << "Effective conductivity: \n" << conductivity << "\n";
	std::cout << "Result Potential: " << VP(7) << "; " << VP(6) << "; " << VP(5) << "; " << VP(4) << "\n";

	mesh[0].readNodalElectroPhysics(VP);
	fe_vtuWrite(0, 0, mesh[0]);
	fe_vtuWrite(0, 0, mesh[1]);

	nodes_host = NULL;
	elements_host = NULL;
	nodes_embed = NULL;
	elements_embed = NULL;

}
