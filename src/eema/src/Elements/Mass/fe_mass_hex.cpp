#include "functions.h"
using namespace Eigen;

MatrixXd fe_mass_hex(MatrixXd nodes, VectorXi elements_row) {

	int nnel = elements_row.size() - 2;
	int nnode = nodes.rows();
	int edof = nnel * ndof;

	VectorXi nodes_local = VectorXi::Zero(nnel);
	VectorXd xcoord = VectorXd::Zero(nnel);
	VectorXd ycoord = VectorXd::Zero(nnel);
	VectorXd zcoord = VectorXd::Zero(nnel);

	double rho = fe_get_mats(elements_row(1), 0, "mechanical"); // material density

	for (int j = 0; j < nnel; j++) {
		int g = -1;
		for (int f = 0; f < nnode; f++) {
			if (elements_row(j + 2) == nodes(f, 0)) {
				g = f;
				break;
			}
		}
		nodes_local(j) = elements_row(j + 2);
		xcoord(j) = nodes(g, 1);
		ycoord(j) = nodes(g, 2);
		zcoord(j) = nodes(g, 3);
	}

	int nglx = 2;
	int ngly = 2;
	int nglz = 2;

	MatrixXd m(edof, edof); // element mass matrix
	m = MatrixXd::Zero(edof, edof);

	MatrixXd shape_mat(3, edof);

	VectorXd shapes_3d(nnel);
	// r,s,t are the isoparametric coordinates.
	VectorXd dndr(nnel);
	VectorXd dnds(nnel);
	VectorXd dndt(nnel);
	// x, y, z are the actual coordinates.
	VectorXd dndx(nnel);
	VectorXd dndy(nnel);
	VectorXd dndz(nnel);
	MatrixXd jacobian(ndof, ndof);
	MatrixXd invJacobian(ndof, ndof);

	MatrixXd points_3d = guass_points_3d(nglx, ngly, nglz);
	MatrixXd weights_3d = guass_weights_3d(ndof, nglx, ngly, nglz);

	for (int intx = 0; intx < nglx; intx++) {
		double x = points_3d(0, intx);
		double wtx = weights_3d(0, intx);
		for (int inty = 0; inty < ngly; inty++) {
			double y = points_3d(1, inty);
			double wty = weights_3d(1, inty);
			for (int intz = 0; intz < nglz; intz++) {
				double z = points_3d(2, intz);
				double wtz = weights_3d(2, intz);

				// std::cout<<"I am here\n";

				//VectorXd shapes_3d(edof);
				shapes_3d = fe_shapes_8(x, y, z);
				//VectorXd dndr(edof);
				fe_dniso_8(dndr, dnds, dndt, x, y, z);

				//MatrixXd jacobian(ndof,ndof);
				jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
				double detJacobian = jacobian.determinant();
				//MatrixXd invJacobian(ndof,ndof);
				//invJacobian = jacobian.inverse();

				//VectorXd dndx(edof);
				//dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
				//VectorXd dndy(edof);
				//dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
				//VectorXd dndz(edof);
				//dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);


				//disp_mat = fe_strDispMatrix(edof,nnel,dndx,dndy,dndz);
				//k = k + (((disp_mat.transpose())*matl_mat*disp_mat)*wtx*wty*wtz*detJacobian);
				shape_mat = fe_shapeMatrix(edof, nnel, shapes_3d);
				m = m + (rho * (shape_mat.transpose()) * shape_mat * detJacobian * wtx * wty * wtz);

				//std::cout<<k<<std::endl;
			}
		}
	}

	return m;
}

MatrixXd fe_mass_truss(double rho, double A_truss, int edof, MatrixXd nodes, MatrixXd elements) {

	/* Direct Lumped mass being used for truss elements */
	MatrixXd m = MatrixXd::Zero(edof, edof);
	double mass = rho * A_truss * 1;

	for (int i = 0; i < edof; i++) {
		for (int j = 0; j < edof; j++) {
			m(i, i) = (mass / 8);
		}
	}

	return m;
}
