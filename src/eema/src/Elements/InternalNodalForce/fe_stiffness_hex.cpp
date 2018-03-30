

#include"functions.h"

using namespace Eigen;

MatrixXd fe_stiffness_hex(double E, double nu, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]){

	int nglx = 2;
	int ngly = 2;
	int nglz = 2;

	MatrixXd k(edof,edof); // element stiffness matrix
	k = MatrixXd::Zero(edof,edof);

	MatrixXd disp_mat(6,edof);

	MatrixXd matl_mat(6,6);
	matl_mat = fe_isoElastic(4,E,nu);

	VectorXd shapes_3d(edof);
	VectorXd dndr(edof);
	VectorXd dnds(edof);
	VectorXd dndt(edof);
	VectorXd dndx(edof);
	VectorXd dndy(edof);
	VectorXd dndz(edof);
	MatrixXd jacobian(edof,edof);
	MatrixXd invJacobian(edof,edof);

	MatrixXd points_3d = guass_points_3d(nglx,ngly,nglz);
	MatrixXd weights_3d = guass_weights_3d(ndof,nglx,ngly,nglz);

	for(int intx=0;intx<nglx;intx++){
		double x = points_3d(0,intx);
		double wtx = weights_3d(0,intx);
		for(int inty=0;inty<ngly;inty++){
			double y = points_3d(1,inty);
			double wty = weights_3d(1,inty);
			for(int intz=0;intz<nglz;intz++){
				double z = points_3d(2,intz);
				double wtz = weights_3d(2,intz);

				std::cout<<"I am here\n";

				//VectorXd shapes_3d(edof);
				shapes_3d = fe_shapes_8(x,y,z);
				//VectorXd dndr(edof);
				dndr = fe_dndr_8(x,y,z);
				//VectorXd dnds(edof);
				dnds = fe_dnds_8(x,y,z);
				//VectorXd dndt(edof);
				dndt = fe_dndt_8(x,y,z);

				//MatrixXd jacobian(ndof,ndof);
				jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
				double detJacobian = jacobian.determinant();
				//MatrixXd invJacobian(ndof,ndof);
				invJacobian = jacobian.inverse();

				//VectorXd dndx(edof);
				fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
				//VectorXd dndy(edof);
				fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
				//VectorXd dndz(edof);
				fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);


				disp_mat = fe_strDispMatrix(edof,nnel,dndx,dndy,dndz);
				k = k + (((disp_mat.transpose())*matl_mat*disp_mat)*wtx*wty*wtz*detJacobian);
				//std::cout<<k<<std::endl;
			}
		}
	}

	return k;
}
