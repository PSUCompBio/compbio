
#include"functions.h"

using namespace Eigen;

MatrixXd fe_stiffness_embed_truss(MatrixXd nodes_truss, MatrixXd elements_truss, double E_truss, double A_truss, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]){

	MatrixXd disp_mat(6,edof);

	VectorXd shapes_3d(edof);
	VectorXd dndr(edof);
	VectorXd dnds(edof);
	VectorXd dndt(edof);
	VectorXd dndx(edof);
	VectorXd dndy(edof);
	VectorXd dndz(edof);
	MatrixXd jacobian(edof,edof);
	MatrixXd invJacobian(edof,edof);

	MatrixXd strain_transform_mat = fe_calTransformation(nodes_truss,2);

	MatrixXd k_truss(edof,edof);
	k_truss = MatrixXd::Zero(edof,edof);

	MatrixXd k_correction(edof,edof);
	k_correction = MatrixXd::Zero(edof,edof);

	double length = 1;
	int truss_intg = 2;
	MatrixXd nodes_intg_global(2,3);
	VectorXd nodes_intg_local(3);
	//VectorXd truss_origin(3); // for local to global transformation
	//truss_origin << nodes_truss(0,1), nodes_truss(0,2), nodes_truss(0,3);
	VectorXd truss_intg_points = guass_points(2);
	VectorXd truss_intg_weights = guass_weights(2);


/// @image html Slide1.png
/// @image latex Slide1.eps "Image caption!" width=\textwidth

	for(int i=0;i<truss_intg;i++){
		double t = truss_intg_points(i);
		//std::cout<<t<<"is t\n";
		double wtt = truss_intg_weights(i);
		//std::cout<<wtt<<"is the weight\n";
		nodes_intg_local << t, 0, 0;
		MatrixXd tmp(3,3);
		tmp = fe_calSimpTransformation(nodes_truss);
		//nodes_intg_global.row(i) = truss_origin + (tmp.transpose()*nodes_intg_local);
		nodes_intg_global.row(i) = (tmp.transpose()*nodes_intg_local);
		//nodes_intg_global.row(i) =

		// Calculating the Parent Strain Displacement Matrix at embedded node positions
		//dndr = fe_dndr_8(nodes_intg_global(i,0),nodes_intg_global(i,1),nodes_intg_global(i,2));
		//dnds = fe_dnds_8(nodes_intg_global(i,0),nodes_intg_global(i,1),nodes_intg_global(i,2));
		//dndt = fe_dndt_8(nodes_intg_global(i,0),nodes_intg_global(i,1),nodes_intg_global(i,2));
		dndr = fe_dndr_8(0,0,t);
		//std::cout<<"DNDR is: "<<dndr<<"\n";
		dnds = fe_dnds_8(0,0,t);
		dndt = fe_dndt_8(0,0,t);
		jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
		//std::cout<<"JACOBIAN: "<<jacobian<<"\n";
                invJacobian = jacobian.inverse();
		//std::cout<<"Inverse Jacobian: "<<invJacobian<<"\n";
				fe_dndx_8_pbr(dndx, nnel, dndr, dnds, dndt, invJacobian);
                fe_dndy_8_pbr(dndy, nnel, dndr, dnds, dndt, invJacobian);
                fe_dndz_8_pbr(dndz, nnel, dndr, dnds, dndt, invJacobian);
                disp_mat = fe_strDispMatrix(edof,nnel,dndx,dndy,dndz);
		//std::cout<<"Displacement Matrix is: "<<disp_mat.transpose()<<"\n";
		//std::cout<<"Size of disp_mat is: "<<disp_mat.transpose().rows()<<"&"<<disp_mat.transpose().cols()<<"\n";
		//std::cout<<"Size of strain_mat is: "<<strain_transform_mat.transpose().rows()<<"&"<<strain_transform_mat.transpose().cols()<<"\n";
		//std::cout<<"Size of strain_mat is: "<<strain_transform_mat.rows()<<"&"<<strain_transform_mat.cols()<<"\n";

		k_truss = k_truss + ((A_truss)*(((disp_mat.transpose()*strain_transform_mat.transpose())*E_truss*strain_transform_mat)*disp_mat)*(length/2)*(wtt));
		//k_correction = k_correction + (A_truss*E*(((disp_mat.transpose()*strain_transform_mat.transpose())*strain_transform_mat)*disp_mat)*(length/2)*wtt);
	}
	std::string truss_stiffness_file = "result_trussStiffness.txt";
	matrix2text(truss_stiffness_file.c_str(),k_truss,0);

	return k_truss;
}
