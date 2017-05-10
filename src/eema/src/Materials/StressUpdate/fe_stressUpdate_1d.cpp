#include"functions.h"

using namespace Eigen;

/** This function calculates the updated stress for 1d elements - hyperelastic material model was implemented so far. */
VectorXd fe_stressUpdate_1d(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u_e, int opt, MatrixXd nodes){

	MatrixXd F = MatrixXd::Zero(3,3); // deformation gradient
	F = fe_calDefGrad(dndx,dndy,dndz,u_e);
	double defJacobian = F.determinant();

	MatrixXd C = F.transpose()*F;

	VectorXd dir_truss = VectorXd::Zero(3);
	dir_truss << (nodes(1,1)-nodes(0,1)) , (nodes(1,2)-nodes(0,2)), (nodes(1,3)-nodes(0,3));
	VectorXd deformed_truss = C*dir_truss;

	double lambda_tmp = dir_truss.dot(deformed_truss);
	double lambda = sqrt(lambda_tmp);
	// Until this point, deformation along the truss direction was calculated. now its time to calculate the stress in this dirn.

	VectorXd sigma = VectorXd::Zero(6); //piola-kirchhoff stress-vector
	VectorXd cauchy_sigma = VectorXd::Zero(6);

	std::string model;
	model = fe_get_model(opt);

	if(model=="simple_elastic"){ // Isotropic truss element
		double E = fe_get_mats(opt,1);
		cauchy_sigma(0) = E*log(lambda);
		// std::cout<<"Stretch is: "<<lambda<<"\n";
	}

	if(model=="mooney-rivlin_hyperelastic"){
		double c_1 = fe_get_mats(opt,3);
		double c_2 = fe_get_mats(opt,4);
		cauchy_sigma(0) = ((2*c_1) + ((2*c_2)/(lambda)))*((pow(lambda,2))-(1/lambda));
	}

	MatrixXd stress_transform_mat = fe_calTransformation(nodes,1); // Stress transformation from truss CSYS to element CSYS.


	VectorXd tmp_cauchy_sigma = stress_transform_mat.transpose()*cauchy_sigma; // Cauchy stress vector in element CSYS.

	MatrixXd element_cauchy_sigma = MatrixXd::Zero(3,3);
	element_cauchy_sigma(0,0) = tmp_cauchy_sigma(0);
	element_cauchy_sigma(1,1) = tmp_cauchy_sigma(1);
	element_cauchy_sigma(2,2) = tmp_cauchy_sigma(2);
	element_cauchy_sigma(0,1) = tmp_cauchy_sigma(3);
	element_cauchy_sigma(1,0) = tmp_cauchy_sigma(3);
	element_cauchy_sigma(1,2) = tmp_cauchy_sigma(4);
	element_cauchy_sigma(2,1) = tmp_cauchy_sigma(4);
	element_cauchy_sigma(0,2) = tmp_cauchy_sigma(5);
	element_cauchy_sigma(2,0) = tmp_cauchy_sigma(5);

	// fe_display_matrix(element_cauchy_sigma);


	MatrixXd F_inv = F.inverse();
	MatrixXd pk_S = MatrixXd::Zero(3,3);
	pk_S = defJacobian*F_inv*element_cauchy_sigma*F_inv.transpose();

	sigma(0) = pk_S(0,0);
	sigma(1) = pk_S(1,1);
	sigma(2) = pk_S(2,2);
	sigma(3) = pk_S(0,1);
	sigma(4) = pk_S(1,2);
	sigma(5) = pk_S(0,2);

	/*sigma(0) = element_cauchy_sigma(0,0);
	sigma(1) = element_cauchy_sigma(1,1);
	sigma(2) = element_cauchy_sigma(2,2);
	sigma(3) = element_cauchy_sigma(0,1);
	sigma(4) = element_cauchy_sigma(1,2);
	sigma(5) = element_cauchy_sigma(0,2);*/

	return sigma;
}
