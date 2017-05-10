#include "functions.h"

using namespace Eigen;

/*! \brief Calculates the resultant nodal force after each time step.
*
* This function represents the 'getforce' step in Belytschko (Box 6.1 - Explicit FEM Algorithm).
*	For each hex element, this function calculates the internal nodal force vector
* and the resultant nodal force vector. Once, this is calculated for each element,
* the resultant vectors are scattered into global vectors.
*/

VectorXd fe_getforce(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext, int size_counter, MatrixXd nodes_truss, MatrixXi elements_truss){

	int nel = elements.rows();
	int nnel = (elements.cols()-2);
	int nnode = nodes.rows();
	int sdof = nnode * ndof;
	int edof = nnel * ndof;

	MatrixXd element_stress_host = MatrixXd::Zero(nel,9);
	MatrixXd element_strain_host = MatrixXd::Zero(nel,9);

	// node_IntForce_system[size_counter] = MatrixXd::Zero(nel,nnel);

	VectorXi nodes_local = VectorXi::Zero(nnel);

	VectorXd xcoord = VectorXd::Zero(nnel);
	VectorXd ycoord = VectorXd::Zero(nnel);
	VectorXd zcoord = VectorXd::Zero(nnel);

	VectorXd f_tot(sdof);
	f_tot = VectorXd::Zero(sdof);

	MatrixXd sigma_all = MatrixXd::Zero(nel,8); // Six Stress Values and 1 element number and 1 time

	for(int i=0;i<nel;i++){

		for(int j=0;j<nnel;j++){
			int g = -1;
			for(int f=0;f<nnode;f++){
				if(elements(i,j+2)==nodes(f,0)){
					g = f;
					break;
				}
			}
			nodes_local(j) = elements(i,j+2);
			xcoord(j) = nodes(g,1);
			ycoord(j) = nodes(g,2);
			zcoord(j) = nodes(g,3);
		}


    //std::cout << "Nodes Local: " << nodes_local << "\n";

	VectorXd u_e = VectorXd::Zero(edof); // element displacements
	u_e = fe_gather(u,u_e,nodes_local,sdof);

	/*int disp_counter = 0;
	for(int j=0;j<xcoord.size();j++){
		xcoord(j) = xcoord(j) + u_e(disp_counter);
		ycoord(j) = ycoord(j) + u_e(disp_counter+1);
		zcoord(j) = zcoord(j) + u_e(disp_counter+2);
		disp_counter = j*3;
	}*/

	VectorXd f_ext_e = VectorXd::Zero(edof);
	f_ext_e = fe_gather(fext,f_ext_e,nodes_local,sdof); // element external nodal forces

	VectorXd f_int_e = VectorXd::Zero(edof);
	f_int_e = VectorXd::Zero(edof); // element internal nodal force

	VectorXd f_tot_e = VectorXd::Zero(edof);
	f_tot_e = VectorXd::Zero(edof); // element total nodal force

	int nglx = 2;
	int ngly = 2;
	int nglz = 2;

	MatrixXd disp_mat(6,edof);

	VectorXd dndr(nnel);
	VectorXd dnds(nnel);
	VectorXd dndt(nnel);
	VectorXd dndx(nnel);
	VectorXd dndy(nnel);
	VectorXd dndz(nnel);
	MatrixXd jacobian(ndof,ndof);
	MatrixXd invJacobian(ndof,ndof);

	MatrixXd points_3d = guass_points_3d(nglx,ngly,nglz);
	MatrixXd weights_3d = guass_weights_3d(ndof,nglx,ngly,nglz);

    if(size_counter!=0){ // if this is not the first time step the go into the loop

	for(int intx=0;intx<nglx;intx++){
		double x = points_3d(0,intx);
		double wtx = weights_3d(0,intx);
		for(int inty=0;inty<ngly;inty++){
			double y = points_3d(1,inty);
			double wty = weights_3d(1,inty);
			for(int intz=0;intz<nglz;intz++){
				double z = points_3d(2,intz);
				double wtz = weights_3d(2,intz);

				//VectorXd dndr(edof);
				dndr = fe_dndr_8(x,y,z);
				//VectorXd dnds(edof);
				dnds = fe_dnds_8(x,y,z);
				//VectorXd dndt(edof);
				dndt = fe_dndt_8(x,y,z);

				jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
				double detJacobian = jacobian.determinant();
				invJacobian = jacobian.inverse();

				dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
				dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
				dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);

				disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);
				// disp_mat = fe_strDispMatrix(edof,nnel,dndx,dndy,dndz);

				VectorXd sigma_e = VectorXd::Zero(6);
				sigma_e = fe_stressUpdate(dndx,dndy,dndz,disp_mat,u_e,elements(i,1),0);

				// f_int_e = f_int_e + ((disp_mat.transpose())*sigma_e*wtx*wty*wtz*detJacobian); (previous correct)
				//std::cout<<k<<std::endl;

				f_int_e = f_int_e + ((disp_mat.transpose())*sigma_e*wtx*wty*wtz*detJacobian);
			}
		}
	}

	/* Calculating Host Element Stresses and Strains at Centroids */
		dndr = fe_dndr_8(0,0,0);
		dnds = fe_dnds_8(0,0,0);
		dndt = fe_dndt_8(0,0,0);
		jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
		invJacobian = jacobian.inverse();
		dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
		dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
		dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
		disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);

		// Stress- Calculation
		VectorXd sigma_centroid = VectorXd::Zero(6);
		sigma_centroid = fe_stressUpdate(dndx,dndy,dndz,disp_mat,u_e,elements(i,1),1);
		element_stress_host(i,0) = sigma_centroid(0);
		element_stress_host(i,1) = sigma_centroid(3);
		element_stress_host(i,2) = sigma_centroid(5);
		element_stress_host(i,3) = sigma_centroid(3);
		element_stress_host(i,4) = sigma_centroid(1);
		element_stress_host(i,5) = sigma_centroid(4);
		element_stress_host(i,6) = sigma_centroid(5);
		element_stress_host(i,7) = sigma_centroid(4);
		element_stress_host(i,8) = sigma_centroid(2);

		// Strain - Calculation
		MatrixXd E = MatrixXd::Zero(3,3);
		MatrixXd I = MatrixXd::Ones(3,3);
		MatrixXd F = fe_calDefGrad(dndx,dndy,dndz,u_e);
		/* Polar Deomposition of F */
		JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV); // F = USV
		MatrixXd Usvd = svd.matrixU();
		MatrixXd Vsvd = svd.matrixV();
		VectorXd Ssvd = svd.singularValues();
		MatrixXd Ssvdmat = MatrixXd::Zero(3,3);
		Ssvdmat(0,0) = Ssvd(0);
		Ssvdmat(1,1) = Ssvd(1);
		Ssvdmat(2,2) = Ssvd(2);
		MatrixXd Rpd = Usvd*Vsvd.adjoint();
		MatrixXd Upd = Vsvd*Ssvdmat*Vsvd.adjoint();
		// E = 0.5*(F.transpose()*F-I);
		E = Upd.log();
		element_strain_host(i,0) = E(0,0);
		element_strain_host(i,1) = E(0,1);
		element_strain_host(i,2) = E(0,2);
		element_strain_host(i,3) = E(1,0);
		element_strain_host(i,4) = E(1,1);
		element_strain_host(i,5) = E(1,2);
		element_strain_host(i,6) = E(2,0);
		element_strain_host(i,7) = E(2,1);
		element_strain_host(i,8) = E(2,2);

		// Host element internal nodal forces

		VectorXd nat_coord = VectorXd::Zero(3);
		nat_coord(0) = 0.005;
		nat_coord(1) = 0.005;
		nat_coord(2) = 0.005;

		VectorXd iso_coord = fe_newtonRhapson(nat_coord,xcoord,ycoord,zcoord);

	/* TRUSS ANALYSIS - WORKS FOR ONLY SINGLE ELEMENT PROBLEMS
	// truss element nodal force

	int truss_intg = 2;
	VectorXd truss_intg_points = guass_points(2);
	VectorXd truss_intg_weights = guass_weights(2);

	for(int q=0;q<truss_intg;q++){
		double t = truss_intg_points(q);
		double wtt = truss_intg_weights(q);
		dndr = fe_dndr_8(0,0,t);
		dnds = fe_dnds_8(0,0,t);
		dndt = fe_dndt_8(0,0,t);
		jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
		invJacobian = jacobian.inverse();
		dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
		dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
		dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
		disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);
		VectorXd sigma_truss = VectorXd::Zero(6);
		sigma_truss = fe_stressUpdate_1d(dndx,dndy,dndz,u_e,elements_truss(0,1),nodes_truss);
		VectorXd f_int_truss = (disp_mat.transpose()*sigma_truss*wtt*0.5*0.5);

		VectorXd sigma_correction = fe_stressUpdate_1d(dndx,dndy,dndz,u_e,elements(i,1),nodes_truss);
		VectorXd f_int_correction = (disp_mat.transpose()*sigma_correction*wtt*0.5*0.5);

		f_int_e = f_int_e + f_int_truss - f_int_correction; // no Volume Redundancy
		// f_int_e = (0.5*f_int_e) + f_int_truss; // no Volume Redundancy using volume fractions

		// f_int_e = f_int_e + f_int_truss ; // With Volume Redundancy
	}

	// Calculating Truss Element Stress and Strains

	dndr = fe_dndr_8(0,0,0);
	dnds = fe_dnds_8(0,0,0);
	dndt = fe_dndt_8(0,0,0);
	jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
	invJacobian = jacobian.inverse();
	dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
	dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
	dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
	disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);

	// Strain - Calculation
		F = fe_calDefGrad(dndx,dndy,dndz,u_e);
		MatrixXd C = F.transpose()*F;
		VectorXd dir_truss = VectorXd::Zero(3);
		dir_truss << (nodes_truss(1,1)-nodes_truss(0,1)) , (nodes_truss(1,2)-nodes_truss(0,2)), (nodes_truss(1,3)-nodes_truss(0,3));
		VectorXd deformed_truss = C*dir_truss;
		double lambda_tmp = dir_truss.dot(deformed_truss);
		double lambda = sqrt(lambda_tmp);
		element_strain_truss(i,0) = lambda-1;
		element_strain_truss(i,1) = 0;
		element_strain_truss(i,2) = 0;
		element_strain_truss(i,3) = 0;
		element_strain_truss(i,4) = 0;
		element_strain_truss(i,5) = 0;
		element_strain_truss(i,6) = 0;
		element_strain_truss(i,7) = 0;
		element_strain_truss(i,8) = 0;

	// Stress- Calculation
		std::string model;
		model = fe_get_model(elements_truss(0,1));
		double sigma = 0;
		if(model=="simple_elastic"){
			double E = fe_get_mats(elements_truss(0,1),1);
			sigma = E*log(lambda);
		}
		element_stress_truss(i,0) = sigma;
		element_stress_truss(i,1) = 0;
		element_stress_truss(i,2) = 0;
		element_stress_truss(i,3) = 0;
		element_stress_truss(i,4) = 0;
		element_stress_truss(i,5) = 0;
		element_stress_truss(i,6) = 0;
		element_stress_truss(i,7) = 0;
		element_stress_truss(i,8) = 0;
	*/

    }

	f_tot_e = f_ext_e - f_int_e;
	f_tot = fe_scatter(f_tot,f_tot_e,nodes_local,sdof);
	}

	mesh[0].readStressStrain(element_stress_host,element_strain_host);

	return f_tot;
}
