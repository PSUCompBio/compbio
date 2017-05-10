#include "functions.h"

using namespace Eigen;

/** For all elements -- this function calculates the minimum critical timestep */
double fe_getTimeStep(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext){

	double deltaT_crit = 0;

	int nel = elements.rows();
	int nnel = (elements.cols()-2);
	int nnode = nodes.rows();
	int sdof = nnode * ndof;
	int edof = nnel * ndof;
	VectorXd xcoord = VectorXd::Zero(nnel);
	VectorXd ycoord = VectorXd::Zero(nnel);
	VectorXd zcoord = VectorXd::Zero(nnel);

	VectorXd deltaT_element = VectorXd::Zero(nel);

	for(int i=0;i<nel;i++){


		VectorXi nodes_local = VectorXi::Zero(nnel);

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

		VectorXd u_e(edof);
		u_e = fe_gather(u,u_e,nodes_local,sdof);

		MatrixXd nodes_new(elements.rows(),ndof);

		int counter;

		for(int j=0;j<nnel;j++){
			counter = j*3;
			xcoord(j) = xcoord(j)+u_e(counter);
			ycoord(j) = ycoord(j)+u_e(counter+1);
			zcoord(j) = zcoord(j)+u_e(counter+2);
		}

		deltaT_element(i) = fe_calTimeStep(xcoord,ycoord,zcoord,elements(i,1)); // reduction factor for time step added.
	}

	deltaT_crit = deltaT_element.minCoeff();

	return deltaT_crit;
}
