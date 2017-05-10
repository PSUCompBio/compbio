#include<iostream>

#include "functions.h"
using namespace Eigen;

MatrixXd fe_transformMass(MatrixXd m, int opt){

	// opt = 1: No transformation
	// opt = 2: Converted to Direct lumped mass
	// opt = 3: Converted to Variational lumped mass

	MatrixXd m_new = MatrixXd::Zero(m.rows(),m.cols());

	if(opt==1){
		m_new = m;
	}

	if (opt==2){
	double sum = 0;
	for(int i = 0; i<m.rows() ; i++){
		for(int j = 0; j< m.cols(); j++){
			sum = sum + m(i,j);
		}
		m_new(i,i) = sum;
		sum = 0;
	}
	}

	return m_new;

}
