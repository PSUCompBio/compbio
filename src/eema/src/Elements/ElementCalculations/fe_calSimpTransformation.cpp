

#include"functions.h"

using namespace Eigen;

MatrixXd fe_calSimpTransformation(MatrixXd truss_nodes){
// choice - 1: Stress Transformation Matrix
// choice - 2: Strain Transformation Matrix

	MatrixXd T(3,3);
	double l_1,l_2,l_3,m_1,m_2,m_3,n_1,n_2,n_3;

	Vector3d dir_truss_x(3);
	dir_truss_x(0) = truss_nodes(1,1) - truss_nodes(0,1);
	dir_truss_x(1) = truss_nodes(1,2) - truss_nodes(0,2);
	dir_truss_x(2) = truss_nodes(1,3) - truss_nodes(0,3);

	Vector3d dir_truss_y(3);
	Vector3d dir_truss_z(3);
	
	if(dir_truss_x(0)!=0){
		dir_truss_y(1) = 1;	
		dir_truss_y(2) = 1;
		dir_truss_y(0) = -((dir_truss_x(1)+dir_truss_x(2))/dir_truss_x(0));
	}
	else{
		if((dir_truss_x(0)==0) && (dir_truss_x(1)!=0)){
			dir_truss_y(0) = 1;
			dir_truss_y(2) = 1;
			dir_truss_y(1) = -((dir_truss_x(2)+dir_truss_x(0))/dir_truss_x(1));		
		}
		else if((dir_truss_x(0)==0) && (dir_truss_x(2)!=0)){
			dir_truss_y(0) = 1;
			dir_truss_y(1) = 1;
			dir_truss_y(2) = -((dir_truss_x(0)+dir_truss_x(1))/dir_truss_x(2)); 
		}
		else{
			std::cout<<"Truss element is not practically possible"<<"\n";
		}
	}
	
	dir_truss_z = dir_truss_x.cross(dir_truss_y);

	VectorXd dir_global_x(3);
	dir_global_x << 1, 0, 0;
	
	VectorXd dir_global_y(3);
	dir_global_y << 0, 1, 0;

	VectorXd dir_global_z(3);
	dir_global_z << 0, 0, 1;

	l_1 = (dir_global_x.dot(dir_truss_x))/(dir_global_x.norm()*dir_truss_x.norm());
	l_2 = (dir_global_x.dot(dir_truss_y))/(dir_global_x.norm()*dir_truss_y.norm());
	l_3 = (dir_global_x.dot(dir_truss_z))/(dir_global_x.norm()*dir_truss_z.norm()); 
	m_1 = (dir_global_y.dot(dir_truss_x))/(dir_global_y.norm()*dir_truss_x.norm());	
	m_2 = (dir_global_y.dot(dir_truss_y))/(dir_global_y.norm()*dir_truss_y.norm());
	m_3 = (dir_global_y.dot(dir_truss_z))/(dir_global_y.norm()*dir_truss_z.norm());
	n_1 = (dir_global_z.dot(dir_truss_x))/(dir_global_z.norm()*dir_truss_x.norm());
	n_2 = (dir_global_z.dot(dir_truss_y))/(dir_global_z.norm()*dir_truss_y.norm());
	n_3 = (dir_global_z.dot(dir_truss_z))/(dir_global_z.norm()*dir_truss_z.norm());

	T.row(0) << l_1, m_1, n_1;
	T.row(1) << l_2, m_2, n_2;
	T.row(2) << l_3, m_3, n_3;

	return T;
}
