

#include"functions.h"

using namespace Eigen;

MatrixXd fe_calTransformation(VectorXd x, VectorXd y, VectorXd z, int choice){
// choice - 1: Stress Transformation Matrix
// choice - 2: Strain Transformation Matrix

	MatrixXd T(6,6);
	double l_1,l_2,l_3,m_1,m_2,m_3,n_1,n_2,n_3;

	Vector3d dir_truss_x(3);
	dir_truss_x(0) = x(1) - x(0);
	dir_truss_x(1) = y(1) - y(0);
	dir_truss_x(2) = z(1) - z(0);

	// std::cout<<"Direction of truss element is: "<<dir_truss_x(0)<<":"<<dir_truss_x(1)<<":"<<dir_truss_x(2)<<"\n";

	Vector3d dir_truss_y(3);
	Vector3d dir_truss_z(3);
	
	if(dir_truss_x(0)!=0){
		dir_truss_y(1) = 1;	
		dir_truss_y(2) = 0;
		dir_truss_y(0) = -(((dir_truss_y(1)*dir_truss_x(1))+(dir_truss_y(2)*dir_truss_x(2)))/dir_truss_x(0));
	}
	else{
		if((dir_truss_x(0)==0) && (dir_truss_x(1)!=0)){
			dir_truss_y(0) = 1;
			dir_truss_y(2) = 1;
			dir_truss_y(1) = -(((dir_truss_y(2)*dir_truss_x(2))+(dir_truss_y(0)*dir_truss_x(0)))/dir_truss_x(1));		
		}
		else if((dir_truss_x(0)==0) && (dir_truss_x(2)!=0)){
			dir_truss_y(0) = 0;
			dir_truss_y(1) = 1;
			dir_truss_y(2) = -(((dir_truss_y(0)*dir_truss_x(0))+(dir_truss_y(1)*dir_truss_x(1)))/dir_truss_x(2)); 
		}
		else{
			std::cout<<"Truss element is not practically possible"<<"\n";
		}
	}
	
	dir_truss_z = dir_truss_x.cross(dir_truss_y);

	// Proof of Right Handed Coordinate System
	//Vector3d tmp1 = dir_truss_y.cross(dir_truss_z);
	//std::cout<<"Y x Z is: "<<tmp1<<"\n";
	//std::cout<<"X is: "<<dir_truss_x<<"\n";

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

	if(choice==2){
	T.row(0) << (l_1*l_1), (m_1*m_1), (n_1*n_1), (l_1*m_1), (m_1*n_1), (n_1*l_1);
	T.row(1) << (l_2*l_2), (m_2*m_2), (n_2*n_2), (l_2*m_2), (m_2*n_2), (n_2*l_2);
	T.row(2) << (l_3*l_3), (m_3*m_3), (n_3*n_3), (l_3*m_3), (m_3*n_3), (n_3*l_3);
	T.row(3) << (2*l_1*l_2), (2*m_1*m_2), (2*n_1*n_2), ((l_1*m_2)+(l_2*m_1)), ((m_1*n_2)+(m_2*n_1)), ((l_1*n_2)+(l_2*n_1));
	T.row(4) << (2*l_2*l_3), (2*m_2*m_3), (2*n_2*n_3), ((l_2*m_3)+(l_3*m_2)), ((m_2*n_3)+(m_3*n_2)), ((l_2*n_3)+(l_3*n_2));
	T.row(5) << (2*l_1*l_3), (2*m_1*m_3), (2*n_1*n_3), ((l_1*m_3)+(l_3*m_1)), ((m_1*n_3)+(m_3*n_1)), ((l_1*n_3)+(l_3*n_1));
	}
	else if(choice==1) {
	T.row(0) << pow(l_1,2), pow(m_1,2), pow(n_1,2), (2*l_1*m_1), (2*m_1*n_1), (2*n_1*l_1);
        T.row(1) << pow(l_2,2), pow(m_2,2), pow(n_2,2), (2*l_2*m_2), (2*m_2*n_2), (2*n_2*l_2);
        T.row(2) << pow(l_3,2), pow(m_3,2), pow(n_3,2), (2*l_3*m_3), (2*m_3*n_3), (2*n_3*l_3);
        T.row(3) << (1*l_1*l_2), (1*m_1*m_2), (1*n_1*n_2), ((l_1*m_2)+(l_2*m_1)), ((m_1*n_2)+(m_2*n_1)), ((l_1*n_2)+(l_2*n_1));
        T.row(4) << (1*l_2*l_3), (1*m_2*m_3), (1*n_2*n_3), ((l_2*m_3)+(l_3*m_2)), ((m_2*n_3)+(m_3*n_2)), ((l_2*n_3)+(l_3*n_2));
        T.row(5) << (1*l_1*l_3), (1*m_1*m_3), (1*n_1*n_3), ((l_1*m_3)+(l_3*m_1)), ((m_1*n_3)+(m_3*n_1)), ((l_1*n_3)+(l_3*n_1));
	}
	else{
		std::cout<<"WRONG CHOICE OF TRANSFORMATION MATRIX\n";
	}
	
	//std::cout<<"Transformation Matrix is: "<<T<<"\n";
	return T;
}
