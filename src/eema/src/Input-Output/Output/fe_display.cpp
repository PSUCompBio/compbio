#include "functions.h"

using namespace Eigen;

void fe_display_matrix(MatrixXd A){

	if(A.rows()==0){
		std::cout<<"The matrix is empty"<<"\n";
	}
	else{

		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.cols();j++){
				std::cout<<A(i,j)<<"\t";
			}
			std::cout<<"\n";
		}

	}

}

void fe_display_matrix(MatrixXi A){

	if(A.rows()==0){
		std::cout<<"The matrix is empty"<<"\n";
	}
	else{

		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.cols();j++){
				std::cout<<A(i,j)<<"\t";
			}
			std::cout<<"\n";
		}

	}

}

void fe_display_vector(VectorXd A){

	if(A.size()==0){
		std::cout<<"Vector is empty"<<"\n";
	}
	else{
		for(int i=0;i<A.size();i++){
			std::cout<<A(i)<<"\t";
		}
		std::cout<<"\n";
	}

}

void fe_display_vector(VectorXi A){

	if(A.size()==0){
		std::cout<<"Vector is empty"<<"\n";
	}
	else{
		for(int i=0;i<A.size();i++){
			std::cout<<A(i)<<"\t";
		}
		std::cout<<"\n";
	}

}
