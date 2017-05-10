#include<iomanip>
#include"functions.h"
using namespace Eigen;

void matrix2text(std::string name, MatrixXd new_slave_master,int width){

	//std::cout<<"**Entering the matrix2text function**\n";

	if(width==0){
		width=3;
	}

	std::ofstream myfile(name.c_str());

        for(int i=0;i<new_slave_master.rows();i++){
        	for(int j=0;j<new_slave_master.cols();j++){
        		//if(j==0){
			//	int tmp;
			//	tmp = int(new_slave_master(i,j));
			//	myfile<<tmp;
			//}
			//else{
				myfile<<std::setw(width)<<new_slave_master(i,j)<<"\t";
			//}
        	}
        	myfile<<"\n";
        }
	//std::cout<<"**Exiting the matrix2text function**\n";

        myfile.close();
}
