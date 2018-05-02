#include<iomanip>
#include"functions.h"


using namespace Eigen;

//typedef vector<double> Vec;
//typedef vector<Vec> Mat;

MatrixXd text2matrix(std::string name, int cols){

	// std::cout<<"**Entering the text2matrix function**"<<std::endl;
        std::ifstream myfile(name.c_str());
        int tmp=0;
        int rows=0;
        // double **A;
        double data;

        while(myfile>>data){
                tmp = tmp+1;
        }
        myfile.close();

        rows = tmp/cols;

        // std::cout<<"The size of the matrix is: "<<rows<<" rows & "<<cols<<" columns. \n";

        // A = new double*[rows];
        MatrixXd A(rows,cols);

        std::ifstream myfile1(name.c_str());

        for(int i=0;i<rows;i++){
                for(int j=0;j<cols;j++){
                        myfile1>>A(i,j);
                }
        }

        myfile1.close();
	// std::cout<<"Exiting the text2matrix function. \n";
        return A;
}

VectorXd text2vector(std::string name){
        std::ifstream myfile(name.c_str());

        int tmp = 0;
        double data;
        
        while(myfile>>data){
          tmp = tmp + 1;
        }

        VectorXd A(tmp);

        std::ifstream myfile1(name.c_str());
        for(int i=0;i<tmp;i++){
          myfile1 >> A(i);
        }

        return A;
}
