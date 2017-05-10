#include "functions.h"
#include "Materials.h"

using namespace Eigen;

void Materials::readMats(int a, std::string b, VectorXd c){

	mat_id = a;
	mat_model = b;

	mat_properties = VectorXd::Zero(100);
	mat_properties = c;
}

void Materials::printMats(){
	std::cout<<"Material ID is: "<<mat_id<<"\n";
	std::cout<<"Material model is: "<<mat_model<<"\n";
	std::cout<<"Material Properties are: "<<mat_properties<<std::endl;
}
