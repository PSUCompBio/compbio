#include<iomanip>
#include"functions.h"
using namespace Eigen;

void new_vector2text(std::string name, VectorXd vector, int width){
	if(width==0){
			width=3;
		}
	std::ofstream myfile(name.c_str());

	for(int i=0;i<vector.size();i++){
		myfile << std::setw(width) << vector(i) << "\n";
	}
	myfile.close();
}

void append_double2text(std::string name, double a){
	std::ofstream myfile(name.c_str(),std::ios_base::app);
	if(myfile.is_open()){
		myfile << std::setw(5)<<std::scientific<<std::setprecision(5) << a << "\n";
	}
	else{
		std::cout << "No Such Text File With This Name Exists - " << name << "\n";
		std::exit(-1);
	}
}

void append_double2textWithTime(std::string name, int frame, double time, double value){
	std::ofstream myfile(name.c_str(),std::ios_base::app);
	if(myfile.is_open()){
		
		myfile << std::setw(5) << frame << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << time << "\t";
		myfile << value << std::endl;
	}
	else{
		std::cout << "No Such Text File With This Name Exists - " << name << "\n";
		std::exit(-1);
	}
}

void new_double2text(std::string name, double a){
	std::ofstream myfile(name.c_str());
	if(myfile.is_open()){
		myfile << std::setw(5)<<std::scientific<<std::setprecision(5) << a << "\n";
	}
	else{
		std::cout << "No Such Text File Is Created - " << name << "\n";
		std::exit(-1);
	}
}

void new_double2textWithTime(std::string name, int frame, double time, double value){
	std::ofstream myfile(name.c_str());
	if(myfile.is_open()){
		myfile << std::setw(5) << frame << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << time << "\t";
		myfile << value << std::endl;
	}
	else{
		std::cout << "No Such Text File Is Created - " << name << "\n";
		std::exit(-1);
	}
}
