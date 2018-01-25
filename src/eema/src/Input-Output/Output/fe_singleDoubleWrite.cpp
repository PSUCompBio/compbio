#include "functions.h"

using namespace Eigen;

void fe_singleDoubleWrite_new(std::string& filename, int plot_state_counter, double& t, double& variable) {

	std::ofstream myfile(filename.c_str());
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << variable << std::endl;
	}
	else{
		std::cout << "No Such Text File Is Created - " << filename << "\n";
		std::exit(-1);
	}

}

void fe_singleDoubleWrite_append(std::string& filename, int plot_state_counter, double& t, double& variable) {

	std::ofstream myfile(filename.c_str(),std::ios_base::app);
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << variable << std::endl;
	}
	else{
		std::cout << "No Such Text File With This Name Exists - " << filename << "\n";
		std::exit(-1);
	}

}
