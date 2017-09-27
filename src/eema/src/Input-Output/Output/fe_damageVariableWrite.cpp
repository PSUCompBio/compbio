#include "functions.h"

using namespace Eigen;

void fe_damageVariableWrite_new(std::string& damage_variables, int plot_state_counter, double& t, double& d_curr, double& delta_d_curr, double& d_tot) {

	std::ofstream myfile(damage_variables.c_str());
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << d_curr << "\t" << delta_d_curr << "\t" << d_tot << std::endl;
	}
	else{
		std::cout << "No Such Text File Is Created - " << damage_variables << "\n";
		std::exit(-1);
	}

}

void fe_damageVariableWrite_append(std::string& damage_variables, int plot_state_counter, double& t, double& d_curr, double& delta_d_curr, double& d_tot) {

	std::ofstream myfile(damage_variables.c_str(),std::ios_base::app);
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << d_curr << "\t" << delta_d_curr << "\t" << d_tot << std::endl;
	}
	else{
		std::cout << "No Such Text File With This Name Exists - " << damage_variables << "\n";
		std::exit(-1);
	}

}
