#include "functions.h"

using namespace Eigen;

void fe_reactionForceWrite_new(std::string& reaction_forces, int plot_state_counter, double& t, double& fr_curr_1, double& fr_curr_2, double& fr_curr_3, double& fr_curr_4) {

	std::ofstream myfile(reaction_forces.c_str());
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << fr_curr_1 << "\t" << fr_curr_2 << "\t" << fr_curr_3 << "\t" << fr_curr_4 << std::endl;
	}
	else{
		std::cout << "No Such Text File Is Created - " << reaction_forces << "\n";
		std::exit(-1);
	}

}

void fe_reactionForceWrite_append(std::string& reaction_forces, int plot_state_counter, double& t, double& fr_curr_1, double& fr_curr_2, double& fr_curr_3, double& fr_curr_4) {

	std::ofstream myfile(reaction_forces.c_str(),std::ios_base::app);
	if(myfile.is_open()){
		myfile << std::setw(5) << (plot_state_counter - 1) << "\t";
		myfile << std::setw(5) << std::scientific << std::setprecision(5) << t << "\t";
		myfile << fr_curr_1 << "\t" << fr_curr_2 << "\t" << fr_curr_3 << "\t" << fr_curr_4 << std::endl;
	}
	else{
		std::cout << "No Such Text File With This Name Exists - " << reaction_forces << "\n";
		std::exit(-1);
	}

}
