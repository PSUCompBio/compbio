#include "functions.h"

using namespace Eigen;

void fe_damageVariableExport(std::string& damage_variables_export, VectorXd& d, VectorXd& delta_d, VectorXd& d_tot) {

	int num_elements_embed = d.size();

	std::ofstream myfile(damage_variables_export.c_str());

	if(myfile.is_open()){

		myfile << num_elements_embed << std::endl;

		for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {

			double tmp_d = d(embed_row);
			double tmp_delta_d = delta_d(embed_row);
			double tmp_d_tot = d_tot(embed_row);

			myfile << std::setprecision(12) << tmp_d << "\t" << tmp_delta_d << "\t" << tmp_d_tot << std::endl;

		}

		myfile.close();

	}

	else{

		std::cout << "No Such Text File Is Created - " << damage_variables_export << "\n";
		std::exit(-1);

	}

}

