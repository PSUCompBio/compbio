#include "functions.h"

using namespace Eigen;

void fe_damageVariableExport(std::string& damage_variables_export, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max) {

	int num_elements_embed = d_static.size();

	std::ofstream myfile(damage_variables_export.c_str());

	if(myfile.is_open()){
		myfile << num_elements_embed << std::endl;
		for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {
			double tmp_d_static = d_static(embed_row);
			double tmp_d_fatigue = d_fatigue(embed_row);
			double tmp_d_tot = d_tot(embed_row);
			double tmp_lambda_min = lambda_min(embed_row);
			double tmp_lambda_max = lambda_max(embed_row);
			myfile << std::setprecision(12) << tmp_d_static << "\t" << tmp_d_fatigue << "\t" << tmp_d_tot << "\t" << tmp_lambda_min << "\t" << tmp_lambda_max << std::endl;
		}
		myfile.close();
	}
	else{
		std::cout << "No Such Text File Is Created - " << damage_variables_export << "\n";
		std::exit(-1);
	}

}

