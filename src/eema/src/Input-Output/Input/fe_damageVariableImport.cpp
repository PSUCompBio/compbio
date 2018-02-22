#include "functions.h"

using namespace Eigen;

void fe_damageVariableImport(std::string& damage_variables_import, VectorXd& d, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max) {

	int num_elements_embed = d.size();

	std::ifstream myfile(damage_variables_import.c_str());

	if(myfile.is_open()){
		int num_imported_damage_values = 0;
		myfile >> num_imported_damage_values;
		if (num_imported_damage_values == num_elements_embed) {
			for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {
				myfile >> d(embed_row);
				myfile >> d_fatigue(embed_row);
				myfile >> d_tot(embed_row);
				myfile >> lambda_min(embed_row);
				myfile >> lambda_max(embed_row);
			}
		}
		else {
			std::cout << '\n' << "Input Fiber Damage File has Incorrect Number of Fibers." << '\n' << '\n';
			std::exit(-1);
		}
		myfile.close();
	}
	else{
		std::cout << "No Such Text File Exists - " << damage_variables_import << "\n";
		std::exit(-1);
	}

}

