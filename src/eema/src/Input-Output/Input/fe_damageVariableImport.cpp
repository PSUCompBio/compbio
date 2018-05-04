#include "functions.h"

using namespace Eigen;

void fe_damageVariableImport(std::string& damage_variables_import, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, double t_healing) {

	int num_elements_embed = d_static.size();

	std::ifstream myfile(damage_variables_import.c_str());

	if(myfile.is_open()){
		int num_imported_damage_values = 0;
		myfile >> num_imported_damage_values;
		if (num_imported_damage_values == num_elements_embed) {
			for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {
				myfile >> d_static(embed_row);
				myfile >> d_fatigue(embed_row);
				myfile >> d_tot(embed_row);
				myfile >> lambda_min(embed_row);
				myfile >> lambda_max(embed_row);
				if (t_healing > 0) {
					double healing_rate = 0.01; // constant healing rate per day (i.e., 0.01 = 1 %/day)
					lambda_max(embed_row) = (lambda_max(embed_row) - 1)*pow((1 - healing_rate), t_healing) + 1;
					lambda_min(embed_row) = 1 - (1 - lambda_min(embed_row))*pow((1 - healing_rate), t_healing);
					d_static(embed_row) = d_static(embed_row)*pow((1 - healing_rate), t_healing);
					d_fatigue(embed_row) = d_fatigue(embed_row)*pow((1 - healing_rate), t_healing);
					if (lambda_max(embed_row) < 1) {
						lambda_max(embed_row) = 1;
					}
					if (lambda_min(embed_row) > 1) {
						lambda_min(embed_row) = 1;
					}
					if (d_static(embed_row) < 0) {
						d_static(embed_row) = 0;
					}
					if (d_fatigue(embed_row) < 0) {
						d_fatigue(embed_row) = 0;
					}
					d_tot(embed_row) = d_static(embed_row) + d_fatigue(embed_row);
					if (d_tot(embed_row) > 1) {
						d_tot(embed_row) = 1;
					}
				}
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

