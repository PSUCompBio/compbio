#include "functions.h"

using namespace Eigen;

/** This function calculates the wave speed for an element based on its material properties */
double fe_calWaveSpeed(int material_id){

	std::string model = fe_get_model(material_id);
	double c_wave = 0;

	if(model=="simple_elastic" || model == "mooney-rivlin_hyperelastic"){

		double E = fe_get_mats(material_id,1);
		double nu = fe_get_mats(material_id,2);
		double rho = fe_get_mats(material_id,0);

		if(nu==0.5){
			c_wave = sqrt(E/rho);
		} else{
			c_wave = sqrt((E*(1-nu))/(rho*(1+nu)*(1-(2*nu))));
		}

	}

	return c_wave;

}
