/* Developed by Harsha Teja Garimella */
#include "functions.h"

using namespace Eigen;

/** \brief This function calculates the updated stress for 3d elements - elastic, hyperelastic material models were implemented so far. */
VectorXd fe_stressUpdate(VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt) {

	/** This function is  for 2d or 3d elements. For 1d elements, different stress update function was included */

	std::string model;
	model = fe_get_model(opt, "mechanical");

	VectorXd sigma = VectorXd::Zero(ndof * 2); //stress-vector

	/** This block develops outputs the updated stress for a 3d elastic material */
	if (model == "simple_elastic") {
		sigma = fe_simple_elastic(dndx, dndy, dndz, disp_mat, u, opt, return_opt);
	}

	if (model == "saint_venant") {
		sigma = fe_saintvenant_elastic(dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block outputs the updated stress for a mooney-rivlin hyperelastic material model */
	if (model == "mooney-rivlin_hyperelastic") {
		sigma = fe_mooneyrivlin_hyperelastic(dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block develops outputs the updated stress for a 3d ogden hyperelastic material */
	if (model == "ogden_hyperelastic") {
		/* TODO: Ogden - this material model needs work */
		// n = number of coefficients - used a double type for consistency
		sigma = fe_ogden_hyperelastic(dndx, dndy, dndz, u, opt, return_opt);
	}

	return sigma;
}

/******************************************/

void fe_stressUpdate_pbr(VectorXd& sigma, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, MatrixXd& disp_mat, VectorXd& u, int opt, int return_opt) {

	/** This function is  for 2d or 3d elements. For 1d elements, different stress update function was included */

	std::string model;
	model = fe_get_model(opt, "mechanical");

	/** This block develops outputs the updated stress for a 3d elastic material */
	if (model == "simple_elastic") {
		fe_simple_elastic_pbr(sigma, dndx, dndy, dndz, disp_mat, u, opt, return_opt);
	}

	if (model == "saint_venant") {
		fe_saintvenant_elastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block outputs the updated stress for a mooney-rivlin hyperelastic material model */
	if (model == "mooney-rivlin_hyperelastic") {
		fe_mooneyrivlin_hyperelastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

	/** This block develops outputs the updated stress for a 3d ogden hyperelastic material */
	if (model == "ogden_hyperelastic") {
		/* TODO: Ogden - this material model needs work */
		// n = number of coefficients - used a double type for consistency
		fe_ogden_hyperelastic_pbr(sigma, dndx, dndy, dndz, u, opt, return_opt);
	}

}
