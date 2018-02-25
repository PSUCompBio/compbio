#include "functions.h"

using namespace Eigen;

void fe_checkEnergies(VectorXd& U_prev, VectorXd& U, VectorXd& fi_prev, VectorXd& fi_curr, VectorXd& f_damp_prev, VectorXd& f_damp_curr, VectorXd& fe_prev, VectorXd& fe, VectorXd& fr_prev, VectorXd& fr_curr, VectorXd& m_system, VectorXd& V, double& energy_int_old, double& energy_int_new, double& energy_vd_old, double& energy_vd_new, double& energy_ext_old, double& energy_ext_new, double& energy_kin, double& energy_total, double& energy_max) {
	/** Calculating the internal energy terms */
	fe_calculateIE(U_prev, U, fi_prev, fi_curr, energy_int_old, energy_int_new);
	/** Calculating the viscous dissipation energy terms */
	fe_calculateVDE(U_prev, U, f_damp_prev, f_damp_curr, energy_vd_old, energy_vd_new);
	/** Calculating the external energy terms */
	fe_calculateEW(U_prev, U, fe_prev, fe, fr_prev, fr_curr, energy_ext_old, energy_ext_new);
	/** Calculating the kinetic energy */
	fe_calculateKE(energy_kin, m_system, V);
	/* Calculating the total energy of the system */
	energy_total = std::abs(energy_kin + energy_int_new - energy_ext_new);
	energy_max   = std::max(std::max(energy_kin, energy_int_new), energy_ext_new);
	/*
	if (energy_total > (eps_energy * (energy_max))) {
		std::cout << "**********************************************" << std::endl;
		std::cout << "ALERT: INSTABILITIES IN THE SYSTEM DETECTED \n BASED ON THE ENERGY BALANCE CHECK \n";
		std::cout << "**********************************************" << std::endl;
	}
	*/
}
