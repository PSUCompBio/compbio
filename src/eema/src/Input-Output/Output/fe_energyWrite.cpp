#include "functions.h"

using namespace Eigen;

void fe_energyWrite_new(std::string& internal_energy, std::string& viscous_dissipation_energy, std::string& external_energy, std::string& kinetic_energy, std::string& total_energy, int plot_state_counter, double& t, double& energy_int_new, double& energy_vd_new, double& energy_ext_new, double& energy_kin, double& energy_total) {

	new_double2textWithTime(internal_energy, plot_state_counter - 1, t, energy_int_new);
	new_double2textWithTime(viscous_dissipation_energy, plot_state_counter - 1, t, energy_vd_new);
	new_double2textWithTime(external_energy, plot_state_counter - 1, t, energy_ext_new);
	new_double2textWithTime(kinetic_energy, plot_state_counter - 1, t, energy_kin);
	new_double2textWithTime(total_energy, plot_state_counter - 1, t, energy_total);

}

void fe_energyWrite_append(std::string& internal_energy, std::string& viscous_dissipation_energy, std::string& external_energy, std::string& kinetic_energy, std::string& total_energy, int plot_state_counter, double& t, double& energy_int_new, double& energy_vd_new, double& energy_ext_new, double& energy_kin, double& energy_total) {
	append_double2textWithTime(internal_energy, plot_state_counter - 1, t, energy_int_new);
	append_double2textWithTime(viscous_dissipation_energy, plot_state_counter - 1, t, energy_vd_new);
	append_double2textWithTime(external_energy, plot_state_counter - 1, t, energy_ext_new);
	append_double2textWithTime(kinetic_energy, plot_state_counter - 1, t, energy_kin);
	append_double2textWithTime(total_energy, plot_state_counter - 1, t, energy_total);
}
