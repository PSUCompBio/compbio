#include "functions.h"

using namespace Eigen;

/** This function calculates the wave speed for an element based on its material properties */
double
fe_calWaveSpeed(int material_id, double volume_initial, double volume_current)
{
    std::string model = fe_get_model(material_id, "mechanical");
    double c_wave     = 0;

    if (model == "simple_elastic") {
        double E   = fe_get_mats(material_id, 1, "mechanical");
        double nu  = fe_get_mats(material_id, 2, "mechanical");
        double rho = fe_get_mats(material_id, 0, "mechanical");

        if (nu == 0.5) {
            c_wave = sqrt(E / rho);
        } else {
            c_wave = sqrt((E * (1 - nu)) / (rho * (1 + nu) * (1 - (2 * nu))));
        }
    }

    if (model == "mooney-rivlin_hyperelastic") {
        double K   = fe_get_mats(material_id, 1, "mechanical");
        double rho = fe_get_mats(material_id, 0, "mechanical");
        c_wave = sqrt(volume_current / volume_initial) * sqrt(K / rho) * 1.3;
    }

    if (model == "ogden_hyperelastic") {
        double K   = fe_get_mats(material_id, 1, "mechanical");
        double rho = fe_get_mats(material_id, 0, "mechanical");
        c_wave = sqrt(volume_current / volume_initial) * sqrt(K / rho) * 1.3;
    }

    return c_wave;
}
