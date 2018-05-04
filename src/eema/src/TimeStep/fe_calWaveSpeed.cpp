#include "functions.h"

using namespace Eigen;

/** This function calculates the wave speed for an element based on its material properties */

double fe_calWaveSpeed(int material_id, double volume_initial, double volume_current)
{
    std::string model_cws = fe_get_model(material_id, "mechanical");
    double c_wave = 0;
    double E   = fe_get_mats(material_id, 1, "mechanical");
    double nu  = fe_get_mats(material_id, 2, "mechanical");
    double rho = fe_get_mats(material_id, 0, "mechanical");

    if (model_cws == "simple_elastic") {
        if (nu == 0.5)
            c_wave = sqrt(E / rho);
        else
            c_wave = sqrt((E * (1 - nu)) / (rho * (1 + nu) * (1 - (2 * nu))));
    }
    else if (model_cws == "mooney-rivlin_hyperelastic" || model_cws == "ogden_hyperelastic")
        c_wave = sqrt(volume_current / volume_initial) * sqrt(E / rho) * 1.3;

    return c_wave;
}
