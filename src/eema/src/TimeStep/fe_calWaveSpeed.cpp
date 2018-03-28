#include "functions.h"

using namespace Eigen;

/** This function calculates the wave speed for an element based on its material properties */
double
fe_calWaveSpeed(int material_id, double volume_initial, double volume_current)
{
    model_cws = fe_get_model(material_id, "mechanical");
    c_wave_cws = 0;
    E_cws   = fe_get_mats(material_id, 1, "mechanical");
    nu_cws  = fe_get_mats(material_id, 2, "mechanical");
    rho_cws = fe_get_mats(material_id, 0, "mechanical");

    if (model_cws == "simple_elastic") {
        if (nu_cws == 0.5)
            c_wave_cws = sqrt(E_cws / rho_cws);
        else
            c_wave_cws = sqrt((E_cws * (1 - nu_cws)) / (rho_cws * (1 + nu_cws) * (1 - (2 * nu_cws))));
    }
    else if (model_cws == "mooney-rivlin_hyperelastic" || model_cws == "ogden_hyperelastic")
        c_wave_cws = sqrt(volume_current / volume_initial) * sqrt(E_cws / rho_cws) * 1.3;

    return c_wave_cws;
}
