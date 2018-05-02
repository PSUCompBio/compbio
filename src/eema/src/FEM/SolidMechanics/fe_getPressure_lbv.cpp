#include "functions.h"

using namespace Eigen;

void fe_getPressure_lbv_pbr(VectorXd& pressure, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id)
{
    // Reference: "Getting Started with Abaqus," Section 9.5.1 "Bulk Viscosity"
    // Linear bulk viscosity is always included to damp "ringing" in the highest element frequency.

    fe_calDefGrad_pbr(F_curr_lbv, dndx, dndy, dndz, u);

    fe_invMatrix_pbr(F_inv_lbv, F_curr_lbv);

    F_invT_lbv = F_inv_lbv.transpose();

    fe_calDefGrad_pbr(F_prev_lbv, dndx, dndy, dndz, u_prev);

    F_dot_lbv = (F_curr_lbv - F_prev_lbv) / dT;

    F_dotT_lbv = F_dot_lbv.transpose();

    D_lbv = 0.5 * ((F_dot_lbv * F_inv_lbv) + (F_invT_lbv * F_dotT_lbv));

    vol_strain_rate_lbv = D_lbv.trace(); // volumetric strain rate

    volume_initial_lbv = fe_calVolume(xcoord, ycoord, zcoord); // initial volume of element

    // calculate current nodal coordinates
    nnel_lbv = xcoord.size();

    VectorXd xcoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd ycoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd zcoord_curr = VectorXd::Zero(nnel_lbv);

    for (i_lbv = 0; i_lbv < nnel_lbv; i_lbv++) {
        xcoord_curr(i_lbv) = xcoord(i_lbv) + u(i_lbv * 3);
        ycoord_curr(i_lbv) = ycoord(i_lbv) + u((i_lbv * 3) + 1);
        zcoord_curr(i_lbv) = zcoord(i_lbv) + u((i_lbv * 3) + 2);
    }

    volume_current_lbv = fe_calVolume(xcoord_curr, ycoord_curr, zcoord_curr); // current volume of element

    lc_lbv = fe_minElementLength(xcoord_curr, ycoord_curr, zcoord_curr); // current element characteristic length

    c_wave_lbv = fe_calWaveSpeed(material_id, volume_initial_lbv, volume_current_lbv); // current dilatational wave speed

    rho_initial_lbv = fe_get_mats(material_id, 0, "mechanical"); // initial material density

    rho_current_lbv = rho_initial_lbv * (volume_initial_lbv / volume_current_lbv); // current material density

    pressure_scalar_lbv = 0.06 * rho_current_lbv * c_wave_lbv * lc_lbv * vol_strain_rate_lbv; // bulk viscosity pressure, scalar format

    pressure_matrix_lbv = pressure_scalar_lbv * I;

    pressure = fe_tensor2voigt(pressure_matrix_lbv); // bulk viscosity pressure, voigt format

}


void fe_getPressure_lbv_pbr_array(VectorXd& pressure, int i, int x, int y, int z, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id)
{
    // Reference: "Getting Started with Abaqus," Section 9.5.1 "Bulk Viscosity"
    // Linear bulk viscosity is always included to damp "ringing" in the highest element frequency.

    fe_calDefGrad_pbr_array(F_curr_lbv, i, x, y, z, u);

    fe_invMatrix_pbr(F_inv_lbv, F_curr_lbv);

    F_invT_lbv = F_inv_lbv.transpose();

    fe_calDefGrad_pbr_array(F_prev_lbv, i, x, y, z, u_prev);

    F_dot_lbv = (F_curr_lbv - F_prev_lbv) / dT;

    F_dotT_lbv = F_dot_lbv.transpose();

    D_lbv = 0.5 * ((F_dot_lbv * F_inv_lbv) + (F_invT_lbv * F_dotT_lbv));

    vol_strain_rate_lbv = D_lbv.trace(); // volumetric strain rate

    volume_initial_lbv = fe_calVolume(xcoord, ycoord, zcoord); // initial volume of element

    // calculate current nodal coordinates
    nnel_lbv = xcoord.size();

    VectorXd xcoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd ycoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd zcoord_curr = VectorXd::Zero(nnel_lbv);

    for (i_lbv = 0; i_lbv < nnel_lbv; i_lbv++) {
        xcoord_curr(i_lbv) = xcoord(i_lbv) + u(i_lbv * 3);
        ycoord_curr(i_lbv) = ycoord(i_lbv) + u((i_lbv * 3) + 1);
        zcoord_curr(i_lbv) = zcoord(i_lbv) + u((i_lbv * 3) + 2);
    }

    volume_current_lbv = fe_calVolume(xcoord_curr, ycoord_curr, zcoord_curr); // current volume of element

    lc_lbv = fe_minElementLength(xcoord_curr, ycoord_curr, zcoord_curr); // current element characteristic length

    c_wave_lbv = fe_calWaveSpeed(material_id, volume_initial_lbv, volume_current_lbv); // current dilatational wave speed

    rho_initial_lbv = fe_get_mats(material_id, 0, "mechanical"); // initial material density

    rho_current_lbv = rho_initial_lbv * (volume_initial_lbv / volume_current_lbv); // current material density

    pressure_scalar_lbv = 0.06 * rho_current_lbv * c_wave_lbv * lc_lbv * vol_strain_rate_lbv; // bulk viscosity pressure, scalar format

    pressure_matrix_lbv = pressure_scalar_lbv * I;

    pressure = fe_tensor2voigt(pressure_matrix_lbv); // bulk viscosity pressure, voigt format

}
