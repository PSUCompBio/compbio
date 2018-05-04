#include "functions.h"

using namespace Eigen;

void fe_getPressure_lbv_pbr(VectorXd& pressure, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id)
{
    // Reference: "Getting Started with Abaqus," Section 9.5.1 "Bulk Viscosity"
    // Linear bulk viscosity is always included to damp "ringing" in the highest element frequency.

    int nnel_lbv;
    double vol_strain_rate, volume_initial, volume_current, lc, c_wave, rho_initial, rho_current, pressure_scalar;

    MatrixXd F_curr = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_invT = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_prev = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_dot = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_dotT = MatrixXd::Zero(ndof, ndof);
    MatrixXd D = MatrixXd::Zero(ndof, ndof);
    MatrixXd pressure_matrix = MatrixXd::Zero(ndof, ndof);

    fe_calDefGrad_pbr(F_curr, dndx, dndy, dndz, u);

    fe_invMatrix_pbr(F_inv, F_curr);

    F_invT = F_inv.transpose();

    fe_calDefGrad_pbr(F_prev, dndx, dndy, dndz, u_prev);

    F_dot = (F_curr - F_prev) / dT;

    F_dotT = F_dot.transpose();

    D = 0.5 * ((F_dot * F_inv) + (F_invT * F_dotT));

    vol_strain_rate = D.trace(); // volumetric strain rate

    volume_initial = fe_calVolume(xcoord, ycoord, zcoord); // initial volume of element

    // calculate current nodal coordinates
    nnel_lbv = xcoord.size();

    VectorXd xcoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd ycoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd zcoord_curr = VectorXd::Zero(nnel_lbv);

    for (int iterator = 0; iterator < nnel_lbv; iterator++) {
        xcoord_curr(iterator) = xcoord(iterator) + u(iterator * 3);
        ycoord_curr(iterator) = ycoord(iterator) + u((iterator * 3) + 1);
        zcoord_curr(iterator) = zcoord(iterator) + u((iterator * 3) + 2);
    }

    volume_current = fe_calVolume(xcoord_curr, ycoord_curr, zcoord_curr); // current volume of element

    lc = fe_minElementLength(xcoord_curr, ycoord_curr, zcoord_curr); // current element characteristic length

    c_wave = fe_calWaveSpeed(material_id, volume_initial, volume_current); // current dilatational wave speed

    rho_initial = fe_get_mats(material_id, 0, "mechanical"); // initial material density

    rho_current = rho_initial * (volume_initial / volume_current); // current material density

    pressure_scalar = 0.06 * rho_current * c_wave * lc * vol_strain_rate; // bulk viscosity pressure, scalar format

    pressure_matrix = pressure_scalar * I;

    pressure = fe_tensor2voigt(pressure_matrix); // bulk viscosity pressure, voigt format

}


void fe_getPressure_lbv_pbr_array(VectorXd& pressure, int i, int x, int y, int z, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id)
{
    // Reference: "Getting Started with Abaqus," Section 9.5.1 "Bulk Viscosity"
    // Linear bulk viscosity is always included to damp "ringing" in the highest element frequency.

    int nnel_lbv;
    double vol_strain_rate, volume_initial, volume_current, lc, c_wave, rho_initial, rho_current, pressure_scalar;

    MatrixXd F_curr = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_inv = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_invT = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_prev = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_dot = MatrixXd::Zero(ndof, ndof);
    MatrixXd F_dotT = MatrixXd::Zero(ndof, ndof);
    MatrixXd D = MatrixXd::Zero(ndof, ndof);
    MatrixXd pressure_matrix = MatrixXd::Zero(ndof, ndof);

    fe_calDefGrad_pbr_array(F_curr, i, x, y, z, u);

    fe_invMatrix_pbr(F_inv, F_curr);

    F_invT = F_inv.transpose();

    fe_calDefGrad_pbr_array(F_prev, i, x, y, z, u_prev);

    F_dot = (F_curr - F_prev) / dT;

    F_dotT = F_dot.transpose();

    D = 0.5 * ((F_dot * F_inv) + (F_invT * F_dotT));

    vol_strain_rate = D.trace(); // volumetric strain rate

    volume_initial = fe_calVolume(xcoord, ycoord, zcoord); // initial volume of element

    // calculate current nodal coordinates
    nnel_lbv = xcoord.size();

    VectorXd xcoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd ycoord_curr = VectorXd::Zero(nnel_lbv);
    VectorXd zcoord_curr = VectorXd::Zero(nnel_lbv);

    for (int iterator = 0; iterator < nnel_lbv; iterator++) {
        xcoord_curr(iterator) = xcoord(iterator) + u(iterator * 3);
        ycoord_curr(iterator) = ycoord(iterator) + u((iterator * 3) + 1);
        zcoord_curr(iterator) = zcoord(iterator) + u((iterator * 3) + 2);
    }

    volume_current = fe_calVolume(xcoord_curr, ycoord_curr, zcoord_curr); // current volume of element

    lc = fe_minElementLength(xcoord_curr, ycoord_curr, zcoord_curr); // current element characteristic length

    c_wave = fe_calWaveSpeed(material_id, volume_initial, volume_current); // current dilatational wave speed

    rho_initial = fe_get_mats(material_id, 0, "mechanical"); // initial material density

    rho_current = rho_initial * (volume_initial / volume_current); // current material density

    pressure_scalar = 0.06 * rho_current * c_wave * lc * vol_strain_rate; // bulk viscosity pressure, scalar format

    pressure_matrix = pressure_scalar * I;

    pressure = fe_tensor2voigt(pressure_matrix); // bulk viscosity pressure, voigt format

}
