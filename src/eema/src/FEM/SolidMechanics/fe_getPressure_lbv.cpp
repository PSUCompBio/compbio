#include "functions.h"

using namespace Eigen;

void fe_getPressure_lbv_pbr(VectorXd& pressure, VectorXd& dndx, VectorXd& dndy, VectorXd& dndz, VectorXd& u, VectorXd& u_prev, double dT, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id)
{
    // Reference: "Getting Started with Abaqus," Section 9.5.1 "Bulk Viscosity"
    // Linear bulk viscosity is always included to damp "ringing" in the highest element frequency.

    double b1 = 0.06; // damping coefficient, Abaqus default value = 0.06

    MatrixXd F_curr = MatrixXd::Zero(ndof, ndof); // current deformation gradient
    fe_calDefGrad_pbr(F_curr, dndx, dndy, dndz, u);

    MatrixXd F_inv = MatrixXd::Zero(ndof, ndof); // inverse deformation gradient
    F_inv = F_curr.inverse();

    MatrixXd F_invT = MatrixXd::Zero(ndof, ndof); // inverse deformation gradient, transposed
    F_invT = F_inv.transpose();

    MatrixXd F_prev = MatrixXd::Zero(ndof, ndof); // previous deformation gradient
    fe_calDefGrad_pbr(F_prev, dndx, dndy, dndz, u_prev);

    MatrixXd F_dot = MatrixXd::Zero(ndof, ndof); // material time derivative of deformation gradient
    F_dot = (F_curr - F_prev) / dT;

    MatrixXd F_dotT = MatrixXd::Zero(ndof, ndof); // material time derivative of deformation gradient,transposed
    F_dotT = F_dot.transpose();

    double one_half   = 1.0 / 2.0;

    MatrixXd D = MatrixXd::Zero(ndof, ndof); // rate-of-deformation tensor
    D = one_half * ((F_dot * F_inv) + (F_invT * F_dotT));

    double vol_strain_rate = D.trace(); // volumetric strain rate

    double volume_initial = fe_calVolume(xcoord, ycoord, zcoord); // initial volume of element

    // calculate current nodal coordinates
    int nnel = xcoord.size();

    VectorXd xcoord_curr = VectorXd::Zero(nnel);
    VectorXd ycoord_curr = VectorXd::Zero(nnel);
    VectorXd zcoord_curr = VectorXd::Zero(nnel);

    int counter;

    for (int i = 0; i < nnel; i++) {
        counter   = i * 3;
        xcoord_curr(i) = xcoord(i) + u(counter);
        ycoord_curr(i) = ycoord(i) + u(counter + 1);
        zcoord_curr(i) = zcoord(i) + u(counter + 2);
    }

    double volume_current = fe_calVolume(xcoord_curr, ycoord_curr, zcoord_curr); // current volume of element

    double lc = fe_minElementLength(xcoord_curr, ycoord_curr, zcoord_curr); // current element characteristic length

    double c_wave = fe_calWaveSpeed(material_id, volume_initial, volume_current); // current dilatational wave speed

    double rho_initial = fe_get_mats(material_id, 0, "mechanical"); // initial material density

    double rho_current = rho_initial * (volume_initial / volume_current); // current material density

    double pressure_scalar = b1 * rho_current * c_wave * lc * vol_strain_rate; // bulk viscosity pressure, scalar format

    MatrixXd I = MatrixXd::Identity(ndof, ndof);

    MatrixXd pressure_matrix = MatrixXd::Zero(ndof, ndof); // bulk viscosity pressure, matrix format

    pressure_matrix = pressure_scalar * I;

    pressure = fe_tensor2voigt(pressure_matrix); // bulk viscosity pressure, voigt format

}
