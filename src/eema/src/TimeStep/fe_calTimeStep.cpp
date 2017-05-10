#include "functions.h"
using namespace Eigen;

/** For a single element - this function calculates the volume of the element and calculates the critical time step based on the wave speed.*/
double fe_calTimeStep(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord, int material_id, VectorXd& u_e, double volume_initial)
{
    double deltaT;

    int counter;

    for (int j = 0; j < xcoord.size(); j++) {
        counter   = j * 3;
        xcoord(j) = xcoord(j) + u_e(counter);
        ycoord(j) = ycoord(j) + u_e(counter + 1);
        zcoord(j) = zcoord(j) + u_e(counter + 2);
    }
    double volume_current = fe_calVolume(xcoord, ycoord, zcoord);

    double lc = fe_minElementLength(xcoord, ycoord, zcoord);

    double c_wave = fe_calWaveSpeed(material_id, volume_initial, volume_current);

    deltaT = lc / c_wave;

    return deltaT;
}
