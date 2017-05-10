#include "functions.h"

using namespace Eigen;

/** For all elements -- this function calculates the minimum critical timestep */
double fe_getTimeStep(void)
{

    //std::cout << "Started Calculating Timestep..." << "\n";

    double deltaT_crit = 0;

    MatrixXd* nodes = mesh[0].getNewNodesPointer();
    MatrixXi* elements = mesh[0].getNewElementsPointer();

    VectorXd u = mesh[0].getNodalDisp();

    int nel   = mesh[0].getNumElements();       /*! number of elements */
    int nnel  = mesh[0].getNumNodesPerElement(); // number of nodes per element
    int nnode = mesh[0].getNumNodes();          // number of nodes
    int sdof  = nnode * ndof;          // system degrees of freedom
    int edof  = nnel * ndof;           // element degrees of freedom
    VectorXd* element_volumes = mesh[0].getElementCharacteristicPointer();

    VectorXd xcoord = VectorXd::Zero(nnel);
    VectorXd ycoord = VectorXd::Zero(nnel);
    VectorXd zcoord = VectorXd::Zero(nnel);

    VectorXd deltaT_element = VectorXd::Zero(nel);

    for (int i = 0; i < nel; i++) {

        //std::cout << "Here..." << i << "\n";

        for (int j = 0; j < nnel; j++) {
            int g = (*elements)(i, j + 2);
            xcoord(j) = (*nodes)(g, 1);
            ycoord(j) = (*nodes)(g, 2);
            zcoord(j) = (*nodes)(g, 3);
        }

        VectorXd u_e = VectorXd::Zero(edof);
        fe_gather_pbr(u, u_e, (*elements).block<1, 8>(i, 2), sdof);
        deltaT_element(i) = fe_calTimeStep(xcoord, ycoord, zcoord, (*elements)(i, 1), u_e, (*element_volumes)(i));
    }

    deltaT_crit = deltaT_element.minCoeff();

    deltaT_crit = deltaT_crit * reduction;

    if (deltaT_crit < failure_time_step) {
        std::cout << "Simulation Failed - Timestep too small" << "\n";
        std::cout << "Timestep is: " << deltaT_crit << "\n";
        std::exit(-1);
    }

    if (deltaT_crit > dt_min) {
        deltaT_crit = dt_min;
    }

    if (deltaT_crit > (t_end / ((double)output_frequency))) {
        deltaT_crit = (t_end / ((double)output_frequency));
    }

    nodes = NULL;
    elements = NULL;

    //std::cout << "Completed Calculating Timestep !!!" << "\n";

    return deltaT_crit;
} // fe_getTimeStep
