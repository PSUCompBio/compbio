#include "functions.h"

using namespace Eigen;

double eps_nr = 1e-6; /** Convergence criteria used in the newton rhapson method */
double iterations_nr = 100;

/**n \brief This functions calculates the isoparametric coordinates of a set of coordinates in global system */

VectorXd fe_newtonRhapson(VectorXd& nat_coord, VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord) {

    VectorXd coord = VectorXd::Zero((nat_coord.size()) * (xcoord.size()));
    int nnel = xcoord.size();
    double eps = 1e6;
    double iterations = 1;

    for (int i = 0; i < nnel; i++) {
        coord((i * 3)) = xcoord(i);
        coord((i * 3) + 1) = ycoord(i);
        coord((i * 3) + 2) = zcoord(i);
    }

    double edof = coord.size();
    VectorXd iso_coord = VectorXd::Zero(3); /** Vector showing the coordinates in iso-parametric system */
    VectorXd iso_coord_new = VectorXd::Zero(3);
    VectorXd shapes;

    if (nnel == 8) {
        shapes = fe_shapes_8(iso_coord(0), iso_coord(1), iso_coord(2));
    }

    MatrixXd shape_mat = fe_shapeMatrix(edof, nnel, shapes);
    VectorXd fx = VectorXd::Zero(3);

    VectorXd dndr = VectorXd::Zero(xcoord.size());
    VectorXd dnds = VectorXd::Zero(xcoord.size());
    VectorXd dndt = VectorXd::Zero(xcoord.size());
    MatrixXd jacobian;

    fe_dniso_8(dndr, dnds, dndt, iso_coord(0), iso_coord(1), iso_coord(2));

    jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
    MatrixXd transJacobian = jacobian.transpose();
    MatrixXd it_Jacobian = transJacobian.inverse();

    fx = (nat_coord - (shape_mat * coord));

    while ((eps >= eps_nr) || (iterations >= iterations_nr)) {

        iso_coord_new = iso_coord + (it_Jacobian * (fx));
        iso_coord = iso_coord_new;

        shapes = fe_shapes_8(iso_coord(0), iso_coord(1), iso_coord(2));
        shape_mat = fe_shapeMatrix(edof, nnel, shapes);
        fx = (nat_coord - (shape_mat * coord));

        eps = fx.norm();

        fe_dniso_8(dndr, dnds, dndt, iso_coord(0), iso_coord(1), iso_coord(2));
        jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
        transJacobian = jacobian.transpose();
        it_Jacobian = transJacobian.inverse();

        iterations = iterations + 1;

    }

    return iso_coord;
}