#include "functions.h"

using namespace Eigen;

double fe_calVolume(VectorXd& xcoord, VectorXd& ycoord, VectorXd& zcoord) {

    double volume = 0;

    if (xcoord.size() == 8) {
        // Area of a side is the sum of area of triangles
        // area of triangle is half the cross-product of two edges

        // TET-1
        Vector3d a1, a2, a3;
        a1 << (xcoord(0) - xcoord(7)), (ycoord(0) - ycoord(7)), (zcoord(0) - zcoord(7));
        a2 << (xcoord(2) - xcoord(7)), (ycoord(2) - ycoord(7)), (zcoord(2) - zcoord(7));
        a3 << (xcoord(3) - xcoord(7)), (ycoord(3) - ycoord(7)), (zcoord(3) - zcoord(7));
        double tet1_vol = a1.dot(a2.cross(a3));

        // TET-2
        Vector3d b1, b2, b3;
        b1 << (xcoord(0) - xcoord(4)), (ycoord(0) - ycoord(4)), (zcoord(0) - zcoord(4));
        b2 << (xcoord(5) - xcoord(4)), (ycoord(5) - ycoord(4)), (zcoord(5) - zcoord(4));
        b3 << (xcoord(7) - xcoord(4)), (ycoord(7) - ycoord(4)), (zcoord(7) - zcoord(4));
        double tet2_vol = b1.dot(b2.cross(b3));

        // TET-3
        Vector3d c1, c2, c3;
        c1 << (xcoord(0) - xcoord(5)), (ycoord(0) - ycoord(5)), (zcoord(0) - zcoord(5));
        c2 << (xcoord(2) - xcoord(5)), (ycoord(2) - ycoord(5)), (zcoord(2) - zcoord(5));
        c3 << (xcoord(1) - xcoord(5)), (ycoord(1) - ycoord(5)), (zcoord(1) - zcoord(5));
        double tet3_vol = c1.dot(c2.cross(c3));

        // TET-4
        Vector3d d1, d2, d3;
        d1 << (xcoord(7) - xcoord(5)), (ycoord(7) - ycoord(5)), (zcoord(7) - zcoord(5));
        d2 << (xcoord(2) - xcoord(5)), (ycoord(2) - ycoord(5)), (zcoord(2) - zcoord(5));
        d3 << (xcoord(6) - xcoord(5)), (ycoord(6) - ycoord(5)), (zcoord(6) - zcoord(5));
        double tet4_vol = d1.dot(d2.cross(d3));

        // TET-5
        Vector3d e1, e2, e3;
        e1 << (xcoord(7) - xcoord(5)), (ycoord(7) - ycoord(5)), (zcoord(7) - zcoord(5));
        e2 << (xcoord(2) - xcoord(5)), (ycoord(2) - ycoord(5)), (zcoord(2) - zcoord(5));
        e3 << (xcoord(0) - xcoord(5)), (ycoord(0) - ycoord(5)), (zcoord(0) - zcoord(5));
        double tet5_vol = e1.dot(e2.cross(e3));

        volume =
            (std::abs(tet1_vol) + std::abs(tet2_vol) + std::abs(tet3_vol) + std::abs(tet4_vol) + std::abs(tet5_vol)) / 6;
    }

    if (xcoord.size() == 2) { // Return length for a 1D element
        Vector3d vec;
        vec << (xcoord(1) - xcoord(0)), (ycoord(1) - ycoord(0)), (zcoord(1) - zcoord(0));

        volume = vec.norm();
    }

    return volume;
}
