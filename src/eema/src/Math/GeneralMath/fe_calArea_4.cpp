#include"functions.h"

using namespace Eigen;

double fe_calArea_4(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4) {

	//Area of a side is the sum of area of triangles
	// area of triangle is half the cross-product of two edges
	double area = 0;
	Vector3d tr1_side1, tr1_side2, tr2_side1, tr2_side2, area_tr1, area_tr2;

	tr1_side1 = Vector3d(3);
    tr1_side2 = Vector3d(3);
    tr2_side1 = Vector3d(3);
    tr2_side2 = Vector3d(3);
    area_tr1 = Vector3d(3);
    area_tr2 = Vector3d(3);

	tr1_side1(0) = x1 - x2;
	tr1_side1(1) = y1 - y2;
	tr1_side1(2) = z1 - z2;

	tr1_side2(0) = x1 - x3;
	tr1_side2(1) = y1 - y3;
	tr1_side2(2) = z1 - z3;

	tr2_side1(0) = x1 - x3;
	tr2_side1(1) = y1 - y3;
	tr2_side1(2) = z1 - z3;

	tr2_side2(0) = x1 - x4;
	tr2_side2(1) = y1 - y4;
	tr2_side2(2) = z1 - z4;

	area_tr1 = tr1_side1.cross(tr1_side2);
	area_tr2 = tr2_side1.cross(tr2_side2);

	area = 0.5 * (area_tr1.norm() + area_tr2.norm());

	return area;
}
