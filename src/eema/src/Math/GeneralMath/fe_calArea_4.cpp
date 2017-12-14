#include"functions.h"

using namespace Eigen;

double fe_calArea_4(double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double c1, double c2, double c3, double c4) {


	//Area of a side is the sum of area of triangles
	// area of triangle is half the cross-product of two edges

	Vector3d tr1_side1(3);
	tr1_side1(0) = a1 - a2;
	tr1_side1(1) = b1 - b2;
	tr1_side1(2) = c1 - c2;

	Vector3d tr1_side2(3);
	tr1_side2(0) = a1 - a3;
	tr1_side2(1) = b1 - b3;
	tr1_side2(2) = c1 - c3;

	Vector3d tr2_side1(3);
	tr2_side1(0) = a1 - a3;
	tr2_side1(1) = b1 - b3;
	tr2_side1(2) = c1 - c3;

	Vector3d tr2_side2(3);
	tr2_side2(0) = a1 - a4;
	tr2_side2(1) = b1 - b4;
	tr2_side2(2) = c1 - c4;

	Vector3d area_tr1(3);
	area_tr1 = tr1_side1.cross(tr1_side2);
	Vector3d area_tr2(3);
	area_tr2 = tr2_side1.cross(tr2_side2);

	double area = 0.5 * (area_tr1.norm() + area_tr2.norm());

	return area;
}
