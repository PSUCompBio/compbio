#include "functions.h"
using namespace std;

void fe_crossProduct(VectorXd& result, VectorXd& A, VectorXd& B) {
	if (A.size() == 3 && B.size() == 3) {

		Vector3d tmp1(A(0), A(1), A(2));
		Vector3d tmp2(B(0), B(1), B(2));
		Vector3d normal = tmp1.cross(tmp2);
		for (int i = 0; i < 3; i++) {
			result(i) = normal(i);
		}

	}
}