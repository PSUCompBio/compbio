#include "functions.h"

using namespace Eigen;

double fe_detMatrix_pbr(MatrixXd& A)
{

  int size = A.rows();
  double detA;

  if (size == 2) {
    detA = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    return detA;
  }

  if (size == 3) {
    double tmp1 = A(0,0)*(A(1,1)*A(2,2) - A(1,2)*A(2,1));
    double tmp2 = A(0,1)*(A(1,0)*A(2,2) - A(1,2)*A(2,0));
    double tmp3 = A(0,2)*(A(1,0)*A(2,1) - A(1,1)*A(2,0));
    detA = tmp1 - tmp2 + tmp3;
    return detA;
  }

  else {
    std::cout << "ALERT: PROBLEM CALCULATING DETERMINANT OF MATRIX." << '\n';
    std::cout << "MATRIX SIZE IS NOT 2X2 OR 3X3. SIMULATION CANCELLED." << '\n';
    std::exit(1);
  }

}

double fe_detMatrix_pbr_array(double** A)
{
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])) - (A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])) + (A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]));
}
