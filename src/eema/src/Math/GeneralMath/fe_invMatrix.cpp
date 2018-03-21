#include "functions.h"

using namespace Eigen;

void fe_invMatrix_pbr(MatrixXd& A_inv, MatrixXd& A)
{

  int size_A = A.rows();
  int size_A_inv = A_inv.rows();
  double detA = fe_detMatrix_pbr(A);

  if (size_A == 2) {
    if (size_A_inv == 2) {
      A_inv(0,0) = 1/detA*(A(1,1));
      A_inv(0,1) = -1/detA*(A(0,1));

      A_inv(1,0) = -1/detA*(A(1,0));
      A_inv(1,1) = 1/detA*(A(0,0));
    }
  }

  if (size_A == 3) {
    if (size_A_inv == 3) {
      A_inv(0,0) = 1/detA*(A(1,1)*A(2,2) - A(1,2)*A(2,1));
      A_inv(0,1) = 1/detA*(A(0,2)*A(2,1) - A(0,1)*A(2,2));
      A_inv(0,2) = 1/detA*(A(0,1)*A(1,2) - A(0,2)*A(1,1));

      A_inv(1,0) = 1/detA*(A(1,2)*A(2,0) - A(1,0)*A(2,2));
      A_inv(1,1) = 1/detA*(A(0,0)*A(2,2) - A(0,2)*A(2,0));
      A_inv(1,2) = 1/detA*(A(0,2)*A(1,0) - A(0,0)*A(1,2));

      A_inv(2,0) = 1/detA*(A(1,0)*A(2,1) - A(1,1)*A(2,0));
      A_inv(2,1) = 1/detA*(A(0,1)*A(2,0) - A(0,0)*A(2,1));
      A_inv(2,2) = 1/detA*(A(0,0)*A(1,1) - A(0,1)*A(1,0));
    }
  }

  else {
    std::cout << "ALERT: PROBLEM CALCULATING INVERSE OF MATRIX." << '\n';
    std::cout << "MATRIX SIZE IS NOT 2X2 OR 3X3. SIMULATION CANCELLED." << '\n';
    std::exit(1);
  }

}

void fe_invMatrix_pbr_array(double** A_inv, double** A, double det)
{
    A_inv[0][0] = (1/det) * (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
    A_inv[0][1] = (1/det) * (A[0][2] * A[2][1] - A[0][1] * A[2][2]);
    A_inv[0][2] = (1/det) * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);

    A_inv[1][0] = (1/det) * (A[1][2] * A[2][0] - A[1][0] * A[2][2]);
    A_inv[1][1] = (1/det) * (A[0][0] * A[2][2] - A[0][2] * A[2][0]);
    A_inv[1][2] = (1/det) * (A[0][2] * A[1][0] - A[0][0] * A[1][2]);

    A_inv[2][0] = (1/det) * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    A_inv[2][1] = (1/det) * (A[0][1] * A[2][0] - A[0][0] * A[2][1]);
    A_inv[2][2] = (1/det) * (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
}
