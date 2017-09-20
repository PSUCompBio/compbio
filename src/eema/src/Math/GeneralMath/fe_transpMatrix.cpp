#include "functions.h"

using namespace Eigen;

void fe_transpMatrix_pbr(MatrixXd& A_T, MatrixXd& A)
{

  int size_A = A.rows();
  int size_A_T = A_T.rows();

  if (size_A == 2) {
    if (size_A_T == 2) {
      A_T(0,0) = A(0,0);
      A_T(0,1) = A(1,0);

      A_T(1,0) = A(0,1);
      A_T(1,1) = A(1,1);
    }
  }

  if (size_A == 3) {
    if (size_A_T == 3) {
      A_T(0,0) = A(0,0);
      A_T(0,1) = A(1,0);
      A_T(0,2) = A(2,0);

      A_T(1,0) = A(0,1);
      A_T(1,1) = A(1,1);
      A_T(1,2) = A(2,1);

      A_T(2,0) = A(0,2);
      A_T(2,1) = A(1,2);
      A_T(2,2) = A(2,2);
    }
  }

  else {
    std::cout << "ALERT: PROBLEM CALCULATING TRANSPOSE OF MATRIX." << '\n';
    std::cout << "MATRIX SIZE IS NOT 2X2 OR 3X3. SIMULATION CANCELLED." << '\n';
    std::exit(1);
  }

}
