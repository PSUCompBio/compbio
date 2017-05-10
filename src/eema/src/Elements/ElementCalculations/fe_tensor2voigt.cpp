#include"functions.h"

using namespace Eigen;

/** This function converts tensor into voigt's vector notation
    The tensor should be either 2X2 or 3X3.
*/
VectorXd fe_tensor2voigt(MatrixXd A){

  VectorXd B; /** The tensor should be symmetric for its transformation into Voigt Vector */

  if(A.rows()==2){
    B = VectorXd::Zero(3);
    B(0) = A(0,0);
    B(1) = A(1,1);
    B(2) = A(0,1);
  }

  if(A.rows()==3){
    B = VectorXd::Zero(6);
    B(0) = A(0,0);
    B(1) = A(1,1);
    B(2) = A(2,2);
    B(3) = A(0,1);
    B(4) = A(1,2);
    B(5) = A(0,2);
  }

  return B;

}
