#include"functions.h"

using namespace Eigen;

/** This function converts vector in voigt's vector notation into a tensor
    The tensor will be either 2X2 or 3X3.
*/
MatrixXd fe_voigt2tensor(VectorXd B){
  MatrixXd A;

  if(B.size()==3){
    A = MatrixXd::Zero(2,2);
    A(0,0) = B(0);
    A(1,1) = B(1);
    A(0,1) = B(2);
    A(1,0) = B(2);
  }

  if(B.size()==6){
    A = MatrixXd::Zero(3,3);
    A(0,0) = B(0);
    A(1,1) = B(1);
    A(2,2) = B(2);
    A(0,1) = B(3);
    A(1,0) = B(3);
    A(1,2) = B(4);
    A(2,1) = B(5);
    A(0,2) = B(5);
    A(2,0) = B(5);
  }

  return A;
}
