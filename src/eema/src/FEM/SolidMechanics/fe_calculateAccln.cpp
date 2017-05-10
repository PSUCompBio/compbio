#include "functions.h"

using namespace Eigen;

VectorXd fe_calculateAccln(VectorXd mm, VectorXd F_net){

  VectorXd A = VectorXd::Zero(F_net.size());

  for(int i=0;i<F_net.size();i++){

    A(i) = (F_net(i)/mm(i));

  }

  return A;
}
