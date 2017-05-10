#include "functions.h"

using namespace Eigen;

double fe_calculateKE(VectorXd mm, VectorXd V){
  double KE;
  VectorXd momentum = VectorXd::Zero(V.size());
  for(int i=0;i<V.size();i++){
    momentum(i) = mm(i)*V(i);
  }
  KE = ( 0.5 * V.transpose() * momentum );
  return KE;
}

double fe_calculateKE(MatrixXd mm, VectorXd V){
  double KE;
  KE = ( 0.5 * V.transpose() * mm * V );
  return KE;
}
