#include "functions.h"
using namespace Eigen;

int fe_find(VectorXd A, double a){

  int number = 100000000;

  for(int i=0;i<A.size();i++){
    if(A(i)==a){
      number = i;
      break;
    }
  }

  return number;
}

int fe_find(VectorXd A, int a){

  int number = 100000000;

  for(int i=0;i<A.size();i++){
    if(A(i)==(double)a){
      number = i;
      break;
    }
  }

  return number;
}
