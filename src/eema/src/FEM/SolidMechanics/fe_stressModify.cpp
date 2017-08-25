#include "functions.h"

using namespace Eigen;

void fe_stressModify(VectorXd& sigma_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed, int choice)
{

	// Eliminate all stress components, except normal stress along fiber axis.

  MatrixXd pk_S  = MatrixXd::Zero(ndof, ndof);
  pk_S = fe_voigt2tensor(sigma_embed);

  MatrixXd T = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 3);

  MatrixXd temp_stress  = MatrixXd::Zero(ndof, ndof);
  temp_stress = T.transpose() * pk_S * T;

  temp_stress(0,1) = 0;
  temp_stress(0,2) = 0;
  temp_stress(1,0) = 0;
  temp_stress(1,1) = 0;
  temp_stress(1,2) = 0;
  temp_stress(2,0) = 0;
  temp_stress(2,1) = 0;
  temp_stress(2,2) = 0;

  pk_S = T * temp_stress * T.transpose();

  sigma_embed = fe_tensor2voigt(pk_S);

}
