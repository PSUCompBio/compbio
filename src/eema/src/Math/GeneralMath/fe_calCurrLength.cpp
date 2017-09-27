#include "functions.h"

using namespace Eigen;

double fe_calCurrLength_pbr(VectorXd& u_embed, VectorXd& xcoord_embed, VectorXd& ycoord_embed, VectorXd& zcoord_embed)
{

  int size = xcoord_embed.size();

  if (size == 2) {

    VectorXd xcoord_embed_curr(size);
    VectorXd ycoord_embed_curr(size);
    VectorXd zcoord_embed_curr(size);

    int counter = 0;

    for (int i = 0; i < size; i++) {
        xcoord_embed_curr(i) = xcoord_embed(i) + u_embed(counter);
        ycoord_embed_curr(i) = xcoord_embed(i) + u_embed(counter + 1);
        zcoord_embed_curr(i) = zcoord_embed(i) + u_embed(counter + 2);
        counter = counter + 3;
    }

    double length_curr = fe_calVolume(xcoord_embed_curr, ycoord_embed_curr, zcoord_embed_curr);

    return length_curr;

  }

}
