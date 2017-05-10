#include "functions.h"

using namespace Eigen;

void fe_calculateKE(double& KE, VectorXd& mm, VectorXd& V) {

	VectorXd momentum = VectorXd::Zero(V.size());
	for (int i = 0; i < V.size(); i++) {
		momentum(i) = mm(i) * V(i);
	}
	KE = ( 0.5 * V.transpose() * momentum );

}

double fe_calculateKE(MatrixXd mm, VectorXd V) {
	double KE;
	KE = ( 0.5 * V.transpose() * mm * V );
	return KE;
}
