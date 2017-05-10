/*
 * fe_write.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: vsg111
 */

#include "functions.h"

using namespace Eigen;

void fe_writeElementStress(MatrixXd sigma_all, double time){
	std::ostringstream ss;
	ss << time;
	std::string name = "./results/elements/stress/time_"+ss.str()+".txt";
	matrix2text(name,sigma_all,0);
}




