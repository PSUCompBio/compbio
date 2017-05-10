#include "functions.h"
using namespace Eigen;

/** For a single element - this function calculates the volume of the element and calculates the critical time step based on the wave speed.*/
double fe_calTimeStep(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord,int material_id){

	double deltaT;

	double lc = fe_minElementLength(xcoord,ycoord,zcoord);

	// std::cout<<"Volume_element: "<<volume_element<<"\n";
	// std::cout<<"Largest Face Area: "<<largest_face<<"\n";
	// std::cout<<"characteristic Length: "<<lc<<"\n";

	double c_wave = fe_calWaveSpeed(material_id);

	// std::cout<<"wave speed: "<<c_wave<<"\n";

	deltaT = lc/c_wave;

	//std::cout<<"in fe_calTimeStep\n"<<"lc = "<<lc<<"\n denominator: "<<sqrt(E/rho)<<"\n";

	return deltaT;
}
