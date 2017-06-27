#include "functions.h"
using namespace Eigen;
/**
 * This is the main file. If you want to submit a new job -- this is where you do it
 */
int main(int argc, char **argv){
	clock_t t;
	t = clock();

	/** Enter the path address for your job folder */
	home_path = ".";
	controls_file = argv[1];
	mesh_file = argv[2];

	fe_mainRead(home_path+"/"+controls_file, home_path+"/"+mesh_file);
	fe_mainEXPLICIT();

	t = clock()-t;
	std::cout<<"--------------------------------------"<<"\n";
	std::cout<<"Total Simulation CPU Time: "<<(((float)t)/CLOCKS_PER_SEC)<<"s \n";
	std::cout<<"Simulation Completed."<<"\n";
	std::cout<<"--------------------------------------"<<"\n";
	return 0;
}
