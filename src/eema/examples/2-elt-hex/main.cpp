#include "functions.h"
#include <mpi.h>
using namespace Eigen;
/**
 * This is the main file. If you want to submit a new job -- this is where you do it
 */
int main(int argc, char **argv){
	clock_t t;
	t = clock();
	int rank, size, nnel, nnelp;	
	/** Enter the path address for your job folder */
	home_path = argv[1];
	job_file = argv[2];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
        /** Read Mesh File **/
	fe_mainRead(home_path+"/"+job_file, size, rank);
	//VectorXi pe = mesh[0].getElementsProcessorID();
	//std::cout<<pe<<" from main";
	/** Main Explicit Loop **/
	fe_mainEXPLICIT(size, rank);
	t = clock()-t;	
	if(rank == 0)
	{
	std::cout<<"--------------------------------------"<<"\n";
	std::cout<<"Total Simulation CPU Time: "<<(((float)t)/CLOCKS_PER_SEC)<<"s \n";
	std::cout<<"Simulation Completed."<<"\n";
	std::cout<<"--------------------------------------"<<"\n";}
	MPI_Finalize();
	return 0;
}
