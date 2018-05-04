#include "functions.h"
using namespace Eigen;

//This is the main file. If you want to submit a new job -- this is where you do it

int main(int argc, char **argv){

	clock_t t;
	t = clock();

	// Enter the path address for your job folder
	home_path = argv[1];
	job_file = argv[2];

	fe_mainRead(home_path+"/"+job_file);
	fe_mainEXPLICIT();

	t = clock()-t;
	std::cout<<"\n--------------------------------------"<<"\n";
	std::cout<<"Total Simulation CPU Time: "<<(((double)t)/CLOCKS_PER_SEC)<<"s \n";
	std::cout<<"Simulation Completed."<<"\n";
	std::cout<<"--------------------------------------"<<"\n";
	return 0;
}

/*
//#pragma omp parallel reduction(+:vec)
//#pragma omp for ordered schedule(static)

for (int intx_normal = 0; intx_normal < 4000; intx_normal++) {
	for (int inty_normal = 0; inty_normal < 4000; inty_normal++) {
		for (int intz_normal = 0; intz_normal < 4000; intz_normal++) {
			vec += vec1*vec1*2;
		}
	}
}

#include "functions.h"
#include <thread>
#include <mutex>

using namespace Eigen;

#pragma omp declare reduction (+: Eigen::VectorXd: omp_out=omp_out+omp_in) initializer(omp_priv=VectorXd::Zero(omp_orig.size()))

VectorXd vec = VectorXd::Zero(3);
VectorXd vec1 = VectorXd::Ones(3);
std::mutex m;

void adder(int start, int total) {

	VectorXd vec_par = VectorXd::Zero(3);

	for (int intx_normal = start; intx_normal < total + start; intx_normal++) {
		for (int inty_normal = 0; inty_normal < 4000; inty_normal++) {
			for (int intz_normal = 0; intz_normal < 4000; intz_normal++) {
				vec_par += vec1*vec1*2;
			}
		}
	}

	m.lock();
	vec += vec_par;
	m.unlock();
}

int main(int argc, char **argv){

	home_path = argv[1];
	job_file = argv[2];
	vec1 = 2*vec1;

	std::thread t[4];

	for (int intx_normal = 0; intx_normal < 4; intx_normal++) {
		t[intx_normal] = std::thread(adder, intx_normal*1000, 1000);
	}

	for (int intx_normal = 0; intx_normal < 4; intx_normal++) {
		t[intx_normal].join();
	}
}

	std::cout << "\n\n vec = \n" << vec;

	return 0;
}
*/
