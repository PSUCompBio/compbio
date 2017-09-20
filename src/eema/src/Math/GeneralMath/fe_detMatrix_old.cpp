

#include"functions.h"

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

double d = 0;

double fe_detMatrix(Mat matrix){

	int rows = matrix.size();

	if(rows==2){
		return ((matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0]));
	}
	else
	{
		Mat submat;
		Vec Row(rows-1);
		for(int i=0;i<rows;i++){
			int subi = 0;
			for(int j=1;j<rows;j++){
				int subj = 0;
				for (int k=0;k<rows;k++){
					if(k==i){continue;}
					Row[subj] = matrix[j][k];
					subj++;
				}
			submat.push_back(Row);
			subi++;
			}
		d = d + (pow(-1,i)*matrix[0][i]*fe_detMatrix(submat));
		//cout<<"submat size and d are: "<<submat.size()<<"\t"<<d<<endl;
		Mat empty;
		submat.swap(empty);
		//cout<<"i is: "<<i<<endl;
		//cout<<"submat size is: "<<submat.size()<<endl;
		}
	}
	
	return d;

}
