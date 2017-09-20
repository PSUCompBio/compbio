

#include"functions.h"

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

Mat fe_invMatrix(Mat matrix){

    double ratio,a;
    int i, j, k, n;
    //printf("Enter order of matrix: ");
    //scanf("%d", &n);
    //Vec Row(n);
    //printf("Enter the matrix: \n");
    //for(i = 0; i < n; i++){
    //    for(j = 0; j < n; j++){
    //		printf("I am here: \n");
    //        	scanf("%lf",&Row[j]);
    //    }
    //	matrix.push_back(Row);
    //}

	n = matrix.size(); // order of the matrix

    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
            if(i==(j-n))
                matrix[i][j] = 1.0;
            else
                matrix[i][j] = 0.0;
        }
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i!=j){
                ratio = matrix[j][i]/matrix[i][i];
                for(k = 0; k < 2*n; k++){
                    matrix[j][k] -= ratio * matrix[i][k];
                }
            }
        }
    }
    for(i = 0; i < n; i++){
        a = matrix[i][i];
        for(j = 0; j < 2*n; j++){
            matrix[i][j] /= a;
        }
    }
    printf("The inverse matrix is: \n");

    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
		if(isnan(matrix[i][j])){
			goto stop;
		}
        }
    }

    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
            printf("%.2f", matrix[i][j]);
            printf("\t");
        }
        printf("\n");
    }

	stop: 
	{
	cout<<"Inverse doesn't exist \n";
	Mat empty;
	matrix.swap(empty);
	}

	return matrix;
}
