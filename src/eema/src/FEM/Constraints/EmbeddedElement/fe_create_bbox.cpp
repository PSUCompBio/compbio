#include "functions.h"
using namespace Eigen;

MatrixXd fe_create_bbox(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host, double length){

    MatrixXi elements_host_tmp;
    
    double tol = (3.0)*length;

    double xmin = A(1) - tol;
    double xmax = A(1) + tol;
    double ymin = A(2) - tol;
    double ymax = A(2) + tol;
    double zmin = A(3) - tol;
    double zmax = A(3) + tol;

    double confirm;
    int index;

    for(int i=0;i<elements_host.rows();i++){
        for(int j=0;j<elements_host.cols()-2;j++){
            confirm = 0;

            index = fe_find(nodes_host.col(1),elements_host(i,j+2));
            x = nodes_host(x_i,1);
            y = nodes_host(x_i,2);
            z = nodes_host(x_i,3);

            if((x > xmin) && (x < xmax)){
                if((y > ymin) && (y < ymax)){
                    if((z > zmin) && (z<zmax)){
                        confirm = 1;
                    }
                }
            }

            if (confirm==1) break;

        }

        if(confirm==1){
            elements_host_tmp = fe_concatenate_vector2matrix(elements_host_tmp,elements_host.row(i),1);
        }

    }

    return elements_host_tmp;

}