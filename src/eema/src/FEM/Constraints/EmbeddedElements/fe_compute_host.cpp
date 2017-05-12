#include "functions.h"
using namespace Eigen;

int fe_compute_host(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host_tmp){

    int host_id;

    MatrixXi faces = MatrixXi::Zero(12,3);
    faces << 0, 1, 2,
            2, 3, 0,
            4, 5, 6,
            6, 7, 4,
            4, 5, 1,
            1, 0, 4,
            2, 6, 7,
            7, 3, 2,
            1, 5, 6,
            6, 2, 1,
            4, 0, 3,
            3, 7, 4;

    for(int i=0;i<elements_host_tmp.rows();i++){
        int counter = 0;
        for(int j=0;j<faces.rows();j++){
            int i_1, i_2, i_3;
            i_1 = fe_find(nodes_host,element_host_tmp.rows((faces(j,0))+2));
            i_2 = fe_find(nodes_host,element_host_tmp.rows((faces(j,1))+2));
            i_3 = fe_find(nodes_host,element_host_tmp.rows((faces(j,2))+2));

            VectorXd center;
            center << ((nodes_host(i_1,1)+nodes_host(i_2,1)+nodes_host(i_3,1))/3), ((nodes_host(i_1,2)+nodes_host(i_2,2)+nodes_host(i_3,2))/3), ((nodes_host(i_1,3)+nodes_host(i_2,3)+nodes_host(i_3,3))/3) ;

            VectorXd edge_1;
            VectorXd edge_2;
            edge_1 << (nodes_host(i_1,1)-nodes_host(i_2,1)), (nodes_host(i_1,2)-nodes_host(i_2,2)), (nodes_host(i_1,3)-nodes_host(i_2,3));
            edge_2 << (nodes_host(i_3,1)-nodes_host(i_2,1)), (nodes_host(i_3,2)-nodes_host(i_2,2)), (nodes_host(i_3,3)-nodes_host(i_2,3)); 

            VectorXd normal;
            normal = edge_1.cross(edge_2);

            VectorXd center_normal = (center(0) - A(1)), (center(1) - A(2)), (center(2) - A(3));

            double angle = center_normal.dot(normal);

            if(angle < 0) break;

            counter = counter + 1;
        }

        if(counter==12){
            host_id = elements_host_tmp(i,0);
            break;
        }
    }

    return host_id;

}