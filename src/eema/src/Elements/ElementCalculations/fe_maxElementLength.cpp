#include "functions.h"

using namespace Eigen;

double fe_maxElementLength(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord){

    double lc = 0;
	int nnel = xcoord.size();

    if(nnel==8){
        VectorXd face_area(6);
	    // Face-1 (0, 1, 2, 3)
	    face_area(0) = fe_calArea_4(xcoord(0),xcoord(1),xcoord(2),xcoord(3),ycoord(0),ycoord(1),ycoord(2),ycoord(3),zcoord(0),zcoord(1),zcoord(2),zcoord(3)); //  Enter the coordinates in CCW direction
	    // Face-2 (4, 5, 6, 7)
	    face_area(1) = fe_calArea_4(xcoord(4),xcoord(5),xcoord(6),xcoord(7),ycoord(4),ycoord(5),ycoord(6),ycoord(7),zcoord(4),zcoord(5),zcoord(6),zcoord(7));
	    // Face-3 (0, 1, 5, 4)
	    face_area(2) = fe_calArea_4(xcoord(0),xcoord(1),xcoord(5),xcoord(4),ycoord(0),ycoord(1),ycoord(5),ycoord(4),zcoord(0),zcoord(1),zcoord(5),zcoord(4));
	    // Face-4 (3, 2, 6, 7)
	    face_area(3) = fe_calArea_4(xcoord(3),xcoord(2),xcoord(6),xcoord(7),ycoord(3),ycoord(2),ycoord(6),ycoord(7),zcoord(3),zcoord(2),zcoord(6),zcoord(7));
	    // Face-5 (0, 4, 7, 3)
	    face_area(4) = fe_calArea_4(xcoord(0),xcoord(4),xcoord(7),xcoord(3),ycoord(0),ycoord(4),ycoord(7),ycoord(3),zcoord(0),zcoord(4),zcoord(7),zcoord(3));
	    // Face-6 (1, 5, 6, 2)
	    face_area(5) = fe_calArea_4(xcoord(1),xcoord(5),xcoord(6),xcoord(2),ycoord(1),ycoord(5),ycoord(6),ycoord(2),zcoord(1),zcoord(5),zcoord(6),zcoord(2));

	    double volume_element = fe_calVolume(xcoord,ycoord,zcoord);
	    double largest_face = face_area.minCoeff();
	    lc = ((volume_element)/(largest_face));
    }

    if(nnel==2){
        /* A 2 noded element - truss, beam etc. */
        lc = sqrt( pow((xcoord(1)-xcoord(0)),2) + pow((ycoord(1)-ycoord(0)),2) + pow((zcoord(1)-zcoord(0)),2));
    }

        return lc;
}