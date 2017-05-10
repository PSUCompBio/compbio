#include "functions.h"
#include "Mesh.h"

using namespace Eigen;

void Mesh::readMesh(MatrixXd n, MatrixXi e){
	nodes = n;
	elements = e;
}

MatrixXd Mesh::getNodes(void){
	return nodes;
}

MatrixXi Mesh::getElements(void){
	return elements;
}

MatrixXd Mesh::getNewNodes(void){
	return nodes_new;
}

MatrixXi Mesh::getNewElements(void){
	return elements_new;
}

void Mesh::readNodalKinematics(VectorXd Usystem, VectorXd Vsystem, VectorXd Asystem){
	U = Usystem;
	V = Vsystem;
	A = Asystem;
}

VectorXd Mesh::getNodalDisp(void){
	return U;
}

VectorXd Mesh::getNodalVel(void){
	return V;
}

VectorXd Mesh::getNodalAcc(void){
	return A;
}

void Mesh::readStressStrain(MatrixXd stress_tmp, MatrixXd strain_tmp){
	stress = stress_tmp;
	strain = strain_tmp;
}

MatrixXd Mesh::getCellStress(){
	return stress;
}

MatrixXd Mesh::getCellStrain(){
	return strain;
}

void Mesh::preprocessMesh(void){

	nodes_new = nodes;
	elements_new = elements;
	/** Nodes Preprocessing - Putting the numbering in order */
	/* for(int i=0;i<nodes.rows();i++){
		nodes_new(i,0) = i;
	} */
	/** Elements Preprocessing - Correcting the element definitions */
	/* for(int i=0;i<elements.rows();i++){
		elements_new(i,0) = i;
		elements_new(i,1) = elements(i,1);
		for(int j=2;j<elements.cols();j++){
			elements_new(i,j) = fe_find(nodes.col(0),elements(i,j));
		}
	}*/
}

void Mesh::replaceNodes(MatrixXd new_nodes, std::string choice){
	if(choice=="old"){
		nodes = new_nodes;
	}
	else
	{
		nodes_new = new_nodes;
	}
}

void Mesh::replaceElements(MatrixXi new_elements, std::string choice){
	if(choice=="old"){
		elements = new_elements;
	}
	else
	{
		elements_new = new_elements;
	}
}

/** \brief Calculates the minimum charateristic length of the mesh */
VectorXd Mesh::getMinCharLength(std::string choice){

	VectorXd min_details;

	MatrixXd nodes_local;
	MatrixXi elements_local;

	if(choice=="old"){
		nodes_local = nodes;
		elements_local = elements;
	}
	else
	{
		nodes_local = nodes_new;
		elements_local = elements_new;
	}

	VectorXd xcoord;
    VectorXd ycoord;
    VectorXd zcoord;
    double min_length = 10000000;
	double id = 0;
    double lc;

    xcoord = VectorXd::Zero(elements_local.cols()-2);
    ycoord = VectorXd::Zero(elements_local.cols()-2);
    zcoord = VectorXd::Zero(elements_local.cols()-2);
    for(int i=0;i<elements_local.rows();i++){
        for(int j=0;j<(elements_local.cols()-2);j++){
			int g = -1;
			for(int f=0;f<nodes_local.rows();f++){
				if(elements_local(i,j+2)==nodes_local(f,0)){
					g = f;
					break;
				}
			}
			xcoord[j] = nodes_local(g,1);
			ycoord[j] = nodes_local(g,2);
			zcoord[j] = nodes_local(g,3);
		}

        lc = fe_minElementLength(xcoord,ycoord,zcoord);

        if(lc<min_length){
            min_length = lc;
			id = i;
        }

    }

	min_details << min_length,id;

	return min_details;
}

/** \brief Calculates the maximum charateristic length of the mesh */
VectorXd Mesh::getMaxCharLength(std::string choice){

	VectorXd max_details;
	double id;

	MatrixXd nodes_local;
	MatrixXi elements_local;

	if(choice=="old"){
		nodes_local = nodes;
		elements_local = elements;
	}
	else
	{
		nodes_local = nodes_new;
		elements_local = elements_new;
	}

	VectorXd xcoord;
    VectorXd ycoord;
    VectorXd zcoord;
    double max_length = 0.0;
    double lc;

    xcoord = VectorXd::Zero(elements_local.cols()-2);
    ycoord = VectorXd::Zero(elements_local.cols()-2);
    zcoord = VectorXd::Zero(elements_local.cols()-2);
    for(int i=0;i<elements_local.rows();i++){
        for(int j=0;j<(elements_local.cols()-2);j++){
			int g = -1;
			for(int f=0;f<nodes_local.rows();f++){
				if(elements_local(i,j+2)==nodes_local(f,0)){
					g = f;
					break;
				}
			}
			xcoord[j] = nodes_local(g,1);
			ycoord[j] = nodes_local(g,2);
			zcoord[j] = nodes_local(g,3);
		}

        lc = fe_maxElementLength(xcoord,ycoord,zcoord);

        if(lc>max_length){
            max_length = lc;
			id = i;
        }

    }

	max_details << max_length, id;
	return max_details;
}

/** \brief Check the mesh for zero charateristic lengths */
void Mesh::checkMesh(){
	VectorXd lc;

	lc = getMinCharLength("new");

	if(lc(0) == 0){
    	std::cout << "**********************************************" << std::endl;
    	std::cout << "ZERO LENGTH ELEMENT IN THE SYSTEM \n FIBER LENGTH PREPROCESSING NOT POSSIBLE" << std::endl;
        std::cout << "**********************************************" << std::endl;
        std::exit(-1);
    }

}
