#include"functions.h"


using namespace Eigen;

void fe_vtkWrite_host(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements){

	std::string name;

	std::ostringstream ss;
	ss << time_step;

	switch(mesh_choice){
	case 1:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 2:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 3:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 4:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 5:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	default:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	}

	std::ofstream myfile(name.c_str());

	myfile << "# vtk DataFile Version 3. 0.\n"; // 1st line - VTK file version
	myfile << output << "\n"; // Title of the result

	/* File Format */
	switch(format_choice){
	case 1:
		myfile << "ASCII" << "\n";
		break;
	case 2:
		myfile << "BINARY" << "\n";
		break;
	default:
		myfile << "ASCII" << "\n";
		break;
	}

	/* Mesh Type */
	switch(mesh_choice){
	case 1:
		myfile << "DATASET STRUCTURED_POINTS \n";
		break;
	case 2:
		myfile << "DATASET STRUCTURED_GRID \n";
		break;
	case 3:
		myfile << "DATASET RECTILLINEAR_GRID \n";
		break;
	case 4:
		myfile << "DATASET POLYDATA \n";
		break;
	case 5:
		myfile << "DATASET UNSTRUCTURED_GRID \n";
		break;
	default:
		myfile << "DATASET UNSTRUCTURED_GRID \n";
		break;
	}

	/* Unstructured Grid */
	if(mesh_choice==5){
		/* POINTS */
		myfile << "POINTS" << " " << nodes.rows() << " " << "double" << "\n";
        for(int i=0;i<nodes.rows();i++){
        	for(int j=1;j<nodes.cols();j++){
				myfile<<std::setw(3)<<nodes(i,j)<<" ";
        	}
        	myfile<<"\n";
        }

        /* CELLS */
        int size = 0;
        for(int k=0;k<elements.rows();k++){
        		size = size + 1 + (elements.row(k).cols()-2);
        }
        myfile << "CELLS" << " " << elements.rows() << " " << size << " " << "\n";
        for(int i=0;i<elements.rows();i++){
        	myfile<<elements.row(i).cols()-2<<" ";
        	for(int j=2;j<elements.cols();j++){
				myfile<<std::setw(3)<<elements(i,j)<<" ";
        	}
        	myfile<<"\n";
        }

        /* CELL TYPES */
        myfile << "CELL_TYPES" << " " << elements.rows() << "\n";
        for(int i=0;i<elements.rows();i++){
        			myfile << "12" << "\n";
        }
        myfile << "\n" ;

        /* NODAL DISPLACEMENTS */
        myfile << "POINT_DATA " << nodes.rows() << "\n";
        myfile << "VECTORS " << "displacement " << "double\n";
        int num = 0;
        for(int i=0;i<nodes.rows();i++){
        		myfile << U_host(num) << " " << U_host(num+1) << " " << U_host(num+2) << " \n";
        		num = num+3;
        }

        /* CELL DATA - STRESSES */
        myfile << "CELL_DATA " << elements.rows() << "\n";
        	myfile << "TENSORS " << "Stress " << "double\n";
        for(int i=0;i<element_stress_host.rows();i++){
        		myfile << element_stress_host(i,0) << " " << element_stress_host(i,1) << " " << element_stress_host(i,2) << "\n";
        		myfile << element_stress_host(i,3) << " " << element_stress_host(i,4) << " " << element_stress_host(i,5) << "\n";
        		myfile << element_stress_host(i,6) << " " << element_stress_host(i,7) << " " << element_stress_host(i,8) << "\n";
        		myfile << "\n";
        }

        /* CELL DATA - STRAINS */
        //myfile << "CELL_DATA " << elements.rows() << "\n";
        	myfile << "TENSORS " << "Strain " << "double\n";
        for(int i=0;i<element_strain_host.rows();i++){
        		myfile << element_strain_host(i,0) << " " << element_strain_host(i,1) << " " << element_strain_host(i,2) << "\n";
        		myfile << element_strain_host(i,3) << " " << element_strain_host(i,4) << " " << element_strain_host(i,5) << "\n";
        		myfile << element_strain_host(i,6) << " " << element_strain_host(i,7) << " " << element_strain_host(i,8) << "\n";
        		myfile << "\n";
        }

	}
}

void fe_vtkWrite_truss(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements){

	std::string name;

	std::ostringstream ss;
	ss << time_step;

	switch(mesh_choice){
	case 1:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 2:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 3:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 4:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	case 5:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	default:
		name = home_path + "results/vtk/" + output + "_" + ss.str() + ".vtk";
		break;
	}

	std::ofstream myfile(name.c_str());

	myfile << "# vtk DataFile Version 3. 0.\n"; // 1st line - VTK file version
	myfile << output << "\n"; // Title of the result

	/* File Format */
	switch(format_choice){
	case 1:
		myfile << "ASCII" << "\n";
		break;
	case 2:
		myfile << "BINARY" << "\n";
		break;
	default:
		myfile << "ASCII" << "\n";
		break;
	}

	/* Mesh Type */
	switch(mesh_choice){
	case 1:
		myfile << "DATASET STRUCTURED_POINTS \n";
		break;
	case 2:
		myfile << "DATASET STRUCTURED_GRID \n";
		break;
	case 3:
		myfile << "DATASET RECTILLINEAR_GRID \n";
		break;
	case 4:
		myfile << "DATASET POLYDATA \n";
		break;
	case 5:
		myfile << "DATASET UNSTRUCTURED_GRID \n";
		break;
	default:
		myfile << "DATASET UNSTRUCTURED_GRID \n";
		break;
	}

	/* Unstructured Grid */
	if(mesh_choice==5){
		/* POINTS */
		myfile << "POINTS" << " " << nodes.rows() << " " << "double" << "\n";
        for(int i=0;i<nodes.rows();i++){
        	for(int j=1;j<nodes.cols();j++){
				myfile<<std::setw(3)<<nodes(i,j)<<" ";
        	}
        	myfile<<"\n";
        }

        /* CELLS */
        int size = 0;
        for(int k=0;k<elements.rows();k++){
        		size = size + 1 + (elements.row(k).cols()-2);
        }
        myfile << "CELLS" << " " << elements.rows() << " " << size << " " << "\n";
        for(int i=0;i<elements.rows();i++){
        	myfile<<elements.row(i).cols()-2<<" ";
        	for(int j=2;j<elements.cols();j++){
				myfile<<std::setw(3)<<elements(i,j)<<" ";
        	}
        	myfile<<"\n";
        }

        /* CELL TYPES */
        myfile << "CELL_TYPES" << " " << elements.rows() << "\n";
        for(int i=0;i<elements.rows();i++){
        			myfile << "3" << "\n";
        }
        myfile << "\n" ;

        /* CELL DATA - STRESSES */
        myfile << "CELL_DATA " << elements.rows() << "\n";
        	myfile << "TENSORS " << "Stress " << "double\n";
        for(int i=0;i<element_stress_truss.rows();i++){
        		myfile << element_stress_truss(i,0) << " " << element_stress_truss(i,1) << " " << element_stress_truss(i,2) << "\n";
        		myfile << element_stress_truss(i,3) << " " << element_stress_truss(i,4) << " " << element_stress_truss(i,5) << "\n";
        		myfile << element_stress_truss(i,6) << " " << element_stress_truss(i,7) << " " << element_stress_truss(i,8) << "\n";
        		myfile << "\n";
        }

        /* CELL DATA - STRAINS */
        //myfile << "CELL_DATA " << elements.rows() << "\n";
        	myfile << "TENSORS " << "Strain " << "double\n";
        for(int i=0;i<element_strain_truss.rows();i++){
        		myfile << element_strain_truss(i,0) << " " << element_strain_truss(i,1) << " " << element_strain_truss(i,2) << "\n";
        		myfile << element_strain_truss(i,3) << " " << element_strain_truss(i,4) << " " << element_strain_truss(i,5) << "\n";
        		myfile << element_strain_truss(i,6) << " " << element_strain_truss(i,7) << " " << element_strain_truss(i,8) << "\n";
        		myfile << "\n";
        }

	}
}


