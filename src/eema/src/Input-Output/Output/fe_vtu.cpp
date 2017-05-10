#include "functions.h"

using namespace Eigen;

void fe_vtuWrite(std::string output, int time_step, Mesh mesh1){

    /** Points Info */
    MatrixXd nodes = mesh1.getNewNodes();

    /** Cells Info */
    MatrixXi elements = mesh1.getNewElements();

    /** Point Data - Displacements */
    VectorXd U_mesh = mesh1.getNodalDisp();

    /** Points Data - Velocities */
    VectorXd V_mesh = mesh1.getNodalVel();

    /** Points Data - Accelerations */
    VectorXd A_mesh = mesh1.getNodalAcc();

    /** Mesh Data - Stresses */
    MatrixXd element_stress = mesh1.getCellStress();

    /** Mesh Data - Strains */
    MatrixXd element_strain = mesh1.getCellStrain();

    /** Output File Name */
    std::string name;
    std::ostringstream ss;
    ss << time_step;
    name = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
    std::ofstream myfile(name.c_str());

    myfile << "<?xml version=\"1.0\"?>\n";
    myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    myfile << "\t<UnstructuredGrid>\n";
    myfile << "\t\t<Piece NumberOfPoints=\"" << nodes.rows() << "\" NumberOfCells=\"" << elements.rows() << "\">\n";

    /** Points Data */
    myfile << "\t\t\t<Points>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0;i<nodes.rows();i++){
	myfile<<"\t\t\t\t\t"<<std::setw(10)<<std::scientific<<std::setprecision(8)
		<<nodes(i,1)<<" "
		<<nodes(i,2)<<" "
		<<nodes(i,3)<<"\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</Points>\n";

    /** Cell Data */
    myfile << "\t\t\t<Cells>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int i=0;i<elements.rows();i++){
                myfile << "\t\t\t\t\t" ;
        	    for(int j=2;j<elements.cols();j++){
				          myfile<<elements(i,j)<<"\t";
                }
    myfile << "\n";
    }
    myfile<<"\n";
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int tmp = elements.row(0).cols()-2;
    for(int i=0;i<elements.rows();i++){
            myfile << "\t\t\t\t\t" << tmp << "\n";
            tmp = tmp + (elements.row(i).cols()-2);
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for(int i=0;i<elements.rows();i++){
        if((elements.row(i).cols()-2)==8) myfile << "\t\t\t\t\t" << "12" << "\n";
        if((elements.row(i).cols()-2)==2) myfile << "\t\t\t\t\t" << "3" << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t</Cells>\n";

    /** Point Vector Data - Displacements */
    myfile << "\t\t\t<PointData>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Displacement\" "
		<<"NumberOfComponents=\"3\" ComponentName0=\"X\" "
		<<"ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n" ;
    int num = 0;
    for(int i=0;i<nodes.rows();i++){
        myfile << "\t\t\t\t\t" << std::setw(10)<<std::scientific << std::setprecision(10)
		<< U_mesh(num) << " "
		<< U_mesh(num+1) << " "
		<< U_mesh(num+2) << " \n";
        	num = num+3;
    }
    myfile << "\t\t\t\t</DataArray>\n";

    /** Point Vector Data - Velocities */
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" "
		<<"NumberOfComponents=\"3\" ComponentName0=\"X\" "
		<<"ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n" ;
    num = 0;
    for(int i=0;i<nodes.rows();i++){
        myfile << "\t\t\t\t\t" << std::setw(10)<<std::scientific << std::setprecision(10)
		<< V_mesh(num) << " "
		<< V_mesh(num+1) << " "
		<< V_mesh(num+2) << " \n";
        	num = num+3;
    }
    myfile << "\t\t\t\t</DataArray>\n";

    /** Point Vector Data - Accelerations */
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Acceleration\" "
		<<"NumberOfComponents=\"3\" ComponentName0=\"X\" "
		<<"ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n" ;
    num = 0;
    for(int i=0;i<nodes.rows();i++){
        myfile << "\t\t\t\t\t" << std::setw(10)<<std::scientific << std::setprecision(10)
		<< A_mesh(num) << " "
		<< A_mesh(num+1) << " "
		<< A_mesh(num+2) << " \n";
        	num = num+3;
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</PointData>\n";

    /** Cell Data - Stresses and Strains */
    myfile << "\t\t\t<CellData>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Stress (Cauchy)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
    for(int i=0;i<element_stress.rows();i++){
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress(i,0) << " " << element_stress(i,1) << " " << element_stress(i,2) << "\n";
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress(i,3) << " " << element_stress(i,4) << " " << element_stress(i,5) << "\n";
        	myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress(i,6) << " " << element_stress(i,7) << " " << element_stress(i,8) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Strain (Logarithmic)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
    for(int i=0;i<element_strain.rows();i++){
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain(i,0) << " " << element_strain(i,1) << " " << element_strain(i,2) << "\n";
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain(i,3) << " " << element_strain(i,4) << " " << element_strain(i,5) << "\n";
        	myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain(i,6) << " " << element_strain(i,7) << " " << element_strain(i,8) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</CellData>\n";

    myfile << "\t\t</Piece>\n";
    myfile << "\t</UnstructuredGrid>\n";
    myfile << "</VTKFile>\n";
}

void fe_pvdNew(std::string output, int time_step, double time){

    /** Output File Name */
    std::string name;
    name = home_path + "/" + "results/vtu/" + output + ".pvd";
    std::ofstream myfile(name.c_str());

    /** Write Header */
    myfile << "<?xml version=\"1.0\"?>\n";
    myfile << "<VTKFile type=\"Collection\" version=\"0.1\"\n";
    myfile << "\t byte_order=\"LittleEndian\"\n";
    myfile << "\t compressor=\"vtkZLibDataCompressor\">\n";
    myfile << "   <Collection>\n";

    /** Write First Line of Data */
    myfile << "     <DataSet timestep= \"" << std::scientific << std::setprecision(5) << time << "\" group=\"\" part=\"0\"\n";
    std::ostringstream ss;
    ss << time_step;
    std::string vtuFileName;
    //vtuFileName = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
    vtuFileName = output + "_" + ss.str() + ".vtu";
    myfile << "\t      file=\"" << vtuFileName << "\"/>\n";
}

void fe_pvdAppend(std::string output, int time_step, double time){

    /** Output File Name */
    std::string name;
    name = home_path + "/" + "results/vtu/" + output + ".pvd";
    std::ofstream myfile(name.c_str(), std::ios_base::app);

    /** Append New Line of Data */
    if(myfile.is_open()){
    	myfile << "     <DataSet timestep= \"" << std::scientific << std::setprecision(5) << time << "\" group=\"\" part=\"0\"\n";
    	std::ostringstream ss;
    	ss << time_step;
    	std::string vtuFileName;
    	//vtuFileName = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
	vtuFileName = output + "_" + ss.str() + ".vtu";
    	myfile << "\t      file=\"" << vtuFileName << "\"/>\n";

    	/** Write Footer */
    	if(time == t_end){
	    myfile << "   </Collection>\n";
	    myfile << "</VTKFile>\n";
    	}
    }

    else{
    	std::cout << "No Such Text File With This Name Exists - " << name << "\n";
    	std::exit(-1);
    }
}
