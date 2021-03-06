#include "functions.h"

using namespace Eigen;

void fe_vtuWrite(int time_step, double time, Mesh& mesh1)
{

    std::string output = mesh1.getName();

    /** Points Info */
    MatrixXd* nodes = mesh1.getNewNodesPointer();

    /** Cells Info */
    MatrixXi* elements = mesh1.getNewElementsPointer();

    /** Point Data - Displacements */
    VectorXd* U_mesh = mesh1.getNodalDispPointer();

    /** Points Data - Velocities */
    VectorXd* V_mesh = mesh1.getNodalVelPointer();

    /** Points Data - Accelerations */
    VectorXd* A_mesh = mesh1.getNodalAccPointer();

    /** Mesh Data - Stresses */
    VectorXd* element_stress = mesh1.getCellStressPointer();

    /** Mesh Data - Strains */
    VectorXd* element_strain = mesh1.getCellStrainPointer();

    VectorXd* electric_potential = mesh1.getNodalEPotentialPointer();

    /** determine if damage should be printed */
    bool print_fiber_damage = 0;
    bool print_host_damage = 0;
    std::string mesh_name = mesh1.getName();

    if (embedded_constraint == 1) {
      for (int i = 0; i < num_constraints; i++) {
        std::string constraint_name = cons[i].getName();
        if (constraint_name == "embedded") {
          std::string constraint_slave = cons[i].get_EmbedSlave();
          if (constraint_slave == mesh_name) {
            print_fiber_damage = 1;
          }
          std::string constraint_master = cons[i].get_EmbedMaster();
          if (constraint_master == mesh_name) {
            print_host_damage = 1;
          }
        }
      }
    }

    /** Output File Name */
    std::string name;
    std::ostringstream ss;
    ss << time_step;
    name = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
    std::ofstream myfile(name.c_str());

    myfile << "<?xml version=\"1.0\"?>\n";
    myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    myfile << "\t<UnstructuredGrid>\n";
    myfile << "\t\t<Piece NumberOfPoints=\"" << (*nodes).rows() << "\" NumberOfCells=\"" << (*elements).rows() << "\">\n";

    /** Points Data */
    myfile << "\t\t\t<Points>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < (*nodes).rows(); i++)
    {
        myfile << "\t\t\t\t\t" << std::setw(10) << std::scientific << std::setprecision(8)
               << (*nodes)(i, 1) << " "
               << (*nodes)(i, 2) << " "
               << (*nodes)(i, 3) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</Points>\n";

    /** Cell Data */
    myfile << "\t\t\t<Cells>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < (*elements).rows(); i++)
    {
        myfile << "\t\t\t\t\t";
        for (int j = 2; j < (*elements).cols(); j++)
        {
            myfile << (*elements)(i, j) << "\t";
        }
        myfile << "\n";
    }
    myfile << "\n";
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int tmp = (*elements).row(0).cols() - 2;
    for (int i = 0; i < (*elements).rows(); i++)
    {
        myfile << "\t\t\t\t\t" << tmp << "\n";
        tmp = tmp + ((*elements).row(i).cols() - 2);
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < (*elements).rows(); i++)
    {
        if (((*elements).row(i).cols() - 2) == 8)
            myfile << "\t\t\t\t\t"
                   << "12"
                   << "\n";
        if (((*elements).row(i).cols() - 2) == 2)
            myfile << "\t\t\t\t\t"
                   << "3"
                   << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t</Cells>\n";

    /** Point Vector Data - Displacements */
    myfile << "\t\t\t<PointData>\n";

    if ((*U_mesh).size() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Displacement\" "
               << "NumberOfComponents=\"3\" ComponentName0=\"X\" "
               << "ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n";
        int num = 0;
        for (int i = 0; i < (*nodes).rows(); i++)
        {
            myfile << "\t\t\t\t\t" << std::setw(10) << std::scientific << std::setprecision(10)
                   << (*U_mesh)(num) << " "
                   << (*U_mesh)(num + 1) << " "
                   << (*U_mesh)(num + 2) << " \n";
            num = num + 3;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }

    /** Point Vector Data - Velocities */
    if ((*V_mesh).size() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Velocity\" "
               << "NumberOfComponents=\"3\" ComponentName0=\"X\" "
               << "ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n";
        int num = 0;
        for (int i = 0; i < (*nodes).rows(); i++)
        {
            myfile << "\t\t\t\t\t" << std::setw(10) << std::scientific << std::setprecision(10)
                   << (*V_mesh)(num) << " "
                   << (*V_mesh)(num + 1) << " "
                   << (*V_mesh)(num + 2) << " \n";
            num = num + 3;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }

    /** Point Vector Data - Accelerations */
    if ((*A_mesh).size() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Acceleration\" "
               << "NumberOfComponents=\"3\" ComponentName0=\"X\" "
               << "ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n";
        int num = 0;
        for (int i = 0; i < (*nodes).rows(); i++)
        {
            myfile << "\t\t\t\t\t" << std::setw(10) << std::scientific << std::setprecision(10)
                   << (*A_mesh)(num) << " "
                   << (*A_mesh)(num + 1) << " "
                   << (*A_mesh)(num + 2) << " \n";
            num = num + 3;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }


    /** Point Scalar Data - Electric Potential */
    if ((*electric_potential).size() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Electric-Potential\" "
               << "format=\"ascii\">\n";
        int num = 0;
        for (int i = 0; i < (*nodes).rows(); i++)
        {
            myfile << "\t\t\t\t\t" << std::setw(10) << std::scientific << std::setprecision(10)
                   << (*electric_potential)(num) << "\n";
            num = num + 1;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }

    myfile << "\t\t\t</PointData>\n";

    /** Cell Data - Stresses, Strains, and Damage */
    myfile << "\t\t\t<CellData>\n";

    if ((*element_stress).rows() != 0 && (*element_stress).cols() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Stress (Cauchy)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
        for (int i = 0; i < (*element_stress).size();)
        {
            VectorXd tmp(9);
            tmp << (*element_stress)(i + 0), (*element_stress)(i + 1), (*element_stress)(i + 2), (*element_stress)(i + 3), (*element_stress)(i + 4), (*element_stress)(i + 5), (*element_stress)(i + 6), (*element_stress)(i + 7), (*element_stress)(i + 8);

            for (int j = 0; j < 9; j++) {
              if (std::abs(tmp(j)) < 1e-30) {
                tmp(j) = 0;
              }
            }

            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(0) << " " << tmp(1) << " " << tmp(2) << "\n";
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(3) << " " << tmp(4) << " " << tmp(5) << "\n";
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(6) << " " << tmp(7) << " " << tmp(8) << "\n";

            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_stress)(i + 0) << " " << (*element_stress)(i + 1) << " " << (*element_stress)(i + 2) << "\n";
            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_stress)(i + 3) << " " << (*element_stress)(i + 4) << " " << (*element_stress)(i + 5) << "\n";
            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_stress)(i + 6) << " " << (*element_stress)(i + 7) << " " << (*element_stress)(i + 8) << "\n";

            i = i + 9;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }

    if ((*element_strain).rows() != 0 && (*element_strain).cols() != 0) {
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Strain (Logarithmic)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
        for (int i = 0; i < (*element_strain).size();)
        {
            VectorXd tmp(9);
            tmp << (*element_strain)(i + 0), (*element_strain)(i + 1), (*element_strain)(i + 2), (*element_strain)(i + 3), (*element_strain)(i + 4), (*element_strain)(i + 5), (*element_strain)(i + 6), (*element_strain)(i + 7), (*element_strain)(i + 8);

            for (int j = 0; j < 9; j++) {
              if (std::abs(tmp(j)) < 1e-30) {
                tmp(j) = 0;
              }
            }

            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(0) << " " << tmp(1) << " " << tmp(2) << "\n";
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(3) << " " << tmp(4) << " " << tmp(5) << "\n";
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << tmp(6) << " " << tmp(7) << " " << tmp(8) << "\n";

            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_strain)(i + 0) << " " << (*element_strain)(i + 1) << " " << (*element_strain)(i + 2) << "\n";
            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_strain)(i + 3) << " " << (*element_strain)(i + 4) << " " << (*element_strain)(i + 5) << "\n";
            // myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*element_strain)(i + 6) << " " << (*element_strain)(i + 7) << " " << (*element_strain)(i + 8) << "\n";

            i = i + 9;
        }
        myfile << "\t\t\t\t</DataArray>\n";
    }

    if (print_fiber_damage == 1) {
        VectorXd* damage = mesh1.getCellDamagePointer();
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Fiber Damage\" NumberOfComponents=\"1\" ComponentName0=\"D\" format=\"ascii\">\n";

        for (int i = 0; i < (*damage).size(); i++)
        {
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*damage)(i) << "\n";
        }
        myfile << "\t\t\t\t</DataArray>\n";
        damage = NULL;

        VectorXd* stretch_min = mesh1.getCellStretchMinPointer();
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Minimum Stretch\" NumberOfComponents=\"1\" ComponentName0=\"Minimum Stretch\" format=\"ascii\">\n";

        for (int i = 0; i < (*stretch_min).size(); i++)
        {
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*stretch_min)(i) << "\n";
        }
        myfile << "\t\t\t\t</DataArray>\n";
        damage = NULL;

        VectorXd* stretch_max = mesh1.getCellStretchMaxPointer();
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Maximum Stretch\" NumberOfComponents=\"1\" ComponentName0=\"Maximum Stretch\" format=\"ascii\">\n";

        for (int i = 0; i < (*stretch_max).size(); i++)
        {
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*stretch_max)(i) << "\n";
        }
        myfile << "\t\t\t\t</DataArray>\n";
        damage = NULL;
    }

    if (print_host_damage == 1) {
        VectorXd* damage = mesh1.getCellDamagePointer();
        myfile << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Average Fiber Damage\" NumberOfComponents=\"1\" ComponentName0=\"D\" format=\"ascii\">\n";

        for (int i = 0; i < (*damage).size(); i++)
        {
            myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << (*damage)(i) << "\n";
        }
        myfile << "\t\t\t\t</DataArray>\n";
        damage = NULL;
    }

    myfile << "\t\t\t</CellData>\n";
    myfile << "\t\t</Piece>\n";
    myfile << "\t</UnstructuredGrid>\n";
    myfile << "</VTKFile>\n";

    if (time_step == 0) {
        fe_pvdNew(output, time_step, time);
    }
    else {
        fe_pvdAppend(output, time_step, time);
    }

    nodes = NULL;
    elements = NULL;
    U_mesh = NULL;
    V_mesh = NULL;
    A_mesh = NULL;
    element_stress = NULL;
    element_strain = NULL;
}

void fe_pvdNew(std::string output, int time_step, double time)
{

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

void fe_pvdAppend(std::string output, int time_step, double time)
{

    /** Output File Name */
    std::string name;
    name = home_path + "/" + "results/vtu/" + output + ".pvd";
    std::ofstream myfile(name.c_str(), std::ios_base::app);

    /** Append New Line of Data */
    if (myfile.is_open())
    {
        myfile << "     <DataSet timestep= \"" << std::scientific << std::setprecision(5) << time << "\" group=\"\" part=\"0\"\n";
        std::ostringstream ss;
        ss << time_step;
        std::string vtuFileName;
        //vtuFileName = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
        vtuFileName = output + "_" + ss.str() + ".vtu";
        myfile << "\t      file=\"" << vtuFileName << "\"/>\n";

        /** Write Footer */
        if (time == t_end)
        {
            myfile << "   </Collection>\n";
            myfile << "</VTKFile>\n";
        }
    }

    else
    {
        std::cout << "No Such Text File With This Name Exists - " << name << "\n";
        std::exit(-1);
    }
}
