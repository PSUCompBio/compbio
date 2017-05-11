#include<iomanip>
#include "functions.h"
using namespace Eigen;

/* Home Folder and Job File */
std::string home_path; /** Job folder name */
std::string job_file; /** Job input file name */

/* Dimension of the Problem */
int ndof;

/* Input meshes - number of meshes and type of meshes */
int num_meshes;
Mesh *mesh;

/* Material Properties are stored in here */
int material_types;
Materials *mat;

/* Time Related Initialized Variables Stored Here */
double t_start; /** Start Time */
double t_end; /** End Time */
int output_frequency;
double dt_initial; /** Initial time step */
double dt_min; /** Minimum time step */
double reduction; /** Reduces the timestep value by this amount */

/* Boundary Conditions Stored Here */
int bc_types;
BC *bc;

int num_constraints;
Constraint *cons;
bool embedded_constraint;

void fe_mainRead(std::string file) {

  num_meshes = 0;
  int num_meshes_counter = 0;
  material_types = 0;
  int material_types_counter = 0;
  bc_types = 0;
  int bc_types_counter = 0;
  num_constraints = 0;
  int num_constraints_counter = 0;

  std::ifstream myfile(file.c_str());
  std::string line;
  std::getline(myfile, line);

  while (line != "*END") {
    if (line[0] == '*') {
      if (line == "*MESH") {
        num_meshes = num_meshes + 1;
      }
      if (line == "*MATERIAL") {
        material_types = material_types + 1;
      }
      if (line == "*BC") {
        bc_types = bc_types + 1;
      }
      if (line == "*CONSTRAINT") {
        num_constraints = num_constraints + 1;
      }
    }
    std::getline(myfile, line);
  }

  mesh = new Mesh[num_meshes];
  mat = new Materials[material_types];
  bc = new BC[bc_types];
  cons = new Constraint[num_constraints];

  std::ifstream myfile1(file.c_str());
  std::getline(myfile1, line);

  while (line != "*END") {

    if (line[0] == '*') {

      if (line == "*DIMENSION") {
        while (line != "*END_DIMENSION") {
          myfile1 >> ndof;
          myfile1 >> line;
        }
      }

      if (line == "*MESH") {
        std::getline(myfile1, line);
        MatrixXd nodes;
        MatrixXi elements;
        std::string name;

        while (line != "*END_MESH") {

          if (line == "*NAME") {
            myfile1 >> name;
          }

          if (line == "*NODES") {
            int rows = 0;
            int cols = 0;
            myfile1 >> rows;
            myfile1 >> cols;
            nodes = MatrixXd::Zero(rows, cols + 1);
            for (int i = 0; i < rows; i++) {
              for (int j = 0; j < (cols + 1); j++) {
                myfile1 >> nodes(i, j);
              }
            }
          }

          if (line == "*ELEMENTS") {
            int rows = 0;
            int cols = 0;
            myfile1 >> rows;
            myfile1 >> cols;
            elements = MatrixXi::Zero(rows, cols + 2);
            for (int i = 0; i < rows; i++) {
              for (int j = 0; j < (cols + 2); j++) {
                myfile1 >> elements(i, j);
              }
            }
          }

          if (nodes.rows() != 0 && elements.rows() != 0 && name.length() != 0) {
            mesh[num_meshes_counter].readMesh(name, nodes, elements);
            num_meshes_counter = num_meshes_counter + 1;
          }

          myfile1 >> line;
        }
      }

      if (line == "*MATERIAL") {
        int mat_id;
        std::string mat_model;
        VectorXd mat_properties = VectorXd::Zero(100);
        myfile1 >> mat_id;
        myfile1 >> mat_model;
        int i = 0;
        std::getline(myfile1, line);
        while (line != "*END_MATERIAL") {
          myfile1.seekg(-(line.length()), std::ios::cur);
          myfile1 >> mat_properties(i);
          myfile1 >> line;
          i = i + 1;
        }
        mat[material_types_counter].readMats(mat_id, mat_model, mat_properties);
        material_types_counter = material_types_counter + 1;
      }

      if (line == "*TIME") {
        while (line != "*END_TIME") {
          myfile1 >> t_start;
          myfile1 >> t_end;
          myfile1 >> output_frequency;
          myfile1 >> dt_min;
          dt_initial = dt_min;
          myfile1 >> reduction;
          myfile1 >> line;
        }
      }

      if (line == "*BC") {
        std::string type;
        double amplitude;
        int num_dof;
        VectorXi dof;
        std::string time_curve;
        while (line != "*END_BC") {
          myfile1 >> type;
          myfile1 >> amplitude;
          myfile1 >> num_dof;
          dof = VectorXi::Zero(num_dof);
          for (int i = 0; i < num_dof; i++) {
            int tmp;
            myfile1 >> tmp;
            dof(i) = tmp;
          }
          myfile1 >> time_curve;
          myfile1 >> line;
        }
        bc[bc_types_counter].readBC(type, amplitude, num_dof, dof, time_curve);
        bc_types_counter = bc_types_counter + 1;
      }

      if (line == "*CONSTRAINT")
      {
        std::cout << "Reading Constraint Information" << "\n";
        std::string type;
        int id;
        std::string master;
        std::string slave;
        bool address_vr;

        while (line != "*END_CONSTRAINT")
        {
          myfile1 >> id;
          myfile1 >> type;
          std::cout << "Type: " << type << "\n";
          if (type == "embedded")
          {
            embedded_constraint = true;

            myfile1 >> master;

            myfile1 >> slave;

            myfile1 >> address_vr;

          }
          myfile1 >> line;
        }
        if (num_constraints_counter < num_constraints) {
          cons[num_constraints_counter].readConstraints(type, id, master, slave, address_vr);
          num_constraints_counter = num_constraints_counter + 1;
        }
      }

    }
    std::getline(myfile1, line);
  }

  std::cout << "EEMA: Reading Input File --> Completed !!" << "\n";

  /*std::cout << "*************************" << '\n';
  std::cout << "Mesh Details: \n";
  for (int i = 0; i < num_meshes; i++) {
    mesh[i].printInfo();
  }
  std::cout << "*************************" << '\n';
  std::cout << "Material Details: " << '\n';
  for (int i = 0; i < material_types; i++) {
    mat[i].printInfo();
  }
  std::cout << "*************************" << '\n';
  std::cout << "Boundary Condition Details: " << '\n';
  for (int i = 0; i < bc_types; i++) {
    bc[i].printInfo();
  }
  std::cout << "*************************" << '\n';
  std::cout << "Constraint Details - Num Constraints: " << num_constraints << '\n';
  for (int i = 0; i < num_constraints; i++) {
    cons[i].printInfo();
  }
  std::cout << "*************************" << '\n';*/

  if (num_meshes == 0) {
    std::cout << "No meshes included - Simulation is not possible !! " << "\n";
    std::exit(-1);
  }

  for (int i = 0; i < num_meshes; i++) {
    mesh[i].preprocessMesh();
  }

  std::cout << "EEMA: Mesh Preprocessing --> Completed !!" << "\n";

}
