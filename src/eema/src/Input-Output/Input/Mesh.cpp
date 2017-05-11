#include "functions.h"
#include "Mesh.h"

using namespace Eigen;

void Mesh::readMesh(std::string name, MatrixXd& n, MatrixXi& e)
{
    mesh_name = name;
    nodes = n;
    elements = e;
}

MatrixXd Mesh::getNodes(void)
{
    return nodes;
}

MatrixXi Mesh::getElements(void)
{
    return elements;
}

MatrixXd Mesh::getNewNodes(void)
{
    return nodes_new;
}

MatrixXi Mesh::getNewElements(void)
{
    return elements_new;
}

MatrixXd* Mesh::getNewNodesPointer(void) {
    return nodes_new_pointer;
}

MatrixXi* Mesh::getNewElementsPointer(void) {
    return elements_new_pointer;
}

int Mesh::getNumElements()
{
    return elements_new.rows();
}

int Mesh::getNumNodes()
{
    return nodes_new.rows();
}

int Mesh::getNumNodesPerElement() {
    return (elements_new.cols() - 2) ;
}

std::string Mesh::getName(void)
{
    return mesh_name;
}

void Mesh::readNodalKinematics(VectorXd& Usystem, VectorXd& Vsystem, VectorXd& Asystem)
{
    U = Usystem;
    V = Vsystem;
    A = Asystem;
}

VectorXd Mesh::getNodalDisp(void)
{
    return U;
}

VectorXd* Mesh::getNodalDispPointer(void) {
    return U_pointer;
}

VectorXd Mesh::getNodalVel(void)
{
    return V;
}

VectorXd* Mesh::getNodalVelPointer(void) {
    return V_pointer;
}

VectorXd Mesh::getNodalAcc(void)
{
    return A;
}

VectorXd* Mesh::getNodalAccPointer(void) {
    return A_pointer;
}

void Mesh::readNodalStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp)
{
    nodal_stress = stress_tmp;
    nodal_strain = strain_tmp;
}

void Mesh::readElementStressStrain(VectorXd& stress_tmp, VectorXd& strain_tmp)
{
    element_stress = stress_tmp;
    element_strain = strain_tmp;
}

VectorXd Mesh::getCellStress()
{
    return element_stress;
}

VectorXd* Mesh::getCellStressPointer() {
    return element_stress_pointer;
}

VectorXd Mesh::getCellStrain()
{
    return element_strain;
}

VectorXd* Mesh::getCellStrainPointer() {
    return element_strain_pointer;
}

VectorXd Mesh::getNodalStress()
{
    return nodal_stress;
}

VectorXd Mesh::getNodalStrain()
{
    return nodal_strain;
}

void Mesh::readNodalElectroPhysics(VectorXd& VP) {
    nodal_electric_potential = VP;
}

VectorXd Mesh::getNodalEPotential() {
    return nodal_electric_potential;
}

VectorXd* Mesh::getNodalEPotentialPointer() {
    return nodal_electric_potential_pointer;
}

void Mesh::preprocessMesh(void)
{

    nodes_new = nodes;
    elements_new = elements;

    int tmp = nodes.col(0).minCoeff();

    if ( tmp != 0) {
        for (int i = 0; i < nodes.rows(); i++) {
            nodes_new(i, 0) = nodes(i, 0) - tmp;
        }
        for (int i = 0; i < elements.rows(); i++) {
            for (int j = 0; j < (elements.cols() - 2); j++) {
                elements_new(i, j + 2) = elements(i, j + 2) - tmp;
            }
        }
    }

    calculateElementCharateristic();

    nodes_new_pointer = &nodes_new;
    elements_new_pointer = &elements_new;
    U_pointer = &U;
    V_pointer = &V;
    A_pointer = &A;
    element_stress_pointer = &element_stress;
    element_strain_pointer = &element_strain;
    element_charateristic_pointer = &element_charateristic;
    nodal_electric_potential_pointer = &nodal_electric_potential;
}

void Mesh::calculateElementCharateristic() {

    element_charateristic = VectorXd::Zero(elements_new.rows());

    VectorXd xcoord = VectorXd::Zero(elements_new.cols() - 2);
    VectorXd ycoord = VectorXd::Zero(elements_new.cols() - 2);
    VectorXd zcoord = VectorXd::Zero(elements_new.cols() - 2);
    double volume;

    for (int i = 0; i < elements_new.rows(); i++) {
        for (int j = 0; j < elements_new.cols() - 2; j++) {
            xcoord(j) = nodes_new(elements_new(i, j + 2), 1);
            ycoord(j) = nodes_new(elements_new(i, j + 2), 2);
            zcoord(j) = nodes_new(elements_new(i, j + 2), 3);
        }
        element_charateristic(i) = fe_calVolume(xcoord, ycoord, zcoord);
    }

}

VectorXd* Mesh::getElementCharacteristicPointer() {
    return element_charateristic_pointer;
}

void Mesh::replaceNodes(MatrixXd new_nodes, std::string choice)
{
    if (choice == "old")
    {
        nodes = new_nodes;
    }
    else
    {
        nodes_new = new_nodes;
    }
}

void Mesh::replaceElements(MatrixXi new_elements, std::string choice)
{
    if (choice == "old")
    {
        elements = new_elements;
    }
    else
    {
        elements_new = new_elements;
    }
}

/** \brief Calculates the minimum charateristic length of the mesh */
VectorXd Mesh::getMinCharLength(std::string choice)
{

    VectorXd min_details;

    MatrixXd nodes_local;
    MatrixXi elements_local;

    if (choice == "old")
    {
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

    xcoord = VectorXd::Zero(elements_local.cols() - 2);
    ycoord = VectorXd::Zero(elements_local.cols() - 2);
    zcoord = VectorXd::Zero(elements_local.cols() - 2);
    for (int i = 0; i < elements_local.rows(); i++)
    {
        for (int j = 0; j < (elements_local.cols() - 2); j++)
        {
            int g = -1;
            for (int f = 0; f < nodes_local.rows(); f++)
            {
                if (elements_local(i, j + 2) == nodes_local(f, 0))
                {
                    g = f;
                    break;
                }
            }
            xcoord[j] = nodes_local(g, 1);
            ycoord[j] = nodes_local(g, 2);
            zcoord[j] = nodes_local(g, 3);
        }

        lc = fe_minElementLength(xcoord, ycoord, zcoord);

        if (lc < min_length)
        {
            min_length = lc;
            id = i;
        }
    }

    min_details << min_length, id;

    return min_details;
}

/** \brief Calculates the maximum charateristic length of the mesh */
VectorXd Mesh::getMaxCharLength(std::string choice)
{

    VectorXd max_details;
    double id;

    MatrixXd nodes_local;
    MatrixXi elements_local;

    if (choice == "old")
    {
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

    xcoord = VectorXd::Zero(elements_local.cols() - 2);
    ycoord = VectorXd::Zero(elements_local.cols() - 2);
    zcoord = VectorXd::Zero(elements_local.cols() - 2);
    for (int i = 0; i < elements_local.rows(); i++)
    {
        for (int j = 0; j < (elements_local.cols() - 2); j++)
        {
            int g = -1;
            for (int f = 0; f < nodes_local.rows(); f++)
            {
                if (elements_local(i, j + 2) == nodes_local(f, 0))
                {
                    g = f;
                    break;
                }
            }
            xcoord[j] = nodes_local(g, 1);
            ycoord[j] = nodes_local(g, 2);
            zcoord[j] = nodes_local(g, 3);
        }

        lc = fe_maxElementLength(xcoord, ycoord, zcoord);

        if (lc > max_length)
        {
            max_length = lc;
            id = i;
        }
    }

    max_details << max_length, id;
    return max_details;
}

/** \brief Check the mesh for zero charateristic lengths */
void Mesh::checkMesh()
{
    VectorXd lc;

    lc = getMinCharLength("new");

    if (lc(0) == 0)
    {
        std::cout << "**********************************************" << std::endl;
        std::cout << "ZERO LENGTH ELEMENT IN THE SYSTEM \n FIBER LENGTH PREPROCESSING NOT POSSIBLE" << std::endl;
        std::cout << "**********************************************" << std::endl;
        std::exit(-1);
    }
}

void Mesh::printInfo()
{
    std::cout << "Mesh: " << mesh_name << "\n";
    std::cout << "Number of elements: " << elements.rows() << "\n";
    std::cout << "Number of nodes: " << nodes.rows() << "\n";
}
