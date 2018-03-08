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

void Mesh::readDamage(VectorXd& d_tot)
{
    damage = d_tot;
}

void Mesh::readStretchMin(VectorXd& lambda_min)
{
    stretch_min = lambda_min;
}

void Mesh::readStretchMax(VectorXd& lambda_max)
{
    stretch_max = lambda_max;
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

VectorXd Mesh::getCellDamage()
{
    return damage;
}

VectorXd* Mesh::getCellDamagePointer() {
    return damage_pointer;
}

VectorXd Mesh::getCellStretchMin()
{
    return stretch_min;
}

VectorXd* Mesh::getCellStretchMinPointer() {
    return stretch_min_pointer;
}

VectorXd Mesh::getCellStretchMax()
{
    return stretch_max;
}

VectorXd* Mesh::getCellStretchMaxPointer() {
    return stretch_max_pointer;
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

void Mesh::preprocessMesh(std::string choice)
{

    nodes_new = nodes;
    elements_new = elements;

    if (choice == "re-number" || choice == "default" || choice.length() == 0) {
        int tmp = nodes.col(0).minCoeff();
        if ( tmp != 0) {
            for (int i = 0; i < nodes_new.rows(); i++) {
                nodes_new(i, 0) = nodes_new(i, 0) - tmp;
            }
            for (int i = 0; i < elements_new.rows(); i++) {
                for (int j = 0; j < (elements_new.cols() - 2); j++) {
                    elements_new(i, j + 2) = elements_new(i, j + 2) - tmp;
                }
            }
        }
    }

    if (choice == "re-order") {
        /** Nodes Preprocessing - Putting the numbering in order */
        for (int i = 0; i < nodes.rows(); i++) {
            nodes_new(i, 0) = i;
        }
        /** Elements Preprocessing - Correcting the element definitions */
        for (int i = 0; i < elements.rows(); i++) {
            elements_new(i, 0) = i;
            elements_new(i, 1) = elements(i, 1);
            for (int j = 2; j < elements.cols(); j++) {
                VectorXd node_list = nodes.col(0);
                elements_new(i, j) = fe_find(node_list, elements(i, j));
            }
        }
        for (int i = 0; i < elements_new.rows(); i++) {
            elements_new(i, 0) = i + 1;
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
    damage_pointer = &damage;
    stretch_min_pointer = &stretch_min;
    stretch_max_pointer = &stretch_max;
    element_charateristic_pointer = &element_charateristic;
    nodal_electric_potential_pointer = &nodal_electric_potential;

    // fe_dniso_8
    if (counter_test == 0) {

        std::cout << "\n Hurray\n";
        counter_test += 1;

        int nglx = 2, ngly = 2, nglz = 2, nnel;
        nnel = getNumNodesPerElement();
        nnel_test = nnel;
        edof_test = nnel * 3;
        std::cout << "\n " << nnel_test << "\n";
        std::cout << "\n " << edof_test << "\n";
        VectorXd points  = guass_points(nglx);
        VectorXd weights = guass_weights(nglx);

        VectorXd dndr(nnel);
        VectorXd dnds(nnel);
        VectorXd dndt(nnel);

        for (int intx = 0; intx < nglx; intx++) {
            double x   = points(intx);
            for (int inty = 0; inty < ngly; inty++) {
                double y   = points(inty);
                for (int intz = 0; intz < nglz; intz++) {
                    double z   = points(intz);
                    fe_dniso_8(dndr, dnds, dndt, x, y, z);
                    Map<VectorXd>(dndr_store[intx][inty][intz], dndr.rows()) = dndr;
                    Map<VectorXd>(dnds_store[intx][inty][intz], dnds.rows()) = dnds;
                    Map<VectorXd>(dndt_store[intx][inty][intz], dndt.rows()) = dndt;
                }
            }
        }
    }
}

void Mesh::append2nodes(std::string choice, VectorXd& b) {

    if (choice == "new") {
        MatrixXd tmp = MatrixXd::Zero(nodes_new.rows(), nodes_new.cols());
        tmp = nodes_new;
        nodes_new = MatrixXd::Zero(tmp.rows() + 1, tmp.cols());
        for (int i = 0; i < nodes_new.rows(); i++) {
            if (i != nodes_new.rows() - 1) {
                nodes_new.row(i) = tmp.row(i);
            }
            else
            {
                nodes_new.row(i) = b;
            }

        }
    }

    if (choice == "old") {
        MatrixXd tmp = MatrixXd::Zero(nodes.rows(), nodes.cols());
        tmp = nodes;
        nodes = MatrixXd::Zero(tmp.rows() + 1, tmp.cols());
        for (int i = 0; i < nodes.rows(); i++) {
            if (i != nodes.rows() - 1) {
                nodes_new.row(i) = tmp.row(i);
            }
            else
            {
                nodes_new.row(i) = b;
            }

        }
    }

}

void Mesh::append2elements(std::string choice, VectorXi& b) {

    if (choice == "new") {
        MatrixXi tmp = MatrixXi::Zero(elements_new.rows(), elements_new.cols());
        tmp = elements_new;
        elements_new = MatrixXi::Zero(tmp.rows() + 1, tmp.cols());
        for (int i = 0; i < elements_new.rows(); i++) {
            if (i != elements_new.rows() - 1) {
                elements_new.row(i) = tmp.row(i);
            }
            else
            {
                elements_new.row(i) = b;
            }

        }
    }

    if (choice == "old") {
        MatrixXi tmp = MatrixXi::Zero(elements.rows(), elements.cols());
        tmp = elements;
        elements = MatrixXi::Zero(tmp.rows() + 1, tmp.cols());
        for (int i = 0; i < elements.rows(); i++) {
            if (i != elements.rows() - 1) {
                elements.row(i) = tmp.row(i);
            }
            else
            {
                elements.row(i) = b;
            }

        }
    }

}

void Mesh::reNumber() {

    MatrixXd tmp1 = nodes_new;
    MatrixXi tmp2 = elements_new;

    /** Nodes Preprocessing - Putting the numbering in order */
    for (int i = 0; i < nodes_new.rows(); i++) {
        nodes_new(i, 0) = i;
    }
    /** Elements Preprocessing - Correcting the element definitions */
    for (int i = 0; i < elements_new.rows(); i++) {
        elements_new(i, 0) = i + 1;
        for (int j = 2; j < elements_new.cols(); j++) {
            VectorXd node_list = tmp1.col(0);
            elements_new(i, j) = fe_find(node_list, tmp2(i, j));
        }
    }
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

    Vector2d min_details;

    MatrixXd* nodes_local;
    MatrixXi* elements_local;

    if (choice == "old")
    {
        nodes_local = &nodes;
        elements_local = &elements;
    }
    else
    {
        nodes_local = &nodes_new;
        elements_local = &elements_new;
    }

    VectorXd xcoord;
    VectorXd ycoord;
    VectorXd zcoord;
    double min_length = 10000000;
    double id = 0;
    double lc;

    xcoord = VectorXd::Zero((*elements_local).cols() - 2);
    ycoord = VectorXd::Zero((*elements_local).cols() - 2);
    zcoord = VectorXd::Zero((*elements_local).cols() - 2);
    for (int i = 0; i < (*elements_local).rows(); i++)
    {
        for (int j = 0; j < ((*elements_local).cols() - 2); j++)
        {
            int g = (*elements_local)(i, j + 2);
            xcoord[j] = (*nodes_local)(g, 1);
            ycoord[j] = (*nodes_local)(g, 2);
            zcoord[j] = (*nodes_local)(g, 3);
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

    Vector2d max_details;
    double id;

    MatrixXd* nodes_local;
    MatrixXi* elements_local;

    if (choice == "old")
    {
        nodes_local = &nodes;
        elements_local = &elements;
    }
    else
    {
        nodes_local = &nodes_new;
        elements_local = &elements_new;
    }

    VectorXd xcoord;
    VectorXd ycoord;
    VectorXd zcoord;
    double max_length = 0.0;
    double lc;

    xcoord = VectorXd::Zero((*elements_local).cols() - 2);
    ycoord = VectorXd::Zero((*elements_local).cols() - 2);
    zcoord = VectorXd::Zero((*elements_local).cols() - 2);
    for (int i = 0; i < (*elements_local).rows(); i++)
    {
        for (int j = 0; j < ((*elements_local).cols() - 2); j++)
        {
            int g = -1;
            for (int f = 0; f < (*nodes_local).rows(); f++)
            {
                if ((*elements_local)(i, j + 2) == (*nodes_local)(f, 0))
                {
                    g = f;
                    break;
                }
            }
            xcoord[j] = (*nodes_local)(g, 1);
            ycoord[j] = (*nodes_local)(g, 2);
            zcoord[j] = (*nodes_local)(g, 3);
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
    Vector2d lc;

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
