#include "functions.h"

using namespace Eigen;

/** This function calculates the updated stress for 1d elements - hyperelastic material model was implemented so far. */
VectorXd
fe_stressUpdate_1d(int material_id, VectorXd u_embed_local, VectorXd xcoord_embed, VectorXd ycoord_embed, VectorXd zcoord_embed, double length_old, int choice)
{
    VectorXd xcoord_embed_new(xcoord_embed.size()), ycoord_embed_new(ycoord_embed.size()), zcoord_embed_new(
        zcoord_embed.size());
    int counter = 0;

    for (int i = 0; i < xcoord_embed.size(); i++) {
        xcoord_embed_new(i) = xcoord_embed(i) + u_embed_local(counter);
        ycoord_embed_new(i) = xcoord_embed(i) + u_embed_local(counter + 1);
        zcoord_embed_new(i) = zcoord_embed(i) + u_embed_local(counter + 2);
        counter = counter + 3;
    }
    double length_new = fe_calVolume(xcoord_embed_new, ycoord_embed_new, zcoord_embed_new);
    double lambda     = length_new / length_old;

    VectorXd sigma = VectorXd::Zero(6); // piola-kirchhoff stress-vector
    double cauchy_sigma;
    VectorXd pk_S = VectorXd::Zero(6);

    std::string model;
    model = fe_get_model(material_id, "mechanical");

    if (model == "simple_elastic") { // Isotropic truss element
        double E = fe_get_mats(material_id, 1, "mechanical");
        //std::cout << "Elastic Modulus - Fiber: " << E << "\n";
        cauchy_sigma = E * log(lambda);
        //std::cout << "Fiber Stretch: " << lambda << "\n";
        //std::cout << "Engineering Strain in the fiber: " << lambda - 1 << "\n";
        //std::cout << "Cauchy Stress - Fiber: " << cauchy_sigma << "\n";
    }

    if (model == "mooney-rivlin_hyperelastic") {
        double c_1 = fe_get_mats(material_id, 3, "mechanical");
        double c_2 = fe_get_mats(material_id, 4, "mechanical");
        cauchy_sigma = ((2 * c_1) + ((2 * c_2) / (lambda))) * ((pow(lambda, 2)) - (1 / lambda));
    }

    double second_pk = cauchy_sigma / lambda;
    // double second_pk = cauchy_sigma;
    pk_S(0) = second_pk;

    sigma(0) = pk_S(0);

    if (choice == 1) {
        sigma(0) = cauchy_sigma;
    }

    return sigma;
} // fe_stressUpdate_1d

/***************************************************/

VectorXd
fe_stressUpdate_1d(int material_id, VectorXd u_embed_local, VectorXd xcoord_embed, VectorXd ycoord_embed, VectorXd zcoord_embed, double length_old, VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u_e, int choice)
{
    VectorXd xcoord_embed_new(xcoord_embed.size()), ycoord_embed_new(ycoord_embed.size()), zcoord_embed_new(zcoord_embed.size());
    int counter = 0;

    for (int i = 0; i < xcoord_embed.size(); i++) {
        xcoord_embed_new(i) = xcoord_embed(i) + u_embed_local(counter);
        ycoord_embed_new(i) = xcoord_embed(i) + u_embed_local(counter + 1);
        zcoord_embed_new(i) = zcoord_embed(i) + u_embed_local(counter + 2);
        counter = counter + 3;
    }

    double length_new = fe_calVolume(xcoord_embed_new, ycoord_embed_new, zcoord_embed_new);
    double lambda     = length_new / length_old;

    VectorXd sigma        = VectorXd::Zero(6); // piola-kirchhoff stress-vector
    VectorXd cauchy_sigma = VectorXd::Zero(6);
    // VectorXd pk_S         = VectorXd::Zero(6);

    std::string model;
    model = fe_get_model(material_id, "mechanical");

    if (model == "simple_elastic") { // Isotropic truss element
        double E = fe_get_mats(material_id, 1, "mechanical");
        cauchy_sigma(0) = E * log(lambda);
    }

    if (model == "mooney-rivlin_hyperelastic") {
        double c_1 = fe_get_mats(material_id, 3, "mechanical");
        double c_2 = fe_get_mats(material_id, 4, "mechanical");
        cauchy_sigma(0) = ((2 * c_1) + ((2 * c_2) / (lambda))) * ((pow(lambda, 2)) - (1 / lambda));
    }

    MatrixXd stress_transformation_mat = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 1);

    VectorXd tmp_cauchy_sigma = stress_transformation_mat.transpose() * cauchy_sigma; // Cauchy stress vector in element CSYS.

    MatrixXd element_cauchy_sigma = MatrixXd::Zero(3, 3);
    element_cauchy_sigma(0, 0) = tmp_cauchy_sigma(0);
    element_cauchy_sigma(1, 1) = tmp_cauchy_sigma(1);
    element_cauchy_sigma(2, 2) = tmp_cauchy_sigma(2);
    element_cauchy_sigma(0, 1) = tmp_cauchy_sigma(3);
    element_cauchy_sigma(1, 0) = tmp_cauchy_sigma(3);
    element_cauchy_sigma(1, 2) = tmp_cauchy_sigma(4);
    element_cauchy_sigma(2, 1) = tmp_cauchy_sigma(4);
    element_cauchy_sigma(0, 2) = tmp_cauchy_sigma(5);
    element_cauchy_sigma(2, 0) = tmp_cauchy_sigma(5);

    // fe_display_matrix(element_cauchy_sigma);

    MatrixXd F = MatrixXd::Zero(3, 3); // deformation gradient
    fe_calDefGrad_pbr(F, dndx, dndy, dndz, u_e);
    double defJacobian = F.determinant();
    MatrixXd F_inv     = F.inverse();
    MatrixXd pk_S      = MatrixXd::Zero(3, 3);
    pk_S = defJacobian * F_inv * element_cauchy_sigma * F_inv.transpose();

    sigma(0) = pk_S(0, 0);
    sigma(1) = pk_S(1, 1);
    sigma(2) = pk_S(2, 2);
    sigma(3) = pk_S(0, 1);
    sigma(4) = pk_S(1, 2);
    sigma(5) = pk_S(0, 2);

    /*sigma(0) = element_cauchy_sigma(0,0);
     *  sigma(1) = element_cauchy_sigma(1,1);
     *  sigma(2) = element_cauchy_sigma(2,2);
     *  sigma(3) = element_cauchy_sigma(0,1);
     *  sigma(4) = element_cauchy_sigma(1,2);
     *  sigma(5) = element_cauchy_sigma(0,2);*/

    return sigma;
} // fe_stressUpdate_1d
