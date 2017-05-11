#include "functions.h"
#include "Materials.h"

using namespace Eigen;

void Materials::readMatId(int a) {
	mat_id = a;
}

void Materials::readMats(std::string b, VectorXd c, std::string type) {
	if (type == "mechanical") {
		mechanical_mat_model = b;
		mechanical_mat_properties = VectorXd::Zero(c.size());
		mechanical_mat_properties = c;
	}
	if (type == "electrical") {
		electrical_mat_model = b;
		electrical_mat_properties = VectorXd::Zero(c.size());
		electrical_mat_properties = c;
	}
}

void Materials::printInfo(std::string type) {
	std::cout << "Material ID is: " << mat_id << "\n";

	if (type == "mechanical") {
		std::cout << "Material Model (Mechanical): " << mechanical_mat_model << "\n";
		std::cout << "Material Properties (Mechanical): \n" << mechanical_mat_properties << "\n";
	}

	if (type == "electrical") {
		std::cout << "Material Model (Electrical): " << electrical_mat_model << "\n";
		std::cout << "Material Properties (Electrical): \n" << electrical_mat_properties << "\n";
	}

	if (type == "all") {
		if (mechanical_mat_properties.size() != 0) {
			std::cout << "Mechanical: " << "\n";
			std::cout << "Material Model (Mechanical): " << mechanical_mat_model << "\n";
			std::cout << "Material Properties (Mechanical): \n" << mechanical_mat_properties << "\n";
		}
		std::cout << "----------------------------------" << "\n";
		if (electrical_mat_properties.size() != 0) {
			std::cout << "Electrical: " << "\n";
			std::cout << "Material Model (Electrical): " << electrical_mat_model << "\n";
			std::cout << "Material Properties (Electrical): \n" << electrical_mat_properties << "\n";
		}
		std::cout << "----------------------------------" << "\n";
	}
}
