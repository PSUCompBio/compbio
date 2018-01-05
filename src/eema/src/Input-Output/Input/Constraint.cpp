#include "functions.h"
#include "Constraint.h"

using namespace Eigen;

void Constraint::readConstraints(std::string name, int id, std::string master, std::string slave, bool address_VR) {
	constraint_name = name;
	constraint_id = id;
	master_name = master;
	slave_name = slave;
	address_volume_redundancy = address_VR;
}

void Constraint::printInfo() {
	std::cout << "Constraint Type: " << constraint_name << "\n";
	std::cout << "Constraint ID: " << constraint_id << "\n";
	std::cout << "According to this constraint definition" << "\n" ;
	std::cout << "\t Master Mesh - " << master_name << "\n";
	std::cout << "\t Slave Mesh - " << slave_name << "\n";
}

std::string Constraint::getName() {
	return constraint_name;
}

std::string Constraint::get_EmbedMaster() {
	return master_name;
}

std::string Constraint::get_EmbedSlave() {
	return slave_name;
}

bool Constraint::get_EmbedAddressVR() {
	return address_volume_redundancy;
}

void Constraint::preprocess() {

	if (constraint_name == "embedded") {
		int host_id, embed_id;
		for (int i = 0; i < num_meshes; ++i) {
			std::string name = mesh[i].getName();
			if (name == master_name) {
				host_id = i;
			}
			if (name == slave_name) {
				embed_id = i;
			}
		}

		embed_map = fe_embed_preprocessing(mesh[host_id], mesh[embed_id]);

		std::cout << "embed_map = " << '\n' << embed_map << '\n';

		embed_map_pointer = &embed_map;
	}

}

VectorXi* Constraint::get_EmbedMapPointer() {
	return embed_map_pointer;
}
