#include "functions.h"
using namespace Eigen;

void fe_electroStatics(double time) {
	if (ndof == 3 && embedded_constraint == true) {
		for (int i = 0; i < num_constraints; ++i) {
			if (cons[i].getName() == "embedded") {
				std::string host = cons[i].get_EmbedMaster();
				std::string slave = cons[i].get_EmbedSlave();
				bool correct_vr = cons[i].get_EmbedAddressVR();
				int host_id, embed_id;

				for (int i = 0; i < num_meshes; ++i) {
					std::string name = mesh[i].getName();

					if (name == host) {
						host_id = i;
					}
					if (name == slave) {
						embed_id = i;
					}
				}

				VectorXi* embed_map = cons[i].get_EmbedMapPointer();

				fe_electroStatics_embed(time, host_id, embed_id, (*embed_map));
			}
		}
	}
	else if (ndof == 3 && embedded_constraint != true) {
		fe_electroStatics_normal(time, 0);
	} else {
		fe_electroStatics_normal(time, 0);
	}
}