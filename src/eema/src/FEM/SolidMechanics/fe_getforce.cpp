#include "functions.h"

using namespace Eigen;

/*! \brief Calculates the resultant nodal force after each time step.
 *
 * This function represents the 'getforce' step in Belytschko (Box 6.1 - Explicit FEM Algorithm).
 *  For each hex element, this function calculates the internal nodal force vector
 * and the resultant nodal force vector. Once, this is calculated for each element,
 * the resultant vectors are scattered into global vectors.
 */

void fe_getforce(VectorXd& f_tot, int ndof, VectorXd& u, VectorXd& fext, int time_step_counter)
{

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
                f_tot = VectorXd::Zero((mesh[host_id].getNumNodes() * ndof));
                fe_getForce_3d_embed(f_tot, u, fext, time_step_counter, host_id, embed_id, correct_vr);
            }
        }
    }

    /*else if(ndof==2 && embedded_constraint==true)
     * {
     *
     * }
     * else if(ndof==2 && embedded_constraint!=true)
     * {
     *
     * }*/
    else if (ndof == 3 && embedded_constraint != true) {
        f_tot = VectorXd::Zero((mesh[0].getNumNodes() * ndof));
        fe_getForce_3d_normal(f_tot, u, fext, time_step_counter, 0);
    } else {
        f_tot = VectorXd::Zero((mesh[0].getNumNodes() * ndof));
        fe_getForce_3d_normal(f_tot, u, fext, time_step_counter, 0);
    }


} // fe_getforce
