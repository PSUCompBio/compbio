#include "functions.h"

using namespace Eigen;

/*! \brief Calculates the resultant nodal force after each time step.
 *
 * This function represents the 'getforce' step in Belytschko (Box 6.1 - Explicit FEM Algorithm).
 *  For each hex element, this function calculates the internal nodal force vector
 * and the resultant nodal force vector. Once, this is calculated for each element,
 * the resultant vectors are scattered into global vectors.
 */

void fe_getforce(VectorXd& f_tot, int ndof, VectorXd& u, VectorXd& fext, int time_step_counter, VectorXd& u_prev, double dT, VectorXd& f_damp, VectorXd& d_static, VectorXd& d_fatigue, VectorXd& d_tot, VectorXd& lambda_min, VectorXd& lambda_max, VectorXd& lambda_min_cycle, VectorXd& lambda_max_cycle, VectorXd& d_avg, VectorXi& n_load_cycle_full, VectorXi& n_load_cycle_partial, double t, int t_plot)
{

    if (ndof == 3 && embedded_constraint == true) {
        for (int i = 0; i < num_constraints; ++i) {
            if (cons[i].getName() == "embedded") {
                std::string host = cons[i].get_EmbedMaster();
                std::string slave = cons[i].get_EmbedSlave();
                bool correct_vr = cons[i].get_EmbedAddressVR();
                bool include_d = cons[i].get_EmbedIncludeDamage();
                int host_id, embed_id;

                for (int j = 0; j < num_meshes; ++j) {
                    std::string name = mesh[j].getName();

                    if (name == host) {
                        host_id = j;
                    }
                    if (name == slave) {
                        embed_id = j;
                    }
                }

                VectorXi* embed_map = cons[i].get_EmbedMapPointer();

                f_tot = VectorXd::Zero((mesh[host_id].getNumNodes() * ndof));
                fe_getForce_3d_embed(f_tot, u, fext, time_step_counter, host_id, embed_id, correct_vr, include_d, (*embed_map), u_prev, dT, f_damp, d_static, d_fatigue, d_tot, lambda_min, lambda_max, lambda_min_cycle, lambda_max_cycle, d_avg, n_load_cycle_full, n_load_cycle_partial, t, t_plot);
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
        fe_getForce_3d_normal(f_tot, u, fext, time_step_counter, 0, u_prev, dT, f_damp, t, t_plot);
    } else {
        f_tot = VectorXd::Zero((mesh[0].getNumNodes() * ndof));
        fe_getForce_3d_normal(f_tot, u, fext, time_step_counter, 0, u_prev, dT, f_damp, t, t_plot);
    }


} // fe_getforce
