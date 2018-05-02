#include "functions.h"

using namespace Eigen;

/** This function carries out all the required pre-processing steps for embedded element constraint */

VectorXi fe_embed_preprocessing(Mesh& host, Mesh& embed) {

    // Retrieve number of embedded elements.
    int num_elements_embed = embed.getNumElements();

    // Define array to store host id for each embedded element.
    VectorXi embed_map = VectorXi::Zero(num_elements_embed);

    // Run sub-routine to assign host ids to embedded elements.
    // During process, embedded elements are split at host interface and new embedded elements are generated.
    // As a result, some new embedded elements are not assigned host ids on first pass of sub-routine.
    // Therefore, sub-routine may need to run multiple times.
    embed_map = fe_embed_preprocessing_host_map(embed_map, host, embed);

    int embed_map_size = embed_map.size();
    int host_id_check = 0;

    // Check if any embedded elements have not been assigned a host id.
    for (int i = 0; i < embed_map_size; i++) {
      if (embed_map(i) == 0) {
        host_id_check = 1;
        break;
      }
    }

    int counter = 0;

    // Run sub-routine until all embedded elements have been assigned host ids.
    while (host_id_check == 1) {

      std::cout << "counter = " << counter << '\n';

      embed_map = fe_embed_preprocessing_host_map(embed_map, host, embed);
      embed_map_size = embed_map.size();

      host_id_check = 0;
      for (int i = 0; i < embed_map_size; i++) {
        if (embed_map(i) == 0) {
          host_id_check = 1;

          break;
        }
      }

      counter = counter + 1;

    }

    return embed_map;

}
