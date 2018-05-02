#include "functions.h"

using namespace Eigen;

void fe_checkFiberVolumeFraction(Mesh& host, Mesh& embed) {

    int num_elements_host = host.getNumElements();
    int num_elements_embed = embed.getNumElements();
    MatrixXi* elements_host = host.getNewElementsPointer();
    VectorXd* element_characteristic_host = host.getElementCharacteristicPointer();
    VectorXd* element_characteristic_embed = embed.getElementCharacteristicPointer();
    VectorXi* embed_map = cons[0].get_EmbedMapPointer();
    VectorXd fiber_volume_fraction = VectorXd::Zero(num_elements_host);
    int host_id = 0;
    double volume_host = 0;
    double sum_fiber_volume_fraction = 0;

    for (int host_row = 0; host_row < num_elements_host; host_row++) {
      host_id = (*elements_host)(host_row, 0);
      volume_host = (*element_characteristic_host)(host_row);
      double total_length_embed = 0;
      for (int embed_row = 0; embed_row < num_elements_embed; embed_row++) {
        int embed_host_id = (*embed_map)(embed_row);
        double length_embed = 0;
        if (embed_host_id == host_id) {
          length_embed = (*element_characteristic_embed)(embed_row);
          total_length_embed = total_length_embed + length_embed;
        }
      }
      fiber_volume_fraction(host_row) = ( total_length_embed * area_truss ) / volume_host;
      sum_fiber_volume_fraction = sum_fiber_volume_fraction + fiber_volume_fraction(host_row);
    }

    double max_fiber_volume_fraction = fiber_volume_fraction.maxCoeff();
    double avg_fiber_volume_fraction = sum_fiber_volume_fraction / num_elements_host;

    // 0.785 is the maximum volume fraction for square packing of circles.
    if (max_fiber_volume_fraction > 0.785) {
      std::cout << "\n" << "Simulation Failed - Max Fiber Volume Fraction is Too Large" << "\n";
      std::cout << "Max Fiber Volume Fraction is: " << max_fiber_volume_fraction << "\n";
      std::cout << "Need to Reduce Truss Area. " << "\n" << "\n";
      std::exit(1);
    }

    // std::cout << "fiber_volume_fraction = " << '\n' << fiber_volume_fraction << '\n';
    // std::cout << "Max Fiber Volume Fraction is: " << max_fiber_volume_fraction << "\n";
    // std::cout << "Avg Fiber Volume Fraction is: " << avg_fiber_volume_fraction << "\n" << "\n";
    // std::exit(1);

}
