#include "titan_path/control.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// Arithmetic-only adapter used by the optional downloaded-evidence regression.
// It does not load a molecular system, select a candidate, or execute MD.
int main() {
  try {
    std::cout << std::setprecision(17);
    std::string cell;
    while (std::cin >> cell) {
      std::array<double, 3> requested{}, old{};
      double angle = 0.0, probe = 0.0;
      if (!(std::cin >> requested[0] >> requested[1] >> requested[2] >> old[0] >>
            old[1] >> old[2] >> angle >> probe))
        throw std::runtime_error("incomplete turn/probe input row");
      const auto projected =
          titan_path::project_vector_to_turn_cone(requested, old, angle);
      if (!projected.valid)
        throw std::runtime_error("turn projection rejected " + cell);
      const auto next = titan_path::scale_vector_to_length(projected.vector, probe);
      if (!next)
        throw std::runtime_error("probe scaling rejected " + cell);
      const double length = std::hypot((*next)[0], (*next)[1], (*next)[2]);
      std::cout << cell << '\t' << projected.capped << '\t'
                << projected.requested_angle_deg << '\t'
                << projected.applied_angle_deg << '\t' << (*next)[0] << '\t'
                << (*next)[1] << '\t' << (*next)[2] << '\t' << length << '\n';
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "B1_NAV_R2_TURN_PROBE_ERROR: " << error.what() << '\n';
    return 1;
  }
}
