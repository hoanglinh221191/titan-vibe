#include "titan_path/route.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

void check(bool value, const char *message) {
  if (!value) throw std::runtime_error(message);
}

int main(int argc, char **argv) {
  try {
    check(argc == 5, "expected original/runtime route pairs for two routes");
    for (int i = 1; i < argc; i += 2) {
      auto before = titan_path::load_route_file(argv[i]);
      auto after = titan_path::load_route_file(argv[i + 1]);
      check(before.ok, before.error.c_str());
      check(after.ok, after.error.c_str());
      const auto &a = before.route;
      const auto &b = after.route;
      check(b.route_id == a.route_id + "-approach2a-v1", "runtime route identity");
      check(a.reference_structure_sha256 == b.reference_structure_sha256,
            "reference identity changed");
      check(a.anchors.size() == b.anchors.size() && a.anchors.size() == 81,
            "anchor count changed");
      for (std::size_t j = 0; j < a.anchors.size(); ++j) {
        const auto &x = a.anchors[j]; const auto &y = b.anchors[j];
        check(x.original_atom_id == y.original_atom_id && x.chain == y.chain &&
              x.atom_name == y.atom_name && x.residue_position == y.residue_position &&
              x.reference_position_angstrom == y.reference_position_angstrom,
              "reference anchor changed");
      }
      check(a.waypoints.size() == 5 && b.waypoints.size() == 6 &&
            b.start_index == 0 && b.terminal_index == 5, "waypoint layout");
      for (std::size_t j = 0; j < a.waypoints.size(); ++j)
        check(a.waypoints[j].position_angstrom == b.waypoints[j+1].position_angstrom &&
              a.waypoints[j].radius_angstrom == b.waypoints[j+1].radius_angstrom,
              "original landmark or radius changed");
      check(std::abs(b.total_length_angstrom-a.total_length_angstrom-2.) < 1e-10,
            "approach length is not 2 A");
      double dot = 0., norm = 0.;
      for (int k = 0; k < 3; ++k) {
        const auto d = a.waypoints[1].position_angstrom[k]-a.waypoints[0].position_angstrom[k];
        const auto e = b.waypoints[1].position_angstrom[k]-b.waypoints[0].position_angstrom[k];
        dot += d*e; norm += d*d;
      }
      check(std::abs(dot-2.*std::sqrt(norm)) < 1e-10, "approach is not the initial tangent");
    }
    std::cout << "ROUTE_APPROACH_NATIVE_PARSER_PASS routes=2 original_landmarks=10 anchors_per_route=81 prefix_A=2 no_md=true\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
