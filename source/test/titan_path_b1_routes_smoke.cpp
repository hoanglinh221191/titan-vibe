#include "titan_path/route.h"

#include <filesystem>
#include <iostream>
#include <set>
#include <string>

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: titan_path_b1_routes_smoke <four .tpath files>\n";
    return 2;
  }
  const std::string reference_sha =
      "d1fd3cd5ab105ff799d84b8b08fa11ae7c94562ddbc73ac8af65730accd0774d";
  std::set<std::string> route_ids;
  for (int index = 1; index < argc; ++index) {
    const auto loaded =
        titan_path::load_route_file(std::filesystem::path(argv[index]));
    if (!loaded.ok) {
      std::cerr << argv[index] << ": " << loaded.error << '\n';
      return 1;
    }
    if (loaded.route.reference_structure_sha256 != reference_sha ||
        loaded.route.anchors.size() != 3 ||
        loaded.route.waypoints.size() != 5 ||
        loaded.route.start_index != 0 || loaded.route.terminal_index != 4 ||
        loaded.route.total_length_angstrom <= 0.0 ||
        loaded.source_file_sha256.size() != 64) {
      std::cerr << argv[index] << ": B1 route contract mismatch\n";
      return 1;
    }
    route_ids.insert(loaded.route.route_id);
  }
  if (route_ids.size() != 4) {
    std::cerr << "B1 route IDs are not unique\n";
    return 1;
  }
  std::cout << "TITAN_PATH_B1_ROUTES_SMOKE_PASS routes=4 anchors_per_route=3 "
               "waypoints_per_route=5\n";
  return 0;
}
