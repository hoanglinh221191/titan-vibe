#include "titan_path/route.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

void write_text(const std::filesystem::path &path, const std::string &text) {
  std::ofstream output(path, std::ios::binary);
  output << text;
}

bool close(double lhs, double rhs, double tolerance = 1.0e-12) {
  return std::fabs(lhs - rhs) <= tolerance;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: titan_path_route_smoke <workdir>\n";
    return 2;
  }
  const std::filesystem::path workdir(argv[1]);
  std::filesystem::create_directories(workdir);

  const std::string sha(64, 'a');
  const auto good = workdir / "good.tpath";
  write_text(good,
             "TITAN_PATH 1\n"
             "route_id = route-A\n"
             "branch_id = front\n"
             "source = synthetic smoke route\n"
             "reference_structure_sha256 = " + sha + "\n"
             "registration = rigid\n"
             "units = nm\n"
             "start_index = 0\n"
             "terminal_index = 2\n"
             "anchor = 1 A 10 CA 0 0 0\n"
             "anchor = 2 A 11 CA 0.1 0 0\n"
             "anchor = 3 A 12 CA 0 0.1 0\n"
             "waypoint = 0 0 0 0 0.2\n"
             "waypoint = 1 0.1 0 0 0.25\n"
             "waypoint = 2 0.1 0.1 0 0.3\n");

  const auto loaded = titan_path::load_route_file(good);
  if (!loaded.ok) {
    std::cerr << loaded.error << '\n';
    return 3;
  }
  if (loaded.route.route_id != "route-A" ||
      loaded.route.registration != titan_path::RegistrationMode::Rigid ||
      loaded.route.anchors.size() != 3 || loaded.route.waypoints.size() != 3 ||
      loaded.source_file_sha256.size() != 64)
    return 4;
  if (!close(loaded.route.waypoints[1].position_angstrom[0], 1.0) ||
      !close(loaded.route.waypoints[2].radius_angstrom, 3.0) ||
      !close(loaded.route.total_length_angstrom, 2.0))
    return 5;

  const auto canonical_path = workdir / "canonical.tpath";
  write_text(canonical_path,
             titan_path::serialize_route_canonical(loaded.route));
  const auto canonical = titan_path::load_route_file(canonical_path);
  if (!canonical.ok || canonical.route.route_id != loaded.route.route_id ||
      !close(canonical.route.total_length_angstrom,
             loaded.route.total_length_angstrom))
    return 6;

  const auto collinear = workdir / "collinear.tpath";
  std::string invalid = titan_path::serialize_route_canonical(loaded.route);
  const std::string needle = "anchor = 3 A 12 CA 0 1 0";
  const std::size_t position = invalid.find(needle);
  if (position == std::string::npos)
    return 7;
  invalid.replace(position, needle.size(), "anchor = 3 A 12 CA 2 0 0");
  write_text(collinear, invalid);
  const auto rejected = titan_path::load_route_file(collinear);
  if (rejected.ok ||
      rejected.error.find("non-collinear") == std::string::npos)
    return 8;

  std::cout << "PASS titan_path_route_smoke\n";
  return 0;
}
