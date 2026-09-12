#include "titan_energy_validation.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

MOL make_frame(double x, double y = 0.0, double z = 0.0) {
  MOL frame;
  frame.pos.x = {x};
  frame.pos.y = {y};
  frame.pos.z = {z};
  return frame;
}

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL titan_energy_validation_smoke: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

} // namespace

int main() {
  MOL receptor = make_frame(0.0);
  MOL ligand = make_frame(0.3); // 0.03 nm: finite, but invalid for GBSA.
  TitanExecutor serial(1);
  TitanExecutor parallel(4);

  require(!titan_energy_validation::find_overlapping_atoms(
               receptor, ligand, false, serial),
          "finite close contact must not be classified as exact overlap");
  const auto serial_close = titan_energy_validation::find_overlapping_atoms(
      receptor, ligand, false, serial,
      titan_energy_validation::kMinimumGbsaPairDistanceNm);
  const auto parallel_close = titan_energy_validation::find_overlapping_atoms(
      receptor, ligand, false, parallel,
      titan_energy_validation::kMinimumGbsaPairDistanceNm);
  require(serial_close.has_value() && parallel_close.has_value(),
          "GBSA close contact must fail closed for nt=1 and nt=4");
  require(std::abs(serial_close->distance_nm - 0.03) < 1.0e-15 &&
              serial_close->mol1_atom == parallel_close->mol1_atom &&
              serial_close->mol2_atom == parallel_close->mol2_atom,
          "close-contact diagnostic must be deterministic");

  ligand.pos.y[0] = std::numeric_limits<double>::quiet_NaN();
  const auto invalid =
      titan_energy_validation::find_nonfinite_coordinate(receptor, ligand);
  require(invalid.has_value() && invalid->molecule == 2 &&
              invalid->atom == 0 && invalid->axis == 'y',
          "NaN coordinate must be rejected with molecule/atom/axis provenance");

  ligand = make_frame(1.0, 2.0, std::numeric_limits<double>::infinity());
  const auto infinite =
      titan_energy_validation::find_nonfinite_coordinate(receptor, ligand);
  require(infinite.has_value() && infinite->axis == 'z',
          "infinite coordinate must be rejected");

  std::cout << "PASS titan_energy_validation_smoke\n";
  return EXIT_SUCCESS;
}
