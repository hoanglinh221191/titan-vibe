#include "titan_sasa.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <span>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double binding_delta_sasa(int point_count) {
  constexpr double radius = 1.7;
  constexpr double separation = 4.0;
  const std::vector<TitanSASA::Atom> complex = {
      {0.0, 0.0, 0.0, radius},
      {separation, 0.0, 0.0, radius},
  };
  TitanSASA::Options options;
  options.probe_radius = TitanSASA::kDefaultProbeRadiusA;
  options.point_count = point_count;
  const std::vector<double> sasa = TitanSASA::compute_per_atom_sasa(
      std::span<const TitanSASA::Atom>(complex.data(), complex.size()),
      options);
  const double expanded_radius = radius + options.probe_radius;
  const double isolated_total = 8.0 * kPi * expanded_radius * expanded_radius;
  return TitanSASA::sum_sasa(sasa) - isolated_total;
}

} // namespace

int main() {
  constexpr double expanded_radius = 1.7 + TitanSASA::kDefaultProbeRadiusA;
  constexpr double cap_height = expanded_radius - 2.0;
  const double oracle = -4.0 * kPi * expanded_radius * cap_height;

  const double delta_96 = binding_delta_sasa(96);
  const double delta_384 = binding_delta_sasa(384);
  const double delta_1536 = binding_delta_sasa(1536);
  const double repeated_1536 = binding_delta_sasa(1536);

  const double error_384 = std::abs(delta_384 - oracle);
  const double error_1536 = std::abs(delta_1536 - oracle);
  const double nonpolar =
      TitanSASA::nonpolar_energy_from_delta_sasa_kj(delta_1536);

  if (!std::isfinite(delta_96) || !std::isfinite(delta_384) ||
      !std::isfinite(delta_1536) || error_384 > 3.0 || error_1536 > 1.5 ||
      std::abs(delta_1536 - repeated_1536) > 1.0e-12 ||
      !(nonpolar < 0.0)) {
    std::cerr << "FAIL titan_gbsa_sasa_convergence_smoke oracle=" << oracle
              << " delta96=" << delta_96 << " delta384=" << delta_384
              << " delta1536=" << delta_1536 << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "PASS titan_gbsa_sasa_convergence_smoke oracle=" << oracle
            << " delta96=" << delta_96 << " delta384=" << delta_384
            << " delta1536=" << delta_1536 << '\n';
  return EXIT_SUCCESS;
}
