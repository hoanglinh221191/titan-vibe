#include "pbe/grid/ses_mask.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

std::size_t index3(std::size_t i, std::size_t j, std::size_t k,
                   const TitanPBE::SesMaskShape &shape) {
  return (k * shape.ny + j) * shape.nx + i;
}

std::vector<std::uint8_t>
sphere_mask(const TitanPBE::SesMaskShape &shape, double radius) {
  std::vector<std::uint8_t> mask(shape.points(), 0U);
  const double cx = 0.5 * static_cast<double>(shape.nx - 1) * shape.spacing;
  const double cy = 0.5 * static_cast<double>(shape.ny - 1) * shape.spacing;
  const double cz = 0.5 * static_cast<double>(shape.nz - 1) * shape.spacing;
  const double radius2 = radius * radius;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    const double z = static_cast<double>(k) * shape.spacing - cz;
    for (std::size_t j = 0; j < shape.ny; ++j) {
      const double y = static_cast<double>(j) * shape.spacing - cy;
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const double x = static_cast<double>(i) * shape.spacing - cx;
        mask[index3(i, j, k, shape)] =
            x * x + y * y + z * z <= radius2 ? 1U : 0U;
      }
    }
  }
  return mask;
}

std::size_t count_inside(const std::vector<std::uint8_t> &mask) {
  return static_cast<std::size_t>(std::count(mask.begin(), mask.end(), 1U));
}

} // namespace

int main() {
  using TitanPBE::SesMaskShape;

  const SesMaskShape shape{65, 65, 65, 0.25};
  constexpr double atom_radius = 2.0;
  constexpr double probe_radius = 1.4;
  const auto vdw = sphere_mask(shape, atom_radius);
  const auto sas = sphere_mask(shape, atom_radius + probe_radius);
  std::vector<std::uint8_t> ses(shape.points(), 0U);

  TitanPBE::erode_expanded_mask_to_ses(shape, sas, probe_radius, ses, nullptr);

  const std::size_t vdw_count = count_inside(vdw);
  const std::size_t sas_count = count_inside(sas);
  const std::size_t ses_count = count_inside(ses);
  assert(vdw_count <= ses_count);
  assert(ses_count < sas_count);

  const double voxel_volume = std::pow(shape.spacing, 3.0);
  const double exact_vdw_volume = 4.0 * std::acos(-1.0) / 3.0 *
                                  atom_radius * atom_radius * atom_radius;
  const double ses_volume = static_cast<double>(ses_count) * voxel_volume;
  assert(std::abs(ses_volume - exact_vdw_volume) / exact_vdw_volume < 0.25);

  std::vector<std::uint8_t> zero_probe(shape.points(), 0U);
  TitanPBE::erode_expanded_mask_to_ses(shape, vdw, 0.0, zero_probe, nullptr);
  assert(zero_probe == vdw);

  return 0;
}
