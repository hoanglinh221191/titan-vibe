#include "titan_path/registration.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool close(double lhs, double rhs, double tolerance = 1.0e-10) {
  return std::fabs(lhs - rhs) <= tolerance;
}

} // namespace

int main() {
  const std::vector<titan_path::Vec3> reference{{0.0, 0.0, 0.0},
                                                 {1.0, 0.0, 0.0},
                                                 {0.0, 2.0, 0.0},
                                                 {0.0, 0.0, 3.0}};
  const std::vector<titan_path::Vec3> current{{4.0, -2.0, 1.0},
                                               {4.0, -1.0, 1.0},
                                               {2.0, -2.0, 1.0},
                                               {4.0, -2.0, 4.0}};
  const auto transform =
      titan_path::fit_rigid_transform(reference, current);
  if (!transform.valid) {
    std::cerr << transform.error << '\n';
    return 1;
  }
  if (!close(transform.determinant, 1.0) ||
      transform.rmsd_angstrom > 1.0e-10)
    return 2;
  for (std::size_t i = 0; i < reference.size(); ++i) {
    const auto mapped = titan_path::apply_transform(transform, reference[i]);
    for (std::size_t axis = 0; axis < 3; ++axis)
      if (!close(mapped[axis], current[i][axis]))
        return 3;
  }

  auto reflected = current;
  for (auto &point : reflected)
    point[2] *= -1.0;
  const auto reflection_fit =
      titan_path::fit_rigid_transform(reference, reflected);
  if (!reflection_fit.valid || reflection_fit.determinant < 0.999999 ||
      reflection_fit.rmsd_angstrom < 0.1)
    return 4;

  std::cout << "PASS titan_path_registration_smoke\n";
  return 0;
}
