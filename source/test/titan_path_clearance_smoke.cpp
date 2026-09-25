#include "titan_path/clearance.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
bool close(double lhs, double rhs, double tolerance = 1.0e-12) {
  return std::fabs(lhs - rhs) <= tolerance;
}

titan_path::CartesianCoordinatesView view(const std::vector<double> &x,
                                          const std::vector<double> &y,
                                          const std::vector<double> &z) {
  return {x, y, z};
}
} // namespace

int main() {
  const auto interior = titan_path::closest_approach_during_translation(
      {0.0, 0.0, 0.0}, {-2.0, 1.0, 0.0}, {4.0, 0.0, 0.0});
  if (!close(interior.distance_squared, 1.0) ||
      !close(interior.fraction, 0.5))
    return 1;
  const auto at_start = titan_path::closest_approach_during_translation(
      {0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {1.0, 0.0, 0.0});
  if (!close(at_start.distance_squared, 0.25) ||
      !close(at_start.fraction, 0.0))
    return 2;
  const auto at_end = titan_path::closest_approach_during_translation(
      {0.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}, {1.5, 0.0, 0.0});
  if (!close(at_end.distance_squared, 0.25) ||
      !close(at_end.fraction, 1.0))
    return 3;
  const auto reverse = titan_path::closest_approach_during_translation(
      {0.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {-4.0, 0.0, 0.0});
  if (!close(reverse.distance_squared, interior.distance_squared) ||
      !close(reverse.fraction, 1.0 - interior.fraction))
    return 4;

  const std::vector<double> rx{0.0}, ry{0.0}, rz{0.0};
  const std::vector<double> lx{-2.0}, ly{1.0}, lz{0.0};
  const auto crossing = titan_path::evaluate_translation_sweep(
      view(rx, ry, rz), view(lx, ly, lz), {4.0, 0.0, 0.0}, 1.1);
  if (crossing.passes_threshold || crossing.global_minimum_complete ||
      !close(crossing.witness_distance_angstrom, 1.0) ||
      !close(crossing.witness_fraction, 0.5) ||
      !(crossing.endpoint_distance_angstrom > 2.0) ||
      !close(crossing.witness_start_distance_angstrom, std::sqrt(5.0)) ||
      !close(crossing.witness_endpoint_distance_angstrom, std::sqrt(5.0)) ||
      !crossing.threshold_violation_hidden_at_endpoints)
    return 5; // both endpoints are safe, the continuous sweep is not
  const auto equality = titan_path::evaluate_translation_sweep(
      view(rx, ry, rz), view(lx, ly, lz), {4.0, 0.0, 0.0}, 1.0);
  if (!equality.passes_threshold || !equality.global_minimum_complete ||
      !close(equality.witness_distance_angstrom, 1.0))
    return 6; // preserve the legacy strict '< threshold' convention

  const std::vector<double> trust_x{-1.0}, trust_y{0.5}, trust_z{0.0};
  const auto long_probe = titan_path::evaluate_translation_sweep(
      view(rx, ry, rz), view(trust_x, trust_y, trust_z), {1.0, 0.0, 0.0},
      0.6);
  const auto short_probe = titan_path::evaluate_translation_sweep(
      view(rx, ry, rz), view(trust_x, trust_y, trust_z), {0.5, 0.0, 0.0},
      0.6);
  if (long_probe.passes_threshold || !short_probe.passes_threshold ||
      !short_probe.global_minimum_complete)
    return 7;

  bool empty_rejected = false;
  try {
    const std::vector<double> empty;
    (void)titan_path::evaluate_translation_sweep(
        view(empty, empty, empty), view(lx, ly, lz), {1.0, 0.0, 0.0}, 0.6);
  } catch (const std::invalid_argument &) {
    empty_rejected = true;
  }
  if (!empty_rejected)
    return 8;

  // A finite input must not silently produce an infinite/NaN intermediate
  // and then report a successful collision-free sweep.
  bool huge_sweep_rejected = false;
  try {
    const std::vector<double> huge_x{1.0e200};
    (void)titan_path::evaluate_translation_sweep(
        view(rx, ry, rz), view(huge_x, ly, lz), {-1.0e200, 0.0, 0.0}, 1.5);
  } catch (const std::overflow_error &) {
    huge_sweep_rejected = true;
  }
  if (!huge_sweep_rejected)
    return 9;
  bool huge_relative_rejected = false;
  try {
    (void)titan_path::closest_approach_during_translation(
        {1.0e308, 0.0, 0.0}, {-1.0e308, 0.0, 0.0}, {1.0, 0.0, 0.0});
  } catch (const std::overflow_error &) {
    huge_relative_rejected = true;
  }
  if (!huge_relative_rejected)
    return 10;
  bool huge_threshold_rejected = false;
  try {
    (void)titan_path::evaluate_translation_sweep(
        view(rx, ry, rz), view(lx, ly, lz), {1.0, 0.0, 0.0}, 1.0e200);
  } catch (const std::overflow_error &) {
    huge_threshold_rejected = true;
  }
  if (!huge_threshold_rejected)
    return 11;
  bool nonfinite_rejected = false;
  try {
    const std::vector<double> invalid_x{
        std::numeric_limits<double>::quiet_NaN()};
    (void)titan_path::evaluate_translation_sweep(
        view(rx, ry, rz), view(invalid_x, ly, lz), {1.0, 0.0, 0.0}, 1.5);
  } catch (const std::invalid_argument &) {
    nonfinite_rejected = true;
  }
  if (!nonfinite_rejected)
    return 12;

  // Exercise the conservative representability boundary, including the
  // largest possible endpoint components (3B) in the unchecked pair kernel.
  const double component_limit =
      titan_path::translation_clearance_component_limit_angstrom();
  const std::vector<double> low{-component_limit}, high{component_limit};
  const auto boundary_sweep = titan_path::evaluate_translation_sweep(
      view(low, low, low), view(high, high, high),
      {component_limit, component_limit, component_limit}, component_limit);
  if (!boundary_sweep.passes_threshold || !boundary_sweep.global_minimum_complete ||
      !std::isfinite(boundary_sweep.witness_distance_angstrom) ||
      !std::isfinite(boundary_sweep.endpoint_distance_angstrom) ||
      !close(boundary_sweep.witness_distance_angstrom / component_limit,
             std::sqrt(12.0)) ||
      !close(boundary_sweep.endpoint_distance_angstrom / component_limit,
             std::sqrt(27.0)))
    return 13;
  const auto boundary_collision = titan_path::evaluate_translation_sweep(
      view(rx, ry, rz), view(low, ry, rz), {component_limit, 0.0, 0.0}, 1.5);
  if (boundary_collision.passes_threshold ||
      boundary_collision.witness_fraction != 1.0 ||
      boundary_collision.witness_distance_angstrom != 0.0)
    return 14;
  const double outside_limit = std::nextafter(
      component_limit, std::numeric_limits<double>::infinity());
  for (int field = 0; field < 3; ++field) {
    bool outside_rejected = false;
    try {
      const std::vector<double> candidate_x{field == 0 ? outside_limit : 2.0};
      (void)titan_path::evaluate_translation_sweep(
          view(rx, ry, rz), view(candidate_x, ry, rz),
          {field == 1 ? outside_limit : 1.0, 0.0, 0.0},
          field == 2 ? outside_limit : 1.5);
    } catch (const std::overflow_error &) {
      outside_rejected = true;
    }
    if (!outside_rejected)
      return 15;
  }
  // Invalid atoms must be rejected before an earlier collision can trigger
  // the pair loop's early exit and hide the invalid input.
  bool late_bad_atom_rejected = false;
  try {
    const std::vector<double> bad_x{0.0, outside_limit}, zeros{0.0, 0.0};
    (void)titan_path::evaluate_translation_sweep(
        view(bad_x, zeros, zeros), view(rx, ry, rz), {1.0, 0.0, 0.0}, 1.5);
  } catch (const std::overflow_error &) {
    late_bad_atom_rejected = true;
  }
  if (!late_bad_atom_rejected)
    return 16;
  const auto distant_short = titan_path::closest_approach_during_translation(
      {0.0, 0.0, 0.0}, {-component_limit, 0.0, 0.0}, {1.0e-11, 0.0, 0.0});
  if (distant_short.fraction != 1.0 || !std::isfinite(distant_short.distance_squared))
    return 17;

  std::cout << "PASS titan_path_clearance_smoke\n";
  return 0;
}
