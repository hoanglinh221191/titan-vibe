#include "titan_path/navigation.h"

#include <cmath>
#include <iostream>
#include <optional>

namespace {
bool close(double lhs, double rhs, double tolerance = 1.0e-10) {
  return std::fabs(lhs - rhs) <= tolerance;
}
} // namespace

int main() {
  const auto inside = titan_path::evaluate_navigation_cost(
      0.4, 0.5, 0.6, 1.0, 2.0, 3.0);
  if (!close(inside.total_cost, 0.0) || !close(inside.tube_cost, 0.0) ||
      !close(inside.backtrack_cost, 0.0))
    return 1;

  const auto outside = titan_path::evaluate_navigation_cost(
      1.5, 0.5, -0.5, 0.5, 2.0, 4.0);
  if (!close(outside.tube_excess_angstrom, 1.0) ||
      !close(outside.backtrack_distance_angstrom, 0.5) ||
      !close(outside.tube_cost, 4.0) ||
      !close(outside.backtrack_cost, 2.0) ||
      !close(outside.total_cost, 6.0))
    return 2;

  const auto positive_a = titan_path::evaluate_navigation_cost(
      1.0, 0.5, 0.1, 1.0, 2.0, 3.0);
  const auto positive_b = titan_path::evaluate_navigation_cost(
      1.0, 0.5, 2.0, 1.0, 2.0, 3.0);
  if (!close(positive_a.total_cost, positive_b.total_cost))
    return 3; // proves there is no hidden forward reward

  titan_path::PathProjection projection;
  projection.valid = true;
  projection.point_angstrom = {0.0, 0.0, 0.0};
  projection.distance_angstrom = 2.0;
  titan_path::CorridorSample frame;
  frame.center_angstrom = projection.point_angstrom;
  frame.tangent = {1.0, 0.0, 0.0};
  frame.normal_1 = {0.0, 1.0, 0.0};
  frame.normal_2 = {0.0, 0.0, 1.0};
  frame.radius_angstrom = 0.5;

  const auto advance = titan_path::build_navigation_guide_direction(
      {0.0, 2.0, 0.0}, projection, frame,
      titan_path::RouteMode::Advance, 1.0);
  if (!(advance[0] > 0.0 && advance[1] < 0.0) ||
      !close(titan_path::norm(advance), 1.0))
    return 4;
  const auto recovery = titan_path::build_navigation_guide_direction(
      {0.0, 2.0, 0.0}, projection, frame,
      titan_path::RouteMode::Recovery, 1.0);
  if (!(recovery[1] < 0.0 && std::fabs(recovery[1]) > recovery[0]) ||
      !close(titan_path::norm(recovery), 1.0))
    return 5;

  const auto seeds = titan_path::build_navigation_seed_directions(
      advance, titan_path::Vec3{1.0, 0.0, 0.0}, frame, 8, 0.35);
  const auto seeds_repeat = titan_path::build_navigation_seed_directions(
      advance, titan_path::Vec3{1.0, 0.0, 0.0}, frame, 8, 0.35);
  if (seeds.size() != seeds_repeat.size() || seeds.size() < 9)
    return 6;
  for (std::size_t index = 0; index < seeds.size(); ++index) {
    if (seeds[index].index != index ||
        seeds[index].direction_unit != seeds_repeat[index].direction_unit ||
        !close(titan_path::norm(seeds[index].direction_unit), 1.0))
      return 7;
  }

  std::cout << "PASS titan_path_navigation_smoke\n";
  return 0;
}
