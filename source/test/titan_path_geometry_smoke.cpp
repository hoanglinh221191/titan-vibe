#include "titan_path/geometry.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

bool close(double lhs, double rhs, double tolerance = 1.0e-10) {
  return std::fabs(lhs - rhs) <= tolerance;
}

bool unit_and_orthogonal(const titan_path::CorridorSample &sample) {
  return close(titan_path::norm(sample.tangent), 1.0) &&
         close(titan_path::norm(sample.normal_1), 1.0) &&
         close(titan_path::norm(sample.normal_2), 1.0) &&
         close(titan_path::dot(sample.tangent, sample.normal_1), 0.0) &&
         close(titan_path::dot(sample.tangent, sample.normal_2), 0.0) &&
         close(titan_path::dot(sample.normal_1, sample.normal_2), 0.0);
}

} // namespace

int main() {
  titan_path::RouteDefinition route;
  route.route_id = "corner";
  route.waypoints = {
      {0, {0.0, 0.0, 0.0}, 2.0},
      {1, {4.0, 0.0, 0.0}, 3.0},
      {2, {4.0, 4.0, 0.0}, 4.0},
  };
  route.total_length_angstrom = 8.0;
  const auto corridor =
      titan_path::build_parallel_transport_corridor(route, 1.0);
  if (!corridor.ok) {
    std::cerr << corridor.error << '\n';
    return 1;
  }
  if (corridor.samples.size() != 9)
    return 2;
  for (std::size_t i = 0; i < corridor.samples.size(); ++i) {
    if (!unit_and_orthogonal(corridor.samples[i]))
      return 3;
    if (i != 0 && titan_path::dot(corridor.samples[i - 1].normal_1,
                                  corridor.samples[i].normal_1) < -1.0e-12)
      return 4;
  }
  if (!close(corridor.samples.front().radius_angstrom, 2.0) ||
      !close(corridor.samples.back().radius_angstrom, 4.0))
    return 5;

  const auto point = titan_path::corridor_point(corridor.samples[2], 1.0, 0.5);
  const titan_path::Vec3 displacement{
      point[0] - corridor.samples[2].center_angstrom[0],
      point[1] - corridor.samples[2].center_angstrom[1],
      point[2] - corridor.samples[2].center_angstrom[2]};
  if (!close(titan_path::norm(displacement), std::sqrt(1.25)))
    return 6;
  const auto middle =
      titan_path::sample_corridor_at_progress(corridor.samples, 3.5);
  if (!close(middle.progress_angstrom, 3.5) ||
      !unit_and_orthogonal(middle))
    return 7;
  try {
    (void)titan_path::corridor_point(corridor.samples[0], 3.0, 0.0);
    return 8;
  } catch (const std::exception &) {
  }

  const std::vector<titan_path::ArcLengthSegment> self_near{
      {0.0, 10.0, {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}},
      {10.0, 10.0, {10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}},
      {20.0, 10.0, {10.0, 10.0, 0.0}, {-10.0, 0.0, 0.0}},
      {30.0, 9.9, {0.0, 10.0, 0.0}, {0.0, -9.9, 0.0}},
  };
  const titan_path::Vec3 near_both_branches{0.01, 0.09, 0.0};
  const auto unrestricted = titan_path::project_point_to_polyline(
      self_near, near_both_branches, 39.0);
  if (!unrestricted.valid || unrestricted.segment_index != 3)
    return 9;
  const auto local = titan_path::project_point_to_polyline(
      self_near, near_both_branches, 1.0, 0.0, 3.0);
  if (!local.valid || local.segment_index != 0 ||
      !close(local.progress_angstrom, 0.01) ||
      !close(local.distance_angstrom, 0.09))
    return 10;
  const auto signed_backward = titan_path::project_point_to_polyline(
      self_near, {0.5, 0.2, 0.0}, 1.0, 0.0, 2.0);
  if (!signed_backward.valid || signed_backward.progress_angstrom >= 1.0 ||
      !close(signed_backward.progress_angstrom, 0.5))
    return 11;
  const auto no_window = titan_path::project_point_to_polyline(
      self_near, near_both_branches, 1.0, 15.0, 16.0);
  if (!no_window.valid || no_window.progress_angstrom < 15.0 ||
      no_window.progress_angstrom > 16.0)
    return 12;

  titan_path::RouteDefinition nearly_exact_spacing;
  nearly_exact_spacing.route_id = "near-exact-spacing";
  nearly_exact_spacing.waypoints = {
      {0, {31.667845918471, 87.357494209188, 49.572293245943}, 0.5},
      {1, {31.805397208101, 85.362229905211, 49.572293245943}, 0.5},
      {2, {31.874172852915, 84.364597753223, 49.572293245943}, 0.5},
  };
  const auto stable_spacing =
      titan_path::build_parallel_transport_corridor(nearly_exact_spacing, 0.25);
  if (!stable_spacing.ok || stable_spacing.samples.size() < 2 ||
      !unit_and_orthogonal(stable_spacing.samples.back()))
    return 13;

  // Reject the count before integer conversion or sample allocation, both
  // for a finite excessive count and for a ratio that overflows.
  const auto excessive_count = titan_path::build_parallel_transport_corridor(
      route, 8.0 / static_cast<double>(titan_path::kMaximumCorridorSamples));
  const auto overflowed_count = titan_path::build_parallel_transport_corridor(
      route, std::numeric_limits<double>::min());
  if (excessive_count.ok || overflowed_count.ok ||
      excessive_count.error.find("sample limit") == std::string::npos ||
      overflowed_count.error.find("sample limit") == std::string::npos)
    return 14;

  auto nonfinite_route = route;
  nonfinite_route.waypoints[1].position_angstrom[0] =
      std::numeric_limits<double>::quiet_NaN();
  if (titan_path::build_parallel_transport_corridor(nonfinite_route, 1.0).ok)
    return 15;
  auto excessive_length = route;
  excessive_length.waypoints[0].position_angstrom = {1.0e308, 0.0, 0.0};
  excessive_length.waypoints[1].position_angstrom = {-1.0e308, 0.0, 0.0};
  if (titan_path::build_parallel_transport_corridor(excessive_length, 1.0).ok)
    return 16;
  if (!std::isfinite(titan_path::norm({1.0e200, 1.0e200, 1.0e200})))
    return 17;
  bool projection_overflow_rejected = false;
  try {
    (void)titan_path::project_point_to_polyline(
        self_near, {1.0e308, 1.0e308, 0.0}, 0.0);
  } catch (const std::invalid_argument &) {
    projection_overflow_rejected = true;
  }
  if (!projection_overflow_rejected)
    return 18;
  bool corridor_offset_overflow_rejected = false;
  try {
    auto large_radius = corridor.samples.front();
    large_radius.radius_angstrom = 1.0e200;
    (void)titan_path::corridor_point(large_radius, 1.0e200, 1.0e200);
  } catch (const std::runtime_error &) {
    corridor_offset_overflow_rejected = true;
  }
  if (!corridor_offset_overflow_rejected)
    return 19;

  std::cout << "PASS titan_path_geometry_smoke\n";
  return 0;
}
