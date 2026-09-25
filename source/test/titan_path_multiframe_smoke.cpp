#include "titan_path/multiframe.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace titan_path;
MultiFrameInput fixture() {
  MultiFrameInput x;
  x.reference_core = {{-2, -2, -2}, {-2, -2, 2}, {-2, 2, -2}, {-2, 2, 2},
                      {2, -2, -2},  {2, -2, 2},  {2, 2, -2},  {2, 2, 2}};
  x.current_core = x.reference_core;
  x.waypoints = {{0, 0, 0}, {5, 0, 0}, {10, 0, 0}};
  x.receptor = {{100, 100, 100}};
  x.ligand = {{3, 1, 0}};
  x.ligand_masses = {12};
  x.old_direction = {1, 0, 0};
  x.displacement = {.25, 0, 0};
  return x;
}
int main() {
  int count = 0;
  auto test = [&](bool ok) {
    ++count;
    if (!ok)
      throw std::runtime_error("case " + std::to_string(count));
  };
  auto throws = [&](auto x) {
    try {
      (void)validate_multiframe_direction(x);
      return false;
    } catch (const std::invalid_argument &) {
      return true;
    }
  };
  try {
    auto x = fixture();
    auto r = validate_multiframe_direction(x);
    test(r.accepted && r.frames.size() == 3);
    test(r.frames[0].end_progress - r.frames[0].start_progress == .25);
    test(validate_multiframe_direction(x).reason == r.reason);
    x.displacement = {1, 0, 0};
    test(validate_multiframe_direction(x).reason ==
         "STOP_PROBE_DISPLACEMENT_MISMATCH");
    x = fixture();
    x.old_direction = {-1, 0, 0};
    test(validate_multiframe_direction(x).reason == "STOP_TURN_CAP");
    x = fixture();
    x.receptor = {{3.125, 1, 0}};
    test(validate_multiframe_direction(x).reason == "STOP_CLEARANCE");
    x = fixture();
    x.displacement = {0, .25, 0};
    x.old_direction = {0, 1, 0};
    test(validate_multiframe_direction(x).frames[0].reason ==
         "STOP_FRAME_FORWARD_PROGRESS");
    x = fixture();
    x.ligand = {{3, 3.99, 0}};
    x.displacement = {.2, .15, 0};
    x.old_direction = x.displacement;
    test(validate_multiframe_direction(x).frames[0].reason ==
         "STOP_FRAME_DISTANCE_BOUND");
    x = fixture();
    x.ligand = {{3, 5, 0}};
    x.native_mode = RouteMode::Recovery;
    x.displacement = {.1, -std::sqrt(.0525), 0};
    x.old_direction = x.displacement;
    test(validate_multiframe_direction(x).accepted);
    x = fixture();
    x.observed_progress = std::array<double,3>{3.,3.,3.};
    test(validate_multiframe_direction(x).accepted);
    x.observed_progress = std::array<double,3>{1.,1.,1.};
    test(!validate_multiframe_direction(x).accepted);
    x.observed_progress = std::array<double,3>{3.,3.2,3.};
    test(validate_multiframe_direction(x).frames[1].reason == "STOP_OBSERVED_PROGRESS_MISMATCH");
    x.observed_progress = std::array<double,3>{3.,-1.,3.};
    test(throws(x));
    x = fixture();
    x.ligand = {{3, 5, 0}};
    x.native_mode = RouteMode::Recovery;
    x.displacement = {.25, 0, 0};
    x.old_direction = x.displacement;
    test(validate_multiframe_direction(x).frames[0].reason ==
         "STOP_FRAME_RECOVERY_DESCENT");
    x.displacement = {-.1, -std::sqrt(.0525), 0};
    x.old_direction = x.displacement;
    test(validate_multiframe_direction(x).frames[0].reason ==
         "STOP_FRAME_BACKTRACK");
    x = fixture();
    x.ligand = {{3, 3.5, 0}};
    x.previous_mode = RouteMode::Recovery;
    x.native_mode = RouteMode::Recovery;
    x.displacement = {.1, -std::sqrt(.0525), 0};
    x.old_direction = x.displacement;
    test(validate_multiframe_direction(x).accepted);
    x = fixture();
    x.ligand = {{3, 3.95, 0}};
    x.native_mode = RouteMode::Recovery;
    x.displacement = {.1, -std::sqrt(.0525), 0};
    x.old_direction = x.displacement;
    for (std::size_t i = 1; i < x.current_core.size(); i += 2)
      x.current_core[i][1] -= .3;
    r = validate_multiframe_direction(x);
    test(!r.accepted && r.frames[1].reason == "STOP_FRAME_MODE_DISAGREEMENT");
    x = fixture();
    x.waypoints = {{0, 0, 0}, {.5, 0, 0}, {.5, .02, 0}, {0, .02, 0}};
    x.ligand = {{.2, .01, 0}};
    test(validate_multiframe_direction(x).frames[0].reason ==
         "STOP_FRAME_AMBIGUOUS");
    x = fixture();
    for (auto &p : x.current_core)
      p = {1, 1, 1};
    test(!validate_multiframe_direction(x).accepted);
    x = fixture();
    x.ligand_masses[0] = 0;
    test(throws(x));
    x = fixture();
    x.old_direction = {0, 0, 0};
    test(throws(x));
    x = fixture();
    x.current_core.pop_back();
    test(throws(x));
    x = fixture();
    x.displacement[0] = std::numeric_limits<double>::quiet_NaN();
    test(throws(x));
    x = fixture();
    x.waypoints[1] = x.waypoints[0];
    test(throws(x));
    x = fixture();
    r = validate_multiframe_direction(x);
    auto move = [](Vec3 p) { return Vec3{-p[1] + 17, p[0] - 4, p[2] + 3}; };
    for (auto *v : {&x.reference_core, &x.current_core, &x.waypoints,
                    &x.receptor, &x.ligand})
      for (auto &p : *v)
        p = move(p);
    x.old_direction = {0, 1, 0};
    x.displacement = {0, .25, 0};
    auto moved = validate_multiframe_direction(x);
    test(moved.accepted && std::abs(moved.clearance - r.clearance) < 1e-10);
    test(std::abs(moved.frames[0].end_progress - r.frames[0].end_progress) <
         1e-10);
    x = fixture();
    std::reverse(x.reference_core.begin(), x.reference_core.end());
    std::reverse(x.current_core.begin(), x.current_core.end());
    test(validate_multiframe_direction(x).accepted);
    std::cout << "TITAN_MULTIFRAME_SMOKE_PASS cases=" << count
              << " no_md=true\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
