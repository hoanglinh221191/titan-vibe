#include "titan_path/continuity.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

template <class F> void rejects(F f) {
  bool failed = false;
  try {
    f();
  } catch (const std::exception &) {
    failed = true;
  }
  assert(failed);
}

std::string replace_field(std::string text, const std::string &key,
                          const std::string &value) {
  const auto begin = text.find(key + "=");
  assert(begin != std::string::npos);
  const auto end = text.find('\n', begin);
  text.replace(begin, end - begin, key + "=" + value);
  return text;
}

int main() {
  using namespace titan_path;
  const std::string good = "SCHEMA=TITAN_PATH_STATUS_V1\nSTATUS=CONTINUE\n"
                           "ROUTE_KIND=EXPLICIT\nROUTE_ID=r\nROUTE_FILE_SHA256="
                           "route\nREFERENCE_SHA256=ref\n"
                           "MODE=Recovery\nPATH_PROGRESS_A=5\nPATH_TOTAL_A=20."
                           "2\nSNAPSHOT_TIME_PS=20\n";
  auto read = [&](const std::string &s) {
    std::istringstream input(s);
    return read_previous_route_state(input, "r", "route", "ref", 20.2, 30);
  };
  auto state = read(good);
  assert(state.mode == RouteMode::Recovery && state.progress_angstrom == 5);
  rejects([&] { read(good + "MODE=Advance\n"); });
  rejects([&] { read(good + "BROKEN\n"); });
  for (const auto &replacement : {std::string("nan"), std::string("inf"),
                                  std::string("21"), std::string("5junk")}) {
    auto bad = good;
    bad.replace(bad.find("PATH_PROGRESS_A=5"), 17,
                "PATH_PROGRESS_A=" + replacement);
    rejects([&] { read(bad); });
  }
  auto terminal = good;
  terminal.replace(terminal.find("STATUS=CONTINUE"), 15, "STATUS=STOP");
  rejects([&] { read(terminal); });
  for (const auto &mismatch :
       {std::pair<std::string, std::string>{"ROUTE_ID", "other"},
        {"ROUTE_FILE_SHA256", "other"},
        {"REFERENCE_SHA256", "other"},
        {"SNAPSHOT_TIME_PS", "31"},
        {"SNAPSHOT_TIME_PS", "-1"},
        {"PATH_TOTAL_A", "20.3"},
        {"STATUS", "COMPLETE"}}) {
    rejects(
        [&] { read(replace_field(good, mismatch.first, mismatch.second)); });
  }
  // Re-evaluating a fixed snapshot at the same timestamp is permitted.
  assert(read(replace_field(good, "SNAPSHOT_TIME_PS", "30")).snapshot_time_ps ==
         30);
  std::vector<ArcLengthSegment> route{{0, 10, {0, 0, 0}, {10, 0, 0}},
                                      {10, .2, {10, 0, 0}, {0, .2, 0}},
                                      {10.2, 10, {10, .2, 0}, {-10, 0, 0}}};
  auto local = project_with_continuity(route, {5, .101, 0}, 5, 1, 20.2);
  auto global = project_point_to_polyline(route, {5, .101, 0}, 0);
  assert(std::fabs(local.progress_angstrom - 5) < 1e-12);
  assert(std::fabs(global.progress_angstrom - 15.2) < 1e-12);
  rejects(
      [&] { (void)project_with_continuity(route, {6.2, 0, 0}, 5, 1, 20.2); });
  auto exact = project_with_continuity(route, {6, 0, 0}, 5, 1, 20.2);
  assert(exact.progress_angstrom == 6);
  auto ingress = project_with_continuity(route, {-1, 0, 0}, 0, 1, 20.2);
  assert(ingress.progress_angstrom == 0); // real route endpoint, not clipping

  const std::vector<ArcLengthSegment> corner{{0, 1, {0, 0, 0}, {1, 0, 0}},
                                             {1, 1, {1, 0, 0}, {0, 1, 0}}};
  // A boundary at a route vertex must also inspect the adjoining segment,
  // not merely the segment selected by the bounded projection.
  rejects(
      [&] { (void)project_with_continuity(corner, {1, .2, 0}, .5, .5, 2); });
  rejects(
      [&] { (void)project_with_continuity(corner, {.8, 0, 0}, 1.5, .5, 2); });
  const auto genuine_corner_high =
      project_with_continuity(corner, {1.2, -.2, 0}, .5, .5, 2);
  const auto genuine_corner_low =
      project_with_continuity(corner, {1.2, -.2, 0}, 1.5, .5, 2);
  assert(genuine_corner_high.progress_angstrom == 1);
  assert(genuine_corner_low.progress_angstrom == 1);
  const auto physical_terminal =
      project_with_continuity(corner, {1, 2, 0}, 1.5, .5, 2);
  assert(physical_terminal.progress_angstrom == 2);
  std::cout << "PASS titan_path_continuity_smoke state_identity_and_window\n";
}
