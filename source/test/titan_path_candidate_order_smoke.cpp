#include "pathway_de_variants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

PathwayCandidateState candidate(double energy, double path_distance = 1.0) {
  PathwayCandidateState result;
  result.valid = true;
  result.fixed_probe = true;
  result.probe_length = 1.0;
  result.pull_distance = 1.0;
  result.delta_progress = 0.5;
  result.applied_path_distance = path_distance;
  result.energy_key = energy;
  result.angle_deg = 20.0;
  return result;
}

void check_strict_weak_order(const std::vector<PathwayCandidateState> &items,
                             PathwayMode mode) {
  const auto less = [&](std::size_t a, std::size_t b) {
    return pathway_candidate_strict_less(items[a], a, items[b], b, mode);
  };
  for (std::size_t a = 0; a < items.size(); ++a) {
    require(!less(a, a), "order is not irreflexive");
    for (std::size_t b = 0; b < items.size(); ++b) {
      require(!(less(a, b) && less(b, a)), "order is not asymmetric");
      for (std::size_t c = 0; c < items.size(); ++c) {
        if (less(a, b) && less(b, c))
          require(less(a, c), "order is not transitive");
        if (!less(a, b) && !less(b, a) && !less(b, c) && !less(c, b))
          require(!less(a, c) && !less(c, a),
                  "order equivalence is not transitive");
      }
    }
  }
}

void test_epsilon_cycle_and_permutations(PathwayMode mode) {
  // IDs oppose the exact key order. Pairwise epsilon ties form A<C<B<A.
  const double eps = mode == PathwayMode::Advance ? 1.0e-9 : 1.0e-6;
  std::vector<PathwayCandidateState> items{
      candidate(0.0), candidate(0.0), candidate(0.0)};
  for (std::size_t index = 0; index < items.size(); ++index) {
    const double value = (2 - index) * 0.75 * eps;
    if (mode == PathwayMode::Advance)
      items[index].energy_key = value;
    else
      items[index].applied_path_distance = value;
  }
  require(pathway_candidate_better(items[2], 2, items[0], 0, mode) &&
              pathway_candidate_better(items[0], 0, items[1], 1, mode) &&
              pathway_candidate_better(items[1], 1, items[2], 2, mode),
          "fixture must retain the old tolerance decision cycle");
  require(pathway_candidate_equivalent(items[1], items[2], mode),
          "tolerance equivalence changed unexpectedly");
  check_strict_weak_order(items, mode);

  std::array<std::size_t, 3> input{0, 1, 2};
  const std::array<std::size_t, 3> expected{2, 1, 0};
  do {
    for (bool stable : {false, true}) {
      auto actual = input;
      const auto less = [&](std::size_t a, std::size_t b) {
        return pathway_candidate_strict_less(items[a], a, items[b], b, mode);
      };
      if (stable)
        std::stable_sort(actual.begin(), actual.end(), less);
      else
        std::sort(actual.begin(), actual.end(), less);
      require(actual == expected, "sort depends on input permutation");
    }
  } while (std::next_permutation(input.begin(), input.end()));

  const auto seeds = select_pathway_pattern_seed_indices(items, 3, 0.0, mode);
  require(seeds == std::vector<std::size_t>({2, 1, 0}),
          "DE pattern ranking did not use exact order");
}

void test_mixed_diagnostic_populations() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<PathwayCandidateState> items;
  for (double value : {-inf, -1.0, -0.0, 0.0, 0.75e-9, 1.5e-9, 1.0, inf, nan}) {
    auto item = candidate(value, value);
    items.push_back(item);
    item.fixed_probe = false;
    items.push_back(item);
    item.valid = false;
    items.push_back(item);
  }
  auto item = candidate(0.0);
  for (double value : {0.5, 0.5 + 0.75e-6, 0.5 + 1.5e-6, inf, nan}) {
    item.probe_length = value;
    items.push_back(item);
  }
  for (PathwayMode mode : {PathwayMode::Advance, PathwayMode::Recovery})
    check_strict_weak_order(items, mode);
}

} // namespace

int main() {
  try {
    test_epsilon_cycle_and_permutations(PathwayMode::Advance);
    test_epsilon_cycle_and_permutations(PathwayMode::Recovery);
    test_mixed_diagnostic_populations();
    std::cout << "PASS titan_path_candidate_order_smoke cases=3 "
                 "sort_permutations=24 tolerance_decisions_preserved=true\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL titan_path_candidate_order_smoke: " << error.what() << '\n';
    return 1;
  }
}
