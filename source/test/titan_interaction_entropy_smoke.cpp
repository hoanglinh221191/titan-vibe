#include "titan_interaction_entropy.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <span>
#include <vector>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL titan_interaction_entropy_smoke: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

} // namespace

int main() {
  TitanEntropy::InteractionEntropyOptions options;
  options.rt_energy = 0.00831446261815324 * 298.15;
  options.filter_sigma = 0.0;
  options.min_samples = 10;
  options.block_size = 3;

  const std::vector<double> constant(10, -42.0);
  const auto stable = TitanEntropy::compute_interaction_entropy(
      std::span<const double>(constant.data(), constant.size()), options);
  require(stable.status == "ok" && stable.samples == 10 &&
              stable.filtered_samples == 10 && stable.ess == 10.0 &&
              std::abs(stable.minus_tds_energy) < 1.0e-14 &&
              stable.block_count == 3,
          "constant ten-frame series should pass with finite block diagnostics");

  const std::vector<double> short_series(9, -42.0);
  const auto short_result = TitanEntropy::compute_interaction_entropy(
      std::span<const double>(short_series.data(), short_series.size()),
      options);
  require(short_result.status == "insufficient_samples",
          "minimum sample gate must fail closed");

  TitanEntropy::InteractionEntropyOptions filtered_options = options;
  filtered_options.filter_sigma = 2.0;
  filtered_options.min_samples = 9;
  filtered_options.block_size = 0;
  const std::vector<double> filtered_constant = {
      -100.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  const auto filtered_result = TitanEntropy::compute_interaction_entropy(
      std::span<const double>(filtered_constant.data(),
                              filtered_constant.size()),
      filtered_options);
  require(filtered_result.status == "ok" &&
              filtered_result.filtered_samples == 9 &&
              filtered_result.ess == 9.0 &&
              std::abs(filtered_result.minus_tds_energy) < 1.0e-14,
          "retained samples must be re-centered after outlier filtering");

  filtered_options.min_samples = 10;
  const auto too_few_after_filter =
      TitanEntropy::compute_interaction_entropy(
          std::span<const double>(filtered_constant.data(),
                                  filtered_constant.size()),
          filtered_options);
  require(too_few_after_filter.status == "insufficient_filtered_samples" &&
              too_few_after_filter.filtered_samples == 9,
          "minimum sample gate must be re-applied after filtering");

  std::vector<double> invalid = constant;
  invalid[4] = std::numeric_limits<double>::quiet_NaN();
  const auto invalid_result = TitanEntropy::compute_interaction_entropy(
      std::span<const double>(invalid.data(), invalid.size()), options);
  require(invalid_result.status == "nonfinite_sample" &&
              !std::isnan(invalid_result.minus_tds_energy),
          "non-finite interaction-energy input must not report status ok");

  std::cout << "PASS titan_interaction_entropy_smoke\n";
  return EXIT_SUCCESS;
}
