// Checks the rgbmax cutoff added to titan_gb::hct_pair_integral_inverse_angstrom.
//
// Three properties matter, and this driver exercises the shipped header rather
// than a copy of it:
//   1. rgbmax = +inf (the TITAN default) must reproduce the uncut integral
//      bit-for-bit, so existing results stay reproducible.
//   2. In the boundary regime (rgbmax - s < r <= rgbmax + s) the clamped upper
//      limit must agree with Amber's expanded egb boundary branch.
//   3. Beyond rgbmax + s the contribution must be exactly zero, and the value
//      must converge monotonically to the uncut integral as rgbmax grows.
#include "../src/titan_gb_model.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>

namespace {

int failures = 0;

void check(bool ok, const char *label) {
  std::printf("%-4s %s\n", ok ? "ok" : "FAIL", label);
  if (!ok)
    ++failures;
}

// Amber egb.F90 boundary branch, valid when offset_i < r - s.
double amber_boundary(double r, double s, double rgbmax) {
  const double uij = 1.0 / (r - s);
  return 0.125 * (1.0 / r) *
         (1.0 + 2.0 * r * uij +
          (r * r - 4.0 * rgbmax * r - s * s) / (rgbmax * rgbmax) +
          2.0 * std::log((r - s) / rgbmax));
}

} // namespace

int main() {
  constexpr double kInf = std::numeric_limits<double>::infinity();
  std::mt19937 rng(20260913u);
  std::uniform_real_distribution<double> offset_dist(0.6, 1.8);
  std::uniform_real_distribution<double> scaled_dist(0.5, 2.2);
  std::uniform_real_distribution<double> dist_dist(0.2, 40.0);

  // 1. Default parameter and explicit +inf must match the three-argument call
  //    exactly, for the ordinary case and for the "i inside j" correction.
  bool default_matches = true;
  bool inf_matches = true;
  for (int trial = 0; trial < 20000; ++trial) {
    const double offset = offset_dist(rng);
    const double scaled = scaled_dist(rng);
    const double distance = dist_dist(rng);
    const double three_arg =
        titan_gb::hct_pair_integral_inverse_angstrom(offset, scaled, distance);
    const double with_inf = titan_gb::hct_pair_integral_inverse_angstrom(
        offset, scaled, distance, kInf);
    if (std::isnan(three_arg) != std::isnan(with_inf))
      default_matches = false;
    else if (!std::isnan(three_arg) && three_arg != with_inf)
      default_matches = false;
    // A huge finite cutoff must also leave these short-range pairs untouched.
    const double with_huge = titan_gb::hct_pair_integral_inverse_angstrom(
        offset, scaled, distance, 1.0e6);
    if (!std::isnan(three_arg) && three_arg != with_huge)
      inf_matches = false;
  }
  check(default_matches, "rgbmax=+inf reproduces the uncut integral exactly");
  check(inf_matches, "a 1e6 A cutoff leaves short-range pairs unchanged");

  // 2. Boundary regime against Amber's expanded form.
  // Tolerance is mixed absolute/relative on purpose. The integral goes to zero
  // as r approaches rgbmax + s, and Amber's expanded polynomial loses accuracy
  // to cancellation exactly there, so a pure relative test would penalise the
  // better-conditioned G(L) - G(U) form for Amber's own rounding error.
  constexpr double kRgbmax = 25.0;
  double worst_rel = 0.0;
  double worst_abs = 0.0;
  bool boundary_ok = true;
  for (int trial = 0; trial < 2000; ++trial) {
    const double scaled = scaled_dist(rng);
    std::uniform_real_distribution<double> boundary(kRgbmax - scaled + 1e-6,
                                                    kRgbmax + scaled);
    const double r = boundary(rng);
    const double offset = 0.9; // safely below r - s
    const double ours = titan_gb::hct_pair_integral_inverse_angstrom(
        offset, scaled, r, kRgbmax);
    const double amber = amber_boundary(r, scaled, kRgbmax);
    const double abs_diff = std::abs(ours - amber);
    worst_abs = std::max(worst_abs, abs_diff);
    worst_rel =
        std::max(worst_rel, abs_diff / std::max(std::abs(amber), 1e-30));
    if (abs_diff > 1.0e-14 + 1.0e-6 * std::abs(amber))
      boundary_ok = false;
  }
  std::printf("     vs Amber boundary form: worst abs %.3e, worst rel %.3e\n",
              worst_abs, worst_rel);
  check(boundary_ok, "clamped upper limit matches Amber's egb boundary");

  // 3. Exact zero beyond rgbmax + s.
  bool zero_beyond = true;
  for (int trial = 0; trial < 2000; ++trial) {
    const double scaled = scaled_dist(rng);
    const double offset = offset_dist(rng);
    std::uniform_real_distribution<double> far(kRgbmax + scaled + 1e-9,
                                               kRgbmax + scaled + 50.0);
    if (titan_gb::hct_pair_integral_inverse_angstrom(offset, scaled, far(rng),
                                                     kRgbmax) != 0.0)
      zero_beyond = false;
  }
  check(zero_beyond, "pairs beyond rgbmax + s contribute exactly zero");

  // 4. Monotone convergence to the uncut value as the cutoff grows.
  const double offset = 1.1, scaled = 1.4, distance = 18.0;
  const double uncut =
      titan_gb::hct_pair_integral_inverse_angstrom(offset, scaled, distance);
  double previous_gap = std::numeric_limits<double>::infinity();
  bool monotone = true;
  for (const double rgbmax : {20.0, 25.0, 40.0, 60.0, 100.0, 400.0}) {
    const double value = titan_gb::hct_pair_integral_inverse_angstrom(
        offset, scaled, distance, rgbmax);
    const double gap = std::abs(value - uncut);
    if (gap > previous_gap)
      monotone = false;
    previous_gap = gap;
  }
  check(monotone, "value converges monotonically to the uncut integral");
  check(previous_gap == 0.0, "a 400 A cutoff reproduces the uncut value");

  std::printf("\n%d check(s) failed\n", failures);
  return failures == 0 ? 0 : 1;
}
