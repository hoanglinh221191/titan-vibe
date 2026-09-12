#include "tabi/tabi_ses_quartic.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

using TitanPBE::Tabi::analytic_ses::QuarticDiagnosticFailure;
using TitanPBE::Tabi::analytic_ses::QuarticDiagnostics;
using TitanPBE::Tabi::analytic_ses::QuarticRefinementTermination;
using TitanPBE::Tabi::analytic_ses::QuarticRoots;
using TitanPBE::Tabi::analytic_ses::solve_real_quartic;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

void check_roots(const QuarticRoots &actual,
                 const std::initializer_list<double> expected,
                 double tolerance = 2.0e-11) {
  if (!actual.certified) {
    std::cerr << "quartic_uncertified count=" << actual.count << " roots=";
    for (std::size_t root = 0; root < actual.count; ++root)
      std::cerr << ' ' << actual.values[root];
    std::cerr << " residual=" << actual.max_normalized_residual << '\n';
    throw std::runtime_error("quartic roots were not certified");
  }
  if (actual.count != expected.size()) {
    std::cerr << "quartic_root_count actual=" << actual.count
              << " expected=" << expected.size() << " roots=";
    for (std::size_t root = 0; root < actual.count; ++root)
      std::cerr << ' ' << actual.values[root];
    std::cerr << " residual=" << actual.max_normalized_residual << '\n';
    throw std::runtime_error("unexpected real root count");
  }
  std::size_t i = 0;
  for (const double value : expected) {
    if (std::abs(actual.values[i] - value) > tolerance) {
      std::cerr << "quartic_root_mismatch index=" << i
                << " actual=" << actual.values[i] << " expected=" << value
                << " residual=" << actual.max_normalized_residual << '\n';
      throw std::runtime_error("quartic real root mismatch");
    }
    ++i;
  }
  require(actual.max_normalized_residual <= 2.0e-12,
          "quartic normalized residual exceeded gate");
}

std::array<double, 5>
quartic_from_real_roots(const std::array<double, 4> &roots) {
  std::array<long double, 5> coefficients{};
  coefficients[0] = 1.0L;
  std::size_t degree = 0U;
  for (const double root : roots) {
    std::array<long double, 5> next{};
    for (std::size_t power = 0U; power <= degree; ++power) {
      next[power] -= static_cast<long double>(root) * coefficients[power];
      next[power + 1U] += coefficients[power];
    }
    coefficients = next;
    ++degree;
  }
  std::array<double, 5> result{};
  for (std::size_t index = 0U; index < result.size(); ++index)
    result[index] = static_cast<double>(coefficients[index]);
  return result;
}

std::array<double, 5> quartic_with_near_complex_pair(double imaginary,
                                                     double center = 1.0,
                                                     double left = -2.0,
                                                     double right = 3.0) {
  const long double l = static_cast<long double>(left);
  const long double r = static_cast<long double>(right);
  const long double c = static_cast<long double>(center);
  const std::array<long double, 3> first = {l * r, -(l + r), 1.0L};
  const long double e = static_cast<long double>(imaginary);
  const std::array<long double, 3> second = {c * c + e * e, -2.0L * c, 1.0L};
  std::array<double, 5> result{};
  for (std::size_t i = 0U; i < first.size(); ++i) {
    for (std::size_t j = 0U; j < second.size(); ++j)
      result[i + j] += static_cast<double>(first[i] * second[j]);
  }
  return result;
}

std::array<double, 5> symmetric_close_real_quartic(double separation) {
  const double outer = 1.0 + separation;
  return {outer * outer, 0.0, -(1.0 + outer * outer), 0.0, 1.0};
}

void check_near_degenerate_certification() {
  for (const double separation :
       {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8}) {
    const std::array<double, 4> expected = {-2.0, 1.0, 1.0 + separation, 3.0};
    const QuarticRoots roots =
        solve_real_quartic(quartic_from_real_roots(expected));
    if (!roots.certified)
      continue;
    if (roots.count != expected.size())
      std::cerr << "close_real_false_certificate separation=" << separation
                << " count=" << roots.count
                << " residual=" << roots.max_normalized_residual << '\n';
    require(roots.count == expected.size(),
            "certified close-real quartic lost a real root");
    for (std::size_t index = 0U; index < expected.size(); ++index) {
      require(std::abs(roots.values[index] - expected[index]) <=
                  std::max(2.0e-9, 0.05 * separation),
              "certified close-real quartic returned a false root");
    }
  }
  for (const double imaginary :
       {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8}) {
    const QuarticRoots roots =
        solve_real_quartic(quartic_with_near_complex_pair(imaginary));
    if (!roots.certified)
      continue;
    require(roots.count == 2U,
            "certified near-complex quartic invented a real-root pair");
    require(std::abs(roots.values[0] + 2.0) <= 2.0e-9 &&
                std::abs(roots.values[1] - 3.0) <= 2.0e-9,
            "certified near-complex quartic returned a false root");
  }

  // Symmetry makes the shallow extrema especially cancellation-prone. These
  // cases previously returned only +/-sqrt(1+delta) with certified=true for
  // delta <= 1e-8 even though the exact-double quartic has four simple roots.
  for (const double separation : {1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10}) {
    const QuarticRoots roots =
        solve_real_quartic(symmetric_close_real_quartic(separation));
    if (!roots.certified)
      continue;
    const std::array<double, 4> expected = {-(1.0 + separation), -1.0, 1.0,
                                            1.0 + separation};
    require(roots.count == expected.size(),
            "symmetric close-real quartic received a false certificate");
    for (std::size_t index = 0U; index < expected.size(); ++index) {
      require(std::abs(roots.values[index] - expected[index]) <=
                  std::max(2.0e-9, 0.10 * separation),
              "symmetric close-real quartic returned a false root");
    }
  }

  // A repeated root elsewhere must not turn the shallow extremum between two
  // close simple roots into a fake tangent root.  The exact discriminant is a
  // global multiplicity predicate; it does not identify which critical point
  // is shared by P and P'.
  for (int exponent = 20; exponent <= 40; ++exponent) {
    const double separation = std::ldexp(1.0, -exponent);
    const std::array<double, 5> coefficients = {
        1.0 + separation, separation, -2.0 - separation, -separation, 1.0};
    const QuarticRoots roots = solve_real_quartic(coefficients);
    if (!roots.certified)
      continue;
    if (roots.count != 3U)
      std::cerr << "localized_repeated_false_certificate exponent=" << exponent
                << " count=" << roots.count
                << " residual=" << roots.max_normalized_residual << '\n';
    require(roots.count == 3U,
            "global multiplicity certified the wrong critical point");
    const std::array<double, 3> expected = {-1.0, 1.0, 1.0 + separation};
    for (std::size_t root = 0U; root < expected.size(); ++root)
      require(std::abs(roots.values[root] - expected[root]) <=
                  std::max(2.0e-11, 0.10 * separation),
              "localized repeated-root regression returned a false root");
  }

  // Deterministic stress over shifted/scaled near-tangent configurations.  A
  // solver may conservatively return certified=false in the uncertainty band,
  // but it must never attach a certificate to a wrong real-root count.
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  for (std::size_t sample = 0U; sample < 65536U; ++sample) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    const double center =
        -1.0 + 2.0 * static_cast<double>((state >> 16U) & 0xffffU) / 65535.0;
    const double separation = sample % 4U == 0U   ? 1.0e-4
                              : sample % 4U == 1U ? 1.0e-5
                              : sample % 4U == 2U ? 1.0e-6
                                                  : 1.0e-7;
    const double left =
        center - 2.0 -
        0.25 * static_cast<double>((state >> 32U) & 0xffU) / 255.0;
    const double right =
        center + 2.0 +
        0.25 * static_cast<double>((state >> 40U) & 0xffU) / 255.0;
    const std::array<double, 4> expected = {left, center, center + separation,
                                            right};
    const QuarticRoots real =
        solve_real_quartic(quartic_from_real_roots(expected));
    if (real.certified) {
      require(real.count == 4U,
              "stress: certified close-real quartic lost a root");
      for (std::size_t root = 0U; root < 4U; ++root)
        require(std::abs(real.values[root] - expected[root]) <=
                    std::max(2.0e-8, 0.10 * separation),
                "stress: certified close-real quartic invented a root");
    }
    const QuarticRoots complex = solve_real_quartic(
        quartic_with_near_complex_pair(separation, center, left, right));
    if (complex.certified) {
      require(complex.count == 2U,
              "stress: certified near-complex quartic invented roots");
      require(std::abs(complex.values[0] - left) <= 2.0e-8 &&
                  std::abs(complex.values[1] - right) <= 2.0e-8,
              "stress: certified near-complex quartic returned false roots");
    }
  }
}

[[nodiscard]] bool bitwise_same_roots(const QuarticRoots &left,
                                      const QuarticRoots &right) noexcept {
  if (left.count != right.count || left.certified != right.certified ||
      std::bit_cast<std::uint64_t>(left.max_normalized_residual) !=
          std::bit_cast<std::uint64_t>(right.max_normalized_residual))
    return false;
  for (std::size_t index = 0U; index < left.values.size(); ++index) {
    if (std::bit_cast<std::uint64_t>(left.values[index]) !=
        std::bit_cast<std::uint64_t>(right.values[index]))
      return false;
  }
  return true;
}

[[nodiscard]] std::uint64_t root_checksum(const QuarticRoots &roots) noexcept {
  std::uint64_t checksum = static_cast<std::uint64_t>(roots.count) ^
                           (static_cast<std::uint64_t>(roots.certified) << 8U);
  checksum ^= std::bit_cast<std::uint64_t>(roots.max_normalized_residual);
  for (const double root : roots.values) {
    checksum ^= std::bit_cast<std::uint64_t>(root) + 0x9e3779b97f4a7c15ULL +
                (checksum << 6U) + (checksum >> 2U);
  }
  return checksum;
}

void check_huge8192_residual_diagnostic() {
  // Exact binary64 coefficients captured from the authoritative huge8192
  // witness. This is a solver regression, not a reconstructed torus formula.
  const std::array<double, 5> coefficients = {-3.1153557361385979e-05,
                                              -42.069168833371904,
                                              19.817147675844499, 0.0, 1.0};
  const std::array<std::uint64_t, 5> coefficient_bits = {
      0xbf00555c1486006bULL, 0xc04508da863a9e0eULL, 0x4033d130970fc129ULL,
      0x0000000000000000ULL, 0x3ff0000000000000ULL};
  for (std::size_t index = 0U; index < coefficients.size(); ++index)
    require(std::bit_cast<std::uint64_t>(coefficients[index]) ==
                coefficient_bits[index],
            "huge8192 coefficient literal changed binary64 identity");

  const QuarticRoots normal = solve_real_quartic(coefficients);
  const QuarticRoots explicit_null = solve_real_quartic(coefficients, nullptr);
  QuarticDiagnostics diagnostics;
  const QuarticRoots audited = solve_real_quartic(coefficients, &diagnostics);
  require(bitwise_same_roots(normal, explicit_null) &&
              bitwise_same_roots(normal, audited),
          "quartic diagnostics changed public result bits");
  require(normal.count == 2U && normal.certified,
          "huge8192 quartic result/certificate changed");
  require(std::bit_cast<std::uint64_t>(normal.values[0]) ==
                  0xbea8d91e15703802ULL &&
              std::bit_cast<std::uint64_t>(normal.values[1]) ==
                  0x3ffd1b0f13aa87d8ULL,
          "huge8192 math-fix root bits changed");
  require(diagnostics.failure_mask == 0U &&
              diagnostics.first_failure == QuarticDiagnosticFailure::none,
          "huge8192 quartic failure reason changed");
  require(diagnostics.output_root_count == normal.count,
          "diagnostic output root count mismatch");
  require(diagnostics.output_max_normalized_residual <= 2.0e-12,
          "huge8192 residual still exceeds the unchanged gate");
  for (std::size_t index = 0U; index < normal.count; ++index) {
    require(std::bit_cast<std::uint64_t>(diagnostics.output_roots[index]) ==
                std::bit_cast<std::uint64_t>(normal.values[index]),
            "diagnostic root bits differ from public result");
    require(diagnostics.bracketed_roots[index] != 0U &&
                diagnostics.repeated_roots[index] == 0U,
            "huge8192 simple-root provenance changed");
  }
  require(diagnostics.refinement_terminations[0] !=
                  QuarticRefinementTermination::width_floor_one &&
              diagnostics.pre_cast_normalized_residuals[0] <= 2.0e-12L &&
              diagnostics.post_cast_normalized_residuals[0] <= 2.0e-12 &&
              diagnostics.pre_cast_normalized_residuals[1] < 2.0e-12L &&
              diagnostics.post_cast_normalized_residuals[1] < 2.0e-12,
          "huge8192 bad-root/refinement attribution changed");

  constexpr std::size_t kTimingCalls = 20000U;
  const auto time_calls = [&](bool explicit_null_path) {
    std::uint64_t checksum = 0U;
    const auto started = std::chrono::steady_clock::now();
    for (std::size_t iteration = 0U; iteration < kTimingCalls; ++iteration) {
      const QuarticRoots roots = explicit_null_path
                                     ? solve_real_quartic(coefficients, nullptr)
                                     : solve_real_quartic(coefficients);
      checksum ^= root_checksum(roots) + iteration;
    }
    const double elapsed_ms = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - started)
                                  .count();
    return std::pair{elapsed_ms, checksum};
  };
  (void)time_calls(false);
  (void)time_calls(true);
  std::array<double, 7> normal_samples{};
  std::array<double, 7> null_samples{};
  for (std::size_t sample = 0U; sample < normal_samples.size(); ++sample) {
    const auto [normal_ms, normal_checksum] = time_calls(false);
    const auto [null_ms, null_checksum] = time_calls(true);
    require(normal_checksum == null_checksum,
            "null-diagnostic timing path changed result checksum");
    normal_samples[sample] = normal_ms;
    null_samples[sample] = null_ms;
  }
  std::sort(normal_samples.begin(), normal_samples.end());
  std::sort(null_samples.begin(), null_samples.end());
  const double normal_ms = normal_samples[normal_samples.size() / 2U];
  const double null_ms = null_samples[null_samples.size() / 2U];

  std::cout << std::setprecision(17)
            << "huge8192_quartic_root0=" << normal.values[0]
            << " root1=" << normal.values[1] << " root0_bits=0x" << std::hex
            << std::bit_cast<std::uint64_t>(normal.values[0])
            << " root1_bits=0x"
            << std::bit_cast<std::uint64_t>(normal.values[1]) << std::dec
            << " pre0=" << diagnostics.pre_cast_normalized_residuals[0]
            << " post0=" << diagnostics.post_cast_normalized_residuals[0]
            << " pre1=" << diagnostics.pre_cast_normalized_residuals[1]
            << " post1=" << diagnostics.post_cast_normalized_residuals[1]
            << " width0=" << diagnostics.final_bracket_widths[0]
            << " bracket0_lower=" << diagnostics.bracket_lowers[0]
            << " bracket0_upper=" << diagnostics.bracket_uppers[0]
            << " bracket0_lower_value=" << diagnostics.bracket_lower_values[0]
            << " bracket0_upper_value=" << diagnostics.bracket_upper_values[0]
            << " final_bracket0_lower=" << diagnostics.final_bracket_lowers[0]
            << " final_bracket0_upper=" << diagnostics.final_bracket_uppers[0]
            << " final_bracket0_lower_value="
            << diagnostics.final_bracket_lower_values[0]
            << " final_bracket0_upper_value="
            << diagnostics.final_bracket_upper_values[0]
            << " derivative0=" << diagnostics.derivative_magnitudes[0]
            << " condition0=" << diagnostics.root_condition_proxies[0]
            << " termination_residual0="
            << diagnostics.termination_normalized_residuals[0]
            << " termination0="
            << static_cast<unsigned int>(diagnostics.refinement_terminations[0])
            << " width1=" << diagnostics.final_bracket_widths[1]
            << " termination1="
            << static_cast<unsigned int>(diagnostics.refinement_terminations[1])
            << " derivative1=" << diagnostics.derivative_magnitudes[1]
            << " condition1=" << diagnostics.root_condition_proxies[1]
            << " normal_ms_median=" << normal_ms
            << " explicit_null_ms_median=" << null_ms
            << " normal_ms_min=" << normal_samples.front()
            << " normal_ms_max=" << normal_samples.back()
            << " explicit_null_ms_min=" << null_samples.front()
            << " explicit_null_ms_max=" << null_samples.back()
            << " null_ratio=" << (normal_ms > 0.0 ? null_ms / normal_ms : 0.0)
            << '\n';
}

void check_scale_separated_simple_roots() {
  for (const std::array<double, 4> expected :
       {std::array<double, 4>{-8.0, -2.0e-4, 0.75, 3.0},
        std::array<double, 4>{-3.0, -7.0e-7, 1.75, 5.0},
        std::array<double, 4>{-12.0, -3.0e-5, 2.0e-3, 2.5}}) {
    check_roots(solve_real_quartic(quartic_from_real_roots(expected)),
                {expected[0], expected[1], expected[2], expected[3]}, 5.0e-10);
  }
}

} // namespace

int main() {
  try {
    // (x + 4)(x + 2)(x - 1)(x - 3)
    check_roots(solve_real_quartic({24.0, -14.0, -13.0, 2.0, 1.0}),
                {-4.0, -2.0, 1.0, 3.0});
    // (x - 1)^2 (x + 2)^2: two distinct tangent roots.
    check_roots(solve_real_quartic({4.0, -4.0, -3.0, 2.0, 1.0}), {-2.0, 1.0});
    check_roots(solve_real_quartic({1.0, 0.0, -2.0, 0.0, 1.0}), {-1.0, 1.0});
    // Repeated-root regressions: certification must count distinct real roots,
    // not nearby Newton iterates, and it must retain an even-multiplicity root
    // that does not change sign.
    check_roots(solve_real_quartic({1.0, -4.0, 6.0, -4.0, 1.0}), {1.0});
    check_roots(solve_real_quartic({4.0, -12.0, 13.0, -6.0, 1.0}), {1.0, 2.0});
    // (x - 1)^2 (x + 2)(x - 3): one tangent root plus two simple roots.
    check_roots(solve_real_quartic({-6.0, 11.0, -3.0, -3.0, 1.0}),
                {-2.0, 1.0, 3.0});
    // (x - 1)^3 (x + 2): a double root of P' may conservatively leave the
    // derivative partition uncertified, but it must never certify a wrong
    // distinct-root count.
    const QuarticRoots triple =
        solve_real_quartic({-2.0, 5.0, -3.0, -1.0, 1.0});
    if (triple.certified) {
      require(triple.count == 2U,
              "certified triple-root quartic lost a distinct root");
      require(std::abs(triple.values[0] + 2.0) <= 2.0e-11 &&
                  std::abs(triple.values[1] - 1.0) <= 2.0e-11,
              "certified triple-root quartic returned a false root");
    }
    // A repeated complex-conjugate pair consumes all four degrees and must not
    // be exposed as a real tangent root.
    const QuarticRoots repeated_complex =
        solve_real_quartic({1.0, 0.0, 2.0, 0.0, 1.0});
    require(repeated_complex.certified && repeated_complex.count == 0U,
            "repeated complex roots were reported as real roots");
    check_roots(solve_real_quartic({0.0, 0.0, -11.2, 0.0, 1.0}),
                {-std::sqrt(11.2), 0.0, std::sqrt(11.2)});
    // x^4 - 5x^2 + 4.
    check_roots(solve_real_quartic({4.0, 0.0, -5.0, 0.0, 1.0}),
                {-2.0, -1.0, 1.0, 2.0});
    const QuarticRoots none =
        solve_real_quartic({1.0e180, 0.0, 0.0, 0.0, 1.0e180});
    require(none.certified && none.count == 0,
            "positive quartic should have no real roots");
    check_near_degenerate_certification();
    check_huge8192_residual_diagnostic();
    check_scale_separated_simple_roots();

    std::cout << "tabi_ses_quartic_smoke=pass max_residual="
              << solve_real_quartic({24.0, -14.0, -13.0, 2.0, 1.0})
                     .max_normalized_residual
              << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_quartic_smoke=fail reason=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
