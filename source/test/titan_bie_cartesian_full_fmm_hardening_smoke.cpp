#include "bie/bie_cartesian_fmm_solver_adapter.h"
#include "bie/fmm/cartesian/bie_cartesian_fmm_checked_size.h"
#include "tabi/tabi_fmm_m2l_fallback.h"
#include "tabi/tabi_near.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if !defined(TITAN_CARTESIAN_FMM_TESTING)
#error "This unit smoke requires TITAN_CARTESIAN_FMM_TESTING."
#endif

namespace {

using titan::bie::CartesianFullFmmPrepareOptions;
using titan::bie::CartesianFullFmmSolverOperator;
using TitanPBE::Tabi::Panel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_bie_cartesian_full_fmm_hardening_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

template <typename Fn>
void expect_throws(Fn &&fn, std::initializer_list<std::string_view> needles,
                   const std::string &context) {
  try {
    fn();
  } catch (const std::exception &error) {
    const std::string diagnostic = error.what();
    for (const std::string_view needle : needles) {
      if (diagnostic.find(needle) == std::string::npos) {
        fail(context + " diagnostic is missing '" + std::string(needle) +
             "': " + diagnostic);
      }
    }
    return;
  }
  fail(context + " did not fail closed");
}

Panel make_panel(std::size_t index) {
  Panel panel;
  const double x = 2.25 * static_cast<double>(index % 4);
  const double y = 2.10 * static_cast<double>(index / 4);
  panel.x = x;
  panel.y = y;
  panel.z = 0.15 * std::sin(static_cast<double>(index));
  panel.nx = 0.0;
  panel.ny = 0.0;
  panel.nz = 1.0;
  panel.area = 0.16;
  panel.diagonal = 1.0;
  panel.atom_index = index;
  panel.has_triangle_geometry = true;
  panel.v0x = x - 0.20;
  panel.v0y = y - 0.12;
  panel.v0z = panel.z;
  panel.v1x = x + 0.20;
  panel.v1y = y - 0.12;
  panel.v1z = panel.z;
  panel.v2x = x;
  panel.v2y = y + 0.28;
  panel.v2z = panel.z;
  panel.v0_index = 3 * index;
  panel.v1_index = 3 * index + 1;
  panel.v2_index = 3 * index + 2;
  return panel;
}

std::vector<Panel> make_panels() {
  std::vector<Panel> panels;
  panels.reserve(12);
  for (std::size_t index = 0; index < 12; ++index)
    panels.push_back(make_panel(index));
  return panels;
}

CartesianFullFmmPrepareOptions base_options() {
  CartesianFullFmmPrepareOptions options;
  options.expansion_order = 1;
  options.leaf_size = 4;
  options.max_depth = 5;
  options.theta = 0.20;
  options.eps_ratio = 80.0;
  options.kappa = 0.12;
  options.near_integration = "analytic_singular";
  options.quad_buffer = "off";
  return options;
}

std::unique_ptr<CartesianFullFmmSolverOperator>
prepare_or_fail(const std::vector<Panel> &panels,
                const CartesianFullFmmPrepareOptions &options,
                TitanExecutor &executor, const std::string &context) {
  auto prepared = CartesianFullFmmSolverOperator::prepare(panels, options,
                                                            executor);
  if (!prepared)
    fail(context + " preparation failed: " + prepared.error());
  return std::move(*prepared);
}

std::vector<double> make_unknowns(std::size_t panel_count) {
  std::vector<double> values(2 * panel_count, 0.0);
  for (std::size_t index = 0; index < panel_count; ++index) {
    const double t = static_cast<double>(index + 1);
    values[index] = 0.25 + std::sin(0.29 * t);
    values[panel_count + index] = -0.10 + std::cos(0.37 * t);
  }
  return values;
}

void expect_complete_coverage(const CartesianFullFmmSolverOperator &op,
                              const std::string &context) {
  const auto &stats = op.stats();
  expect(stats.last_apply.expected_direct_panel_pairs > 0,
         context + " must have direct panel pairs");
  expect(stats.last_apply.direct_panel_pairs ==
             stats.last_apply.expected_direct_panel_pairs,
         context + " applied direct coverage differs from the plan");
  expect(stats.total_apply.direct_panel_pairs ==
             stats.total_apply.expected_direct_panel_pairs,
         context + " total direct coverage differs from expected total");
}

void expect_finite_output(const std::vector<double> &values,
                          const std::string &context) {
  expect(std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); }),
         context + " output must be finite");
}

void test_checked_sizes() {
  namespace checked = titan::bie::fmm::cartesian::detail;
  std::size_t result = 0;
  expect(!checked::checked_add(std::numeric_limits<std::size_t>::max(),
                               std::size_t{1}, result),
         "checked size addition must reject overflow");
  expect(!checked::checked_multiply(std::numeric_limits<std::size_t>::max(),
                                    std::size_t{2}, result),
         "checked size multiplication must reject overflow");
  expect(checked::checked_multiply(std::size_t{17}, std::size_t{19}, result) &&
             result == 323,
         "checked size multiplication must preserve valid products");
}

void test_zero_self_contract() {
  TitanPBE::Tabi::JcpNearPairKernel zero_self;
  zero_self.self_pair = true;
  expect(TitanPBE::Tabi::is_intentional_zero_self_kernel(
             zero_self, false, "analytic_singular"),
         "uncoupled self invalid must be an intentional zero");
  expect(TitanPBE::Tabi::is_intentional_zero_self_kernel(
             zero_self, true, "centroid"),
         "centroid self invalid must be an intentional zero");
  expect(!TitanPBE::Tabi::is_intentional_zero_self_kernel(
             zero_self, true, "analytic_singular"),
         "coupled analytic self invalid must fail closed");

  TitanPBE::Tabi::JcpNearPairKernel invalid_nonself;
  expect(!TitanPBE::Tabi::is_intentional_zero_self_kernel(
             invalid_nonself, false, "analytic_singular"),
         "nonself invalid must never be accepted as a zero-self term");

  TitanPBE::Tabi::JcpNearPairKernel valid_self;
  valid_self.self_pair = true;
  valid_self.valid = true;
  expect(!TitanPBE::Tabi::is_intentional_zero_self_kernel(
             valid_self, false, "analytic_singular"),
         "valid self kernel must remain a numerical contribution");
}

} // namespace

int main() {
  test_checked_sizes();
  test_zero_self_contract();
  const auto panels = make_panels();
  const auto unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);
  std::vector<double> output;

  auto scalar_options = base_options();
  scalar_options.force_scalar_direct_for_testing = true;
  auto scalar = prepare_or_fail(panels, scalar_options, executor, "scalar");
  scalar->apply(unknowns, 0.5, 0.5, true, executor, output);
  expect(scalar->stats().direct_kernel == "scalar",
         "forced scalar test must not select packed direct evaluation");
  expect_complete_coverage(*scalar, "scalar");
  expect_finite_output(output, "scalar");

  auto uncoupled_options = scalar_options;
  auto uncoupled = prepare_or_fail(panels, uncoupled_options, executor,
                                   "uncoupled scalar zero-self");
  uncoupled->apply(unknowns, 0.5, 0.5, false, executor, output);
  expect(uncoupled->stats().direct_kernel == "scalar",
         "uncoupled zero-self test must use scalar direct evaluation");
  expect_complete_coverage(*uncoupled, "uncoupled scalar zero-self");
  expect_finite_output(output, "uncoupled scalar zero-self");

  auto centroid_options = scalar_options;
  centroid_options.near_integration = "centroid";
  auto centroid = prepare_or_fail(panels, centroid_options, executor,
                                  "centroid scalar zero-self");
  centroid->apply(unknowns, 0.5, 0.5, true, executor, output);
  expect(centroid->stats().direct_kernel == "scalar",
         "centroid zero-self test must use scalar direct evaluation");
  expect_complete_coverage(*centroid, "centroid scalar zero-self");
  expect_finite_output(output, "centroid scalar zero-self");

  scalar->hold_apply_guard_for_testing();
  expect_throws(
      [&] { scalar->apply(unknowns, 0.5, 0.5, true, executor, output); },
      {"concurrent or reentrant apply"}, "apply guard");
  scalar->release_apply_guard_for_testing();
  scalar->apply(unknowns, 0.5, 0.5, true, executor, output);
  expect_complete_coverage(*scalar, "scalar after guard release");

  auto invalid_options = scalar_options;
  invalid_options.inject_invalid_scalar_direct_pair_for_testing = true;
  auto invalid = prepare_or_fail(panels, invalid_options, executor,
                                 "invalid scalar injection");
  expect_throws(
      [&] { invalid->apply(unknowns, 0.5, 0.5, true, executor, output); },
      {"scalar direct kernel invalid", "target_panel=", "source_panel=",
       "target_node=", "source_node=", "near_integration=",
       "quad_buffer=", "forced_test=true"},
      "scalar kernel.valid=false");
  // The apply guard must be released by the exception path before a caller
  // can acquire it again.
  invalid->hold_apply_guard_for_testing();
  invalid->release_apply_guard_for_testing();

  auto coverage_options = scalar_options;
  coverage_options.inject_direct_coverage_mismatch_for_testing = true;
  auto coverage = prepare_or_fail(panels, coverage_options, executor,
                                  "coverage injection");
  expect_throws(
      [&] { coverage->apply(unknowns, 0.5, 0.5, true, executor, output); },
      {"direct coverage mismatch", "expected_direct_panel_pairs=",
       "applied_direct_panel_pairs=", "worklist_hash=",
       "node_patch_correction_pairs="},
      "direct coverage invariant");

  if (titan::tabi::fmm_m2l::packed_jcp_fallback_avx2_available()) {
    auto packed = prepare_or_fail(panels, base_options(), executor, "packed");
    packed->apply(unknowns, 0.5, 0.5, true, executor, output);
    expect(packed->stats().direct_kernel.find("packed-avx2") == 0,
           "available AVX2 path must select packed direct evaluation");
    expect_complete_coverage(*packed, "packed");
  } else {
    std::cout << "packed_direct=unavailable\n";
  }

  std::cout << "titan_bie_cartesian_full_fmm_hardening_smoke ok\n";
  return 0;
}
