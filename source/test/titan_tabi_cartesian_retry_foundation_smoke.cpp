#include "tabi/tabi_cartesian_retry_foundation.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

using TitanPBE::FmmCartesianRetryAttemptDisposition;
using TitanPBE::FmmCartesianRetryAttemptRole;
using TitanPBE::FmmCartesianRetryAttemptTelemetry;
using TitanPBE::FmmCartesianRetryTelemetry;
using TitanPBE::FmmCartesianRetryWrapperAction;
using TitanPBE::Tabi::Result;
using TitanPBE::Tabi::NodePatchGeometry;
using TitanPBE::Tabi::NodePatchTriangle;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::cartesian_retry::fingerprint_cartesian_retry_geometry;
using TitanPBE::Tabi::sanitize_cartesian_retry_failure;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_tabi_cartesian_retry_foundation_smoke: " << message
            << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::shared_ptr<const NodePatchGeometry> make_geometry() {
  auto geometry = std::make_shared<NodePatchGeometry>();
  geometry->vertex_xyz = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                          {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
  geometry->triangles = {
      NodePatchTriangle{{0, 1, 2}, 0.0, 0.0, 1.0, 0.5, 3},
      NodePatchTriangle{{0, 2, 3}, 0.0, 0.0, 1.0, 0.5, 4},
  };
  geometry->incident_offsets = {0, 2, 3, 5, 6};
  geometry->incident_triangles = {0, 1, 0, 0, 1, 1};
  geometry->adjacent_offsets = {0, 2, 4, 6, 8};
  geometry->adjacent_nodes = {1, 3, 0, 2, 1, 3, 0, 2};
  geometry->solid_angle_fallbacks = 2;
  return geometry;
}

std::vector<Panel>
make_panels(const std::shared_ptr<const NodePatchGeometry> &geometry) {
  std::vector<Panel> panels;
  panels.reserve(4);
  for (std::size_t index = 0; index < 4; ++index) {
    Panel panel;
    const double x = static_cast<double>(index % 2);
    const double y = static_cast<double>(index / 2);
    panel.x = x;
    panel.y = y;
    panel.z = 0.05 * static_cast<double>(index);
    panel.source_x = x + 0.015;
    panel.source_y = y - 0.010;
    panel.source_z = panel.z + 0.020;
    panel.has_source_quadrature_center = true;
    panel.nx = 0.0;
    panel.ny = 0.0;
    panel.nz = 1.0;
    panel.area = 0.25;
    panel.diagonal = 0.75 + 0.01 * static_cast<double>(index);
    panel.solid_angle = 0.5;
    panel.atom_index = 10 + index;
    panel.v0x = x - 0.1;
    panel.v0y = y - 0.1;
    panel.v0z = panel.z;
    panel.v1x = x + 0.1;
    panel.v1y = y - 0.1;
    panel.v1z = panel.z;
    panel.v2x = x;
    panel.v2y = y + 0.2;
    panel.v2z = panel.z;
    panel.v0_index = 3 * index;
    panel.v1_index = 3 * index + 1;
    panel.v2_index = 3 * index + 2;
    panel.has_triangle_geometry = true;
    panel.node_patch_geometry = geometry;
    panel.node_patch_index = index;
    panels.push_back(panel);
  }
  return panels;
}

void expect_hash_changes(std::vector<Panel> panels,
                         const std::uint64_t baseline,
                         const std::string &label) {
  const auto fingerprint = fingerprint_cartesian_retry_geometry(panels);
  expect(fingerprint.valid, label + " fingerprint unexpectedly invalid");
  expect(fingerprint.value != baseline, label + " did not change the hash");
}

FmmCartesianRetryAttemptTelemetry make_valid_attempt(double theta) {
  FmmCartesianRetryAttemptTelemetry attempt;
  attempt.attempt = theta < 0.18 ? 2 : 1;
  attempt.started = true;
  attempt.completed = true;
  attempt.status = "ok";
  attempt.converged = true;
  attempt.finite = true;
  attempt.theta = theta;
  attempt.salt_molar = 0.15;
  attempt.coupled = true;
  attempt.discretization = "node_patch";
  attempt.node_patch_operator_audit_enabled = true;
  attempt.true_residual_rel = 2.0e-6;
  attempt.cartesian_full_active = true;
  attempt.cartesian_full_finite = true;
  attempt.cartesian_full_requested = "cartesian-full-experimental";
  attempt.cartesian_full_effective = "cartesian-full-experimental";
  attempt.cartesian_full_order = 1;
  attempt.cartesian_full_leaf_size = 64;
  attempt.panel_count = 4;
  attempt.geometry_hash = 101;
  attempt.tree_hash = 202;
  attempt.worklist_hash = 303;
  attempt.policy_hash = theta < 0.18 ? 505 : 404;
  attempt.direct_pairs_expected = 12;
  attempt.direct_pairs_applied_last = 12;
  attempt.direct_pairs_applied_total = 36;
  attempt.pre_audit_sampled_target_hash = 606;
  attempt.solution_audit.enabled = true;
  attempt.solution_audit.finite = true;
  attempt.solution_audit.relative_available = true;
  attempt.solution_audit.status = "relative";
  attempt.solution_audit.relative_error = 5.0e-5;
  attempt.solution_audit.absolute_l2 = 2.5e-4;
  attempt.solution_audit.direct_l2 = 5.0;
  attempt.solution_audit.seconds = 0.025;
  attempt.solution_audit.targets = 3;
  attempt.solution_audit.rows = 6;
  attempt.solution_audit.target_hash = 606;
  attempt.solution_audit.direct_pairs_expected = 12;
  attempt.solution_audit.direct_pairs_applied = 12;
  attempt.preconditioner_requested = "fmm-nearfield-ilu0";
  attempt.preconditioner_effective = "fmm_nearfield_ilu0";
  attempt.preconditioner_quality_status = "pass";
  attempt.preconditioner_coarse_status = "off";
  attempt.preconditioner_ilu_regularization = "adaptive";
  attempt.preconditioner_ilu_regularization_selected = 0.0;
  attempt.preconditioner_ilu_regularization_shifted_rows = 0;
  attempt.preconditioner_ilu_unshifted_failure_stage = "off";
  attempt.preconditioner_valid_panels = 4;
  return attempt;
}

void test_geometry_fingerprint() {
  const auto first_geometry = make_geometry();
  const auto first = make_panels(first_geometry);
  const auto baseline = fingerprint_cartesian_retry_geometry(first);
  expect(baseline.valid && baseline.value != 0,
         "baseline geometry must have a nonzero fingerprint");

  // Independent shared_ptr ownership is deliberately invisible to the value
  // hash when the operator-relevant geometry content is equal.
  const auto equal_content = make_panels(make_geometry());
  const auto equal_hash = fingerprint_cartesian_retry_geometry(equal_content);
  expect(equal_hash.valid && equal_hash.value == baseline.value,
         "content-equal geometry with a different shared_ptr changed hash");

  auto source_mutation = first;
  source_mutation[0].source_x += 0.125;
  expect_hash_changes(std::move(source_mutation), baseline.value,
                      "effective source center");
  auto normal_mutation = first;
  normal_mutation[0].ny = 0.125;
  expect_hash_changes(std::move(normal_mutation), baseline.value, "normal");
  auto area_mutation = first;
  area_mutation[0].area *= 1.25;
  expect_hash_changes(std::move(area_mutation), baseline.value, "area");
  auto diagonal_mutation = first;
  diagonal_mutation[0].diagonal += 0.125;
  expect_hash_changes(std::move(diagonal_mutation), baseline.value,
                      "diagonal");
  auto angle_mutation = first;
  angle_mutation[0].solid_angle += 0.125;
  expect_hash_changes(std::move(angle_mutation), baseline.value,
                      "solid angle");
  auto index_mutation = first;
  index_mutation[0].node_patch_index = 1;
  expect_hash_changes(std::move(index_mutation), baseline.value,
                      "node-patch index");
  auto triangle_mutation = first;
  triangle_mutation[0].v2y += 0.125;
  expect_hash_changes(std::move(triangle_mutation), baseline.value,
                      "triangle field");

  auto topology = std::make_shared<NodePatchGeometry>(*first_geometry);
  topology->adjacent_nodes[0] = 2;
  auto topology_mutation = make_panels(topology);
  expect_hash_changes(std::move(topology_mutation), baseline.value,
                      "node-patch adjacency");

  auto mixed_ownership = first;
  mixed_ownership[1].node_patch_geometry = make_geometry();
  expect(!fingerprint_cartesian_retry_geometry(mixed_ownership).valid,
         "mixed node-patch ownership must fail closed");
  auto missing_topology = first;
  for (Panel &panel : missing_topology)
    panel.node_patch_geometry.reset();
  expect(!fingerprint_cartesian_retry_geometry(missing_topology).valid,
         "missing node-patch topology must fail closed");
}

void test_retry_classifier() {
  using TitanPBE::classify_fmm_cartesian_retry_attempt;
  using TitanPBE::classify_fmm_cartesian_retry_identity;
  using TitanPBE::decide_fmm_cartesian_retry_after_attempt1;
  using TitanPBE::decide_fmm_cartesian_retry_after_attempt2;

  const auto first = make_valid_attempt(0.20);
  const auto accepted = classify_fmm_cartesian_retry_attempt(
      first, FmmCartesianRetryAttemptRole::initial_theta_020);
  expect(accepted.disposition == FmmCartesianRetryAttemptDisposition::accept,
         "valid initial attempt must be accepted");

  auto retry = first;
  retry.solution_audit.relative_error = 2.0e-4;
  const auto retry_result = classify_fmm_cartesian_retry_attempt(
      retry, FmmCartesianRetryAttemptRole::initial_theta_020);
  expect(retry_result.disposition ==
             FmmCartesianRetryAttemptDisposition::retry_trigger &&
             retry_result.reason == "solution_audit_above_gate",
         "only a technically valid high solution audit may trigger retry");
  const auto first_retry_decision =
      decide_fmm_cartesian_retry_after_attempt1(retry_result);
  expect(first_retry_decision.action ==
             FmmCartesianRetryWrapperAction::run_attempt2 &&
             first_retry_decision.accepted_attempt == 0 &&
             !first_retry_decision.sanitize_result,
         "only a valid initial retry trigger may start attempt 2");
  const auto tightened_high = classify_fmm_cartesian_retry_attempt(
      retry, FmmCartesianRetryAttemptRole::tightened_theta_015);
  expect(tightened_high.disposition ==
             FmmCartesianRetryAttemptDisposition::fail_closed &&
             tightened_high.reason == "theta_mismatch",
         "role-specific theta must be enforced before a tightened decision");

  auto second = make_valid_attempt(0.15);
  second.solution_audit.relative_error = 2.0e-4;
  const auto tightened_result = classify_fmm_cartesian_retry_attempt(
      second, FmmCartesianRetryAttemptRole::tightened_theta_015);
  expect(tightened_result.disposition ==
             FmmCartesianRetryAttemptDisposition::fail_closed &&
             tightened_result.reason == "solution_audit_above_gate_after_tighten",
         "tightened high audit must fail closed rather than request a third run");

  auto regularized = first;
  regularized.preconditioner_effective = "fmm_nearfield_ilu0_regularized";
  expect(classify_fmm_cartesian_retry_attempt(
             regularized, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "plain_fmm_nearfield_ilu0_required",
         "regularized ILU0 must be outside the calibrated retry scope");
  auto shifted = first;
  shifted.preconditioner_ilu_regularization_shifted_rows = 1;
  expect(classify_fmm_cartesian_retry_attempt(
             shifted, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "plain_fmm_nearfield_ilu0_required",
         "an applied ILU shift must be outside the calibrated retry scope");
  auto bad_rows = first;
  bad_rows.solution_audit.rows = 5;
  expect(classify_fmm_cartesian_retry_attempt(
             bad_rows, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "solution audit rows must exactly match two rows per target");
  auto bad_target_hash = first;
  bad_target_hash.solution_audit.target_hash = 607;
  expect(classify_fmm_cartesian_retry_attempt(
             bad_target_hash, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "solution audit target hash must match the pre-audit target hash");
  auto negative_residual = first;
  negative_residual.true_residual_rel = -1.0e-6;
  expect(classify_fmm_cartesian_retry_attempt(
             negative_residual, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "true_residual_gate_failed",
         "negative true residual must fail closed");
  auto negative_relative_error = first;
  negative_relative_error.solution_audit.relative_error = -1.0e-6;
  expect(classify_fmm_cartesian_retry_attempt(
             negative_relative_error,
             FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "negative solution relative error must fail closed");
  auto negative_absolute_error = first;
  negative_absolute_error.solution_audit.absolute_l2 = -1.0e-6;
  expect(classify_fmm_cartesian_retry_attempt(
             negative_absolute_error,
             FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "negative solution absolute error must fail closed");
  auto zero_direct_norm = first;
  zero_direct_norm.solution_audit.direct_l2 = 0.0;
  expect(classify_fmm_cartesian_retry_attempt(
             zero_direct_norm, FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "zero relative reference norm must fail closed");
  auto nonfinite_direct_norm = first;
  nonfinite_direct_norm.solution_audit.direct_l2 =
      std::numeric_limits<double>::quiet_NaN();
  expect(classify_fmm_cartesian_retry_attempt(
             nonfinite_direct_norm,
             FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "nonfinite solution reference norm must fail closed");
  auto invalid_audit_seconds = first;
  invalid_audit_seconds.solution_audit.seconds =
      std::numeric_limits<double>::quiet_NaN();
  expect(classify_fmm_cartesian_retry_attempt(
             invalid_audit_seconds,
             FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "nonfinite solution audit seconds must fail closed");
  auto negative_audit_seconds = first;
  negative_audit_seconds.solution_audit.seconds = -1.0e-3;
  expect(classify_fmm_cartesian_retry_attempt(
             negative_audit_seconds,
             FmmCartesianRetryAttemptRole::initial_theta_020)
                 .reason == "solution_audit_unavailable_or_invalid",
         "negative solution audit seconds must fail closed");

  const auto identity = classify_fmm_cartesian_retry_identity(first, second);
  expect(identity.valid, "matched attempts must pass identity classification");
  auto second_invalid = second;
  second_invalid.solution_audit.finite = false;
  const auto second_invalid_result = classify_fmm_cartesian_retry_attempt(
      second_invalid, FmmCartesianRetryAttemptRole::tightened_theta_015);
  expect(second_invalid_result.disposition ==
             FmmCartesianRetryAttemptDisposition::fail_closed,
         "an invalid attempt 2 must fail closed");
  const auto second_fail_closed = decide_fmm_cartesian_retry_after_attempt2(
      second_invalid_result, identity);
  expect(second_fail_closed.action ==
             FmmCartesianRetryWrapperAction::fail_closed &&
             second_fail_closed.accepted_attempt == 0 &&
             second_fail_closed.sanitize_result,
         "an A2 fail-closed result cannot request A3 or acceptance");
  auto changed_geometry = second;
  changed_geometry.geometry_hash = 102;
  expect(classify_fmm_cartesian_retry_identity(first, changed_geometry).reason ==
             "geometry_hash_mismatch",
         "geometry mismatch must fail closed");
  auto same_policy = second;
  same_policy.policy_hash = first.policy_hash;
  expect(classify_fmm_cartesian_retry_identity(first, same_policy).reason ==
             "policy_hash_transition_invalid",
         "theta transition must change policy hash");
  auto missing_worklist = second;
  missing_worklist.worklist_hash = 0;
  expect(classify_fmm_cartesian_retry_identity(first, missing_worklist).reason ==
             "worklist_hash_missing",
         "each retry attempt must record a worklist hash");
}

void test_retry_failure_sanitization() {
  FmmCartesianRetryTelemetry retry;
  retry.mode = "channel-shadow-retry-global-theta015";
  retry.enabled = true;
  retry.attempts_started = 2;
  retry.attempts_completed = 2;
  retry.accepted_attempt = 2;
  retry.final_reason = "solution_audit_unavailable_or_invalid";
  for (const char *origin : {"attempt1", "attempt2"}) {
    Result result;
    result.status = "ok";
    result.converged = true;
    result.reaction_energy_source = origin;
    result.energy_kcal_mol = 11.0;
    result.energy_kj_mol = 22.0;
    result.energy_kt = 33.0;
    result.boundary_energy_kcal_mol = 44.0;
    result.reference_energy_kcal_mol = 55.0;
    result.selected_energy_scale = 66.0;
    result.atom_reaction_energy_closure_kj_mol = 77.0;
    result.atom_reaction_energy_closure_error_kj_mol = 88.0;
    result.atom_reaction_potential_kj_mol_e = {1.0};
    result.atom_reaction_energy_kj_mol = {2.0};
    result.output_panels.emplace_back();
    result.output_solution = {3.0};
    result.output_rhs = {4.0};
    const Result sanitized =
        sanitize_cartesian_retry_failure(std::move(result), retry);
    expect(sanitized.status == "failed" && !sanitized.converged &&
               sanitized.fmm_cartesian_retry.accepted_attempt == 0 &&
               sanitized.reaction_energy_source ==
                   "not_computed_cartesian_retry_failed",
           std::string(origin) + " failure did not sanitize its result state");
    expect(!std::isfinite(sanitized.energy_kcal_mol) &&
               !std::isfinite(sanitized.energy_kj_mol) &&
               !std::isfinite(sanitized.energy_kt) &&
               !std::isfinite(sanitized.boundary_energy_kcal_mol) &&
               !std::isfinite(sanitized.reference_energy_kcal_mol) &&
               !std::isfinite(sanitized.selected_energy_scale) &&
               !std::isfinite(sanitized.atom_reaction_energy_closure_kj_mol) &&
               !std::isfinite(
                   sanitized.atom_reaction_energy_closure_error_kj_mol) &&
               sanitized.atom_reaction_potential_kj_mol_e.empty() &&
               sanitized.atom_reaction_energy_kj_mol.empty() &&
               sanitized.output_panels.empty() &&
               sanitized.output_solution.empty() && sanitized.output_rhs.empty(),
           std::string(origin) +
               " failure leaked an unaccepted energy or output vector");
  }
}

} // namespace

int main() {
  test_geometry_fingerprint();
  test_retry_classifier();
  test_retry_failure_sanitization();
  std::cout << "PASS titan_tabi_cartesian_retry_foundation_smoke\n";
  return 0;
}
