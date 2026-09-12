param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_component_search_optimality_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/component_search.h"
#include "pbctopo/certification.h"
#include "pbctopo/reference_strain.h"
#include "pbctopo/steric_clash.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>
#include <unordered_set>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_component_search_optimality_smoke: "
              << message << "\n";
    std::exit(1);
  }
}

PbctopoComponentCandidate build_seed(
    const PbctopoCutCandidate &,
    const std::vector<std::vector<std::size_t>> &,
    const std::vector<Vec3> &) {
  return {};
}

PbctopoComponentCandidate build_incumbent_seed(
    const PbctopoCutCandidate &cut,
    const std::vector<std::vector<std::size_t>> &,
    const std::vector<Vec3> &) {
  PbctopoComponentCandidate candidate;
  candidate.valid = cut.valid;
  candidate.ultra_clash = cut.ultra_clash;
  candidate.shift_norm2 = cut.shift_norm2;
  candidate.continuity_path = cut.continuity_path;
  candidate.continuity_max_d2 = cut.continuity_max_d2;
  return candidate;
}

double target_x1 = 10.2;
double target_x2 = -5.0;

double objective(const std::vector<Vec3> &points) {
  if (points.size() < 3)
    return std::numeric_limits<double>::infinity();
  const double e1 = points[1].x - target_x1;
  const double e2 = points[2].x - target_x2;
  return e1 * e1 + e2 * e2;
}

PbctopoComponentCandidate score_assignment(
    const std::vector<PbctopoAtom> &, const std::vector<Vec3> &points,
    const std::vector<std::vector<std::size_t>> &, double, std::size_t,
    std::size_t, const std::vector<PbctopoContinuityPair> &,
    const std::vector<std::size_t> &, double, double,
    const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
    const PbctopoAnchorGraphContext *) {
  PbctopoComponentCandidate candidate;
  candidate.valid = true;
  candidate.shift_norm2 = objective(points);
  candidate.continuity_path = 0.0;
  candidate.continuity_max_d2 = 0.0;
  return candidate;
}

PbctopoComponentCandidate score_tied_assignment(
    const std::vector<PbctopoAtom> &, const std::vector<Vec3> &,
    const std::vector<std::vector<std::size_t>> &, double, std::size_t,
    std::size_t, const std::vector<PbctopoContinuityPair> &,
    const std::vector<std::size_t> &, double, double,
    const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
    const PbctopoAnchorGraphContext *) {
  PbctopoComponentCandidate candidate;
  candidate.valid = true;
  candidate.shift_norm2 = 0.0;
  candidate.continuity_path = 0.0;
  candidate.continuity_max_d2 = 0.0;
  return candidate;
}

PbctopoComponentCandidate score_selectively_invalid_assignment(
    const std::vector<PbctopoAtom> &, const std::vector<Vec3> &points,
    const std::vector<std::vector<std::size_t>> &, double, std::size_t,
    std::size_t, const std::vector<PbctopoContinuityPair> &,
    const std::vector<std::size_t> &, double, double,
    const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
    const PbctopoAnchorGraphContext *) {
  PbctopoComponentCandidate candidate;
  candidate.valid = points.size() == 3 && points[1].x > 0.0;
  if (candidate.valid) {
    candidate.shift_norm2 = objective(points);
    candidate.continuity_path = 0.0;
    candidate.continuity_max_d2 = 0.0;
  }
  return candidate;
}

bool better_candidate(const PbctopoComponentCandidate &lhs,
                      const PbctopoComponentCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  return lhs.shift_norm2 + 1.0e-12 < rhs.shift_norm2;
}

PbctopoComponentCandidate score_continuity_assignment(
    const std::vector<PbctopoAtom> &, const std::vector<Vec3> &points,
    const std::vector<std::vector<std::size_t>> &, double, std::size_t,
    std::size_t,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &, double, double shift_norm2,
    const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
    const PbctopoAnchorGraphContext *) {
  PbctopoComponentCandidate candidate;
  candidate.valid = true;
  candidate.shift_norm2 = shift_norm2;
  candidate.continuity_path = 0.0;
  candidate.continuity_max_d2 = 0.0;
  for (const auto &pair : continuity_pairs) {
    const Vec3 delta = points[pair.a] - points[pair.b];
    const double d2 =
        delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
    candidate.continuity_path +=
        pair.weight * pbctopo_symmetric_reference_strain(d2, pair.ref_d2);
  }
  return candidate;
}

bool better_continuity_candidate(const PbctopoComponentCandidate &lhs,
                                 const PbctopoComponentCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.continuity_path + 1.0e-12 < rhs.continuity_path)
    return true;
  if (rhs.continuity_path + 1.0e-12 < lhs.continuity_path)
    return false;
  return lhs.shift_norm2 + 1.0e-12 < rhs.shift_norm2;
}

PbctopoComponentCandidate score_hydrogen_band_assignment(
    const std::vector<PbctopoAtom> &, const std::vector<Vec3> &points,
    const std::vector<std::vector<std::size_t>> &, double, std::size_t,
    std::size_t, const std::vector<PbctopoContinuityPair> &,
    const std::vector<std::size_t> &, double, double,
    const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
    const PbctopoAnchorGraphContext *) {
  PbctopoComponentCandidate candidate;
  candidate.valid = points.size() == 2;
  if (!candidate.valid)
    return candidate;
  const double displacement = points[1].x - 3.0;
  candidate.shift_norm2 = displacement * displacement;
  candidate.continuity_path = 0.0;
  candidate.continuity_max_d2 = 0.0;
  return candidate;
}

std::vector<PbctopoShift> make_shifts(double length) {
  std::vector<PbctopoShift> shifts;
  for (int image : {-1, 0, 1}) {
    PbctopoShift shift;
    shift.shift_index = shifts.size();
    shift.ix = image;
    shift.tx = image * length;
    shifts.push_back(shift);
  }
  return shifts;
}

bool assignment_clashes(const std::vector<PbctopoAtom> &atoms,
                        const std::vector<double> &x,
                        double fallback_cutoff) {
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    for (std::size_t b = a + 1; b < atoms.size(); ++b) {
      const double dx = x[a] - x[b];
      if (pbctopo_atoms_clash(dx * dx, atoms[a], atoms[b], fallback_cutoff))
        return true;
    }
  }
  return false;
}

int main() {
  box Box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
  RuntimeCutoffs cutoffs;
  cutoffs.clash_cutoff = 1.05;
  cutoffs.pair_contact_cutoff = 4.5;
  cutoffs.continuity_cutoff = 20.0;

  std::vector<PbctopoAtom> atoms(3);
  atoms[0].x = 0.0;
  atoms[1].x = 0.2;
  atoms[2].x = 5.0;
  for (auto &atom : atoms) {
    atom.is_heavy = true;
    atom.has_vdw_radius = true;
    atom.vdw_radius = 1.0;
    atom.steric_radius = 0.6;
  }
  const std::vector<std::vector<std::size_t>> components{{0}, {1}, {2}};
  const auto shifts = make_shifts(10.0);
  const PbctopoComponentSearchCallbacks callbacks{
      build_seed, score_assignment, better_candidate, true, true};

  double brute_best = std::numeric_limits<double>::infinity();
  for (std::size_t s1 = 0; s1 < shifts.size(); ++s1) {
    for (std::size_t s2 = 0; s2 < shifts.size(); ++s2) {
      const std::vector<double> x{atoms[0].x,
                                  atoms[1].x + shifts[s1].tx,
                                  atoms[2].x + shifts[s2].tx};
      if (assignment_clashes(atoms, x, cutoffs.clash_cutoff))
        continue;
      const std::vector<Vec3> points{{x[0], 0.0, 0.0},
                                     {x[1], 0.0, 0.0},
                                     {x[2], 0.0, 0.0}};
      brute_best = std::min(brute_best, objective(points));
    }
  }
  expect(std::fabs(brute_best) < 1.0e-12,
         "brute-force fixture has the wrong optimum");

  std::vector<Vec3> exact_points;
  PbctopoTimingStats exact_timing;
  PbctopoComponentSearchDiagnostics exact_diagnostics;
  std::vector<PbctopoComponentSearchSurvivor> exact_survivors;
  const auto exact = search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, exact_points,
      callbacks, nullptr, nullptr, nullptr, &exact_timing,
      AnchorBackend::Auto, nullptr, &exact_diagnostics, &exact_survivors, 8);
  expect(exact.valid && std::fabs(objective(exact_points) - brute_best) < 1.0e-12,
         "exact search differs from brute force");
  expect(exact.search_domain_audit_known &&
             exact.component_search_shell_radius == 1 &&
             exact.selected_assignment_touches_shell_boundary &&
             exact_diagnostics.selected_assignment_boundary_known &&
             exact_diagnostics.selected_assignment_touches_shell_boundary,
         "exact winner did not expose its shell-boundary witness");
  expect(exact_timing.component_pair_table_relative_deltas > 0,
          "exact search did not build the relative pair table");
  expect(!exact_diagnostics.exhaustive &&
             exact_diagnostics.search_completed_under_admissible_bounds &&
             exact_diagnostics.best_objective_proven &&
             exact_diagnostics.equivalent_optima_complete &&
             !exact_diagnostics.feasible_assignment_set_enumerated &&
             exact_diagnostics.feasible_assignments_within_domain == 0 &&
             exact_diagnostics.bounded_domain_feasibility_known &&
             exact_diagnostics.bounded_domain_feasible &&
             exact_diagnostics.optimality_gap_known &&
             exact_diagnostics.optimality_gap == 0.0,
          "normal exact search confused objective proof with set enumeration");
  expect(exact_survivors.size() > 1 &&
             exact_survivors.front().component_offsets.size() ==
                 components.size() &&
             std::fabs(exact_survivors.front().candidate.shift_norm2 -
                       exact.shift_norm2) < 1.0e-12,
         "exact search did not expose its ordered certification band");

  const std::vector<std::pair<double, double>> objective_cases{
      {-9.8, -5.0}, {0.2, 5.0}, {10.2, 15.0}, {-9.8, 15.0},
      {10.2, -5.0}, {0.2, -5.0}};
  for (const auto &[case_target_x1, case_target_x2] : objective_cases) {
    target_x1 = case_target_x1;
    target_x2 = case_target_x2;
    double case_brute_best = std::numeric_limits<double>::infinity();
    for (std::size_t s1 = 0; s1 < shifts.size(); ++s1) {
      for (std::size_t s2 = 0; s2 < shifts.size(); ++s2) {
        const std::vector<double> x{atoms[0].x,
                                    atoms[1].x + shifts[s1].tx,
                                    atoms[2].x + shifts[s2].tx};
        if (assignment_clashes(atoms, x, cutoffs.clash_cutoff))
          continue;
        const std::vector<Vec3> points{{x[0], 0.0, 0.0},
                                       {x[1], 0.0, 0.0},
                                       {x[2], 0.0, 0.0}};
        case_brute_best = std::min(case_brute_best, objective(points));
      }
    }
    std::vector<Vec3> case_points;
    PbctopoComponentSearchDiagnostics case_diagnostics;
    const auto case_exact = search_component_shift_candidate(
        atoms, Box, components, shifts, {}, {}, cutoffs, case_points,
        callbacks, nullptr, nullptr, nullptr, nullptr, AnchorBackend::Auto,
        nullptr, &case_diagnostics);
    expect(case_exact.valid && case_diagnostics.optimality_gap_known &&
               std::fabs(objective(case_points) - case_brute_best) < 1.0e-12,
           "exact branch-and-bound disagrees with deterministic brute force");
  }
  target_x1 = 10.2;
  target_x2 = -5.0;

  std::size_t brute_safe_assignments = 0;
  for (std::size_t s1 = 0; s1 < shifts.size(); ++s1) {
    for (std::size_t s2 = 0; s2 < shifts.size(); ++s2) {
      const std::vector<double> x{atoms[0].x,
                                  atoms[1].x + shifts[s1].tx,
                                  atoms[2].x + shifts[s2].tx};
      brute_safe_assignments +=
          !assignment_clashes(atoms, x, cutoffs.clash_cutoff);
    }
  }

  PbctopoContinuityPair soft_pair;
  soft_pair.a = 0;
  soft_pair.b = 2;
  soft_pair.ref_d2 = 25.0;
  soft_pair.inv_ref_d2 = 1.0 / 25.0;
  soft_pair.kind = PbctopoContinuityKind::MetadataSoft;
  soft_pair.weight = 1.0;
  const std::vector<PbctopoContinuityPair> soft_continuity{soft_pair};
  const PbctopoComponentSearchCallbacks continuity_callbacks{
      build_seed, score_continuity_assignment, better_continuity_candidate,
      true, true};
  std::vector<Vec3> pruned_points;
  PbctopoComponentSearchDiagnostics pruned_diagnostics;
  (void)search_component_shift_candidate(
      atoms, Box, components, shifts, soft_continuity, {}, cutoffs,
      pruned_points, continuity_callbacks, nullptr, nullptr, nullptr, nullptr,
      AnchorBackend::Auto, nullptr, &pruned_diagnostics);
  expect(pruned_diagnostics.incumbent_bound_pruning_used &&
             !pruned_diagnostics.feasible_assignment_set_enumerated,
         "normal exact search did not expose objective-bound pruning");

  PbctopoComponentSearchLimits enumeration_limits;
  enumeration_limits.enumerate_feasible_assignments = true;
  std::vector<Vec3> enumerated_points;
  PbctopoComponentSearchDiagnostics enumerated_diagnostics;
  const auto enumerated = search_component_shift_candidate(
      atoms, Box, components, shifts, soft_continuity, {}, cutoffs,
      enumerated_points, continuity_callbacks, nullptr, nullptr, nullptr,
      nullptr, AnchorBackend::Auto, nullptr, &enumerated_diagnostics, nullptr,
      0, nullptr, &enumeration_limits);
  expect(enumerated.valid &&
             !enumerated_diagnostics.incumbent_bound_pruning_used &&
             enumerated_diagnostics.feasible_assignment_set_enumerated &&
             enumerated_diagnostics.feasible_assignments_within_domain ==
                 brute_safe_assignments &&
             enumerated_diagnostics
                     .soft_score_valid_assignments_within_domain ==
                 brute_safe_assignments,
         "audit enumeration did not recover the full feasible assignment set");

  const PbctopoComponentSearchCallbacks selectively_invalid_callbacks{
      build_seed, score_selectively_invalid_assignment, better_candidate,
      true, true};
  std::vector<Vec3> selectively_invalid_points;
  PbctopoComponentSearchDiagnostics selectively_invalid_diagnostics;
  const auto selectively_invalid = search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs,
      selectively_invalid_points, selectively_invalid_callbacks, nullptr,
      nullptr, nullptr, nullptr, AnchorBackend::Auto, nullptr,
      &selectively_invalid_diagnostics, nullptr, 0, nullptr,
      &enumeration_limits);
  expect(selectively_invalid.valid &&
             selectively_invalid_diagnostics
                 .feasible_assignment_set_enumerated &&
             selectively_invalid_diagnostics
                     .feasible_assignments_within_domain ==
                 brute_safe_assignments &&
             selectively_invalid_diagnostics
                     .soft_score_valid_assignments_within_domain > 0 &&
             selectively_invalid_diagnostics
                     .soft_score_valid_assignments_within_domain <
                 selectively_invalid_diagnostics
                     .feasible_assignments_within_domain,
         "hard-feasible leaf count still depends on soft scorer validity");

  auto incomplete_certificate_callbacks = continuity_callbacks;
  incomplete_certificate_callbacks.hard_certificate_domain_complete = false;
  std::vector<Vec3> incomplete_certificate_points;
  PbctopoComponentSearchDiagnostics incomplete_certificate_diagnostics;
  (void)search_component_shift_candidate(
      atoms, Box, components, shifts, soft_continuity, {}, cutoffs,
      incomplete_certificate_points, incomplete_certificate_callbacks,
      nullptr, nullptr, nullptr, nullptr, AnchorBackend::Auto, nullptr,
      &incomplete_certificate_diagnostics, nullptr, 0, nullptr,
      &enumeration_limits);
  expect(!incomplete_certificate_diagnostics
              .feasible_assignment_set_enumerated &&
             !incomplete_certificate_diagnostics
                  .bounded_domain_feasibility_known &&
             incomplete_certificate_diagnostics
                     .feasible_assignments_within_domain == 0 &&
             incomplete_certificate_diagnostics
                     .soft_score_valid_assignments_within_domain == 0,
         "incomplete certificate domain exposed a hard-feasible count");

  const PbctopoComponentSearchCallbacks tied_callbacks{
      build_seed, score_tied_assignment, better_candidate, true, true};
  std::vector<Vec3> tied_points;
  PbctopoComponentSearchDiagnostics tied_diagnostics;
  const auto tied = search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, tied_points,
      tied_callbacks, nullptr, nullptr, nullptr, nullptr, AnchorBackend::Auto,
      nullptr, &tied_diagnostics);
  expect(tied.valid && tied_diagnostics.optimality_gap_known &&
             tied_diagnostics.equivalent_optima_complete &&
             tied_diagnostics.equivalent_best_assignments_within_domain ==
                 brute_safe_assignments,
         "exact search lost comparator-equivalent tied assignments");

  PbctopoCutCandidate incumbent;
  incumbent.valid = true;
  incumbent.shift_norm2 = -1.0;
  incumbent.continuity_path = 0.0;
  incumbent.continuity_max_d2 = 0.0;
  const PbctopoComponentSearchCallbacks bounded_callbacks{
      build_incumbent_seed, score_assignment, better_candidate, true, true};
  const std::vector<Vec3> out_of_domain_incumbent_points{
      {0.0, 0.0, 0.0}, {123.0, 0.0, 0.0}, {-456.0, 0.0, 0.0}};
  std::vector<Vec3> bounded_points;
  PbctopoTimingStats bounded_timing;
  PbctopoComponentSearchDiagnostics bounded_diagnostics;
  const auto bounded = search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, bounded_points,
      bounded_callbacks, nullptr, &incumbent,
      &out_of_domain_incumbent_points, &bounded_timing, AnchorBackend::Auto,
      nullptr, &bounded_diagnostics);
  expect(bounded.valid &&
             std::fabs(objective(bounded_points) - brute_best) < 1.0e-12 &&
             bounded_points.size() == atoms.size() &&
             std::fabs(bounded_points[1].x -
                       out_of_domain_incumbent_points[1].x) > 1.0,
         "exact search accepted an incumbent outside its bounded shift domain");
  expect(bounded_timing.component_search_seed_graph == 0 &&
             bounded_timing.component_search_seed_none == 1,
         "exact search used an unverified graph incumbent as a proof bound");
  expect(bounded_diagnostics.search_completed_under_admissible_bounds &&
             bounded_diagnostics.best_objective_proven &&
             bounded_diagnostics.equivalent_optima_complete &&
             bounded_diagnostics.bounded_domain_feasibility_known &&
             bounded_diagnostics.bounded_domain_feasible &&
             bounded_diagnostics.optimality_gap_known &&
             bounded_diagnostics.optimality_gap == 0.0,
         "domain-enumerated exact search did not prove its bounded optimum");

  PbctopoComponentSearchLimits state_limits;
  state_limits.max_states = 1;
  std::vector<Vec3> state_limited_points;
  PbctopoComponentSearchDiagnostics state_limited_diagnostics;
  (void)search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, state_limited_points,
      callbacks, nullptr, nullptr, nullptr, nullptr, AnchorBackend::Auto,
      nullptr, &state_limited_diagnostics, nullptr, 0, nullptr,
      &state_limits);
  expect(state_limited_diagnostics.resource_limit_hit &&
             state_limited_diagnostics.stop_reason ==
                 PbctopoSearchStopReason::StateLimit &&
             state_limited_diagnostics.states_visited == 1 &&
             !state_limited_diagnostics.best_objective_proven &&
             !state_limited_diagnostics.equivalent_optima_complete &&
             !state_limited_diagnostics.bounded_domain_feasibility_known &&
             !state_limited_diagnostics
                  .search_completed_under_admissible_bounds &&
             !state_limited_diagnostics.optimality_gap_known,
         "state-limited exact audit was reported as complete");

  PbctopoComponentSearchLimits memory_limits;
  memory_limits.max_memory_bytes = 1;
  std::vector<Vec3> memory_limited_points;
  PbctopoComponentSearchDiagnostics memory_limited_diagnostics;
  (void)search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, memory_limited_points,
      callbacks, nullptr, nullptr, nullptr, nullptr, AnchorBackend::Auto,
      nullptr, &memory_limited_diagnostics, nullptr, 0, nullptr,
      &memory_limits);
  expect(memory_limited_diagnostics.resource_limit_hit &&
             memory_limited_diagnostics.stop_reason ==
                 PbctopoSearchStopReason::MemoryLimit &&
             memory_limited_diagnostics.memory_bytes_estimate > 1 &&
             memory_limited_diagnostics.memory_budget_is_estimate &&
             !memory_limited_diagnostics.optimality_gap_known,
         "memory-limited exact audit was reported as complete");

  const std::vector<PbctopoRequiredContact> required_contacts{
      {1, 2, 5.0, 1.0}};
  std::vector<Vec3> constrained_points;
  PbctopoTimingStats constrained_timing;
  PbctopoComponentSearchDiagnostics constrained_diagnostics;
  std::vector<PbctopoComponentSearchSurvivor> constrained_survivors;
  const auto constrained = search_component_shift_candidate(
      atoms, Box, components, shifts, {}, {}, cutoffs, constrained_points,
      callbacks, nullptr, nullptr, nullptr, &constrained_timing,
      AnchorBackend::Auto, nullptr, &constrained_diagnostics,
      &constrained_survivors, 8, &required_contacts);
  expect(constrained.valid && !constrained_survivors.empty() &&
             constrained_survivors.front().component_offsets[1].x ==
                 constrained_survivors.front().component_offsets[2].x &&
             constrained_survivors.front().component_offsets[1].x != 0,
         "required interface contact was not enforced in exact domains");
  expect(constrained_timing.component_required_contact_constraint_pairs == 1 &&
             constrained_timing
                     .component_required_contact_relative_deltas_checked > 0 &&
             constrained_timing
                     .component_required_contact_relative_deltas_rejected > 0,
         "required-contact domain telemetry was not recorded");

  std::vector<Vec3> beam_points;
  PbctopoTimingStats beam_timing;
  PbctopoComponentSearchDiagnostics beam_diagnostics;
  std::vector<PbctopoComponentSearchSurvivor> beam_survivors;
  const auto beam = search_component_shift_candidate_beam(
      atoms, Box, components, shifts, {}, {}, cutoffs, beam_points, 64,
      callbacks, nullptr, nullptr, nullptr, &beam_timing,
      AnchorBackend::Auto, nullptr, 0, &beam_diagnostics, &beam_survivors, 8);
  expect(beam.valid && std::fabs(objective(beam_points) - brute_best) < 1.0e-12,
         "beam search differs from brute force when all states fit the beam");
  expect(beam_timing.component_large_mode_beam_final_rescore_all_frames == 1,
          "beam search did not exact-rescore every surviving complete state");
  expect(beam_diagnostics.beam_used && !beam_diagnostics.exhaustive &&
             !beam_diagnostics.optimality_gap_known &&
             beam_diagnostics.beam_survivors_exactly_rescored > 0,
          "beam diagnostics incorrectly claim global optimality");
  expect(beam_survivors.size() > 1 &&
             beam_survivors.front().component_offsets.size() ==
                 components.size() &&
             std::fabs(beam_survivors.front().candidate.shift_norm2 -
                       beam.shift_norm2) < 1.0e-12,
         "beam search did not expose its exact-rescored survivor band");

  std::vector<PbctopoAtom> heavy_search_atoms(2);
  heavy_search_atoms[0].x = 0.0;
  heavy_search_atoms[1].x = 3.0;
  for (auto &atom : heavy_search_atoms) {
    atom.is_heavy = true;
    atom.has_vdw_radius = true;
    atom.vdw_radius = 1.0;
    atom.steric_radius = 1.0;
  }
  const std::vector<std::vector<std::size_t>> hydrogen_components{{0}, {1}};
  const PbctopoComponentSearchCallbacks hydrogen_callbacks{
      build_seed, score_hydrogen_band_assignment, better_candidate, true,
      true};
  std::vector<Vec3> hydrogen_search_points;
  std::vector<PbctopoComponentSearchSurvivor> hydrogen_survivors;
  const auto hydrogen_search = search_component_shift_candidate(
      heavy_search_atoms, Box, hydrogen_components, shifts, {}, {}, cutoffs,
      hydrogen_search_points, hydrogen_callbacks, nullptr, nullptr, nullptr,
      nullptr, AnchorBackend::Auto, nullptr, nullptr, &hydrogen_survivors, 3);
  expect(hydrogen_search.valid && hydrogen_survivors.size() >= 2 &&
             hydrogen_survivors.front().component_offsets[1].x == 0,
         "hydrogen-band fixture did not rank the clashing heavy-safe state first");

  std::vector<PbctopoAtom> all_atoms(4);
  for (std::size_t atom_idx : {std::size_t{0}, std::size_t{2}}) {
    all_atoms[atom_idx].is_heavy = true;
    all_atoms[atom_idx].has_vdw_radius = true;
    all_atoms[atom_idx].vdw_radius = 1.0;
    all_atoms[atom_idx].steric_radius = 1.0;
  }
  for (std::size_t atom_idx : {std::size_t{1}, std::size_t{3}}) {
    all_atoms[atom_idx].is_heavy = false;
    all_atoms[atom_idx].has_vdw_radius = true;
    all_atoms[atom_idx].vdw_radius = 0.7;
    all_atoms[atom_idx].steric_radius = 0.7;
  }
  const std::vector<Vec3> base_all_points{{0.0, 0.0, 0.0},
                                           {1.1, 0.0, 0.0},
                                           {3.0, 0.0, 0.0},
                                           {1.9, 0.0, 0.0}};
  const std::vector<std::size_t> all_atom_components{0, 0, 1, 1};
  const std::unordered_set<std::uint64_t> no_exclusions;
  std::size_t first_all_atom_safe = hydrogen_survivors.size();
  for (std::size_t survivor_idx = 0;
       survivor_idx < hydrogen_survivors.size(); ++survivor_idx) {
    std::vector<Vec3> placed = base_all_points;
    const double translation =
        10.0 * hydrogen_survivors[survivor_idx].component_offsets[1].x;
    placed[2].x += translation;
    placed[3].x += translation;
    const auto audit = evaluate_pbctopo_final_sterics(
        all_atoms, placed, all_atom_components, no_exclusions,
        cutoffs.clash_cutoff);
    expect(audit.input_valid,
           "all-atom certification-band audit rejected valid inputs");
    if (audit.hard_clashes == 0) {
      first_all_atom_safe = survivor_idx;
      break;
    }
  }
  expect(first_all_atom_safe == 1,
         "sequential all-atom certification did not recover the second survivor");

  std::vector<Vec3> narrow_beam_points;
  PbctopoTimingStats narrow_beam_timing;
  PbctopoComponentSearchDiagnostics narrow_beam_diagnostics;
  (void)search_component_shift_candidate_beam(
      atoms, Box, components, shifts, {}, {}, cutoffs, narrow_beam_points, 1,
      callbacks, nullptr, nullptr, nullptr, &narrow_beam_timing,
      AnchorBackend::Auto, nullptr, 0, &narrow_beam_diagnostics);
  expect(narrow_beam_timing.component_search_beam_state_copies_avoided > 0,
         "beam expansion materialized children that were outside the beam");

  std::vector<PbctopoAtom> wide_radius_atoms(2);
  wide_radius_atoms[0].x = 0.0;
  wide_radius_atoms[1].x = 2.0;
  for (auto &atom : wide_radius_atoms) {
    atom.is_heavy = true;
    atom.has_vdw_radius = true;
    atom.vdw_radius = 1.2;
    atom.steric_radius = 1.2;
  }
  const std::vector<std::vector<std::size_t>> two_components{{0}, {1}};
  PbctopoShift identity;
  identity.shift_index = 0;
  const std::vector<PbctopoShift> identity_only{identity};
  std::vector<Vec3> clash_points;
  PbctopoComponentSearchDiagnostics clash_diagnostics;
  const auto exact_clash = search_component_shift_candidate(
      wide_radius_atoms, Box, two_components, identity_only, {}, {}, cutoffs,
      clash_points, callbacks, nullptr, nullptr, nullptr, nullptr,
      AnchorBackend::Auto, nullptr, &clash_diagnostics);
  expect(!exact_clash.valid,
         "exact search missed a vdW clash beyond the scalar fallback cutoff");
  expect(clash_diagnostics.bounded_domain_feasibility_known &&
             !clash_diagnostics.bounded_domain_feasible &&
             !clash_diagnostics.best_objective_proven &&
             !clash_diagnostics.optimality_gap_known,
         "bounded-domain infeasibility was reported as a physical optimum");
  const auto beam_clash = search_component_shift_candidate_beam(
      wide_radius_atoms, Box, two_components, identity_only, {}, {}, cutoffs,
      clash_points, 8, callbacks);
  expect(!beam_clash.valid,
         "beam search missed a vdW clash beyond the scalar fallback cutoff");

  std::cout << "PASS pbctopo_component_search_optimality_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\component_search.cpp `
        src\pbctopo\component_search_pair_table.cpp `
        src\pbctopo\component_shift_domain.cpp `
        src\pbctopo\component_shift_geometry.cpp `
        src\pbctopo\component_metric.cpp `
        src\pbctopo\candidate_order.cpp `
        src\pbctopo\required_contacts.cpp `
        src\pbctopo\certification.cpp `
        src\pbctopo\steric_clash.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO component-search optimality smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO component-search optimality smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
