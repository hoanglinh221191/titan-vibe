#include "titan_path/control.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  using titan_path::RouteMode;
  if (titan_path::next_route_mode(std::nullopt, 4.1, 4.0, 3.0) !=
      RouteMode::Recovery)
    return 1;
  if (titan_path::next_route_mode(RouteMode::Recovery, 3.5, 4.0, 3.0) !=
      RouteMode::Recovery)
    return 2;
  if (titan_path::next_route_mode(RouteMode::Recovery, 2.9, 4.0, 3.0) !=
      RouteMode::Advance)
    return 3;

  const std::vector<titan_path::SelectionCandidate> candidates{
      {-10.00, 0.5, 0.5, 2.0, 5.0, true, 0.0, false, 0.2},
      {-9.97, 1.8, 1.8, 2.1, 10.0, true, 0.0, false, 0.3},
      {-9.80, 2.0, 2.0, 2.2, 4.0, true, 0.0, false, 0.4},
  };
  const auto epsilon = titan_path::select_epsilon_progress_candidate(
      candidates, RouteMode::Advance, 0.05, false);
  if (!epsilon.found || epsilon.index != 1 || epsilon.eligible_count != 2)
    return 4;
  const auto geometry = titan_path::select_epsilon_progress_candidate(
      candidates, RouteMode::Advance, 0.0, true);
  if (!geometry.found || geometry.index != 2)
    return 5;
  const auto recovery = titan_path::select_epsilon_progress_candidate(
      candidates, RouteMode::Recovery, 0.05, false);
  if (!recovery.found || recovery.index != 0)
    return 6;

  const std::vector<titan_path::SelectionCandidate> zero_progress_low_energy{
      {-100.0, 0.0, 0.4, 2.0, 1.0, true},
      {-9.0, 0.25, 0.7, 1.8, 2.0, true},
      {-8.9, 0.50, 1.0, 1.9, 3.0, true},
  };
  const auto forward_only = titan_path::select_epsilon_progress_candidate(
      zero_progress_low_energy, RouteMode::Advance, 0.05, false);
  if (!forward_only.found || forward_only.index != 1 ||
      forward_only.valid_count != 3 || forward_only.selectable_count != 2 ||
      forward_only.nonforward_count != 1 || forward_only.eligible_count != 1)
    return 7;
  const std::vector<titan_path::SelectionCandidate> stalled{
      {-100.0, 0.0, 0.4, 2.0, 1.0, true},
  };
  const auto no_forward = titan_path::select_epsilon_progress_candidate(
      stalled, RouteMode::Advance, 0.05, false);
  if (no_forward.found || no_forward.selectable_count != 0 ||
      no_forward.nonforward_count != 1)
    return 8;
  using titan_path::FinalSelectionFailure;
  const auto forward_failure = titan_path::classify_final_selection_failure(
      no_forward, RouteMode::Advance);
  if (forward_failure != FinalSelectionFailure::ForwardProgressExhausted ||
      std::string(titan_path::final_selection_failure_name(forward_failure)) !=
          "FORWARD_PROGRESS_EXHAUSTED")
    return 17;
  if (titan_path::classify_final_selection_failure(forward_only,
                                                    RouteMode::Advance) !=
      FinalSelectionFailure::None)
    return 18;

  const auto capped = titan_path::project_vector_to_turn_cone(
      {0.0, 2.0, 0.0}, {1.0, 0.0, 0.0}, 20.0);
  if (!capped.valid || !capped.capped ||
      std::fabs(capped.requested_angle_deg - 90.0) > 1.0e-10 ||
      std::fabs(capped.applied_angle_deg - 20.0) > 1.0e-10 ||
      std::fabs(std::sqrt(capped.vector[0] * capped.vector[0] +
                         capped.vector[1] * capped.vector[1] +
                         capped.vector[2] * capped.vector[2]) -
                2.0) > 1.0e-10)
    return 9;
  const auto unchanged = titan_path::project_vector_to_turn_cone(
      {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 20.0);
  if (!unchanged.valid || unchanged.capped ||
      std::fabs(unchanged.vector[0] - 1.0) > 1.0e-12)
    return 10;
  const auto antiparallel = titan_path::project_vector_to_turn_cone(
      {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 20.0);
  if (antiparallel.valid)
    return 11;

  // The old slerp calculation declared 20 degrees but produced 20.215...
  // degrees at epsilon=5e-8. Check actual vectors, not reported telemetry.
  for (const double epsilon : {1.0e-3, 1.0e-6, 5.0e-8, 1.0e-8}) {
    const auto near_antiparallel = titan_path::project_vector_to_turn_cone(
        {-std::cos(epsilon), std::sin(epsilon), 0.0},
        {1.0, 0.0, 0.0}, 20.0);
    const double actual_angle =
        std::atan2(std::hypot(near_antiparallel.vector[1],
                              near_antiparallel.vector[2]),
                   near_antiparallel.vector[0]) *
        180.0 / std::acos(-1.0);
    if (!near_antiparallel.valid || !near_antiparallel.capped ||
        std::fabs(actual_angle - 20.0) > 1.0e-10 ||
        std::fabs(actual_angle - near_antiparallel.applied_angle_deg) > 1.0e-10 ||
        std::fabs(std::hypot(near_antiparallel.vector[0],
                             near_antiparallel.vector[1],
                             near_antiparallel.vector[2]) - 1.0) > 1.0e-12)
      return 31;
  }
  const double inverse_sqrt3 = 1.0 / std::sqrt(3.0);
  const double inverse_sqrt2 = 1.0 / std::sqrt(2.0);
  const std::array<double, 3> oblique_old{
      inverse_sqrt3, inverse_sqrt3, inverse_sqrt3};
  const double oblique_epsilon = 5.0e-8;
  const std::array<double, 3> oblique_request{
      -std::cos(oblique_epsilon) * inverse_sqrt3 +
          std::sin(oblique_epsilon) * inverse_sqrt2,
      -std::cos(oblique_epsilon) * inverse_sqrt3 -
          std::sin(oblique_epsilon) * inverse_sqrt2,
      -std::cos(oblique_epsilon) * inverse_sqrt3};
  const auto oblique_turn = titan_path::project_vector_to_turn_cone(
      oblique_request, oblique_old, 60.0);
  const double oblique_cosine =
      (oblique_turn.vector[0] + oblique_turn.vector[1] + oblique_turn.vector[2]) *
      inverse_sqrt3 /
      std::hypot(oblique_turn.vector[0], oblique_turn.vector[1],
                  oblique_turn.vector[2]);
  if (!oblique_turn.valid || !oblique_turn.capped ||
      std::fabs(std::acos(oblique_cosine) * 180.0 / std::acos(-1.0) - 60.0) >
          1.0e-10)
    return 32;

  for (const double length : {1.0, 0.5, 0.25}) {
    const auto probe =
        titan_path::scale_vector_to_length({3.0, 4.0, 0.0}, length);
    if (!probe.has_value() ||
        std::fabs(std::sqrt((*probe)[0] * (*probe)[0] +
                           (*probe)[1] * (*probe)[1] +
                           (*probe)[2] * (*probe)[2]) -
                  length) > 1.0e-12 ||
        std::fabs((*probe)[0] / (*probe)[1] - 0.75) > 1.0e-12)
      return 19;
  }
  if (titan_path::scale_vector_to_length({0.0, 0.0, 0.0}, 1.0).has_value())
    return 20;
  bool rejected_probe_length = false;
  try {
    (void)titan_path::scale_vector_to_length({1.0, 0.0, 0.0}, 0.0);
  } catch (const std::invalid_argument &) {
    rejected_probe_length = true;
  }
  if (!rejected_probe_length)
    return 21;

  const auto large_probe =
      titan_path::scale_vector_to_length({1.0e200, 0.0, 0.0}, 1.0);
  const auto large_target =
      titan_path::scale_vector_to_length({1.0e200, 1.0e200, 0.0}, 1.0e200);
  if (!large_probe.has_value() || (*large_probe)[0] != 1.0 ||
      !large_target.has_value() ||
      std::fabs(std::hypot((*large_target)[0], (*large_target)[1],
                           (*large_target)[2]) / 1.0e200 - 1.0) > 1.0e-12)
    return 33;
  bool rejected_unrepresentable_norm = false;
  try {
    const double maximum = std::numeric_limits<double>::max();
    (void)titan_path::scale_vector_to_length({maximum, maximum, maximum}, 1.0);
  } catch (const std::invalid_argument &) {
    rejected_unrepresentable_norm = true;
  }
  if (!rejected_unrepresentable_norm)
    return 34;
  const auto large_turn = titan_path::project_vector_to_turn_cone(
      {0.0, 1.0e200, 0.0}, {1.0e200, 0.0, 0.0}, 20.0);
  if (!large_turn.valid ||
      std::fabs(large_turn.applied_angle_deg - 20.0) > 1.0e-10 ||
      std::fabs(std::hypot(large_turn.vector[0], large_turn.vector[1],
                           large_turn.vector[2]) / 1.0e200 - 1.0) > 1.0e-12)
    return 35;

  std::vector<titan_path::SelectionCandidate> trust_region{
      {-100.0, 0.20, 0.25, 2.0, 1.0, true},
      {-9.0, 0.15, 0.50, 2.0, 1.0, true},
      {-8.0, 0.10, 1.00, 2.0, 1.0, true},
  };
  for (auto &candidate : trust_region) {
    candidate.probe_length = candidate.pull_distance;
    candidate.fixed_probe = true;
  }
  const auto longest_probe = titan_path::select_epsilon_progress_candidate(
      trust_region, RouteMode::Advance, 0.05, false, true);
  if (!longest_probe.found || longest_probe.index != 2 ||
      longest_probe.valid_count != 3 || longest_probe.selectable_count != 1)
    return 22;

  std::vector<titan_path::SelectionCandidate> tier_after_progress{
      {-100.0, 0.0, 1.00, 2.0, 1.0, true},
      {-9.0, 0.20, 0.50, 2.0, 1.0, true},
      {-8.0, 0.10, 0.25, 2.0, 1.0, true},
  };
  for (auto &candidate : tier_after_progress) {
    candidate.probe_length = candidate.pull_distance;
    candidate.fixed_probe = true;
  }
  const auto shorter_forward = titan_path::select_epsilon_progress_candidate(
      tier_after_progress, RouteMode::Advance, 0.05, false, true);
  if (!shorter_forward.found || shorter_forward.index != 1 ||
      shorter_forward.valid_count != 3 ||
      shorter_forward.selectable_count != 1 ||
      shorter_forward.nonforward_count != 1)
    return 25;

  std::vector<titan_path::SelectionCandidate> recovery_tiers{
      {-100.0, 0.1, 1.00, 2.0, 1.0, true},
      {-9.0, 0.0, 0.50, 2.0, 1.0, true},
      {-8.0, 0.0, 0.25, 2.0, 1.0, true},
  };
  for (auto &candidate : recovery_tiers) {
    candidate.probe_length = candidate.pull_distance;
    candidate.fixed_probe = true;
    candidate.reference_path_distance = 4.5;
  }
  recovery_tiers[0].path_distance = 4.7;
  recovery_tiers[1].path_distance = 4.0;
  recovery_tiers[2].path_distance = 3.9;
  const auto shorter_descent = titan_path::select_epsilon_progress_candidate(
      recovery_tiers, RouteMode::Recovery, 0.05, false, true);
  if (!shorter_descent.found || shorter_descent.index != 1 ||
      shorter_descent.selectable_count != 1)
    return 26;
  recovery_tiers[1].path_distance = 4.5;
  recovery_tiers[2].path_distance = 4.6;
  const auto no_descent = titan_path::select_epsilon_progress_candidate(
      recovery_tiers, RouteMode::Recovery, 0.05, false, true);
  if (no_descent.found || no_descent.selectable_count != 0 ||
      titan_path::classify_final_selection_failure(no_descent,
                                                   RouteMode::Recovery) !=
          FinalSelectionFailure::RecoveryDescentExhausted)
    return 27;

  std::vector<titan_path::SelectionCandidate> navigation_pool{
      {-100.0, 0.90, 0.50, 2.0, 1.0, true}, // shorter tier: excluded first
      {-10.00, 0.20, 1.00, 2.0, 1.0, true}, // energy minimum
      {-9.97, 0.40, 1.00, 2.0, 1.0, true},  // nav minimum
      {-9.96, 0.70, 1.00, 2.0, 1.0, true},  // nav epsilon + max progress
      {-9.80, 1.00, 1.00, 2.0, 1.0, true},  // outside energy epsilon
  };
  const double navigation_cost[] = {0.0, 0.50, 0.10, 0.12, 0.0};
  for (std::size_t index = 0; index < navigation_pool.size(); ++index) {
    auto &candidate = navigation_pool[index];
    candidate.probe_length = candidate.pull_distance;
    candidate.fixed_probe = true;
    candidate.navigation_enabled = true;
    candidate.navigation_cost = navigation_cost[index];
  }
  const auto nested_navigation =
      titan_path::select_epsilon_progress_candidate(
          navigation_pool, RouteMode::Advance, 0.05, false, true, 0.02);
  if (!nested_navigation.found || nested_navigation.index != 3 ||
      nested_navigation.valid_count != 5 ||
      nested_navigation.mode_eligible_count != 5 ||
      nested_navigation.maximum_probe_tier_count != 4 ||
      nested_navigation.eligible_count != 3 ||
      nested_navigation.navigation_eligible_count != 2 ||
      !std::isfinite(nested_navigation.minimum_navigation_cost) ||
      std::fabs(nested_navigation.minimum_navigation_cost - 0.10) > 1.0e-12)
    return 28;
  const auto navigation_disabled =
      titan_path::select_epsilon_progress_candidate(
          navigation_pool, RouteMode::Advance, 0.05, false, true);
  if (!navigation_disabled.found || navigation_disabled.index != 3 ||
      navigation_disabled.navigation_eligible_count !=
          navigation_disabled.eligible_count)
    return 29;
  bool rejected_recovery_navigation = false;
  try {
    (void)titan_path::select_epsilon_progress_candidate(
        navigation_pool, RouteMode::Recovery, 0.05, false, true, 0.02);
  } catch (const std::invalid_argument &) {
    rejected_recovery_navigation = true;
  }
  if (!rejected_recovery_navigation)
    return 30;

  const auto recovery_turn = titan_path::project_vector_to_turn_cone(
      {0.0, -4.5, 0.0}, {1.0, 0.0, 0.0}, 60.0);
  if (!recovery_turn.valid || !recovery_turn.capped)
    return 23;
  const auto recovery_probe =
      titan_path::scale_vector_to_length(recovery_turn.vector, 1.0);
  if (!recovery_probe.has_value() ||
      std::fabs(4.5 + (*recovery_probe)[1]) >= 4.5)
    return 24;

  if (!titan_path::candidate_within_scan_shell(0.5, 0.5, 4.0) ||
      !titan_path::candidate_within_scan_shell(4.0, 0.5, 4.0) ||
      titan_path::candidate_within_scan_shell(0.49, 0.5, 4.0) ||
      titan_path::candidate_within_scan_shell(4.01, 0.5, 4.0))
    return 12;

  using titan_path::CandidateExhaustion;
  const auto scan_exhausted =
      titan_path::classify_candidate_exhaustion(12, 0, 12, 0, 0);
  if (scan_exhausted != CandidateExhaustion::ScanShellExhausted ||
      std::string(titan_path::candidate_exhaustion_name(scan_exhausted)) !=
          "SCAN_SHELL_EXHAUSTED")
    return 13;
  const auto clearance_exhausted =
      titan_path::classify_candidate_exhaustion(10, 6, 4, 6, 0);
  if (clearance_exhausted != CandidateExhaustion::ClearanceExhausted ||
      std::string(titan_path::candidate_exhaustion_name(clearance_exhausted)) !=
          "CLEARANCE_EXHAUSTED")
    return 14;
  if (titan_path::classify_candidate_exhaustion(11, 6, 4, 5, 1) !=
      CandidateExhaustion::None)
    return 15;
  bool rejected_inconsistent_accounting = false;
  try {
    (void)titan_path::classify_candidate_exhaustion(5, 5, 1, 0, 0);
  } catch (const std::invalid_argument &) {
    rejected_inconsistent_accounting = true;
  }
  if (!rejected_inconsistent_accounting)
    return 16;

  std::cout << "PASS titan_path_control_smoke\n";
  return 0;
}
