#include "titan_path/clearance.h"
#include "titan_path/transit.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace titan_path;
namespace {
int cases = 0;
void check(bool ok, const char *msg) {
  if (!ok)
    throw std::runtime_error(msg);
}
void test(const char *name, const std::function<void()> &f) {
  f();
  ++cases;
  std::cout << "PASS " << name << "\n";
}
void rejects(const std::function<void()> &f) {
  bool failed = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  check(failed, "expected invalid input rejection");
}
TransitSnapshot approach() {
  TransitSnapshot s;
  s.route_sha256 = std::string(64, 'a');
  s.snapshot_sha256 = std::string(64, 'b');
  s.geometry_validated = true;
  s.registered_waypoints = {{0, 0, 0}, {5, 0, 0}};
  s.waypoint_radii = {2, 2};
  s.receptor = {{0, 20, 0}};
  s.ligand = {{-1.6, .2, 0}};
  s.ligand_masses = {12};
  s.old_direction = {1, 0, 0};
  return s;
}
TransitSnapshot corner() {
  auto s = approach();
  s.registered_waypoints = {{-3, 0, 0}, {0, 0, 0}, {0, 5, 0}};
  s.waypoint_radii = {3, 3, 3};
  s.ligand = {{.4, -.6, 0}};
  s.old_direction = {.6, .8, 0};
  return s;
}
void validate_trace(const TransitPlan &r, const TransitPolicy &p) {
  check(r.counters.evaluations <= p.candidate_budget,
        "candidate budget overrun");
  check(r.model_steps.size() <= p.max_steps + 1, "step overrun");
  for (std::size_t i = 1; i < r.model_steps.size(); ++i) {
    const auto &a = r.model_steps[i - 1];
    const auto &b = r.model_steps[i];
    check(b.merit_angstrom2 < a.merit_angstrom2 - p.merit_decrease_angstrom2,
          "merit did not decrease");
    check(b.progress_angstrom >= a.progress_angstrom - 1e-8, "backtrack");
    check(b.turn_deg <= p.turn_cap_deg + 1e-8, "turn cap");
    check(b.travel_angstrom <= p.max_total_travel_angstrom + 1e-10,
          "travel cap");
    check(b.sweep_clearance_angstrom >= p.clearance_angstrom,
          "sweep threshold");
  }
  if (r.found)
    check(r.model_steps.back().phase == TransitPhase::HandoffAdvance,
          "found without exit");
}
} // namespace
int main() {
  try {
    test("approach_flat_progress_then_entry", [] {
      auto s = approach();
      auto r = plan_bounded_transit(s);
      check(r.initial_phase == TransitPhase::Approach && r.found,
            "approach failed");
      check(r.model_steps.size() >= 3 &&
                r.model_steps[1].progress_angstrom == 0,
            "no ingress plateau");
      validate_trace(r, {});
    });
    test("corner_plateau_then_exit", [] {
      auto r = plan_bounded_transit(corner());
      check(r.initial_phase == TransitPhase::CornerTransit && r.found,
            "corner failed");
      validate_trace(r, {});
    });
    test("ordinary_advance_not_hijacked", [] {
      auto s = approach();
      s.ligand = {{1, .2, 0}};
      auto r = plan_bounded_transit(s);
      check(r.initial_phase == TransitPhase::NotApplicable && !r.found,
            "ordinary route changed");
    });
    test("terminal_not_corner", [] {
      auto s = approach();
      s.ligand = {{6, 0, 0}};
      check(plan_bounded_transit(s).reason == "TRANSIT_NOT_APPLICABLE",
            "terminal intercepted");
    });
    test("straight_vertex_not_corner", [] {
      auto s = approach();
      s.registered_waypoints = {{-3, 0, 0}, {0, 0, 0}, {3, 0, 0}};
      s.waypoint_radii = {2, 2, 2};
      s.ligand = {{0, 1, 0}};
      check(plan_bounded_transit(s).initial_phase ==
                TransitPhase::NotApplicable,
            "straight corner invented");
    });
    test("initial_overlap_retained", [] {
      auto s = approach();
      s.receptor = s.ligand;
      auto r = plan_bounded_transit(s);
      check(r.reason == "INITIAL_CLEARANCE_VIOLATION" && !r.found,
            "initial collision bypassed");
    });
    test("candidate_budget_global", [] {
      auto p = TransitPolicy{};
      p.candidate_budget = 1;
      auto r = plan_bounded_transit(approach(), p);
      check(r.counters.evaluations == 1 &&
                r.reason == "TRANSIT_CANDIDATE_BUDGET_EXHAUSTED",
            "budget semantics");
      validate_trace(r, p);
    });
    test("valid_exit_retained_at_budget", [] {
      auto s = approach();
      s.ligand = {{-.1, 0, 0}};
      auto p = TransitPolicy{};
      p.candidate_budget = 1;
      auto r = plan_bounded_transit(s, p);
      check(r.found && r.reason == "GEOMETRY_PLAN_FOUND_BUDGET_LIMITED" &&
                r.counters.evaluations == 1,
            "discarded valid budget-limited exit");
      validate_trace(r, p);
    });
    test("clearance_specific_exhaustion", [] {
      auto s = approach();
      s.receptor = {{0, .2, 0}};
      auto r = plan_bounded_transit(s);
      check(!r.found && r.reason == "TRANSIT_CLEARANCE_EXHAUSTED" &&
                r.counters.clearance_rejections > 0,
            "clearance reason");
    });
    test("corridor_specific_exhaustion", [] {
      auto s = corner();
      s.waypoint_radii = {.1, .1, .1};
      s.old_direction = {1, 0, 0};
      auto r = plan_bounded_transit(s);
      check(!r.found && r.reason == "TRANSIT_CORRIDOR_EXHAUSTED" &&
                r.counters.clearance_rejections == 0,
            "corridor reason");
    });
    test("step_limit", [] {
      auto p = TransitPolicy{};
      p.max_steps = 1;
      auto r = plan_bounded_transit(approach(), p);
      check(!r.found && r.reason == "TRANSIT_STEP_LIMIT", "step limit ignored");
      validate_trace(r, p);
    });
    test("travel_limit", [] {
      auto p = TransitPolicy{};
      p.max_total_travel_angstrom = .5;
      auto r = plan_bounded_transit(approach(), p);
      check(!r.found && r.counters.bound_rejections > 0,
            "travel limit ignored");
      validate_trace(r, p);
    });
    test("seed_trust_region", [] {
      auto p = TransitPolicy{};
      p.max_seed_displacement_angstrom = .3;
      auto r = plan_bounded_transit(approach(), p);
      check(!r.found && r.counters.bound_rejections > 0, "seed bound ignored");
      validate_trace(r, p);
    });
    test("no_merit_relaxation", [] {
      auto s = approach();
      s.old_direction = {-1, 0, 0};
      auto r = plan_bounded_transit(s);
      check(!r.found && r.model_steps.size() == 1, "away direction accepted");
    });
    test("deterministic_repeat", [] {
      auto a = plan_bounded_transit(corner());
      auto b = plan_bounded_transit(corner());
      check(a.reason == b.reason &&
                a.counters.evaluations == b.counters.evaluations &&
                a.model_steps.size() == b.model_steps.size(),
            "repeat mismatch");
      for (std::size_t i = 0; i < a.model_steps.size(); ++i)
        check(a.model_steps[i].com == b.model_steps[i].com,
              "coordinate nondeterminism");
    });
    test("rigid_transform_invariance", [] {
      auto s = approach();
      auto a = plan_bounded_transit(s);
      auto transform = [](Vec3 v) {
        return Vec3{-v[1] + 7, v[0] - 4, v[2] + 3};
      };
      for (auto &v : s.registered_waypoints)
        v = transform(v);
      for (auto &v : s.receptor)
        v = transform(v);
      for (auto &v : s.ligand)
        v = transform(v);
      s.old_direction = {0, 1, 0};
      auto b = plan_bounded_transit(s);
      check(a.found == b.found, "rigid motion changed feasibility");
      check(std::fabs(a.model_steps.back().progress_angstrom -
                      b.model_steps.back().progress_angstrom) < 1e-8,
            "rigid motion changed progress");
    });
    test("interior_swept_collision", [] {
      std::vector<double> rx{0}, ry{0}, rz{0}, mx{-2}, my{0}, mz{0};
      auto r = evaluate_translation_sweep({rx, ry, rz}, {mx, my, mz}, {4, 0, 0},
                                          1.5);
      check(!r.passes_threshold && r.threshold_violation_hidden_at_endpoints,
            "endpoint-only collision check");
    });
    test("counter_partition", [] {
      auto r = plan_bounded_transit(corner());
      auto c = r.counters;
      check(c.evaluations == c.bound_rejections + c.continuity_rejections +
                                 c.ambiguity_rejections +
                                 c.backtrack_rejections +
                                 c.corridor_rejections + c.merit_rejections +
                                 c.clearance_rejections + c.accepted_edges,
            "counter partition mismatch");
    });
    test("reject_unvalidated", [] {
      auto s = approach();
      s.geometry_validated = false;
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_bad_identity", [] {
      auto s = approach();
      s.route_sha256 = "bad";
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_nan", [] {
      auto s = approach();
      s.ligand[0][0] = std::numeric_limits<double>::quiet_NaN();
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_zero_direction", [] {
      auto s = approach();
      s.old_direction = {0, 0, 0};
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_zero_mass", [] {
      auto s = approach();
      s.ligand_masses[0] = 0;
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_duplicate_waypoint", [] {
      auto s = approach();
      s.registered_waypoints[1] = s.registered_waypoints[0];
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_radius_count", [] {
      auto s = approach();
      s.waypoint_radii.clear();
      rejects([&] { (void)plan_bounded_transit(s); });
    });
    test("reject_zero_budget", [] {
      auto p = TransitPolicy{};
      p.candidate_budget = 0;
      rejects([&] { (void)plan_bounded_transit(approach(), p); });
    });
    test("reject_nondecreasing_ladder", [] {
      auto p = TransitPolicy{};
      p.probe_lengths_angstrom = {1, 1, .25};
      rejects([&] { (void)plan_bounded_transit(approach(), p); });
    });
    std::cout << "NAV6B_TRANSIT_SMOKE_PASS cases=" << cases
              << " md_executed=false\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "NAV6B_TRANSIT_SMOKE_FAIL: " << e.what() << "\n";
    return 1;
  }
}
