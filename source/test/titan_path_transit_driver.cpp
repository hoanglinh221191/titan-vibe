// Simple text-in/JSON-out bridge to the production C++ geometry-only planner.
#include "titan_path/transit.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace titan_path;
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::invalid_argument("usage: transit_driver input.txt");
    std::ifstream in(argv[1]);
    std::string schema;
    in >> schema;
    if (schema != "NAV6B_INPUT_V1")
      throw std::invalid_argument("input schema");
    TransitSnapshot s;
    TransitPolicy p;
    std::size_t nw, nr, nl;
    in >> s.route_sha256 >> s.snapshot_sha256 >> nw >> nr >> nl;
    if (nw < 2 || nw > 10000 || nr < 1 || nr > 1000000 || nl < 1 ||
        nl > 1000000)
      throw std::invalid_argument("input array limits");
    s.geometry_validated = true;
    s.registered_waypoints.resize(nw);
    s.waypoint_radii.resize(nw);
    s.receptor.resize(nr);
    s.ligand.resize(nl);
    s.ligand_masses.resize(nl);
    for (std::size_t i = 0; i < nw; ++i) {
      for (double &x : s.registered_waypoints[i])
        in >> x;
      in >> s.waypoint_radii[i];
    }
    for (auto &v : s.receptor)
      for (double &x : v)
        in >> x;
    for (std::size_t i = 0; i < nl; ++i) {
      for (double &x : s.ligand[i])
        in >> x;
      in >> s.ligand_masses[i];
    }
    for (double &x : s.old_direction)
      in >> x;
    in >> p.clearance_angstrom >> p.turn_cap_deg;
    for (double &l : p.probe_lengths_angstrom)
      in >> l;
    in >> p.lookahead_angstrom >> p.exit_progress_angstrom >>
        p.merit_decrease_angstrom2 >> p.continuity_half_window_angstrom >>
        p.ambiguity_distance_angstrom >>
        p.ambiguity_progress_separation_angstrom >>
        p.max_total_travel_angstrom >> p.max_seed_displacement_angstrom >>
        p.max_steps >> p.beam_width >> p.candidate_budget;
    if (!in)
      throw std::invalid_argument("truncated input");
    std::string extra;
    if (in >> extra)
      throw std::invalid_argument("trailing input");
    const auto r = plan_bounded_transit(s, p);
    const auto &c = r.counters;
    auto vector = [](const Vec3 &v) {
      std::cout << "[" << v[0] << "," << v[1] << "," << v[2] << "]";
    };
    std::cout << std::setprecision(17)
              << "{\"schema\":\"nav6b-model-plan-v1\",\"model_only\":true,\"md_"
                 "authorized\":false,\"found\":"
              << (r.found ? "true" : "false") << ",\"reason\":\"" << r.reason
              << "\",\"initial_phase\":\""
              << transit_phase_name(r.initial_phase)
              << "\",\"vertex\":" << r.vertex << ",\"goal\":";
    vector(r.fixed_goal);
    std::cout << ",\"counters\":{\"proposals\":" << c.proposals
              << ",\"evaluations\":" << c.evaluations
              << ",\"bound_rejections\":" << c.bound_rejections
              << ",\"continuity_rejections\":" << c.continuity_rejections
              << ",\"ambiguity_rejections\":" << c.ambiguity_rejections
              << ",\"backtrack_rejections\":" << c.backtrack_rejections
              << ",\"corridor_rejections\":" << c.corridor_rejections
              << ",\"merit_rejections\":" << c.merit_rejections
              << ",\"clearance_rejections\":" << c.clearance_rejections
              << ",\"accepted_edges\":" << c.accepted_edges
              << ",\"swept_pair_checks\":" << c.swept_pair_checks
              << "},\"steps\":[";
    for (std::size_t i = 0; i < r.model_steps.size(); ++i) {
      if (i)
        std::cout << ",";
      const auto &step = r.model_steps[i];
      std::cout << "{\"phase\":\"" << transit_phase_name(step.phase)
                << "\",\"com\":";
      vector(step.com);
      std::cout << ",\"displacement\":";
      vector(step.displacement);
      std::cout << ",\"direction\":";
      vector(step.direction);
      std::cout << ",\"s\":" << step.progress_angstrom
                << ",\"d\":" << step.distance_angstrom
                << ",\"merit\":" << step.merit_angstrom2
                << ",\"travel\":" << step.travel_angstrom
                << ",\"turn_deg\":" << step.turn_deg
                << ",\"clearance\":" << step.sweep_clearance_angstrom << "}";
    }
    std::cout << "]}\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "NAV6B_INPUT_OR_GEOMETRY_ERROR: " << e.what() << "\n";
    return 2;
  }
}
