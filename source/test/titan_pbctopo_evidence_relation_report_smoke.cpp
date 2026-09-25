#include "pbctopo/evidence_relation_report.h"

#include "pbctopo/evidence_support_policy.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <vector>

using namespace titan_pbctopo;

namespace {

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL pbctopo_evidence_relation_report_smoke: " << message
            << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

PbctopoImageAssignment make_assignment(
    std::uint64_t layout_signature, std::uint64_t signature,
    std::vector<std::uint64_t> component_ids,
    std::vector<Int3> component_offsets) {
  PbctopoImageAssignment assignment;
  assignment.layout_signature = layout_signature;
  assignment.signature = signature;
  assignment.component_ids = std::move(component_ids);
  assignment.component_offsets = std::move(component_offsets);
  assignment.valid = true;
  return assignment;
}

PbctopoCutCandidate make_candidate(PbctopoImageAssignment assignment) {
  PbctopoCutCandidate candidate;
  candidate.assignment = std::move(assignment);
  candidate.search_equivalence_known = true;
  candidate.equivalent_to_search_best = true;
  candidate.component_search_shell_radius = 2;
  candidate.search_domain_audit_known = true;
  candidate.search_domain_evidence_complete = true;
  candidate.domain_expansion_attempted = true;
  candidate.domain_expansion_exhausted = true;
  candidate.domain_completeness = PbctopoDomainCompleteness::GlobalCertified;
  candidate.hard_feasibility.evaluated = true;
  candidate.hard_feasibility.assignment_available = true;
  candidate.hard_feasibility.steric_audit_available = true;
  candidate.hard_feasibility.steric_decision_complete = true;
  candidate.hard_feasibility.steric_pair_enumeration_complete = true;
  candidate.hard_feasibility.steric_scan_complete = true;
  candidate.valid = true;
  return candidate;
}

PbctopoSelectedEvidenceRelationAudit make_relation(
    std::size_t component_a, std::size_t component_b,
    std::uint64_t component_id_a, std::uint64_t component_id_b,
    bool explicit_relation, std::span<const Int3> observed_shifts,
    const Int3 &selected_shift, std::size_t compatible_contacts = 0) {
  PbctopoSelectedEvidenceRelationAudit relation;
  relation.component_a = component_a;
  relation.component_b = component_b;
  relation.component_id_a = component_id_a;
  relation.component_id_b = component_id_b;
  relation.evaluated = true;
  relation.observed_hypothesis_relation = !observed_shifts.empty();
  relation.compatible_contact_count = compatible_contacts;
  relation.support = evaluate_pbctopo_selected_compatibility_support(
      explicit_relation, observed_shifts, selected_shift);
  return relation;
}

PbctopoSelectedEvidenceAudit make_audit(
    std::size_t component_count,
    std::vector<PbctopoSelectedEvidenceRelationAudit> relations,
    PbctopoEvidenceState state) {
  PbctopoSelectedEvidenceAudit audit;
  audit.state = state;
  audit.evaluated = true;
  audit.consistent = true;
  audit.evidence_component_count = component_count;
  audit.relation_audits = std::move(relations);
  for (const auto &relation : audit.relation_audits) {
    if (relation.support.active_at_threshold)
      ++audit.evidence_edge_count;
    if (relation.support.explicit_relation)
      ++audit.explicit_contact_evidence_edges;
    if (relation.support.compatibility_match)
      ++audit.scored_contact_evidence_edges;
  }

  std::vector<std::size_t> parent(component_count);
  for (std::size_t index = 0; index < component_count; ++index)
    parent[index] = index;
  const auto find_root = [&](std::size_t value) {
    while (parent[value] != value) {
      parent[value] = parent[parent[value]];
      value = parent[value];
    }
    return value;
  };
  for (const auto &relation : audit.relation_audits) {
    if (!relation.support.active_at_threshold)
      continue;
    const std::size_t left = find_root(relation.component_a);
    const std::size_t right = find_root(relation.component_b);
    if (left != right)
      parent[right] = left;
  }
  for (std::size_t index = 0; index < component_count; ++index) {
    if (find_root(index) == index)
      ++audit.evidence_connected_component_count;
  }
  audit.evidence_graph_connected =
      audit.evidence_connected_component_count == 1;
  audit.unconstrained_relative_components =
      audit.evidence_connected_component_count - 1;
  audit.relative_lattice_degrees_of_freedom =
      3 * audit.unconstrained_relative_components;
  return audit;
}

PbctopoFrameReport make_frame(std::size_t frame,
                              std::size_t trajectory_frame_index,
                              double time_ps, std::uint64_t topology_epoch_id,
                              std::uint64_t layout_signature,
                              std::uint64_t assignment_signature,
                              std::uint64_t hash_seed) {
  PbctopoFrameReport report;
  report.frame = frame;
  report.trajectory_frame_index = trajectory_frame_index;
  report.time_ps = time_ps;
  report.topology_epoch_id = topology_epoch_id;
  report.component_layout_signature = layout_signature;
  report.assignment_signature = assignment_signature;
  report.box_hash = {hash_seed + 1, hash_seed + 2};
  report.atom_selection_hash = {hash_seed + 3, hash_seed + 4};
  report.canonical_atom_universe_hash = {hash_seed + 5, hash_seed + 6};
  report.wrapped_coordinate_hash = {hash_seed + 7, hash_seed + 8};
  report.frame_binding_valid = true;
  return report;
}

std::string hash128_text(const PbctopoHash128 &hash) {
  constexpr char digits[] = "0123456789abcdef";
  std::string output(32, '0');
  const std::uint64_t words[] = {hash.hi, hash.lo};
  std::size_t cursor = 0;
  for (const auto word : words) {
    for (int shift = 60; shift >= 0; shift -= 4)
      output[cursor++] = digits[(word >> shift) & 0xfU];
  }
  return output;
}

void write_frame_report_row(std::ostream &out, const PbctopoFrameReport &frame,
                            const PbctopoSelectedEvidenceAudit &audit) {
  out << frame.frame << ',' << frame.trajectory_frame_index << ','
      << frame.time_ps << ',' << frame.topology_epoch_id << ','
      << frame.component_layout_signature << ',' << frame.assignment_signature
      << ',' << hash128_text(frame.box_hash) << ','
      << hash128_text(frame.atom_selection_hash) << ','
      << hash128_text(frame.canonical_atom_universe_hash) << ','
      << hash128_text(frame.wrapped_coordinate_hash) << ','
      << (audit.consistent ? 1 : 0) << ',' << audit.evidence_component_count
      << ',' << audit.evidence_edge_count << ','
      << audit.explicit_contact_evidence_edges << ','
      << audit.scored_contact_evidence_edges << ','
      << audit.evidence_connected_component_count << ','
      << (audit.evidence_graph_connected ? 1 : 0) << ','
      << audit.relative_lattice_degrees_of_freedom << '\n';
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    fail("expected one output-directory argument");
  const std::filesystem::path output_dir = argv[1];
  std::filesystem::create_directories(output_dir);

  const auto weak_assignment = make_assignment(
      9001, 101, {10, 20, 30}, {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
  const auto connected_assignment = make_assignment(
      9001, 202, {10, 20, 30}, {{0, 0, 0}, {0, 0, 0}, {1, 0, 0}});
  const auto singleton_assignment =
      make_assignment(9002, 303, {40}, {{0, 0, 0}});

  const std::vector<Int3> compatible_shift{{0, 0, 0}};
  const std::vector<Int3> incompatible_shift{{1, 0, 0}};
  const Int3 selected_shift{0, 0, 0};
  const auto weak_audit = make_audit(
      3,
      {make_relation(0, 1, 10, 20, false, compatible_shift, selected_shift,
                     3),
       make_relation(1, 2, 20, 30, false, incompatible_shift, selected_shift)},
      PbctopoEvidenceState::WeakAmbiguous);
  const auto connected_audit = make_audit(
      3,
      {make_relation(0, 1, 10, 20, false, compatible_shift, selected_shift,
                     3),
       make_relation(1, 2, 20, 30, true, std::span<const Int3>{},
                     selected_shift)},
      PbctopoEvidenceState::SupportedConnected);
  const auto singleton_audit = make_audit(
      1, {}, PbctopoEvidenceState::SupportedConnected);

  expect(!weak_audit.evidence_graph_connected &&
             weak_audit.evidence_connected_component_count == 2 &&
             weak_audit.relative_lattice_degrees_of_freedom == 3,
         "weak fixture aggregate is wrong");
  expect(connected_audit.evidence_graph_connected &&
             connected_audit.evidence_edge_count == 2,
         "connected fixture aggregate is wrong");

  std::vector<PbctopoCandidateEvidenceAuditSnapshot> snapshots;
  auto inserted = upsert_pbctopo_candidate_evidence_audit_snapshot(
      snapshots, weak_assignment, weak_audit);
  expect(inserted && *inserted, "weak snapshot was not inserted");
  inserted = upsert_pbctopo_candidate_evidence_audit_snapshot(
      snapshots, weak_assignment, weak_audit);
  expect(inserted && !*inserted, "identical snapshot was not deduplicated");
  inserted = upsert_pbctopo_candidate_evidence_audit_snapshot(
      snapshots, connected_assignment, connected_audit);
  expect(inserted && *inserted, "connected snapshot was not inserted");
  inserted = upsert_pbctopo_candidate_evidence_audit_snapshot(
      snapshots, singleton_assignment, singleton_audit);
  expect(inserted && *inserted, "singleton snapshot was not inserted");

  auto collision_assignment = weak_assignment;
  collision_assignment.component_offsets[2] = {2, 0, 0};
  const auto collision = upsert_pbctopo_candidate_evidence_audit_snapshot(
      snapshots, collision_assignment, weak_audit);
  expect(!collision, "assignment-signature collision was accepted");

  const auto weak_candidate = make_candidate(weak_assignment);
  const auto connected_candidate = make_candidate(connected_assignment);
  const auto singleton_candidate = make_candidate(singleton_assignment);

  std::ofstream summary(output_dir / "pbctopo_candidate_evidence_summary.csv",
                        std::ios::binary | std::ios::trunc);
  std::ofstream relations(
      output_dir / "pbctopo_candidate_evidence_relations.csv",
      std::ios::binary | std::ios::trunc);
  std::ofstream frames(output_dir / "pbctopo_frame_report.csv",
                       std::ios::binary | std::ios::trunc);
  expect(summary.is_open() && relations.is_open() && frames.is_open(),
         "could not open fixture outputs");
  write_pbctopo_candidate_evidence_summary_header(summary);
  write_pbctopo_candidate_evidence_relation_header(relations);
  frames << "frame,trajectory_frame_index,time_ps,topology_epoch_id,"
            "component_layout_signature,assignment_signature,box_hash,"
            "atom_selection_hash,canonical_atom_universe_hash,"
            "wrapped_coordinate_hash,evidence_consistent,"
            "evidence_component_count,evidence_edge_count,"
            "explicit_contact_evidence_edges,scored_contact_evidence_edges,"
            "evidence_connected_component_count,evidence_graph_connected,"
            "relative_lattice_degrees_of_freedom\n";

  const auto frame1 = make_frame(7, 10, 1.5, 77, 9001, 202, 100);
  const std::vector<PbctopoCutCandidate> retained1{
      connected_candidate, weak_candidate, weak_candidate};
  const auto stats1 = write_pbctopo_candidate_evidence_frame(
      summary, relations, frame1, retained1, &connected_candidate, true,
      snapshots);
  expect(stats1 && stats1->candidate_groups == 2 &&
             stats1->relation_rows == 4 &&
             stats1->duplicate_candidates_skipped == 1,
         stats1 ? "unexpected frame-1 write statistics" : stats1.error());
  write_frame_report_row(frames, frame1, connected_audit);

  const auto frame2 = make_frame(8, 11, 2.5, 77, 9002, 303, 200);
  const std::vector<PbctopoCutCandidate> retained2{singleton_candidate};
  const auto stats2 = write_pbctopo_candidate_evidence_frame(
      summary, relations, frame2, retained2, &singleton_candidate, false,
      snapshots);
  expect(stats2 && stats2->candidate_groups == 1 &&
             stats2->relation_rows == 0 &&
             stats2->duplicate_candidates_skipped == 0,
         stats2 ? "unexpected frame-2 write statistics" : stats2.error());
  write_frame_report_row(frames, frame2, singleton_audit);

  std::ostringstream atomic_summary;
  std::ostringstream atomic_relations;
  auto missing_candidate = connected_candidate;
  missing_candidate.assignment.signature = 404;
  missing_candidate.assignment.component_offsets[2] = {4, 0, 0};
  const std::vector<PbctopoCutCandidate> missing{weak_candidate,
                                                 missing_candidate};
  const auto rejected = write_pbctopo_candidate_evidence_frame(
      atomic_summary, atomic_relations, frame1, missing, nullptr, false,
      snapshots);
  expect(!rejected && atomic_summary.str().empty() &&
             atomic_relations.str().empty(),
         "failed frame export was not atomic");

  summary.close();
  relations.close();
  frames.close();
  expect(summary.good() && relations.good() && frames.good(),
         "fixture output flush failed");
  std::cout << "PASS pbctopo_evidence_relation_report_smoke "
               "candidate_groups=3 relation_rows=4\n";
}
