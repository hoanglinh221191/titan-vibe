param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_image_assignment_smoke')
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
#include "pbctopo/certification.h"
#include "pbctopo/image_assignment.h"
#include "pbctopo/report.h"
#include "pbctopo/steric_clash.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_image_assignment_smoke: " << message << "\n";
    std::exit(1);
  }
}

std::vector<std::string> split_csv_simple(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ','))
    fields.push_back(field);
  return fields;
}

} // namespace

int main() {
  using namespace titan_pbctopo;

  expect(std::string_view(pbctopo_certificate_scope_name(
             PbctopoCertificateScope::ImageIntercomponentSteric)) ==
             "image_intercomponent_steric" &&
             std::string_view(pbctopo_certificate_scope_name(
                 PbctopoCertificateScope::
                     ImageIntercomponentStericInterface)) ==
                 "image_intercomponent_steric_interface",
         "inter-component certificate scope names are unstable");
  expect(std::string_view(pbctopo_certificate_graph_source_name(
             PbctopoCertificateGraphSource::None)) == "none" &&
             std::string_view(pbctopo_certificate_graph_source_name(
                 PbctopoCertificateGraphSource::GromacsTopology)) ==
                 "gromacs_topology" &&
             std::string_view(pbctopo_certificate_graph_source_name(
                 PbctopoCertificateGraphSource::MetadataChainResidue)) ==
                 "metadata_chain_residue" &&
             std::string_view(pbctopo_certificate_graph_source_name(
                 PbctopoCertificateGraphSource::GeometryCutoff)) ==
                 "geometry_cutoff" &&
             std::string_view(pbctopo_certificate_graph_source_name(
                 PbctopoCertificateGraphSource::Mixed)) == "mixed",
         "certificate graph source names are unstable");

  PbctopoFrameReport coverage_report;
  coverage_report.hard_graph_source =
      PbctopoHardGraphSource::ExplicitTopology;
  coverage_report.certificate_graph_source =
      PbctopoCertificateGraphSource::GeometryCutoff;
  coverage_report.hard_feasibility.hard_connectivity_edges_available = 10;
  coverage_report.hard_feasibility.topology_edges_available = 4;
  coverage_report.hard_feasibility.topology_edges_checked = 3;
  coverage_report.hard_feasibility.metadata_edges_available = 6;
  coverage_report.hard_feasibility.metadata_edges_checked = 5;
  std::ostringstream header_stream;
  std::ostringstream row_stream;
  write_pbctopo_frame_report_header(header_stream);
  write_pbctopo_frame_report_row(row_stream, coverage_report);
  const auto header_fields = split_csv_simple(header_stream.str());
  const auto row_fields = split_csv_simple(row_stream.str());
  const auto report_field = [&](std::string_view name) -> std::string {
    const auto found = std::find(header_fields.begin(), header_fields.end(),
                                 std::string(name));
    expect(found != header_fields.end(), "report field is missing");
    const std::size_t index =
        static_cast<std::size_t>(found - header_fields.begin());
    expect(index < row_fields.size(), "report row is shorter than its header");
    return row_fields[index];
  };
  expect(report_field("certificate_graph_source") == "gromacs_topology",
         "geometry search evidence leaked into certificate graph source");
  expect(std::fabs(std::stod(report_field("hard_topology_edge_fraction")) -
                       0.75) <
             1.0e-12,
         "topology coverage used the combined hard-edge denominator");
  expect(std::fabs(std::stod(report_field("hard_metadata_edge_fraction")) -
                       (5.0 / 6.0)) <
             1.0e-6,
         "metadata coverage fraction used the wrong denominator");

  const double matrix[3][3] = {{10.0, 0.0, 0.0},
                               {0.0, 10.0, 0.0},
                               {0.0, 0.0, 10.0}};
  const double inverse[3][3] = {{0.1, 0.0, 0.0},
                                {0.0, 0.1, 0.0},
                                {0.0, 0.0, 0.1}};
  const PbctopoLatticeMetric lattice(matrix, inverse);

  pbctopo_reset_lattice_query_audit();
  const auto optional_failure = lattice.nearest_image_scaled(
      std::numeric_limits<double>::infinity(), 0.0, 0.0);
  expect(!optional_failure.valid,
         "invalid optional lattice query did not fail closed");
  const auto valid_query_checkpoint = pbctopo_lattice_query_audit();
  const auto valid_after_optional = lattice.nearest_image_scaled(0.2, 0.0, 0.0);
  const auto valid_query_delta = pbctopo_lattice_query_audit_delta(
      valid_query_checkpoint, pbctopo_lattice_query_audit());
  expect(valid_after_optional.valid && valid_query_delta.queries == 1 &&
             valid_query_delta.failures == 0,
         "prior optional lattice failure poisoned a later valid query scope");

  std::int64_t rounded_image = 0;
  expect(pbctopo_checked_round_image(1099511627776.25, rounded_image) &&
             rounded_image == INT64_C(1099511627776),
         "large int64 image offset was not rounded safely");
  expect(!pbctopo_checked_round_image(1.0e30, rounded_image),
         "out-of-range image offset was accepted");

  std::vector<PbctopoAtom> wrapped_atoms(2);
  wrapped_atoms[0].x = 9.9;
  wrapped_atoms[0].sx = 0.99;
  wrapped_atoms[0].owner = 0;
  wrapped_atoms[1].x = 0.1;
  wrapped_atoms[1].sx = 0.01;
  wrapped_atoms[1].owner = 1;
  const std::vector<PbctopoLocalEdge> topology_edges{
      {0, 1, 0.04, 25.0, 1.0, PbctopoLocalEdgeKind::TopologyBond},
      {0, 1, 0.04, 25.0, 1.0, PbctopoLocalEdgeKind::GeometryHeuristic}};
  const std::vector<Vec3> intact_points{{9.9, 0.0, 0.0},
                                        {10.1, 0.0, 0.0}};
  const auto intact_equations = evaluate_pbctopo_image_equations(
      wrapped_atoms, intact_points, topology_edges, lattice);
  const auto compact_geometry =
      classify_pbctopo_final_geometry(intact_points, lattice);
  expect(compact_geometry.evaluated &&
             compact_geometry.geometry_class ==
                 PbctopoGeometryClass::CompactUnitBox,
         "compact final geometry was not classified as unit-box compact");
  const std::vector<Vec3> winding_points{{0.0, 0.0, 0.0},
                                          {10.5, 0.0, 0.0}};
  const auto winding_geometry =
      classify_pbctopo_final_geometry(winding_points, lattice);
  expect(winding_geometry.evaluated &&
             winding_geometry.geometry_class ==
                 PbctopoGeometryClass::Extended &&
             !winding_geometry.fits_one_fundamental_cell,
         "extended final geometry was not classified explicitly");
  std::vector<PbctopoAtom> cycle_atoms(3);
  cycle_atoms[0].x = 1.0;
  cycle_atoms[0].sx = 0.1;
  cycle_atoms[1].x = 4.0;
  cycle_atoms[1].sx = 0.4;
  cycle_atoms[2].x = 8.0;
  cycle_atoms[2].sx = 0.8;
  const std::vector<Vec3> cycle_points{{1.0, 0.0, 0.0},
                                        {4.0, 0.0, 0.0},
                                        {8.0, 0.0, 0.0}};
  const std::vector<PbctopoLocalEdge> cycle_edges{
      {0, 1, 9.0, 1.0 / 9.0, 100.0,
       PbctopoLocalEdgeKind::TopologyBond},
      {1, 2, 16.0, 1.0 / 16.0, 100.0,
       PbctopoLocalEdgeKind::TopologyBond},
      {2, 0, 9.0, 1.0 / 9.0, 100.0,
       PbctopoLocalEdgeKind::TopologyBond}};
  const auto cycle_equations = evaluate_pbctopo_image_equations(
      cycle_atoms, cycle_points, cycle_edges, lattice);
  const auto topological_winding_geometry =
      classify_pbctopo_final_geometry(cycle_points, lattice, 1.0e-6,
                                      &cycle_equations);
  expect(cycle_equations.topological_winding_rank == 1 &&
             topological_winding_geometry.geometry_class ==
                 PbctopoGeometryClass::Winding &&
             topological_winding_geometry.fits_one_fundamental_cell,
         "topological winding was conflated with physical extension");
  expect(intact_equations.input_valid,
         "intact image equations were not evaluated");
  expect(intact_equations.topology_edges_checked == 1,
         "non-topology edge entered the hard equation set");
  expect(intact_equations.topology_edge_residuals == 0,
         "valid seam bond has a topology image residual");
  expect(intact_equations.atom_image_residuals == 0,
         "integer atom images were rejected");
  const auto intact_required = evaluate_pbctopo_required_contacts(
      wrapped_atoms, intact_points, topology_edges);
  expect(intact_required.input_valid,
         "hard-required contacts were not evaluated");
  expect(intact_required.contacts_checked == 1,
         "soft-observed geometry edge entered the hard contact set");
  expect(intact_required.contacts_lost == 0,
         "intact topology contact was reported lost");

  std::vector<PbctopoAtom> half_cell_atoms(2);
  half_cell_atoms[0].x = 0.0;
  half_cell_atoms[0].sx = 0.0;
  half_cell_atoms[1].x = 5.0;
  half_cell_atoms[1].sx = 0.5;
  const std::vector<PbctopoLocalEdge> half_cell_edge{
      {0, 1, 25.0, 1.0 / 25.0, 100.0,
       PbctopoLocalEdgeKind::TopologyBond}};
  const std::vector<Vec3> half_cell_image_zero{{0.0, 0.0, 0.0},
                                                {5.0, 0.0, 0.0}};
  const std::vector<Vec3> half_cell_image_minus_one{{0.0, 0.0, 0.0},
                                                     {-5.0, 0.0, 0.0}};
  const auto half_zero = evaluate_pbctopo_image_equations(
      half_cell_atoms, half_cell_image_zero, half_cell_edge, lattice);
  const auto half_minus = evaluate_pbctopo_image_equations(
      half_cell_atoms, half_cell_image_minus_one, half_cell_edge, lattice);
  expect(half_zero.input_valid && half_minus.input_valid &&
             half_zero.ambiguous_topology_edges == 1 &&
             half_minus.ambiguous_topology_edges == 1,
         "half-cell topology ambiguity was not reported");
  expect(half_zero.topology_edge_residuals == 0 &&
             half_minus.topology_edge_residuals == 0,
         "an equivalent half-cell nearest image was rejected");

  const std::vector<std::vector<std::size_t>> components{{0, 1, 2}, {3, 4}};
  const std::vector<Int3> offsets_a{{3, -2, 1}, {4, -2, 1}};
  const std::vector<Int3> offsets_global_shifted{{-5, 7, 9}, {-4, 7, 9}};
  const std::vector<Int3> offsets_relative_change{{3, -2, 1}, {5, -2, 1}};

  const auto assignment_a =
      make_pbctopo_image_assignment(offsets_a, components);
  const auto assignment_global =
      make_pbctopo_image_assignment(offsets_global_shifted, components);
  const auto assignment_relative =
      make_pbctopo_image_assignment(offsets_relative_change, components);

  expect(assignment_a.valid, "baseline assignment is invalid");
  expect(assignment_a.signature != 0, "baseline signature is zero");
  expect(same_pbctopo_image_assignment(assignment_a, assignment_global),
         "global lattice translation changed canonical assignment");
  expect(!same_pbctopo_image_assignment(assignment_a, assignment_relative),
         "relative component shift did not change assignment");
  const std::vector<std::vector<std::size_t>> reordered_components{
      {3, 4}, {0, 1, 2}};
  const std::vector<Int3> reordered_offsets{{4, -2, 1}, {3, -2, 1}};
  const auto reordered_assignment = make_pbctopo_image_assignment(
      reordered_offsets, reordered_components);
  expect(same_pbctopo_image_assignment(assignment_a, reordered_assignment) &&
             assignment_a.layout_signature ==
                 reordered_assignment.layout_signature,
         "component reorder changed semantic assignment identity");

  std::vector<PbctopoAtom> source_atoms(5);
  for (std::size_t atom_idx = 0; atom_idx < source_atoms.size(); ++atom_idx) {
    source_atoms[atom_idx].owner = atom_idx < 3 ? 0 : 1;
    source_atoms[atom_idx].source_atom_id =
        atom_idx < 3 ? 100 + atom_idx : 200 + atom_idx - 3;
  }
  const auto stable_assignment = make_pbctopo_image_assignment(
      offsets_a, components, source_atoms);
  const std::vector<PbctopoAtom> permuted_source_atoms{
      source_atoms[3], source_atoms[4], source_atoms[0], source_atoms[1],
      source_atoms[2]};
  const std::vector<std::vector<std::size_t>> permuted_components{
      {0, 1}, {2, 3, 4}};
  const std::vector<Int3> permuted_offsets{{4, -2, 1}, {3, -2, 1}};
  const auto stable_permuted_assignment = make_pbctopo_image_assignment(
      permuted_offsets, permuted_components, permuted_source_atoms);
  expect(stable_assignment.valid && stable_permuted_assignment.valid &&
             same_pbctopo_image_assignment(stable_assignment,
                                            stable_permuted_assignment),
         "local atom reindexing changed owner/source-based assignment identity");
  auto changed_source_atoms = source_atoms;
  ++changed_source_atoms[0].source_atom_id;
  const auto changed_source_assignment = make_pbctopo_image_assignment(
      offsets_a, components, changed_source_atoms);
  expect(changed_source_assignment.valid &&
             changed_source_assignment.layout_signature !=
                 stable_assignment.layout_signature,
         "source atom identity change reused a component layout signature");

  const std::vector<std::vector<std::size_t>> different_membership{
      {0, 1, 3}, {2, 4}};
  const auto different_layout =
      make_pbctopo_image_assignment(offsets_a, different_membership);
  expect(!same_pbctopo_image_assignment(assignment_a, different_layout) &&
             assignment_a.layout_signature != different_layout.layout_signature,
         "different component membership reused an assignment signature");

  const std::vector<Int3> offsets_large{
      {INT64_C(1099511627776), INT64_C(-1099511627776), 3},
      {INT64_C(1099511627777), INT64_C(-1099511627776), 3}};
  const auto assignment_large =
      make_pbctopo_image_assignment(offsets_large, components);
  expect(assignment_large.valid &&
             same_pbctopo_image_assignment(assignment_a, assignment_large),
         "large int64 gauge changed a valid relative assignment");
  const std::vector<Int3> offsets_overflow{
      {std::numeric_limits<std::int64_t>::min(), 0, 0},
      {std::numeric_limits<std::int64_t>::max(), 0, 0}};
  expect(!make_pbctopo_image_assignment(offsets_overflow, components).valid,
         "overflowing canonical assignment was accepted");
  const Int3 max_offset{std::numeric_limits<std::int64_t>::max(), 0, 0};
  const Int3 min_offset{std::numeric_limits<std::int64_t>::min(), 0, 0};
  Int3 checked_result{};
  expect(!checked_int3_add(max_offset, Int3{1, 0, 0}, checked_result) &&
             !checked_int3_subtract(min_offset, Int3{1, 0, 0},
                                    checked_result) &&
             !checked_int3_negate(min_offset, checked_result),
         "overflowing Int3 arithmetic did not fail closed");

  const std::vector<Int3> wrong_size{{0, 0, 0}};
  expect(!make_pbctopo_image_assignment(wrong_size, components).valid,
         "size mismatch produced a valid assignment");

  std::vector<PbctopoAtom> component_reference_atoms = wrapped_atoms;
  component_reference_atoms[1].x = 10.1;
  const std::vector<std::size_t> intact_atom_components{0, 0};
  const std::vector<std::vector<std::size_t>> intact_components{{0, 1}};
  const std::vector<Int3> intact_offsets{{0, 0, 0}};
  const auto intact_assignment =
      make_pbctopo_image_assignment(intact_offsets, intact_components);
  const auto intact_assignment_points =
      evaluate_pbctopo_assignment_point_consistency(
          component_reference_atoms, intact_points, intact_atom_components,
          intact_assignment, lattice);
  expect(intact_assignment_points.input_valid &&
             intact_assignment_points.point_residuals == 0,
         "intact component assignment did not reproduce every atom");
  const std::unordered_set<std::uint64_t> no_exclusions;
  const auto intact_sterics = evaluate_pbctopo_final_sterics(
      wrapped_atoms, intact_points, intact_atom_components, no_exclusions,
      1.0);
  expect(intact_sterics.input_valid && intact_sterics.hard_clashes == 0,
         "same-component geometry entered final inter-component sterics");

  PbctopoCutCandidate safe;
  safe.valid = true;
  safe.broken_edges = 99;
  safe.ultra_clash = false;
  safe.assignment = intact_assignment;
  expect(certify_pbctopo_candidate(safe, intact_equations,
                                    intact_assignment_points, intact_sterics,
                                    intact_required),
         "hard-feasible assignment was not certified");
  expect(safe.hard_feasibility.topology_edges_checked == 1,
         "hard certificate did not report checked topology edges");
  expect(safe.hard_feasibility.hard_connectivity_edges_available == 1 &&
             safe.hard_feasibility.topology_edges_available == 1 &&
             safe.hard_feasibility.metadata_edges_available == 0,
         "hard certificate did not report topology coverage");
  expect(safe.hard_feasibility.certified(),
         "certified tuple does not report certified state");

  const std::vector<PbctopoLocalEdge> metadata_edges{
      {0, 1, 0.04, 25.0, 1.0,
       PbctopoLocalEdgeKind::MetadataSkeleton}};
  const auto metadata_equations = evaluate_pbctopo_image_equations(
      wrapped_atoms, intact_points, metadata_edges, lattice);
  const auto metadata_required = evaluate_pbctopo_required_contacts(
      wrapped_atoms, intact_points, metadata_edges);
  PbctopoCutCandidate metadata_only = safe;
  metadata_only.hard_feasibility.hard_connectivity_edges_available = 1;
  metadata_only.hard_feasibility.topology_edges_available = 0;
  metadata_only.hard_feasibility.metadata_edges_available = 1;
  expect(certify_pbctopo_candidate(
             metadata_only, metadata_equations, intact_assignment_points,
             intact_sterics, metadata_required),
         "metadata-only hard connectivity was made a topology requirement");
  expect(metadata_only.hard_feasibility.topology_edges_checked == 0 &&
             metadata_only.hard_feasibility.topology_edges_available == 0 &&
             metadata_only.hard_feasibility.hard_connectivity_edges_available ==
                 1 &&
             metadata_only.hard_feasibility.metadata_edges_available == 1,
         "metadata-only certificate coverage was not retained");

  PbctopoCutCandidate inferred_metadata_availability = safe;
  inferred_metadata_availability.hard_feasibility
      .hard_connectivity_edges_available = 0;
  inferred_metadata_availability.hard_feasibility.topology_edges_available = 0;
  inferred_metadata_availability.hard_feasibility.metadata_edges_available = 0;
  expect(certify_pbctopo_candidate(
             inferred_metadata_availability, metadata_equations,
             intact_assignment_points, intact_sterics, metadata_required),
         "audited metadata edges did not infer their availability");
  expect(inferred_metadata_availability.hard_feasibility
                 .metadata_edges_available == 1 &&
             inferred_metadata_availability.hard_feasibility
                 .hard_connectivity_edges_available == 1,
         "metadata availability was not recovered from the audit");

  const std::vector<Vec3> split_points{{9.9, 0.0, 0.0},
                                       {0.1, 0.0, 0.0}};
  const auto split_equations = evaluate_pbctopo_image_equations(
      wrapped_atoms, split_points, topology_edges, lattice);
  const auto split_required = evaluate_pbctopo_required_contacts(
      wrapped_atoms, split_points, topology_edges);
  expect(split_equations.topology_edge_residuals == 1,
         "split seam bond did not produce a direct image residual");
  expect(split_required.contacts_lost == 1,
         "split hard-required topology contact was not reported lost");
  PbctopoCutCandidate broken = safe;
  const auto split_assignment_points =
      evaluate_pbctopo_assignment_point_consistency(
          component_reference_atoms, split_points, intact_atom_components,
          intact_assignment, lattice);
  expect(!certify_pbctopo_candidate(broken, split_equations,
                                    split_assignment_points, intact_sterics,
                                    split_required),
         "topology image residual was certified");

  const std::vector<Vec3> noninteger_points{{9.9, 0.0, 0.0},
                                            {12.6, 0.0, 0.0}};
  const auto noninteger_equations = evaluate_pbctopo_image_equations(
      wrapped_atoms, noninteger_points, topology_edges, lattice);
  expect(noninteger_equations.atom_image_residuals == 1,
         "non-lattice atom displacement was not rejected");
  const std::vector<PbctopoAtom> huge_wrapped_atoms(1);
  const std::vector<Vec3> huge_points{{1.0e30, 0.0, 0.0}};
  const auto huge_image_equations = evaluate_pbctopo_image_equations(
      huge_wrapped_atoms, huge_points, {}, lattice);
  expect(huge_image_equations.atom_image_residuals == 1,
         "out-of-range atom image conversion was not rejected safely");

  std::vector<PbctopoAtom> clash_atoms(2);
  clash_atoms[0].x = 1.0;
  clash_atoms[1].x = 1.5;
  for (auto &atom : clash_atoms) {
    atom.is_heavy = true;
    atom.has_vdw_radius = true;
    atom.steric_radius_source =
        PbctopoStericRadiusSource::ForceFieldParameter;
    atom.vdw_radius = 1.0;
    atom.steric_radius = 1.0;
  }
  const std::vector<Vec3> clash_points{{1.0, 0.0, 0.0},
                                       {1.5, 0.0, 0.0}};
  const std::vector<std::size_t> split_components_by_atom{0, 1};
  const auto clash_sterics = evaluate_pbctopo_final_sterics(
      clash_atoms, clash_points, split_components_by_atom, no_exclusions, 1.0);
  expect(clash_sterics.input_valid && clash_sterics.hard_clashes == 1,
         "final heavy-atom steric audit missed an overlap");
  const std::vector<std::vector<std::size_t>> clash_components{{0}, {1}};
  const auto clash_assignment = make_pbctopo_image_assignment(
      std::vector<Int3>{{0, 0, 0}, {0, 0, 0}}, clash_components);
  const auto clash_assignment_points =
      evaluate_pbctopo_assignment_point_consistency(
          clash_atoms, clash_points, split_components_by_atom,
          clash_assignment, lattice);
  const auto clash_equations = evaluate_pbctopo_image_equations(
      clash_atoms, clash_points, {}, lattice);
  const auto clash_required = evaluate_pbctopo_required_contacts(
      clash_atoms, clash_points, {});
  PbctopoCutCandidate clashing;
  clashing.valid = true;
  clashing.ultra_clash = false;
  clashing.assignment = clash_assignment;
  expect(!certify_pbctopo_candidate(
             clashing, clash_equations, clash_assignment_points,
             clash_sterics, clash_required),
         "independent final steric clash was certified");

  clash_atoms[1].is_heavy = false;
  const auto hydrogen_sterics = evaluate_pbctopo_final_sterics(
      clash_atoms, clash_points, split_components_by_atom, no_exclusions, 1.0);
  expect(hydrogen_sterics.input_valid &&
             hydrogen_sterics.hard_clashes == 1,
         "all-atom certificate missed a hydrogen-only overlap");
  expect(hydrogen_sterics.force_field_radius_atoms == 2,
         "steric-radius provenance was not retained");

  clash_atoms[1].steric_radius_source =
      PbctopoStericRadiusSource::ExplicitZeroLennardJones;
  clash_atoms[1].vdw_radius = 0.0;
  clash_atoms[1].steric_radius = 0.0;
  const auto zero_lj_sterics = evaluate_pbctopo_final_sterics(
      clash_atoms, clash_points, split_components_by_atom, no_exclusions, 1.0);
  expect(zero_lj_sterics.input_valid &&
             zero_lj_sterics.hard_clashes == 0 &&
             zero_lj_sterics.explicit_zero_lj_radius_atoms == 1,
         "explicit zero-LJ site did not follow no-hard-core policy");
  expect(pbctopo_pair_soft_overlap_penalty(
             0.01, clash_atoms[0], clash_atoms[1], 1.0) == 0.0,
         "explicit zero-LJ site contributed a soft LJ overlap penalty");

  std::vector<PbctopoAtom> metadata_atoms(2);
  metadata_atoms[0] = clash_atoms[0];
  metadata_atoms[1] = clash_atoms[0];
  metadata_atoms[0].x = 1.0;
  metadata_atoms[1].x = 5.0;
  for (auto &atom : metadata_atoms) {
    atom.has_topology_key = true;
    atom.has_topology_graph = false;
  }
  const std::vector<Vec3> metadata_bad_points{{1.0, 0.0, 0.0},
                                               {1.5, 0.0, 0.0}};
  const std::vector<std::size_t> metadata_component_by_atom{0, 0};
  const auto metadata_sterics = evaluate_pbctopo_final_sterics(
      metadata_atoms, metadata_bad_points, metadata_component_by_atom,
      no_exclusions, 1.0, nullptr, &lattice);
  expect(metadata_sterics.input_valid &&
             metadata_sterics.metadata_internal_pairs_checked == 1 &&
             metadata_sterics.metadata_internal_clashes == 1 &&
             metadata_sterics.hard_clashes == 1,
         "metadata unwrap-induced internal overlap was not rejected");

  const std::vector<Vec3> inconsistent_points{{1.0, 0.0, 0.0},
                                               {11.5, 0.0, 0.0}};
  const auto inconsistent_assignment_points =
      evaluate_pbctopo_assignment_point_consistency(
          clash_atoms, inconsistent_points, split_components_by_atom,
          clash_assignment, lattice);
  expect(inconsistent_assignment_points.point_residuals == 1,
         "non-root atom assignment mismatch was not detected");
  const auto inconsistent_equations = evaluate_pbctopo_image_equations(
      clash_atoms, inconsistent_points, {}, lattice);
  const auto inconsistent_sterics = evaluate_pbctopo_final_sterics(
      clash_atoms, inconsistent_points, split_components_by_atom,
      no_exclusions, 1.0);
  PbctopoCutCandidate inconsistent = clashing;
  inconsistent.ultra_clash = false;
  expect(!certify_pbctopo_candidate(
             inconsistent, inconsistent_equations,
             inconsistent_assignment_points, inconsistent_sterics,
             clash_required),
         "assignment-to-point mismatch was certified");

  PbctopoCutCandidate cyclic = safe;
  auto cyclic_equations = intact_equations;
  cyclic_equations.topology_cycle_residuals = 1;
  expect(!certify_pbctopo_candidate(cyclic, cyclic_equations,
                                    intact_assignment_points, intact_sterics,
                                    intact_required),
         "topology cycle residual was certified");

  auto unsupported_winding_equations = cycle_equations;
  unsupported_winding_equations.topology_edge_residuals = 0;
  unsupported_winding_equations.topology_cycle_residuals = 0;
  PbctopoCutCandidate unsupported_winding = safe;
  expect(!certify_pbctopo_candidate(
             unsupported_winding, unsupported_winding_equations,
             intact_assignment_points, intact_sterics, intact_required),
         "unsupported VIBE winding was certified");
  expect(unsupported_winding.hard_feasibility.unsupported_winding_rank == 1,
         "topological winding rank was not exposed by hard feasibility");

  PbctopoCutCandidate missing_assignment = safe;
  missing_assignment.assignment = {};
  expect(!certify_pbctopo_candidate(
             missing_assignment, intact_equations, intact_assignment_points,
             intact_sterics, intact_required),
         "missing image assignment was certified");

  std::cout << "PASS pbctopo_image_assignment_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\image_assignment.cpp `
        src\pbctopo\certification.cpp `
        src\pbctopo\observation_binding.cpp `
        src\pbctopo\report.cpp `
        src\pbctopo\steric_clash.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO image assignment smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO image assignment smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
