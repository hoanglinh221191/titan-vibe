param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_observation_schema_smoke')
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'pbctopo_observation_schema_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_observation_schema_smoke.exe'

@'
#include "pbctopo/assignment_set_report.h"
#include "pbctopo/image_assignment.h"
#include "pbctopo/observation_layout.h"
#include "pbctopo/observation_schema.h"
#include "pbctopo/report.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_observation_schema_smoke: " << message << "\n";
    std::exit(1);
  }
}

std::vector<std::string> split_lines(const std::string &text) {
  std::vector<std::string> lines;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line))
    lines.push_back(line);
  return lines;
}

std::vector<std::string> split_csv(const std::string &line) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (std::size_t idx = 0; idx < line.size(); ++idx) {
    const char c = line[idx];
    if (c == '"') {
      if (quoted && idx + 1 < line.size() && line[idx + 1] == '"') {
        field.push_back('"');
        ++idx;
      } else {
        quoted = !quoted;
      }
    } else if (c == ',' && !quoted) {
      fields.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  fields.push_back(field);
  return fields;
}

bool better(const PbctopoCutCandidate &lhs,
            const PbctopoCutCandidate &rhs) {
  return lhs.shift_norm2 < rhs.shift_norm2;
}

int main() {
  expect(VIBE_OBSERVATION_SCHEMA_VERSION == 7,
         "compiled observation schema version is not v7");
  expect(validate_vibe_observation_schema_version(7).has_value(),
         "schema v7 was rejected");
  const auto legacy = validate_vibe_observation_schema_version(6);
  expect(!legacy.has_value() &&
             legacy.error().find("unsupported VIBE observation schema") !=
                 std::string::npos,
         "legacy schema v6 rejection was not explicit");
  const auto future = validate_vibe_observation_schema_version(8);
  expect(!future.has_value() &&
             future.error().find("unsupported VIBE observation schema") !=
                 std::string::npos,
         "future schema rejection was not explicit");

  std::vector<PbctopoAtom> atoms(4);
  atoms[0].owner = 0;
  atoms[0].source_atom_id = 12;
  atoms[0].atomic_mass = 12.0;
  atoms[0].has_atomic_mass = true;
  atoms[1].owner = 0;
  atoms[1].source_atom_id = 10;
  atoms[1].atomic_mass = 10.0;
  atoms[1].has_atomic_mass = true;
  atoms[2].owner = 1;
  atoms[2].source_atom_id = 21;
  atoms[2].atomic_mass = 21.0;
  atoms[2].has_atomic_mass = true;
  atoms[3].owner = 1;
  atoms[3].source_atom_id = 20;
  atoms[3].atomic_mass = 20.0;
  atoms[3].has_atomic_mass = true;
  const std::vector<std::vector<std::size_t>> components{{0, 1}, {2, 3}};
  const std::vector<PbctopoLocalEdge> hard_edges{
      {0, 1, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond},
      {2, 3, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond}};
  const auto layout = build_pbctopo_component_layout_observation(
      components, atoms, hard_edges,
      PbctopoHardGraphSource::ExplicitTopology);
  expect(layout.valid && layout.topology_epoch_id != 0 &&
             layout.layout_signature != 0 && layout.components.size() == 2 &&
              layout.hard_edges.size() == 2 &&
              !layout.canonical_atom_universe_hash.empty(),
         "canonical component layout was not built");
  expect(layout.layout_signature == 16197695885957244026ULL,
         "golden component-layout identity changed");

  const std::vector<PbctopoAtom> reordered_atoms{
      atoms[3], atoms[2], atoms[1], atoms[0]};
  const std::vector<std::vector<std::size_t>> reordered_components{
      {3, 2}, {1, 0}};
  const std::vector<PbctopoLocalEdge> reordered_hard_edges{
      {3, 2, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond},
      {1, 0, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond}};
  const auto reordered_layout = build_pbctopo_component_layout_observation(
      reordered_components, reordered_atoms, reordered_hard_edges,
      PbctopoHardGraphSource::ExplicitTopology);
  expect(same_pbctopo_component_layout_observation(layout, reordered_layout),
         "exact layout membership changed under atom/component reorder");

  const std::vector<PbctopoLocalEdge> changed_hard_edges{
      {0, 1, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond},
      {1, 2, 1.0, 1.0, 4.0, PbctopoLocalEdgeKind::TopologyBond}};
  const auto changed_topology = build_pbctopo_component_layout_observation(
      components, atoms, changed_hard_edges,
      PbctopoHardGraphSource::ExplicitTopology);
  expect(changed_topology.valid &&
             changed_topology.layout_signature == layout.layout_signature &&
             changed_topology.topology_epoch_id != layout.topology_epoch_id,
         "hard-edge change did not create a new topology epoch");

  std::ostringstream layout_out;
  write_pbctopo_component_layout_header(layout_out);
  write_pbctopo_component_layout_rows(layout_out, layout);
  const auto layout_lines = split_lines(layout_out.str());
  expect(layout_lines.size() == 3,
         "component layout did not emit one row per component");
  const std::size_t layout_columns = split_csv(layout_lines.front()).size();
  for (std::size_t idx = 1; idx < layout_lines.size(); ++idx)
    expect(split_csv(layout_lines[idx]).size() == layout_columns,
           "component layout header/row column mismatch");
  expect(layout_out.str().find("0:10|0:12") != std::string::npos &&
             layout_out.str().find("1:20|1:21") != std::string::npos &&
             layout_lines.front().find("exact_atom_masses") !=
                 std::string::npos &&
             layout_out.str().find("10|12") != std::string::npos &&
             layout_out.str().find("20|21") != std::string::npos,
         "exact atom membership or mass is missing from layout sidecar");
  std::ostringstream hard_edge_out;
  write_pbctopo_hard_edge_header(hard_edge_out);
  write_pbctopo_hard_edge_rows(hard_edge_out, layout);
  const auto hard_edge_lines = split_lines(hard_edge_out.str());
  expect(hard_edge_lines.size() == 3 &&
             hard_edge_out.str().find("topology_bond") != std::string::npos,
         "hard-edge sidecar is incomplete");

  const std::vector<Int3> offsets_a{{0, 0, 0}, {1, -1, 0}};
  const std::vector<Int3> offsets_b{{3, 4, -2}, {4, 3, -2}};
  const auto assignment_a =
      make_pbctopo_image_assignment(offsets_a, components, atoms);
  const auto assignment_b =
      make_pbctopo_image_assignment(offsets_b, components, atoms);
  expect(same_pbctopo_image_assignment(assignment_a, assignment_b),
         "global gauge changed the canonical assignment");

  PbctopoFrameReport frame_report;
  frame_report.frame = 3;
  frame_report.time_ps = 1.5;
  frame_report.topology_epoch_id = layout.topology_epoch_id;
  frame_report.temporal_inference_experimental = true;
  std::ostringstream frame_out;
  write_pbctopo_frame_report_header(frame_out);
  write_pbctopo_frame_report_row(frame_out, frame_report);
  const auto frame_lines = split_lines(frame_out.str());
  expect(frame_lines.size() == 2 &&
             split_csv(frame_lines[0]).size() ==
                 split_csv(frame_lines[1]).size(),
         "frame report header/row column mismatch");
  const auto frame_fields = split_csv(frame_lines[1]);
  expect(frame_fields.size() > 13 && frame_fields[0] == "7" &&
               frame_fields[2] == "angstrom" && frame_fields[3] == "ps" &&
               frame_fields[5] == "IEEE754-binary64" &&
               frame_fields[6] == "sha256-128" &&
               frame_fields[7] == "sha256-128" &&
               frame_fields[8] == "sha256-256" &&
               frame_fields[9] == "bit_exact" &&
               frame_fields[13] == "experimental_viterbi_diagnostic",
         "frame report schema metadata is incomplete");

  PbctopoCutCandidate candidate;
  candidate.valid = true;
  candidate.shift_norm2 = 0.0;
  candidate.assignment = assignment_a;
  candidate.search_equivalence_known = true;
  candidate.equivalent_to_search_best = true;
  candidate.component_search_shell_radius = 2;
  candidate.search_domain_audit_known = true;
  candidate.search_domain_evidence_complete = true;
  candidate.bounded_domain_only = true;
  candidate.domain_completeness =
      PbctopoDomainCompleteness::BoundedExhaustive;
  candidate.hard_feasibility.evaluated = true;
  candidate.hard_feasibility.assignment_available = true;
  candidate.hard_feasibility.steric_audit_available = true;
  candidate.hard_feasibility.steric_decision_complete = true;
  candidate.hard_feasibility.steric_pair_enumeration_complete = true;
  candidate.hard_feasibility.steric_scan_complete = true;
  PbctopoIdentifiability identifiability;
  identifiability.identified = true;
  identifiability.uniqueness_search_exhaustive = true;
  identifiability.objective_uniqueness_known = true;
  identifiability.objective_unique_within_search_domain = true;
  identifiability.evidence_uniqueness_known = true;
  identifiability.evidence_unique_within_search_domain = true;
  identifiability.feasible_assignment_set_enumerated = true;
  identifiability.feasible_assignments_within_domain = 1;
  identifiability.equivalent_best_assignments_within_domain = 1;
  const std::vector<PbctopoCutCandidate> candidates{candidate};
  std::ostringstream assignment_out;
  write_pbctopo_assignment_set_header(assignment_out);
  write_pbctopo_assignment_set_rows(
      assignment_out, 3, 1.5, candidates, &candidates[0], identifiability,
      {1.0, 0.0, 0.0, 1.0, 0.0, 1.0}, better,
      layout.topology_epoch_id, 1);
  const auto assignment_lines = split_lines(assignment_out.str());
  expect(assignment_lines.size() == 2 &&
             split_csv(assignment_lines[0]).size() ==
                 split_csv(assignment_lines[1]).size(),
         "assignment report header/row column mismatch");
  expect(assignment_out.str().find("strong_bounded_anchor") !=
             std::string::npos,
         "bounded interior observation was not classified explicitly");

  std::cout << "PASS pbctopo_observation_schema_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\assignment_set_report.cpp `
        src\pbctopo\image_assignment.cpp `
        src\pbctopo\observation_binding.cpp `
        src\pbctopo\observation_layout.cpp `
        src\pbctopo\observation_schema.cpp `
        src\pbctopo\report.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO observation schema smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO observation schema smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
