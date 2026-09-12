param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_temporal_assignment_replay_smoke')
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
#include "pbctopo/temporal_replay.h"
#include "pbctopo/selected_report.h"
#include "pbctopo/image_assignment.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_temporal_assignment_replay_smoke: " << message
              << "\n";
    std::exit(1);
  }
}

int main() {
  using namespace titan_pbctopo;
  PbctopoContinuityPair edge_specific_cutoff;
  edge_specific_cutoff.kind = PbctopoContinuityKind::TopologyHard;
  edge_specific_cutoff.hard_cutoff_d2 = 4.0;
  expect(std::fabs(pbctopo_continuity_hard_margin_A(
                       edge_specific_cutoff, 2.25, 16.0) -
                   0.5) <
             1.0e-12,
         "continuity margin ignored the edge-specific hard cutoff");
  box Box(10.0, 12.0, 14.0, 90.0, 90.0, 90.0);
  const std::vector<Vec3> applied{{1.0, 2.0, 3.0}, {2.0, 2.0, 3.0},
                                  {8.0, 4.0, 5.0}};
  const std::vector<std::uint32_t> layout{0, 0, 1};
  std::vector<PbctopoAtom> identity_atoms(3);
  identity_atoms[0].owner = 0;
  identity_atoms[0].source_atom_id = 100;
  identity_atoms[1].owner = 0;
  identity_atoms[1].source_atom_id = 101;
  identity_atoms[2].owner = 1;
  identity_atoms[2].source_atom_id = 200;
  const std::vector<std::vector<std::size_t>> applied_components{
      {0, 1}, {2}};
  const std::vector<std::vector<std::size_t>> selected_components{
      {2}, {0, 1}};
  const std::vector<Int3> greedy_offsets{{0, 0, 0}, {1, 0, 0}};
  const std::vector<Int3> chosen_offsets_reordered{
      {0, -1, 0}, {0, 0, 0}};
  const PbctopoImageAssignment greedy = make_pbctopo_image_assignment(
      greedy_offsets, applied_components, identity_atoms);
  const PbctopoImageAssignment chosen = make_pbctopo_image_assignment(
      chosen_offsets_reordered, selected_components, identity_atoms);
  expect(greedy.valid && chosen.valid &&
             greedy.layout_signature == chosen.layout_signature,
         "stable component assignments were not constructed");
  std::vector<Vec3> replayed;
  expect(replay_pbctopo_assignment_delta(applied, layout, greedy, chosen, Box,
                                          replayed),
         "valid immutable assignment delta was rejected");
  expect(replayed.size() == applied.size(), "replay size changed");
  expect(std::fabs(replayed[0].x - applied[0].x) < 1.0e-12 &&
             std::fabs(replayed[1].x - applied[1].x) < 1.0e-12,
         "fixed component moved");
  expect(std::fabs(replayed[2].x - (applied[2].x - 10.0)) < 1.0e-12 &&
             std::fabs(replayed[2].y - (applied[2].y - 12.0)) < 1.0e-12,
         "selected assignment delta was not applied in Cartesian space");
  PbctopoImageAssignment wrong_layout = chosen;
  ++wrong_layout.layout_signature;
  expect(!replay_pbctopo_assignment_delta(applied, layout, greedy,
                                           wrong_layout, Box, replayed),
         "layout mismatch was replayed instead of failing closed");
  PbctopoImageAssignment duplicate_ids = chosen;
  duplicate_ids.component_ids[1] = duplicate_ids.component_ids[0];
  expect(!replay_pbctopo_assignment_delta(applied, layout, greedy,
                                           duplicate_ids, Box, replayed),
         "duplicate component IDs were replayed");

  PbctopoCutCandidate temporal_choice;
  temporal_choice.valid = true;
  temporal_choice.assignment = chosen;
  temporal_choice.assignment.signature = 42;
  temporal_choice.assignment.layout_signature = 77;
  temporal_choice.min_pair_d2 = 4.0;
  temporal_choice.continuity_max_d2 = 9.0;
  temporal_choice.continuity_path = 2.0;
  temporal_choice.continuity_min_margin_A = 0.75;
  temporal_choice.continuity_margin_atom_a = 0;
  temporal_choice.continuity_margin_atom_b = 2;
  temporal_choice.anchor_max_d2 = 16.0;
  temporal_choice.anchor_mst2 = 8.0;
  temporal_choice.shift_norm2 = 5.0;
  temporal_choice.split_owner_pairs = 2;
  temporal_choice.max_owner_raw_minus_minimg_delta = 0.25;
  temporal_choice.geometry_classification.evaluated = true;
  temporal_choice.geometry_classification.geometry_class =
      PbctopoGeometryClass::CompactUnitBox;
  temporal_choice.geometry_classification
      .fits_current_basis_fundamental_parallelepiped = true;
  temporal_choice.hard_feasibility.evaluated = true;
  temporal_choice.hard_feasibility.assignment_available = true;
  temporal_choice.hard_feasibility.steric_audit_available = true;
  temporal_choice.hard_feasibility.steric_decision_complete = true;
  temporal_choice.hard_feasibility.steric_pair_enumeration_complete = true;
  temporal_choice.hard_feasibility.steric_scan_complete = true;

  PbctopoFrameReport report;
  report.best_score.gap_score = 123.0;
  std::vector<PbctopoAtom> atoms(3);
  const std::vector<PbctopoCutCandidate> pool{temporal_choice};
  expect(!populate_pbctopo_selected_candidate_report(
             report, temporal_choice, atoms, 1.0, 4.0, pool),
         "single-assignment temporal pool was marked anchor-ambiguous");
  expect(report.hard_feasible && report.assignment_signature == 42 &&
             report.component_layout_signature == 77 &&
             report.assignment_component_count == 2,
         "temporal report did not replace assignment/certificate fields");
  expect(report.geometry_classification.geometry_class ==
             PbctopoGeometryClass::CompactUnitBox &&
             report.current_basis_fundamental_parallelepiped_certified &&
             report.compact_unit_box_certified &&
             report.split_owner_pairs == 2 &&
             std::fabs(report.split_owner_gap_scaled - 0.25) < 1.0e-12,
         "temporal report did not replace geometry fields");
  expect(std::fabs(report.clash_margin_A - 1.0) < 1.0e-12 &&
             std::fabs(report.continuity_margin_A - 0.75) < 1.0e-12 &&
             report.continuity_margin_atom_a == 0 &&
             report.continuity_margin_atom_b == 2 &&
             std::fabs(report.anchor_max_edge_A - 4.0) < 1.0e-12 &&
             report.selected_score.valid &&
             std::fabs(report.selected_score.shift_norm2 - 5.0) < 1.0e-12,
         "temporal report did not replace score/margin fields");
  expect(std::fabs(report.best_score.gap_score - 123.0) < 1.0e-12,
         "selected temporal candidate overwrote the pool-best score");
  std::cout << "PASS pbctopo_temporal_assignment_replay_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\temporal_replay.cpp `
        src\pbctopo\selected_report.cpp `
        src\pbctopo\candidate_scoring.cpp `
        src\pbctopo\candidate_order.cpp `
        src\pbctopo\anchor_graph.cpp `
        src\pbctopo\image_assignment.cpp `
        src\pbctopo\steric_clash.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO temporal assignment replay smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO temporal assignment replay smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
