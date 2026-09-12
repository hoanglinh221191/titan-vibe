param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_identifiability_smoke')
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
#include "pbctopo/identifiability.h"

#include <cstdlib>
#include <iostream>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_identifiability_smoke: " << message << "\n";
    std::exit(1);
  }
}

PbctopoPairShiftCandidateSet make_pair(std::size_t a, std::size_t b,
                                       std::size_t contacts,
                                       double fraction,
                                       double quality) {
  PbctopoPairShiftCandidateSet set;
  set.a = a;
  set.b = b;
  PbctopoPairShiftCandidate best;
  best.score.valid = true;
  best.score.complete = true;
  best.score.contact_pairs = contacts;
  best.score.contact_fraction = fraction;
  best.score.quality_score = quality;
  set.candidates.push_back(best);
  PbctopoPairShiftCandidate second = best;
  second.score.quality_score = quality - 1.0;
  set.candidates.push_back(second);
  return set;
}

int main() {
  const std::vector<std::size_t> atom_components{0, 1, 2};
  const auto disconnected = evaluate_pbctopo_identifiability(
      3, {}, {}, atom_components, true, 1);
  expect(!disconnected.identified &&
             disconnected.evidence_connected_component_count == 3 &&
             disconnected.relative_lattice_degrees_of_freedom == 6,
         "disconnected evidence graph was identified");
  expect(!pbctopo_relative_images_identified(disconnected),
         "disconnected graph was allowed to bypass image search");

  const std::vector<PbctopoPairShiftCandidateSet> bbox_only{
      make_pair(0, 1, 0, 0.0, 0.0)};
  const auto no_false_evidence = evaluate_pbctopo_identifiability(
      3, bbox_only, {}, atom_components, false, 0);
  expect(no_false_evidence.evidence_edge_count == 0,
         "bbox-only pair became an evidence edge");

  const std::vector<PbctopoPairShiftCandidateSet> scored_contact{
      make_pair(0, 1, 4, 0.5, 3.0)};
  const std::vector<PbctopoRequiredContact> explicit_contact{
      {1, 2, 5.0, 1.0}};
  const auto identified = evaluate_pbctopo_identifiability(
      3, scored_contact, explicit_contact, atom_components, true, 1);
  expect(identified.identified && identified.evidence_graph_connected &&
             identified.evidence_edge_count == 2 &&
             identified.explicit_contact_evidence_edges == 1 &&
             identified.scored_contact_evidence_edges == 1 &&
             identified.objective_uniqueness_known &&
             identified.objective_unique_within_search_domain &&
             identified.unique_within_search_domain &&
             !identified.evidence_uniqueness_known &&
             !identified.evidence_unique_within_search_domain,
         "objective uniqueness was not separated from evidence uniqueness");
  expect(pbctopo_relative_images_identified(identified),
         "connected evidence graph was blocked from the fast path");

  const auto ambiguous = evaluate_pbctopo_identifiability(
      3, scored_contact, explicit_contact, atom_components, true, 2, true, 2);
  expect(ambiguous.identified && !ambiguous.unique_within_search_domain &&
             ambiguous.evidence_uniqueness_known &&
             !ambiguous.evidence_unique_within_search_domain &&
             ambiguous.feasible_assignments_within_domain == 2,
         "equivalent best assignments were reported as unique");

  const auto objective_unique_but_evidence_ambiguous =
      evaluate_pbctopo_identifiability(
          3, scored_contact, explicit_contact, atom_components, true, 1, true,
          2);
  expect(objective_unique_but_evidence_ambiguous
             .objective_unique_within_search_domain &&
             objective_unique_but_evidence_ambiguous
                 .evidence_uniqueness_known &&
             !objective_unique_but_evidence_ambiguous
                  .evidence_unique_within_search_domain,
         "a unique optimum was confused with a unique feasible assignment");

  const auto fully_unique = evaluate_pbctopo_identifiability(
      3, scored_contact, explicit_contact, atom_components, true, 1, true, 1);
  expect(fully_unique.objective_unique_within_search_domain &&
             fully_unique.evidence_unique_within_search_domain &&
             fully_unique.feasible_assignment_set_enumerated,
         "fully enumerated singleton evidence was not reported as unique");

  std::cout << "PASS pbctopo_identifiability_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ascii

& $Compiler -std=c++23 -O2 -Isrc $src `
    src/pbctopo/identifiability.cpp `
    src/pbctopo/component_graph_math.cpp `
    src/pbctopo/component_metric.cpp `
    src/pbctopo/component_shift_geometry.cpp `
    src/pbctopo/lattice_math.cpp `
    -o $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile identifiability smoke"
}
& $exe
if ($LASTEXITCODE -ne 0) {
    throw "Identifiability smoke failed"
}
