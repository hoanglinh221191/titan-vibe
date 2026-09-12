param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_assignment_set_report_smoke')
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
#include "pbctopo/assignment_set_report.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_assignment_set_report_smoke: " << message
              << "\n";
    std::exit(1);
  }
}

bool better(const PbctopoCutCandidate &lhs,
            const PbctopoCutCandidate &rhs) {
  return lhs.shift_norm2 < rhs.shift_norm2;
}

PbctopoCutCandidate make_candidate(std::uint64_t signature,
                                   std::int64_t second_offset) {
  PbctopoCutCandidate candidate;
  candidate.valid = true;
  candidate.shift_norm2 = 1.0;
  candidate.assignment.valid = true;
  candidate.assignment.signature = signature;
  candidate.assignment.layout_signature = 99;
  candidate.assignment.component_ids = {20, 10};
  candidate.assignment.component_offsets = {{0, second_offset, 0},
                                              {0, 0, 0}};
  candidate.assignment.global_gauge = {4, 5, 6};
  candidate.search_equivalence_known = true;
  candidate.equivalent_to_search_best = true;
  candidate.component_search_shell_radius = 2;
  candidate.search_domain_audit_known = true;
  candidate.search_domain_evidence_complete = true;
  candidate.selected_assignment_touches_shell_boundary = false;
  candidate.bounded_domain_only = true;
  candidate.domain_completeness =
      PbctopoDomainCompleteness::BoundedExhaustive;
  candidate.hard_feasibility.evaluated = true;
  candidate.hard_feasibility.assignment_available = true;
  candidate.hard_feasibility.steric_audit_available = true;
  candidate.hard_feasibility.steric_decision_complete = true;
  candidate.hard_feasibility.steric_pair_enumeration_complete = true;
  candidate.hard_feasibility.steric_scan_complete = true;
  return candidate;
}

int main() {
  const std::vector<PbctopoCutCandidate> candidates{
      make_candidate(101, 1), make_candidate(102, -1)};
  PbctopoIdentifiability identifiability;
  identifiability.identified = true;
  identifiability.unique_within_search_domain = false;
  identifiability.uniqueness_search_exhaustive = true;
  identifiability.equivalent_best_assignments_within_domain = 2;
  identifiability.objective_uniqueness_known = true;
  identifiability.objective_unique_within_search_domain = false;
  identifiability.evidence_uniqueness_known = true;
  identifiability.evidence_unique_within_search_domain = false;
  identifiability.feasible_assignment_set_enumerated = true;
  identifiability.feasible_assignments_within_domain = 2;

  std::ostringstream out;
  write_pbctopo_assignment_set_header(out);
  write_pbctopo_assignment_set_rows(
      out, 7, 2.5, candidates, &candidates[1], identifiability,
      {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  const std::string text = out.str();
  expect(text.find("component_assignment") != std::string::npos,
         "stable component-assignment column is missing");
  expect(text.find("schema_version,producer_version,coordinate_unit") !=
             std::string::npos,
         "observation schema metadata is missing");
  expect(text.find(
             ",7,2.5,1,0,retained_certified_band,ambiguous_anchor,1,1") !=
             std::string::npos,
         "first equivalent assignment row has incorrect semantics");
  expect(text.find(
             ",7,2.5,2,1,retained_certified_band,ambiguous_anchor,1,1,1,1") !=
             std::string::npos,
         "framewise selected assignment was not marked");
  expect(text.find("10:0:0:0|20:0:1:0") != std::string::npos,
         "component assignments were not serialized by stable component ID");

  auto unknown_candidates = candidates;
  for (auto &candidate : unknown_candidates) {
    candidate.search_equivalence_known = false;
    candidate.equivalent_to_search_best = false;
  }
  std::ostringstream unknown_out;
  write_pbctopo_assignment_set_rows(
      unknown_out, 8, 3.0, unknown_candidates, &unknown_candidates[0],
      identifiability, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(unknown_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             unknown_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             unknown_out.str().find("ambiguous_anchor,0,na") !=
                 std::string::npos,
         "unknown search equivalence was reported as a complete set");
  PbctopoIdentifiability non_exhaustive_unique = identifiability;
  non_exhaustive_unique.identified = true;
  non_exhaustive_unique.unique_within_search_domain = true;
  non_exhaustive_unique.uniqueness_search_exhaustive = false;
  non_exhaustive_unique.objective_uniqueness_known = false;
  non_exhaustive_unique.objective_unique_within_search_domain = false;
  non_exhaustive_unique.evidence_uniqueness_known = false;
  non_exhaustive_unique.evidence_unique_within_search_domain = false;
  non_exhaustive_unique.feasible_assignment_set_enumerated = false;
  non_exhaustive_unique.feasible_assignments_within_domain = 0;
  std::ostringstream non_exhaustive_out;
  write_pbctopo_assignment_set_rows(
      non_exhaustive_out, 9, 4.0, candidates, &candidates[0],
      non_exhaustive_unique, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(non_exhaustive_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             non_exhaustive_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             non_exhaustive_out.str().find("ambiguous_anchor") !=
                 std::string::npos,
         "non-exhaustive uniqueness was promoted to a strong anchor");

  PbctopoIdentifiability exhaustive_unique = identifiability;
  exhaustive_unique.identified = true;
  exhaustive_unique.unique_within_search_domain = true;
  exhaustive_unique.uniqueness_search_exhaustive = true;
  exhaustive_unique.equivalent_best_assignments_within_domain = 1;
  exhaustive_unique.objective_uniqueness_known = true;
  exhaustive_unique.objective_unique_within_search_domain = true;
  exhaustive_unique.evidence_uniqueness_known = true;
  exhaustive_unique.evidence_unique_within_search_domain = true;
  exhaustive_unique.feasible_assignment_set_enumerated = true;
  exhaustive_unique.feasible_assignments_within_domain = 1;

  auto incomplete_candidates = candidates;
  for (auto &candidate : incomplete_candidates) {
    candidate.search_equivalence_known = false;
    candidate.equivalent_to_search_best = false;
  }
  std::ostringstream incomplete_out;
  write_pbctopo_assignment_set_rows(
      incomplete_out, 10, 5.0, incomplete_candidates,
      &incomplete_candidates[0], exhaustive_unique,
      {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(incomplete_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             incomplete_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             incomplete_out.str().find("ambiguous_anchor,0,na") !=
                 std::string::npos,
         "incomplete equivalent export was promoted to a strong anchor");

  auto unique_candidates = candidates;
  unique_candidates[0].shift_norm2 = 0.5;
  unique_candidates[0].search_equivalence_known = true;
  unique_candidates[0].equivalent_to_search_best = true;
  unique_candidates[1].shift_norm2 = 2.0;
  unique_candidates[1].search_equivalence_known = false;
  unique_candidates[1].equivalent_to_search_best = false;
  std::ostringstream unique_out;
  write_pbctopo_assignment_set_rows(
      unique_out, 11, 6.0, unique_candidates, &unique_candidates[0],
      exhaustive_unique, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  const std::string unique_text = unique_out.str();
  expect(unique_text.find("strong_bounded_anchor") == std::string::npos &&
             unique_text.find("strong_global_anchor") == std::string::npos &&
             unique_text.find(
                 ",11,6,1,1,retained_certified_band,ambiguous_anchor,1,1") !=
                 std::string::npos,
         "unknown retained hard-certified assignment did not block a strong anchor");
  expect(unique_text.find(
             ",11,6,2,0,retained_certified_band,ambiguous_anchor,0,na") !=
             std::string::npos,
         "non-best retained candidate was promoted to a strong anchor");

  const std::vector<PbctopoCutCandidate> sole_candidate{
      unique_candidates[0]};
  std::ostringstream sole_out;
  write_pbctopo_assignment_set_rows(
      sole_out, 13, 8.0, sole_candidate, &sole_candidate[0],
      exhaustive_unique, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(sole_out.str().find(
             ",13,8,1,1,retained_certified_band,strong_bounded_anchor,1,1,1,1") !=
             std::string::npos,
         "complete unique selected optimum was not a bounded strong anchor");

  auto global_candidate = sole_candidate;
  global_candidate[0].bounded_domain_only = false;
  global_candidate[0].domain_expansion_attempted = true;
  global_candidate[0].domain_expansion_exhausted = true;
  global_candidate[0].domain_completeness =
      PbctopoDomainCompleteness::GlobalCertified;
  std::ostringstream global_out;
  write_pbctopo_assignment_set_rows(
      global_out, 16, 11.0, global_candidate, &global_candidate[0],
      exhaustive_unique, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(global_out.str().find("strong_global_anchor") == std::string::npos &&
             global_out.str().find(
                 ",16,11,1,1,retained_certified_band,ambiguous_anchor,1,1,1,1") !=
                 std::string::npos,
          "unsupported global claim was promoted without outside-domain proof");

  auto boundary_candidate = sole_candidate;
  boundary_candidate[0].selected_assignment_touches_shell_boundary = true;
  std::ostringstream boundary_out;
  write_pbctopo_assignment_set_rows(
      boundary_out, 14, 9.0, boundary_candidate, &boundary_candidate[0],
      exhaustive_unique, {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(boundary_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             boundary_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             boundary_out.str().find("ambiguous_anchor") !=
                 std::string::npos,
         "shell-boundary assignment was promoted to a strong anchor");

  auto incomplete_domain_candidate = sole_candidate;
  incomplete_domain_candidate[0].search_domain_evidence_complete = false;
  std::ostringstream incomplete_domain_out;
  write_pbctopo_assignment_set_rows(
      incomplete_domain_out, 15, 10.0, incomplete_domain_candidate,
      &incomplete_domain_candidate[0], exhaustive_unique,
      {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(incomplete_domain_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             incomplete_domain_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             incomplete_domain_out.str().find("ambiguous_anchor") !=
                 std::string::npos,
         "incomplete search-domain evidence was promoted to a strong anchor");

  PbctopoIdentifiability objective_only = exhaustive_unique;
  objective_only.evidence_unique_within_search_domain = false;
  objective_only.feasible_assignments_within_domain = 2;
  auto objective_only_candidates = unique_candidates;
  objective_only_candidates[1].search_equivalence_known = true;
  std::ostringstream objective_only_out;
  write_pbctopo_assignment_set_rows(
      objective_only_out, 12, 7.0, objective_only_candidates,
      &objective_only_candidates[0], objective_only,
      {1.0, 0.1, 0.2, 2.0, 0.3, 3.0}, better);
  expect(objective_only_out.str().find("strong_bounded_anchor") ==
                 std::string::npos &&
             objective_only_out.str().find("strong_global_anchor") ==
                 std::string::npos &&
             objective_only_out.str().find("ambiguous_anchor") !=
                 std::string::npos,
         "a unique objective optimum was promoted without unique evidence");
  std::cout << "PASS pbctopo_assignment_set_report_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\assignment_set_report.cpp `
        src\pbctopo\observation_binding.cpp `
        src\pbctopo\image_assignment.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO assignment-set report smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO assignment-set report smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
