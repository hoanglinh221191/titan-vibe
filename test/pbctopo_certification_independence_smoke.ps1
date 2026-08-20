param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_certification_independence_smoke')
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
#include "pbctopo/steric_clash.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_certification_independence_smoke: "
              << message << "\n";
    std::exit(1);
  }
}

std::uint64_t pair_key(std::size_t a, std::size_t b) {
  const std::size_t lo = std::min(a, b);
  const std::size_t hi = std::max(a, b);
  return (static_cast<std::uint64_t>(lo) << 32) |
         static_cast<std::uint64_t>(hi);
}

box make_box() {
  box value;
  value.set_a(10.0);
  value.set_b(10.0);
  value.set_c(10.0);
  return value;
}

int main() {
  expect(pbctopo_edge_is_hard_connectivity(
             PbctopoLocalEdgeKind::TopologyBond),
         "topology bond is not hard connectivity");
  expect(pbctopo_edge_is_hard_connectivity(
             PbctopoLocalEdgeKind::MetadataSkeleton),
         "trusted metadata is not hard connectivity");
  expect(!pbctopo_edge_is_hard_connectivity(
             PbctopoLocalEdgeKind::GeometryHeuristic),
         "geometry heuristic became hard connectivity");
  expect(pbctopo_edge_is_steric_exclusion(
             PbctopoLocalEdgeKind::TopologyBond),
         "topology bond is not a steric exclusion");
  expect(!pbctopo_edge_is_steric_exclusion(
             PbctopoLocalEdgeKind::MetadataSkeleton) &&
             !pbctopo_edge_is_steric_exclusion(
                 PbctopoLocalEdgeKind::GeometryHeuristic),
         "heuristic edge opened a steric-certificate hole");

  std::vector<PbctopoAtom> atoms(2);
  atoms[0].is_heavy = true;
  atoms[0].has_vdw_radius = true;
  atoms[0].vdw_radius = 1.70;
  atoms[0].steric_radius_source =
      PbctopoStericRadiusSource::ForceFieldParameter;
  atoms[1].is_heavy = false;
  atoms[1].has_vdw_radius = true;
  atoms[1].vdw_radius = 1.20;
  atoms[1].steric_radius_source =
      PbctopoStericRadiusSource::ForceFieldParameter;
  prepare_pbctopo_steric_radii(atoms, 0.55);

  const std::vector<Vec3> hydrogen_clash{{0.0, 0.0, 0.0},
                                          {0.5, 0.0, 0.0}};
  const std::vector<std::size_t> components{0, 1};
  const std::unordered_set<std::uint64_t> no_exclusions;
  const auto steric = evaluate_pbctopo_final_sterics(
      atoms, hydrogen_clash, components, no_exclusions, 1.05);
  expect(steric.input_valid && steric.pair_checks == 1 &&
             steric.hard_clashes == 1 && steric.decision_complete &&
             !steric.pair_enumeration_complete &&
             steric.termination_reason ==
                 PbctopoStericTerminationReason::FirstHardClash,
         "all-atom audit missed a hydrogen-only clash");

  PbctopoCutCandidate candidate;
  candidate.valid = true;
  candidate.ultra_clash = false;
  candidate.assignment.valid = true;
  PbctopoImageEquationResiduals image_equations;
  image_equations.input_valid = true;
  PbctopoAssignmentPointResiduals assignment_points;
  assignment_points.input_valid = true;
  PbctopoRequiredContactResiduals contacts;
  contacts.input_valid = true;
  expect(!certify_pbctopo_candidate(candidate, image_equations,
                                     assignment_points, steric, contacts),
         "search clash flag overrode the independent final steric audit");

  std::unordered_set<std::uint64_t> topology_exclusion{pair_key(0, 1)};
  const auto excluded = evaluate_pbctopo_final_sterics(
      atoms, hydrogen_clash, components, topology_exclusion, 1.05);
  expect(excluded.input_valid && excluded.pair_checks == 0 &&
             excluded.hard_clashes == 0 && excluded.decision_complete &&
             excluded.pair_enumeration_complete &&
             excluded.termination_reason ==
                 PbctopoStericTerminationReason::Completed,
         "explicit topology exclusion was not honored");

  const box Box = make_box();
  const PbctopoLatticeMetric lattice(Box);
  atoms[0].has_topology_key = true;
  atoms[1].has_topology_key = true;
  atoms[0].x = 0.0;
  atoms[1].x = 4.0;
  const std::vector<std::size_t> same_component{0, 0};
  const auto fixed_clash = evaluate_pbctopo_fixed_certificate_invariants(
      atoms, hydrogen_clash, same_component, {}, {}, no_exclusions, 1.05,
      lattice);
  expect(fixed_clash.evaluated && fixed_clash.input_valid &&
             fixed_clash.metadata_internal_clashes == 1 &&
             !fixed_clash.passed(),
         "fixed metadata-internal clash was not rejected");

  const std::vector<Vec3> separated{{0.0, 0.0, 0.0},
                                    {3.0, 0.0, 0.0}};
  const std::vector<PbctopoRequiredContact> lost_contact{
      {0, 1, 1.0, 1.0}};
  const auto fixed_contact = evaluate_pbctopo_fixed_certificate_invariants(
      atoms, separated, same_component, {}, lost_contact, no_exclusions, 1.05,
      lattice);
  expect(fixed_contact.input_valid &&
             fixed_contact.required_contact_scan_complete &&
             fixed_contact.required_contacts_checked == 1 &&
             fixed_contact.required_contacts_lost == 1 &&
             !fixed_contact.passed(),
         "fixed same-component required contact was not rejected");

  std::vector<PbctopoAtom> fallback_atoms(2);
  fallback_atoms[0].steric_radius_source =
      PbctopoStericRadiusSource::Unavailable;
  fallback_atoms[1].steric_radius_source =
      PbctopoStericRadiusSource::Unavailable;
  const std::vector<std::size_t> separate_components{0, 1};
  const auto fallback_model = evaluate_pbctopo_fixed_certificate_invariants(
      fallback_atoms, separated, separate_components, {}, {}, no_exclusions,
      1.05, lattice);
  expect(fallback_model.radius_model_valid && fallback_model.passed(),
         "declared generic cutoff fallback was treated as an invalid model");

  fallback_atoms[0].steric_radius_source =
      PbctopoStericRadiusSource::GenericFallback;
  const auto missing_generic_radius =
      evaluate_pbctopo_fixed_certificate_invariants(
          fallback_atoms, separated, separate_components, {}, {},
          no_exclusions, 1.05, lattice);
  expect(!missing_generic_radius.radius_model_valid &&
             !missing_generic_radius.passed(),
         "generic radius source without a radius passed the fixed audit");

  fallback_atoms[0].has_vdw_radius = true;
  fallback_atoms[0].vdw_radius = 1.70;
  fallback_atoms[0].steric_radius = 0.935;
  const auto generic_radius = evaluate_pbctopo_fixed_certificate_invariants(
      fallback_atoms, separated, separate_components, {}, {}, no_exclusions,
      1.05, lattice);
  expect(generic_radius.radius_model_valid && generic_radius.passed(),
         "declared positive generic radius was rejected");

  fallback_atoms[0].vdw_radius = 0.0;
  fallback_atoms[0].steric_radius = 0.0;
  fallback_atoms[0].steric_radius_source =
      PbctopoStericRadiusSource::ExplicitZeroLennardJones;
  const auto explicit_zero = evaluate_pbctopo_fixed_certificate_invariants(
      fallback_atoms, separated, separate_components, {}, {}, no_exclusions,
      1.05, lattice);
  expect(explicit_zero.radius_model_valid && explicit_zero.passed(),
         "consistent explicit zero-LJ state was rejected");

  fallback_atoms[0].vdw_radius = 1.0;
  fallback_atoms[0].steric_radius = 0.55;
  const auto invalid_zero = evaluate_pbctopo_fixed_certificate_invariants(
      fallback_atoms, separated, separate_components, {}, {}, no_exclusions,
      1.05, lattice);
  expect(!invalid_zero.radius_model_valid && !invalid_zero.passed(),
         "positive explicit zero-LJ radius passed the fixed audit");

  fallback_atoms[0].has_vdw_radius = true;
  fallback_atoms[0].vdw_radius = std::numeric_limits<double>::quiet_NaN();
  fallback_atoms[0].steric_radius = 0.55;
  fallback_atoms[0].steric_radius_source =
      PbctopoStericRadiusSource::ForceFieldParameter;
  const auto invalid_radius = evaluate_pbctopo_fixed_certificate_invariants(
      fallback_atoms, separated, separate_components, {}, {}, no_exclusions,
      1.05, lattice);
  expect(!invalid_radius.radius_model_valid && !invalid_radius.passed(),
         "inconsistent force-field radius state passed the fixed audit");

  std::cout << "PASS pbctopo_certification_independence_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\certification.cpp `
        src\pbctopo\steric_clash.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO certification-independence smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO certification-independence smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
