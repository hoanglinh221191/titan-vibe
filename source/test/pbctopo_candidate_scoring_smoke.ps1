param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_candidate_scoring_smoke')
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
#include "pbctopo/candidate_scoring.h"
#include "pbctopo/candidate_order.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace titan_pbctopo {

int compare_anchor_graph_score(const PbctopoAnchorGraphScore &lhs,
                               const PbctopoAnchorGraphScore &rhs) {
  if (lhs.residue_order_violations != rhs.residue_order_violations)
    return lhs.residue_order_violations < rhs.residue_order_violations ? -1 : 1;
  return compare_pbctopo_metric(lhs.order_max_d2, rhs.order_max_d2,
                                {1.0e-6, 0.0});
}

} // namespace titan_pbctopo

namespace {

void require_true(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL pbctopo_candidate_scoring_smoke: " << message << '\n';
    std::exit(1);
  }
}

titan_pbctopo::PbctopoCutCandidate make_candidate(int index) {
  using namespace titan_pbctopo;
  PbctopoCutCandidate candidate;
  candidate.valid = true;
  candidate.ultra_clash = index % 17 == 0;
  candidate.broken_edges = static_cast<std::size_t>(index % 4);
  candidate.continuity_max_d2 = 20.0 + static_cast<double>((index * 7) % 19);
  candidate.continuity_path = 80.0 + static_cast<double>((index * 11) % 31);
  candidate.anchor_graph.residue_order_violations =
      static_cast<std::size_t>((index * 3) % 2);
  candidate.anchor_graph.contact_map_lost =
      static_cast<std::size_t>((index * 5) % 3);
  candidate.anchor_graph.order_max_d2 =
      10.0 + static_cast<double>((index * 13) % 23);
  candidate.anchor_max_d2 = 30.0 + static_cast<double>((index * 17) % 29);
  candidate.anchor_mst2 = 100.0 + static_cast<double>((index * 19) % 37);
  candidate.principal_rg2_max = 200.0 + static_cast<double>((index * 23) % 61);
  candidate.principal_rg2_mid = 100.0 + static_cast<double>((index * 29) % 47);
  candidate.principal_rg2_min = 20.0 + static_cast<double>((index * 31) % 41);
  candidate.shape_kappa2 = 0.01 * static_cast<double>((index * 7) % 25);
  candidate.shape_asphericity = static_cast<double>((index * 11) % 43);
  candidate.shape_acylindricity = static_cast<double>((index * 13) % 37);
  candidate.edge_stretch = 0.01 * static_cast<double>((index * 17) % 31);
  candidate.min_pair_d2 = 2.0 + 0.1 * static_cast<double>((index * 19) % 17);
  candidate.shift_norm2 = static_cast<double>((index * 23) % 29);
  candidate.cut_x = 0.5 + 0.01 * static_cast<double>(index);
  candidate.cut_y = 0.6 + 0.01 * static_cast<double>(index);
  candidate.cut_z = 0.7 + 0.01 * static_cast<double>(index);
  return candidate;
}

} // namespace

int main() {
  using namespace titan_pbctopo;
  std::vector<PbctopoCutCandidate> candidates;
  for (int index = 1; index <= 48; ++index)
    candidates.push_back(make_candidate(index));
  candidates.push_back(candidates[3]);
  candidates.push_back(candidates[11]);

  auto less = [](const auto &lhs, const auto &rhs) {
    return better_cut_candidate(lhs, rhs);
  };
  auto equivalent = [&](const auto &lhs, const auto &rhs) {
    return !less(lhs, rhs) && !less(rhs, lhs);
  };
  for (const auto &a : candidates) {
    require_true(!less(a, a), "comparator is not irreflexive");
    for (const auto &b : candidates) {
      require_true(!(less(a, b) && less(b, a)),
                   "comparator is not asymmetric");
      for (const auto &c : candidates) {
        require_true(!(less(a, b) && less(b, c)) || less(a, c),
                     "comparator is not transitive");
        require_true(!(equivalent(a, b) && equivalent(b, c)) ||
                         equivalent(a, c),
                     "comparator equivalence is not transitive");
      }
    }
  }
  auto sorted_reference = candidates;
  std::stable_sort(sorted_reference.begin(), sorted_reference.end(), less);
  for (int rotation = 1; rotation < 8; ++rotation) {
    auto permuted = candidates;
    std::rotate(permuted.begin(), permuted.begin() + rotation, permuted.end());
    std::stable_sort(permuted.begin(), permuted.end(), less);
    require_true(permuted.size() == sorted_reference.size(),
                 "permutation changed candidate count");
    for (std::size_t idx = 0; idx < permuted.size(); ++idx)
      require_true(equivalent(permuted[idx], sorted_reference[idx]),
                   "sorting depends on input permutation");
  }

  auto feasible = make_candidate(2);
  feasible.ultra_clash = false;
  feasible.broken_edges = 0;
  feasible.principal_rg2_max = 5000.0;
  auto compact_but_broken = feasible;
  compact_but_broken.broken_edges = 1;
  compact_but_broken.principal_rg2_max = 10.0;
  require_true(better_cut_candidate(feasible, compact_but_broken),
               "feasibility must rank before compactness");

  auto compact_incumbent = feasible;
  compact_incumbent.principal_rg2_max = 10.0;
  compact_incumbent.principal_rg2_mid = 5.0;
  compact_incumbent.principal_rg2_min = 1.0;
  compact_incumbent.shape_kappa2 = 0.0;
  compact_incumbent.shape_asphericity = 0.0;
  auto elongated = compact_incumbent;
  elongated.principal_rg2_max = 10000.0;
  elongated.principal_rg2_mid = 5000.0;
  elongated.principal_rg2_min = 1000.0;
  elongated.shape_kappa2 = 0.95;
  elongated.shape_asphericity = 9000.0;
  require_true(!cut_candidate_cannot_beat_incumbent(
                   elongated, &compact_incumbent),
               "elongated geometry was hard-pruned by compactness");
  require_true(compare_candidate_prefix_before_edge_stretch(
                   elongated, compact_incumbent,
                   PbctopoRefineScorePhase::CoarseStage1) == 0,
               "principal radii entered the hard edge-scan prefix");

  auto same_strain_longer_edge = feasible;
  same_strain_longer_edge.continuity_max_d2 *= 4.0;
  require_true(!better_cut_candidate(feasible, same_strain_longer_edge) &&
                   !better_cut_candidate(same_strain_longer_edge, feasible),
               "continuity max distance was used as a soft compactness rank");
  auto incompatible_contact_annotation = feasible;
  incompatible_contact_annotation.anchor_graph.contact_map_lost = 999;
  incompatible_contact_annotation.anchor_graph.contact_map_lost_fraction = 1.0;
  incompatible_contact_annotation.anchor_graph.contact_map_max_d2 = 1.0e9;
  require_true(equivalent(feasible, incompatible_contact_annotation),
               "optional contact annotation entered candidate ordering");
  auto worse_strain = feasible;
  worse_strain.continuity_path += 0.1;
  require_true(better_cut_candidate(feasible, worse_strain),
               "lower reference strain must rank first");

  PbctopoComponentCandidate elongated_component;
  elongated_component.valid = true;
  elongated_component.continuity_max_d2 = 1.0;
  elongated_component.continuity_path = 0.0;
  elongated_component.principal_rg2_max = 10000.0;
  PbctopoComponentCandidate compact_component = elongated_component;
  compact_component.principal_rg2_max = 1.0;
  require_true(!component_candidate_cannot_beat_incumbent_prefix(
                   elongated_component, &compact_component,
                   PbctopoComponentCandidatePrefixStage::Principal),
               "component scorer pruned an elongated candidate");

  auto clashing = feasible;
  clashing.ultra_clash = true;
  require_true(better_cut_candidate(feasible, clashing),
               "safe candidate must rank before ultra clash");

  auto second = feasible;
  second.broken_edges = 2;
  const auto gap = candidate_lexicographic_gap(feasible, second);
  require_true(gap.metric == "broken_edges" && gap.gap == 2.0,
               "score gap must report the first distinguishing metric");
  require_true(candidate_temporal_emission_cost(feasible, feasible) == 0.0,
               "best candidate must have zero temporal emission cost");
  auto cut_only = feasible;
  cut_only.cut_x += 0.1;
  require_true(candidate_temporal_emission_cost(feasible, second) >
                   candidate_temporal_emission_cost(feasible, cut_only),
               "hard geometry gap must outweigh a cut tie-break gap");

  PbctopoComponentCandidate component;
  component.valid = true;
  component.shift_norm2 = 12.5;
  component.min_pair_d2 = 3.25;
  const auto converted = component_candidate_to_cut_candidate(component);
  require_true(converted.valid && converted.shift_norm2 == 12.5 &&
                   converted.min_pair_d2 == 3.25,
               "component-to-cut conversion lost ranking fields");

  std::cout << "PASS pbctopo_candidate_scoring_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\candidate_scoring.cpp `
        src\pbctopo\candidate_order.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO candidate scoring smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO candidate scoring smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
