param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_evidence_support_policy_smoke')
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
#include "pbctopo/evidence_support_policy.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_evidence_support_policy_smoke: " << message
              << "\n";
    std::exit(1);
  }
}

bool legacy_scored_pair_decision(
    const PbctopoPairShiftCandidateSet &pair_set) {
  if (pair_set.candidates.empty())
    return false;
  const auto &best = pair_set.candidates.front();
  const double confidence =
      pair_set.candidates.size() > 1
          ? pbctopo_pair_shift_confidence(
                best.score, pair_set.candidates[1].score)
          : 0.0;
  return best.score.valid && best.score.complete && !best.score.ultra_clash &&
         best.score.contact_pairs > 0 &&
         best.score.contact_fraction > 0.0 &&
         (confidence >= PBCTOPO_SCORED_PAIR_SUPPORT_THRESHOLD ||
          best.score.contact_fraction >=
              PBCTOPO_SCORED_PAIR_SUPPORT_THRESHOLD);
}

PbctopoPairShiftCandidate make_candidate(bool valid, bool complete,
                                          bool ultra_clash,
                                          std::size_t contact_pairs,
                                          double contact_fraction,
                                          double quality_score) {
  PbctopoPairShiftCandidate candidate;
  candidate.score.valid = valid;
  candidate.score.complete = complete;
  candidate.score.ultra_clash = ultra_clash;
  candidate.score.contact_pairs = contact_pairs;
  candidate.score.contact_fraction = contact_fraction;
  candidate.score.quality_score = quality_score;
  return candidate;
}

int main() {
  expect(PBCTOPO_EVIDENCE_SUPPORT_POLICY_VERSION == 1,
         "unexpected support policy version");
  expect(pbctopo_evidence_support_reason_name(
             PbctopoEvidenceSupportReason::NoCompatibleHypothesis) ==
             "no_compatible_hypothesis",
         "support reason vocabulary changed");

  std::size_t scored_cases = 0;
  constexpr double fractions[] = {0.0, 0.01, 0.05, 0.2};
  constexpr double qualities[] = {0.0, 100.0};
  for (const bool valid : {false, true}) {
    for (const bool complete : {false, true}) {
      for (const bool ultra_clash : {false, true}) {
        for (const std::size_t contact_pairs : {std::size_t{0},
                                                std::size_t{1}}) {
          for (const double fraction : fractions) {
            for (const double quality : qualities) {
              for (int second_mode = 0; second_mode < 5; ++second_mode) {
                PbctopoPairShiftCandidateSet pair_set;
                pair_set.a = 0;
                pair_set.b = 1;
                pair_set.candidates.push_back(make_candidate(
                    valid, complete, ultra_clash, contact_pairs, fraction,
                    quality));
                if (second_mode != 0) {
                  bool second_valid = true;
                  bool second_complete = true;
                  double second_quality = quality;
                  if (second_mode == 2)
                    second_quality = quality - 10.0;
                  if (second_mode == 3)
                    second_valid = false;
                  if (second_mode == 4)
                    second_complete = false;
                  pair_set.candidates.push_back(make_candidate(
                      second_valid, second_complete, false, contact_pairs,
                      fraction, second_quality));
                }

                const auto support =
                    evaluate_pbctopo_scored_pair_evidence_support(pair_set);
                const bool legacy = legacy_scored_pair_decision(pair_set);
                expect(support.active_at_threshold == legacy,
                       "finite scored-pair parity mismatch");
                expect(support.active_at_threshold ==
                           (support.support_primary >=
                            support.decision_threshold),
                       "scored-pair threshold does not replay its decision");
                expect(std::isfinite(support.support_primary) &&
                           std::isfinite(support.decision_threshold) &&
                           support.raw_valid,
                       "finite scored-pair input produced invalid support");
                ++scored_cases;
              }
            }
          }
        }
      }
    }
  }
  expect(scored_cases == 640, "unexpected finite scored-pair case count");

  PbctopoPairShiftCandidateSet empty_pair;
  const auto empty_support =
      evaluate_pbctopo_scored_pair_evidence_support(empty_pair);
  expect(!empty_support.active_at_threshold &&
             empty_support.reason ==
                 PbctopoEvidenceSupportReason::NoCandidates &&
             empty_support.support_primary == 0.0,
         "empty pair set did not map to inactive finite support");

  PbctopoPairShiftCandidateSet boundary_pair;
  boundary_pair.candidates.push_back(
      make_candidate(true, true, false, 1, 0.05, 0.0));
  const auto boundary_support =
      evaluate_pbctopo_scored_pair_evidence_support(boundary_pair);
  expect(boundary_support.active_at_threshold &&
             boundary_support.support_primary == 0.05,
         "inclusive scored-pair threshold was not preserved");

  PbctopoPairShiftCandidateSet confidence_pair;
  confidence_pair.candidates.push_back(
      make_candidate(true, true, false, 1, 0.01, 10.0));
  confidence_pair.candidates.push_back(
      make_candidate(true, true, false, 1, 0.01, 0.0));
  const auto confidence_support =
      evaluate_pbctopo_scored_pair_evidence_support(confidence_pair);
  expect(confidence_support.active_at_threshold &&
             confidence_support.rank_confidence > 0.05 &&
             confidence_support.support_primary ==
                 confidence_support.rank_confidence,
         "rank-confidence branch did not activate scored support");

  for (const double invalid_fraction : {
           std::numeric_limits<double>::quiet_NaN(),
           std::numeric_limits<double>::infinity()}) {
    PbctopoPairShiftCandidateSet invalid_pair;
    invalid_pair.candidates.push_back(
        make_candidate(true, true, false, 1, invalid_fraction, 0.0));
    const auto invalid_support =
        evaluate_pbctopo_scored_pair_evidence_support(invalid_pair);
    expect(!invalid_support.raw_valid &&
               !invalid_support.active_at_threshold &&
               invalid_support.support_primary == 0.0 &&
               invalid_support.reason ==
                   PbctopoEvidenceSupportReason::NonfiniteContactFraction,
           "nonfinite contact fraction did not fail closed");
  }

  const std::vector<Int3> shift_catalog{{-1, 0, 0},
                                        {0, 0, 0},
                                        {1, 0, 0}};
  std::size_t selected_cases = 0;
  for (std::uint32_t subset = 0; subset < 8; ++subset) {
    std::vector<Int3> observed;
    for (std::size_t index = 0; index < shift_catalog.size(); ++index) {
      if ((subset & (std::uint32_t{1} << index)) != 0)
        observed.push_back(shift_catalog[index]);
    }
    for (const bool explicit_relation : {false, true}) {
      for (std::int64_t selected_x = -2; selected_x <= 2; ++selected_x) {
        const Int3 selected{selected_x, 0, 0};
        const bool legacy_match = std::any_of(
            observed.begin(), observed.end(),
            [&](const Int3 &shift) { return shift == selected; });
        const auto support =
            evaluate_pbctopo_selected_compatibility_support(
                explicit_relation, observed, selected);
        expect(support.compatibility_match == legacy_match,
               "selected compatibility-match parity mismatch");
        expect(support.active_at_threshold ==
                   (explicit_relation || legacy_match),
               "selected compatibility active-edge parity mismatch");
        expect(support.active_at_threshold ==
                   (support.support_primary >= support.decision_threshold),
               "selected support threshold does not replay its decision");
        expect(support.support_primary ==
                   (support.active_at_threshold ? 1.0 : 0.0),
               "selected decision embedding was not binary");
        expect(support.observed_hypothesis_count == observed.size(),
               "selected observed-hypothesis count mismatch");
        if (observed.empty()) {
          expect(!support.closest_shift_linf_distance_defined,
                 "empty shift set acquired a closest-distance value");
        } else {
          std::uint64_t expected_distance =
              std::numeric_limits<std::uint64_t>::max();
          for (const auto &shift : observed) {
            const auto distance = static_cast<std::uint64_t>(
                std::llabs(selected_x - shift.x));
            expected_distance = std::min(expected_distance, distance);
          }
          expect(support.closest_shift_linf_distance_defined &&
                     support.closest_shift_linf_distance == expected_distance,
                 "selected closest lattice distance mismatch");
        }
        ++selected_cases;
      }
    }
  }
  expect(selected_cases == 80, "unexpected selected-support case count");

  const std::vector<Int3> extreme_observed{{
      std::numeric_limits<std::int64_t>::min(), 0, 0}};
  const Int3 extreme_selected{std::numeric_limits<std::int64_t>::max(), 0, 0};
  const auto extreme = evaluate_pbctopo_selected_compatibility_support(
      false, extreme_observed, extreme_selected);
  expect(extreme.closest_shift_linf_distance_defined &&
             extreme.closest_shift_linf_distance ==
                 std::numeric_limits<std::uint64_t>::max() &&
             !extreme.active_at_threshold,
         "integer-extreme lattice distance overflowed");

  std::cout << "PASS pbctopo_evidence_support_policy_smoke scored_cases="
            << scored_cases << " selected_cases=" << selected_cases << "\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ascii

& $Compiler -std=c++23 -O2 -Wall -Wextra -pedantic -Isrc $src `
    src/pbctopo/component_graph_math.cpp `
    src/pbctopo/component_metric.cpp `
    src/pbctopo/component_shift_geometry.cpp `
    src/pbctopo/lattice_math.cpp `
    -o $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile evidence support-policy smoke"
}
& $exe
if ($LASTEXITCODE -ne 0) {
    throw "Evidence support-policy smoke failed"
}
