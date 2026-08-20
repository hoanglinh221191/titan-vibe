param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_pair_graph_math_smoke')
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
#include "pbctopo/component_graph_math.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_pair_graph_math_smoke: " << message << "\n";
    std::exit(1);
  }
}

PbctopoPairShiftCandidateSet make_edge(std::size_t a, std::size_t b,
                                       Int3 delta, double weight) {
  PbctopoPairShiftCandidateSet edge;
  edge.a = a;
  edge.b = b;
  PbctopoPairShiftCandidate candidate;
  candidate.delta = delta;
  candidate.weight = weight;
  edge.candidates.push_back(candidate);
  return edge;
}

int main() {
  PbctopoPairShiftScore ab;
  ab.valid = true;
  ab.lhs_contact_atoms = 2;
  ab.rhs_contact_atoms = 8;
  ab.contact_pairs = 7;
  ab.soft_penalty = 1.75;
  ab.shift_norm2 = 3.0;
  finalize_normalized_pair_shift_score(ab, 4, 16, 4.5);

  PbctopoPairShiftScore ba;
  ba.valid = true;
  ba.lhs_contact_atoms = ab.rhs_contact_atoms;
  ba.rhs_contact_atoms = ab.lhs_contact_atoms;
  ba.contact_pairs = ab.contact_pairs;
  ba.soft_penalty = ab.soft_penalty;
  ba.shift_norm2 = ab.shift_norm2;
  finalize_normalized_pair_shift_score(ba, 16, 4, 4.5);
  expect(std::fabs(ab.contact_fraction - ba.contact_fraction) < 1.0e-15,
         "contact fraction depends on component ordering");
  expect(std::fabs(ab.soft_penalty_mean - ba.soft_penalty_mean) < 1.0e-15,
         "soft pair penalty depends on component ordering");
  expect(std::fabs(ab.quality_score - ba.quality_score) < 1.0e-12,
         "pair potential is not symmetric under edge reversal");

  PbctopoPairShiftScore second = ab;
  second.quality_score -= 50.0;
  const double confidence = pbctopo_pair_shift_confidence(ab, second);
  const double weight = pbctopo_pair_shift_candidate_weight(ab, confidence);
  expect(confidence > 0.0 && confidence < 1.0,
         "finite score gap produced invalid confidence");
  expect(weight >= 1.0 && weight <= 9.0,
         "cycle-sync pair weight is not bounded");

  PbctopoPairShiftScore incomplete = second;
  incomplete.complete = false;
  incomplete.quality_score = -std::numeric_limits<double>::infinity();
  expect(pbctopo_pair_shift_confidence(ab, incomplete) == 0.0,
         "incomplete runner-up produced artificial confidence");
  expect(pbctopo_pair_shift_candidate_weight(incomplete, 1.0) == 0.0,
         "incomplete score entered the cycle-sync shortlist");

  const box rectangular(10.0, 20.0, 30.0, 90.0, 90.0, 90.0);
  std::vector<PbctopoPairShiftCandidateSet> consistent{
      make_edge(0, 1, {1, 0, 0}, 10.0),
      make_edge(1, 2, {0, 1, 0}, 9.0),
      make_edge(0, 2, {1, 1, 0}, 1.0)};
  const auto zero = compute_fundamental_cycle_residuals(3, consistent,
                                                         rectangular);
  expect(zero.fundamental_cycle_count == 1 &&
             zero.fundamental_cycle_residual_max < 1.0e-12,
         "cycle-consistent graph has a nonzero residual");

  consistent.back() = make_edge(0, 2, {2, 1, 0}, 1.0);
  const auto x_error = compute_fundamental_cycle_residuals(3, consistent,
                                                            rectangular);
  expect(std::fabs(x_error.fundamental_cycle_residual_max - 10.0) < 1.0e-12,
         "cycle residual ignored the x lattice-vector length");

  consistent.back() = make_edge(0, 2, {1, 2, 0}, 1.0);
  const auto y_error = compute_fundamental_cycle_residuals(3, consistent,
                                                            rectangular);
  expect(std::fabs(y_error.fundamental_cycle_residual_max - 20.0) < 1.0e-12,
         "cycle residual ignored the Gram metric anisotropy");

  std::cout << "PASS pbctopo_pair_graph_math_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\component_graph_math.cpp `
        src\pbctopo\component_metric.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO pair-graph math smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO pair-graph math smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
