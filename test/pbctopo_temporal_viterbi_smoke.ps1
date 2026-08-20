param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_temporal_viterbi_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/temporal_path.h"
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using namespace titan_pbctopo;
  std::vector<std::vector<PbctopoTemporalCandidate>> frames(3);
  frames[0] = {{{0.0, 0.0, 0.0}, 0.0}, {{0.5, 0.0, 0.0}, 0.3}};
  frames[1] = {{{0.5, 0.0, 0.0}, 0.0}, {{0.5, 0.0, 0.0}, 0.3}};
  frames[2] = {{{0.0, 0.0, 0.0}, 0.0}, {{0.5, 0.0, 0.0}, 0.3}};
  const auto path = solve_pbctopo_temporal_viterbi(frames, 4.0);
  if (!path.valid || path.selected.size() != 3) return 1;
  if (path.selected[0] != 1 || path.selected[1] != 0 || path.selected[2] != 1)
    return 2;
  if (path.transition_count != 0) return 4;
  if (std::fabs(pbctopo_temporal_center_jump2({2.45, 0.0, 0.0},
                                               {0.55, 0.0, 0.0}) - 0.01) >
      1.0e-12) return 3;

  auto assignment_candidate = [](int relative_x, double emission) {
    PbctopoTemporalCandidate candidate;
    candidate.emission_cost = emission;
    candidate.component_offsets = {{0, 0, 0}, {relative_x, 0, 0}};
    candidate.certified = true;
    return candidate;
  };
  std::vector<std::vector<PbctopoTemporalCandidate>> assignments(3);
  assignments[0] = {assignment_candidate(0, 0.0),
                    assignment_candidate(1, 0.3)};
  assignments[1] = {assignment_candidate(1, 0.0),
                    assignment_candidate(0, 0.3)};
  assignments[2] = {assignment_candidate(0, 0.0),
                    assignment_candidate(1, 0.3)};
  const auto assignment_path =
      solve_pbctopo_temporal_viterbi(assignments, 4.0);
  if (!assignment_path.valid || assignment_path.selected !=
                                    std::vector<std::size_t>{0, 1, 0})
    return 5;
  if (assignment_path.transition_count != 0)
    return 6;

  PbctopoTemporalCandidate gauge_lhs;
  gauge_lhs.component_offsets = {{0, 0, 0}, {1, 0, 0}};
  gauge_lhs.component_ids = {10, 20};
  gauge_lhs.component_weights = {100.0, 50.0};
  gauge_lhs.certified = true;
  PbctopoTemporalCandidate gauge_rhs = gauge_lhs;
  gauge_rhs.component_offsets = {{5, -3, 2}, {6, -3, 2}};
  if (pbctopo_temporal_assignment_jump2(gauge_lhs, gauge_rhs) > 1.0e-12)
    return 7;

  PbctopoTemporalCandidate reordered_rhs = gauge_lhs;
  reordered_rhs.component_offsets = {{1, 0, 0}, {0, 0, 0}};
  reordered_rhs.component_ids = {20, 10};
  reordered_rhs.component_weights = {50.0, 100.0};
  if (pbctopo_temporal_assignment_jump2(gauge_lhs, reordered_rhs) > 1.0e-12)
    return 8;

  const auto anisotropic_gram = pbctopo_temporal_normalized_gram(
      {150.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 30.0});
  PbctopoTemporalCandidate anisotropic_lhs = gauge_lhs;
  anisotropic_lhs.normalized_gram = anisotropic_gram;
  anisotropic_lhs.metric_available = true;
  PbctopoTemporalCandidate x_jump = anisotropic_lhs;
  x_jump.component_offsets = {{0, 0, 0}, {2, 0, 0}};
  PbctopoTemporalCandidate y_jump = anisotropic_lhs;
  y_jump.component_offsets = {{0, 0, 0}, {1, 1, 0}};
  if (!(pbctopo_temporal_assignment_jump2(anisotropic_lhs, x_jump) >
        pbctopo_temporal_assignment_jump2(anisotropic_lhs, y_jump)))
    return 9;
  std::cout << "PASS pbctopo_temporal_viterbi_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\temporal_path.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile temporal Viterbi smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw "Temporal Viterbi smoke failed: $LASTEXITCODE" }
}
finally { Pop-Location }
