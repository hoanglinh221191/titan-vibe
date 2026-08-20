param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_warm_state_smoke')
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
#include "pbctopo/warm_state.h"

#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_warm_state_smoke: " << message << "\n";
    std::exit(1);
  }
}

} // namespace

int main() {
  using namespace titan_pbctopo;

  const std::vector<std::vector<std::size_t>> layout_a{{0, 1}, {2, 3}};
  const std::vector<std::vector<std::size_t>> layout_b{{0, 2}, {1, 3}};
  const auto signature_a = build_component_layout_signature(layout_a);
  expect(signature_a.valid, "valid layout did not produce a signature");
  expect(same_component_layout(layout_a, signature_a),
         "identical component membership was rejected");
  expect(!same_component_layout(layout_b, signature_a),
         "same-sized components with different membership were accepted");

  PbctopoWarmState frame_start;
  frame_start.has_kinematic_latest_virtual_candidate = true;
  frame_start.kinematic_latest_virtual_frame = 7;

  PbctopoWarmState trial = frame_start;
  trial.previous_component_layout = signature_a;
  trial.previous_component_offsets = {{9, 0, 0}, {8, 0, 0}};
  trial.previous_stage1_delta_hints_by_component = {{{1, 0, 0}}, {}};
  trial.previous_stage2_delta_hints_by_component = {{{0, 1, 0}}, {}};
  trial.has_component_offsets = true;
  trial.has_component_trial_hints = true;

  PbctopoCutCandidate selected;
  selected.valid = true;
  selected.assignment.valid = true;
  selected.assignment.component_offsets = {{0, 0, 0}, {1, 0, 0}};
  selected.hard_feasibility.evaluated = true;
  selected.hard_feasibility.assignment_available = true;
  selected.hard_feasibility.steric_audit_available = true;
  selected.hard_feasibility.steric_decision_complete = true;
  selected.hard_feasibility.steric_pair_enumeration_complete = true;
  selected.hard_feasibility.steric_scan_complete = true;

  PbctopoCutCandidate incomplete_steric = selected;
  incomplete_steric.hard_feasibility.steric_scan_complete = false;
  PbctopoWarmState incomplete_rejected_state = trial;
  expect(!commit_pbctopo_certified_warm_state(
             incomplete_rejected_state, frame_start, trial,
             incomplete_steric, signature_a),
         "incomplete steric scan was accepted as certified");

  PbctopoWarmState committed = trial;
  expect(commit_pbctopo_certified_warm_state(
             committed, frame_start, trial, selected, signature_a),
         "certified candidate was not committed");
  expect(committed.previous_component_offsets ==
             selected.assignment.component_offsets,
         "trial offsets leaked into committed warm state");
  expect(!committed.has_component_trial_hints,
         "hints from a different trial assignment were retained");
  expect(committed.has_kinematic_latest_virtual_candidate &&
             committed.kinematic_latest_virtual_frame == 7,
         "committed kinematic history was not preserved transactionally");

  PbctopoWarmState matching_trial = trial;
  matching_trial.previous_component_offsets =
      selected.assignment.component_offsets;
  expect(commit_pbctopo_certified_warm_state(
             committed, frame_start, matching_trial, selected, signature_a),
         "matching certified trial was not committed");
  expect(committed.has_component_trial_hints,
         "matching trial hints were not retained");

  PbctopoCutCandidate rejected = selected;
  rejected.hard_feasibility.hard_clashes = 1;
  const PbctopoWarmState before_reject = committed;
  expect(!commit_pbctopo_certified_warm_state(
             committed, frame_start, trial, rejected, signature_a),
         "uncertified candidate entered warm state");
  expect(committed.previous_component_offsets ==
             before_reject.previous_component_offsets,
         "failed commit mutated persistent warm state");

  expect(pbctopo_warm_state_has_data(committed),
         "populated warm state was reported empty");
  reset_pbctopo_warm_state(committed);
  expect(!pbctopo_warm_state_has_data(committed),
         "warm reset left reusable state behind");

  std::cout << "PASS pbctopo_warm_state_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\warm_state.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO warm-state smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO warm-state smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
