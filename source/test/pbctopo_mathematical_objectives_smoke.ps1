param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_mathematical_objectives_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/component_graph_math.h"
#include "pbctopo/steric_clash.h"
#include <cmath>
#include <iostream>

int main() {
  using namespace titan_pbctopo;
  PbctopoPairShiftScore small;
  small.valid = true;
  small.lhs_contact_atoms = 5;
  small.rhs_contact_atoms = 4;
  small.contact_pairs = 2;
  small.soft_penalty = 10.0;
  small.bbox_d2 = 4.0;
  small.shift_norm2 = 1.0;
  finalize_normalized_pair_shift_score(small, 10, 8, 4.0);

  PbctopoPairShiftScore large = small;
  large.lhs_contact_atoms = 50;
  large.rhs_contact_atoms = 40;
  large.contact_pairs = 20;
  large.soft_penalty = 100.0;
  finalize_normalized_pair_shift_score(large, 100, 80, 4.0);
  if (std::fabs(small.contact_fraction - 0.5) > 1.0e-12) return 1;
  if (std::fabs(small.soft_penalty_mean - 5.0) > 1.0e-12) return 2;
  if (std::fabs(small.quality_score - large.quality_score) > 1.0e-12) return 3;

  auto set = [](std::size_t a, std::size_t b, titan_pbctopo::Int3 delta) {
    PbctopoPairShiftCandidateSet result;
    result.a = a;
    result.b = b;
    result.candidates.push_back({delta, {}, 1.0});
    return result;
  };
  std::vector<PbctopoPairShiftCandidateSet> consistent{
      set(0, 1, {1, 0, 0}), set(1, 2, {1, 0, 0}), set(0, 2, {2, 0, 0})};
  const box unit_box(1.0, 1.0, 1.0, 90.0, 90.0, 90.0);
  const auto good =
      compute_fundamental_cycle_residuals(3, consistent, unit_box);
  if (good.fundamental_cycle_count != 1 ||
      good.fundamental_cycle_residual_max != 0.0) return 4;
  consistent[2].candidates[0].delta = {0, 0, 0};
  const auto bad =
      compute_fundamental_cycle_residuals(3, consistent, unit_box);
  if (bad.fundamental_cycle_count != 1 ||
      std::fabs(bad.fundamental_cycle_residual_rms - 2.0) > 1.0e-12) return 5;

  PbctopoAtom a, b;
  a.has_vdw_radius = b.has_vdw_radius = true;
  a.vdw_radius = b.vdw_radius = 1.5;
  std::vector<PbctopoAtom> atoms{a, b};
  prepare_pbctopo_steric_radii(atoms, 0.60);
  const auto threshold = pbctopo_pair_clash_threshold(atoms[0], atoms[1], 1.05);
  if (!threshold.used_vdw_radii || std::fabs(threshold.cutoff - 1.8) > 1.0e-12)
    return 6;
  if (!pbctopo_atoms_clash(1.5 * 1.5, atoms[0], atoms[1], 1.05)) return 7;
  if (pbctopo_atoms_clash(2.0 * 2.0, atoms[0], atoms[1], 1.05)) return 8;

  PbctopoAtom zero_lj, oxygen;
  zero_lj.has_vdw_radius = oxygen.has_vdw_radius = true;
  zero_lj.vdw_radius = 0.0;
  oxygen.vdw_radius = 1.5;
  std::vector<PbctopoAtom> zero_lj_pair{zero_lj, oxygen};
  prepare_pbctopo_steric_radii(zero_lj_pair, 0.60);
  const auto zero_lj_threshold = pbctopo_pair_clash_threshold(
      zero_lj_pair[0], zero_lj_pair[1], 1.05);
  if (!zero_lj_threshold.used_vdw_radii ||
      std::fabs(zero_lj_threshold.cutoff - 0.9) > 1.0e-12) return 9;
  if (pbctopo_atoms_clash(1.0, zero_lj_pair[0], zero_lj_pair[1], 1.05))
    return 10;

  std::cout << "PASS pbctopo_mathematical_objectives_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\component_graph_math.cpp `
        src\pbctopo\component_metric.cpp `
        src\pbctopo\lattice_math.cpp `
        src\pbctopo\steric_clash.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile mathematical objectives smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw "Mathematical objectives smoke failed: $LASTEXITCODE" }
}
finally { Pop-Location }
