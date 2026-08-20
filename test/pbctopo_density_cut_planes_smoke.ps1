param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_density_cut_planes_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_density_cut_planes_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_density_cut_planes_smoke.exe'

@'
#include "pbctopo/virtual_planes.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void require_true(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL pbctopo_density_cut_planes_smoke: " << message << '\n';
    std::exit(1);
  }
}

bool has_close(const std::vector<double> &values, double target, double eps) {
  return std::any_of(values.begin(), values.end(), [&](double value) {
    return std::fabs(value - target) <= eps;
  });
}

titan_pbctopo::PbctopoAtom make_atom(double sx, double sy, double sz) {
  titan_pbctopo::PbctopoAtom atom;
  atom.sx = sx;
  atom.sy = sy;
  atom.sz = sz;
  return atom;
}

} // namespace

int main() {
  using titan_pbctopo::PbctopoVirtualPlaneBuildOptions;
  using titan_pbctopo::PbctopoVirtualAssignmentSignatureTable;
  using titan_pbctopo::PbctopoVirtualStateTracker;
  using titan_pbctopo::build_dense_virtual_box_event_planes;
  using titan_pbctopo::build_density_cut_plane_shortlist_or_dense_fallback;
  using titan_pbctopo::build_density_virtual_box_event_planes;

  std::vector<titan_pbctopo::PbctopoAtom> clustered_atoms;
  for (double sx : {0.05, 0.06, 0.07, 0.08, 0.09}) {
    clustered_atoms.push_back(make_atom(sx, 0.20, 0.30));
  }
  for (double sx : {0.72, 0.73, 0.74, 0.75, 0.76}) {
    clustered_atoms.push_back(make_atom(sx, 0.20, 0.30));
  }

  const std::vector<double> sparse_x =
      build_density_virtual_box_event_planes(clustered_atoms, 'x', 64, 8);
  require_true(!sparse_x.empty(), "sparse planes should not be empty");
  require_true(sparse_x.size() <= 8, "sparse plane count should honor max");
  require_true(std::is_sorted(sparse_x.begin(), sparse_x.end()),
               "sparse planes should be sorted");
  require_true(std::adjacent_find(sparse_x.begin(), sparse_x.end(),
                                  [](double lhs, double rhs) {
                                    return std::fabs(lhs - rhs) <= 1.0e-9;
                                  }) == sparse_x.end(),
               "sparse planes should be unique");
  require_true(has_close(sparse_x, 0.5, 1.0e-12), "missing lower boundary");
  require_true(has_close(sparse_x, 1.5, 1.0e-12), "missing center boundary");
  require_true(has_close(sparse_x, 2.5, 1.0e-12), "missing upper boundary");
  require_true(has_close(sparse_x, 1.0, 0.20),
               "expected a sparse plane near the void between clusters");

  PbctopoVirtualStateTracker tracker;
  tracker.reset(sparse_x, sparse_x, sparse_x);
  require_true(tracker.try_insert(1.000000001, 1.0, 1.0),
               "first center should insert");
  require_true(tracker.try_insert(1.000000009, 1.0, 1.0),
               "same sparse bucket but distinct quantized center should insert");
  require_true(!tracker.try_insert(1.000000009, 1.0, 1.0),
               "identical quantized center should dedupe");
  require_true(tracker.collision_count() > 0,
               "tracker should count sparse bucket collisions");

  PbctopoVirtualStateTracker signature_tracker;
  signature_tracker.reset(sparse_x, sparse_x, sparse_x);
  require_true(signature_tracker.try_insert(1.05, 1.0, 1.0, 42),
               "first signature center should insert");
  require_true(!signature_tracker.try_insert(1.10, 1.0, 1.0, 42),
               "same sparse bucket and same shift signature should dedupe");
  require_true(signature_tracker.try_insert(1.10, 1.0, 1.0, 43),
               "same sparse bucket but different shift signature should insert");

  PbctopoVirtualAssignmentSignatureTable signature_table;
  signature_table.build(clustered_atoms);
  const auto signature_before_boundary =
      signature_table.state_signature(0.51, 1.0, 1.0);
  const auto signature_same_interval =
      signature_table.state_signature(0.54, 1.0, 1.0);
  const auto signature_after_boundary =
      signature_table.state_signature(0.56, 1.0, 1.0);
  require_true(signature_before_boundary == signature_same_interval,
               "signature should be stable inside a dense-equivalent interval");
  require_true(signature_before_boundary != signature_after_boundary,
               "signature should change after crossing an atom event boundary");

  PbctopoVirtualPlaneBuildOptions fallback_options;
  fallback_options.max_sparse_planes_per_axis = 2;
  const auto fallback_planes =
      build_density_cut_plane_shortlist_or_dense_fallback(
          clustered_atoms, fallback_options);
  const auto dense_x = build_dense_virtual_box_event_planes(clustered_atoms, 'x');
  require_true(fallback_planes.used_dense_fallback,
               "max_sparse_planes_per_axis < 3 should force dense fallback");
  require_true(fallback_planes.x == dense_x,
               "dense fallback x planes should match legacy dense planes");

  std::cout << "PASS pbctopo_density_cut_planes_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src src\pbctopo\virtual_planes.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to compile pbctopo density cut planes smoke"
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "pbctopo density cut planes smoke failed with exit code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
