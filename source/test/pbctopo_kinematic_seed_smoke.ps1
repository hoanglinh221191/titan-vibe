param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_kinematic_seed_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_kinematic_seed_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_kinematic_seed_smoke.exe'

@'
#include "pbctopo/virtual_seed.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void expect_close(double actual, double expected, const char *message) {
  if (std::fabs(actual - expected) > 1.0e-9) {
    std::cerr << message << ": expected " << expected << ", got " << actual
              << '\n';
    std::exit(1);
  }
}

void expect_in_range(double value, const char *message) {
  if (value < 0.5 || value > 2.5) {
    std::cerr << message << ": center outside [0.5,2.5]: " << value << '\n';
    std::exit(1);
  }
}

} // namespace

int main() {
  using titan_pbctopo::PbctopoVirtualSeedCenter;
  using titan_pbctopo::extrapolate_virtual_seed_center;
  using titan_pbctopo::pbctopo_center_delta_mic;
  using titan_pbctopo::pbctopo_normalize_virtual_center;

  const double matrix[3][3] = {{10.0, 0.0, 0.0},
                               {0.0, 10.0, 0.0},
                               {0.0, 0.0, 10.0}};
  const double inverse[3][3] = {{0.1, 0.0, 0.0},
                                {0.0, 0.1, 0.0},
                                {0.0, 0.0, 0.1}};
  const titan_pbctopo::PbctopoLatticeMetric lattice(matrix, inverse);

  const PbctopoVirtualSeedCenter older{1.10, 1.40, 1.70};
  const PbctopoVirtualSeedCenter latest{1.20, 1.55, 1.90};
  const PbctopoVirtualSeedCenter predicted =
      extrapolate_virtual_seed_center(older, latest, lattice);
  expect_close(predicted.cx, 1.30, "linear drift x");
  expect_close(predicted.cy, 1.70, "linear drift y");
  expect_close(predicted.cz, 2.10, "linear drift z");

  const PbctopoVirtualSeedCenter wrap_older{2.45, 1.50, 1.50};
  const PbctopoVirtualSeedCenter wrap_latest{0.55, 1.50, 1.50};
  const PbctopoVirtualSeedCenter wrap_delta =
      pbctopo_center_delta_mic(wrap_latest, wrap_older, lattice);
  const PbctopoVirtualSeedCenter wrap_predicted =
      extrapolate_virtual_seed_center(wrap_older, wrap_latest, lattice);
  expect_close(wrap_delta.cx, 0.10, "wrapped MIC delta x");
  expect_close(wrap_predicted.cx, 0.65, "wrapped extrapolated x");

  const double normalized_low = pbctopo_normalize_virtual_center(-2.75);
  const double normalized_high = pbctopo_normalize_virtual_center(5.75);
  expect_in_range(normalized_low, "normalized low");
  expect_in_range(normalized_high, "normalized high");

  std::cout << "PASS pbctopo_kinematic_seed_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src src\pbctopo\virtual_seed.cpp src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to compile pbctopo kinematic seed smoke"
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "pbctopo kinematic seed smoke failed with exit code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
