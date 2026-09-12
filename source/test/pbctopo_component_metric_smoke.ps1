param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_component_metric_smoke')
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
#include "pbctopo/component_metric.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_component_metric_smoke: " << message << "\n";
    std::exit(1);
  }
}

int main() {
  using titan_pbctopo::Int3;
  using titan_pbctopo::checked_int3_add;
  using titan_pbctopo::pbctopo_component_shift_metric;
  using titan_pbctopo::rebase_component_offsets;

  box orthogonal(10.0, 12.0, 14.0, 90.0, 90.0, 90.0);
  const std::vector<Int3> offsets{{0, 0, 0}, {1, 0, 0}};
  const std::vector<Int3> gauge_shifted{{3, -2, 1}, {4, -2, 1}};
  const double base = pbctopo_component_shift_metric(orthogonal, offsets);
  const double shifted =
      pbctopo_component_shift_metric(orthogonal, gauge_shifted);
  expect(std::fabs(base - 50.0) < 1.0e-9,
         "orthogonal physical metric has wrong units");
  expect(std::fabs(base - shifted) < 1.0e-9,
         "global lattice translation changed metric");

  std::vector<Int3> canonical = gauge_shifted;
  rebase_component_offsets(canonical, orthogonal);
  expect(std::fabs(pbctopo_component_shift_metric(orthogonal, canonical) - base) <
             1.0e-9,
         "physical rebase changed metric");

  box triclinic(10.0, 11.0, 13.0, 74.0, 83.0, 68.0);
  const std::vector<Int3> skew_offsets{{-2, 1, 0}, {1, -1, 1}, {2, 0, -1}};
  std::vector<Int3> skew_shifted = skew_offsets;
  for (Int3 &offset : skew_shifted) {
    Int3 translated;
    expect(checked_int3_add(offset, Int3{5, -4, 3}, translated),
           "checked gauge translation unexpectedly overflowed");
    offset = translated;
  }
  expect(std::fabs(pbctopo_component_shift_metric(triclinic, skew_offsets) -
                   pbctopo_component_shift_metric(triclinic, skew_shifted)) <
             1.0e-9,
         "triclinic gauge invariance failed");

  Int3 overflow_result;
  expect(!checked_int3_add(
             Int3{std::numeric_limits<std::int64_t>::max(), 0, 0},
             Int3{1, 0, 0}, overflow_result),
         "checked Int3 addition accepted overflow");

  titan_pbctopo::PbctopoShift shift;
  shift.tx = 3.0;
  shift.ty = 4.0;
  expect(std::fabs(titan_pbctopo::pbctopo_translation_norm2(shift) - 25.0) <
             1.0e-12,
         "shift ordering still uses coefficient norm");
  std::cout << "PASS pbctopo_component_metric_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\component_metric.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO component metric smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO component metric smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
