param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_lattice_math_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/lattice_math.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

int main() {
  const double matrix[3][3] = {{10.0, 0.0, 0.0},
                               {9.0, 2.0, 0.0},
                               {0.0, 0.0, 10.0}};
  const double inverse[3][3] = {{0.1, 0.0, 0.0},
                                {-0.45, 0.5, 0.0},
                                {0.0, 0.0, 0.1}};
  titan_pbctopo::PbctopoLatticeMetric lattice(matrix, inverse);
  const auto mic = lattice.nearest_image_scaled(0.49, 0.49, 0.0);
  if (!mic.valid || mic.image_x != 0 || mic.image_y != -1 || mic.d2 > 1.2)
    return 1;

  double brute = std::numeric_limits<double>::infinity();
  for (int nx = -4; nx <= 4; ++nx)
    for (int ny = -4; ny <= 4; ++ny)
      for (int nz = -4; nz <= 4; ++nz) {
        const auto cart = lattice.scaled_to_cartesian(0.49 - nx, 0.49 - ny,
                                                       0.0 - nz);
        brute = std::min(brute, cart[0] * cart[0] + cart[1] * cart[1] +
                                   cart[2] * cart[2]);
      }
  if (std::fabs(mic.d2 - brute) > 1.0e-10) return 2;

  const auto cart_mic = lattice.nearest_image_cartesian(9.31, 0.98, 0.0);
  if (!cart_mic.valid || std::fabs(cart_mic.d2 - brute) > 1.0e-10) return 3;
  std::cout << "PASS pbctopo_lattice_math_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile lattice smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw "Lattice smoke failed: $LASTEXITCODE" }
}
finally { Pop-Location }
