param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_lattice_basis_invariance_smoke')
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
#include "pbctopo/lattice_math.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_lattice_basis_invariance_smoke: " << message
              << "\n";
    std::exit(1);
  }
}

void invert(const double m[3][3], double inv[3][3]) {
  const double det =
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
  inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
  inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
  inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
  inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
  inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
  inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
}

double brute_d2(const double m[3][3], const std::array<double, 3> &r) {
  double best = std::numeric_limits<double>::infinity();
  for (int i = -8; i <= 8; ++i) {
    for (int j = -8; j <= 8; ++j) {
      for (int k = -8; k <= 8; ++k) {
        const double x = r[0] - i * m[0][0] - j * m[1][0] - k * m[2][0];
        const double y = r[1] - i * m[0][1] - j * m[1][1] - k * m[2][1];
        const double z = r[2] - i * m[0][2] - j * m[1][2] - k * m[2][2];
        best = std::min(best, x * x + y * y + z * z);
      }
    }
  }
  return best;
}

int main() {
  using titan_pbctopo::PbctopoLatticeMetric;
  const double basis_a[3][3] = {
      {10.0, 0.0, 0.0}, {9.0, 2.0, 0.0}, {3.0, 1.0, 8.0}};
  const double basis_b[3][3] = {
      {9.0, 2.0, 0.0}, {1.0, -2.0, 0.0}, {3.0, 1.0, 8.0}};
  double inverse_a[3][3]{};
  double inverse_b[3][3]{};
  invert(basis_a, inverse_a);
  invert(basis_b, inverse_b);
  const PbctopoLatticeMetric lattice_a(basis_a, inverse_a);
  const PbctopoLatticeMetric lattice_b(basis_b, inverse_b);
  expect(lattice_a.basis_reduced() || lattice_b.basis_reduced(),
         "skew lattice was not reduced");

  const std::vector<std::array<double, 3>> probes{
      {14.2, 3.7, 4.1}, {-7.1, 11.3, -5.2}, {25.4, -8.8, 17.6},
      {0.49, 0.51, 0.25}, {19.9, 4.1, 8.2}};
  for (const auto &probe : probes) {
    const auto a = lattice_a.nearest_image_cartesian(
        probe[0], probe[1], probe[2]);
    const auto b = lattice_b.nearest_image_cartesian(
        probe[0], probe[1], probe[2]);
    expect(a.valid && b.valid, "MIC result is invalid");
    expect(std::fabs(a.d2 - brute_d2(basis_a, probe)) < 1.0e-9,
           "reduced-basis MIC differs from brute-force CVP");
    expect(std::fabs(a.d2 - b.d2) < 1.0e-9,
           "equivalent lattice bases changed MIC distance");
    expect(std::fabs(a.x - b.x) < 1.0e-9 &&
               std::fabs(a.y - b.y) < 1.0e-9 &&
               std::fabs(a.z - b.z) < 1.0e-9,
           "equivalent lattice bases changed Cartesian MIC vector");
  }
  std::cout << "PASS pbctopo_lattice_basis_invariance_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO lattice basis invariance smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO lattice basis invariance smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
