param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_lattice_stability_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'pbctopo_lattice_stability_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_lattice_stability_smoke.exe'

@'
#include "BoxMatrix.h"
#include "pbctopo/lattice_math.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_lattice_stability_smoke: " << message << "\n";
    std::exit(1);
  }
}

void require_close(double actual, double expected, double tolerance,
                   const char *message) {
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << "FAIL pbctopo_lattice_stability_smoke: " << message
              << ": got " << actual << " expected " << expected << "\n";
    std::exit(1);
  }
}

void invert3(const double matrix[3][3], double inverse[3][3]) {
  const double det =
      matrix[0][0] * (matrix[1][1] * matrix[2][2] -
                      matrix[1][2] * matrix[2][1]) -
      matrix[0][1] * (matrix[1][0] * matrix[2][2] -
                      matrix[1][2] * matrix[2][0]) +
      matrix[0][2] * (matrix[1][0] * matrix[2][1] -
                      matrix[1][1] * matrix[2][0]);
  inverse[0][0] = (matrix[1][1] * matrix[2][2] -
                   matrix[1][2] * matrix[2][1]) /
                  det;
  inverse[0][1] = (matrix[0][2] * matrix[2][1] -
                   matrix[0][1] * matrix[2][2]) /
                  det;
  inverse[0][2] = (matrix[0][1] * matrix[1][2] -
                   matrix[0][2] * matrix[1][1]) /
                  det;
  inverse[1][0] = (matrix[1][2] * matrix[2][0] -
                   matrix[1][0] * matrix[2][2]) /
                  det;
  inverse[1][1] = (matrix[0][0] * matrix[2][2] -
                   matrix[0][2] * matrix[2][0]) /
                  det;
  inverse[1][2] = (matrix[0][2] * matrix[1][0] -
                   matrix[0][0] * matrix[1][2]) /
                  det;
  inverse[2][0] = (matrix[1][0] * matrix[2][1] -
                   matrix[1][1] * matrix[2][0]) /
                  det;
  inverse[2][1] = (matrix[0][1] * matrix[2][0] -
                   matrix[0][0] * matrix[2][1]) /
                  det;
  inverse[2][2] = (matrix[0][0] * matrix[1][1] -
                   matrix[0][1] * matrix[1][0]) /
                  det;
}

double brute_force_d2(const titan_pbctopo::PbctopoLatticeMetric &metric,
                      double sx, double sy, double sz) {
  double best = std::numeric_limits<double>::infinity();
  const int seed_x = static_cast<int>(std::floor(sx + 0.5));
  const int seed_y = static_cast<int>(std::floor(sy + 0.5));
  const int seed_z = static_cast<int>(std::floor(sz + 0.5));
  for (int nx = seed_x - 7; nx <= seed_x + 7; ++nx) {
    for (int ny = seed_y - 7; ny <= seed_y + 7; ++ny) {
      for (int nz = seed_z - 7; nz <= seed_z + 7; ++nz) {
        const auto cart = metric.scaled_to_cartesian(
            sx - static_cast<double>(nx), sy - static_cast<double>(ny),
            sz - static_cast<double>(nz));
        const double d2 = cart[0] * cart[0] + cart[1] * cart[1] +
                          cart[2] * cart[2];
        best = std::min(best, d2);
      }
    }
  }
  return best;
}

} // namespace

int main() {
  using namespace titan_pbctopo;

  const box orthogonal(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
  const PbctopoLatticeMetric orthogonal_metric(orthogonal);
  require(orthogonal_metric.valid(), "orthogonal lattice rejected");
  const auto plus_half = orthogonal_metric.nearest_image_scaled(0.5, 0.0, 0.0);
  const auto minus_half =
      orthogonal_metric.nearest_image_scaled(-0.5, 0.0, 0.0);
  require(plus_half.valid && minus_half.valid, "half-cell MIC failed");
  require(plus_half.ambiguous && minus_half.ambiguous,
          "half-cell ambiguity was not reported");
  require(plus_half.equivalent_minima == 2 &&
              minus_half.equivalent_minima == 2,
          "half-cell equivalent minimum count mismatch");
  require(plus_half.image_x == -minus_half.image_x,
          "half-cell image tie-break is not antisymmetric");
  require_close(plus_half.x, -minus_half.x, 1.0e-12,
                "half-cell Cartesian MIC is not antisymmetric");

  const box triclinic(10.0, 11.0, 12.0, 70.0, 80.0, 35.0);
  const PbctopoLatticeMetric metric(triclinic);
  require(metric.valid(), "triclinic lattice rejected");
  for (int idx = 0; idx < 160; ++idx) {
    const double sx = 2.3 * std::sin(0.37 * idx) + 0.013 * idx;
    const double sy = 1.9 * std::cos(0.23 * idx) - 0.009 * idx;
    const double sz = 2.1 * std::sin(0.19 * idx + 0.4);
    const auto mic = metric.nearest_image_scaled(sx, sy, sz);
    const auto reverse = metric.nearest_image_scaled(-sx, -sy, -sz);
    if (!mic.valid || !reverse.valid) {
      std::cerr << "triclinic failure idx=" << idx
                << " forward_status="
                << pbctopo_lattice_status_name(mic.lattice_status)
                << " reverse_status="
                << pbctopo_lattice_status_name(reverse.lattice_status)
                << " forward_budget=" << mic.node_budget_exhausted
                << " reverse_budget=" << reverse.node_budget_exhausted << "\n";
      require(false, "triclinic sphere decoder failed");
    }
    const double brute = brute_force_d2(metric, sx, sy, sz);
    require_close(mic.d2, brute, 1.0e-9 * std::max(1.0, brute),
                  "triclinic MIC differs from brute-force CVP");
    require_close(mic.x, -reverse.x, 1.0e-9,
                  "triclinic MIC x is not antisymmetric");
    require_close(mic.y, -reverse.y, 1.0e-9,
                  "triclinic MIC y is not antisymmetric");
    require_close(mic.z, -reverse.z, 1.0e-9,
                  "triclinic MIC z is not antisymmetric");
  }

  const BoxMatrix source(triclinic);
  const long long transform[3][3]{{1, 1, 0}, {0, 1, 1}, {0, 0, 1}};
  double transformed[3][3]{};
  for (int row = 0; row < 3; ++row)
    for (int col = 0; col < 3; ++col)
      for (int inner = 0; inner < 3; ++inner)
        transformed[row][col] += transform[row][inner] * source.m[inner][col];
  double transformed_inverse[3][3]{};
  invert3(transformed, transformed_inverse);
  const PbctopoLatticeMetric transformed_metric(transformed,
                                                 transformed_inverse);
  require(transformed_metric.valid(), "unimodularly transformed lattice rejected");
  const auto cart = metric.scaled_to_cartesian(1.37, -0.82, 0.44);
  const auto original_mic = metric.nearest_image_cartesian(
      cart[0], cart[1], cart[2]);
  const auto transformed_mic = transformed_metric.nearest_image_cartesian(
      cart[0], cart[1], cart[2]);
  require(original_mic.valid && transformed_mic.valid,
          "unimodular MIC failed");
  require_close(original_mic.d2, transformed_mic.d2, 1.0e-9,
                "MIC norm changed under unimodular basis transform");
  require_close(original_mic.x, transformed_mic.x, 1.0e-9,
                "MIC x changed under unimodular basis transform");
  require_close(original_mic.y, transformed_mic.y, 1.0e-9,
                "MIC y changed under unimodular basis transform");
  require_close(original_mic.z, transformed_mic.z, 1.0e-9,
                "MIC z changed under unimodular basis transform");

  const double ill_matrix[3][3]{{1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0e-12}};
  const double ill_inverse[3][3]{{1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0e12}};
  const PbctopoLatticeMetric ill_metric(ill_matrix, ill_inverse);
  require(!ill_metric.valid(), "ill-conditioned lattice was accepted");
  require(ill_metric.status() == PbctopoLatticeStatus::IllConditionedMatrix,
          "ill-conditioned lattice status mismatch");
  require(!ill_metric.nearest_image_scaled(0.1, 0.1, 0.1).valid,
          "invalid lattice produced a valid MIC");

  std::cout << "PASS pbctopo_lattice_stability_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\pbctopo\lattice_math.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile PBCTOPO lattice stability smoke'
}

& $exe
if ($LASTEXITCODE -ne 0) {
    throw "PBCTOPO lattice stability smoke failed with exit code $LASTEXITCODE"
}
