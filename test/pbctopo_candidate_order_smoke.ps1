param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_candidate_order_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/candidate_order.h"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

int main() {
  using namespace titan_pbctopo;
  const PbctopoMetricResolution resolution{5.0, 0.08};
  std::vector<double> values;
  for (int i = 0; i <= 160; ++i)
    values.push_back(static_cast<double>(i));
  auto less = [&](double a, double b) {
    return compare_pbctopo_metric(a, b, resolution) < 0;
  };
  auto equivalent = [&](double a, double b) {
    return compare_pbctopo_metric(a, b, resolution) == 0;
  };
  for (double a : values) {
    if (less(a, a)) return 1;
    for (double b : values) {
      if (less(a, b) && less(b, a)) return 2;
      for (double c : values) {
        if (less(a, b) && less(b, c) && !less(a, c)) return 3;
        if (equivalent(a, b) && equivalent(b, c) && !equivalent(a, c))
          return 4;
      }
    }
  }
  if (!equivalent(0.0, 2.0) || equivalent(0.0, 8.0)) return 5;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (compare_pbctopo_metric(-inf, 0.0, resolution) >= 0) return 6;
  if (compare_pbctopo_metric(inf, 0.0, resolution) <= 0) return 7;
  if (compare_pbctopo_metric(nan, inf, resolution) <= 0) return 8;
  std::cout << "PASS pbctopo_candidate_order_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src src\pbctopo\candidate_order.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile comparator smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw "Comparator smoke failed: $LASTEXITCODE" }
}
finally { Pop-Location }
