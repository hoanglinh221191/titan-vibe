param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_adaptive_image_shell_smoke')
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
#include "pbctopo/component_shift_domain.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <tuple>

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_adaptive_image_shell_smoke: " << message << "\n";
    std::exit(1);
  }
}

int main() {
  using namespace titan_pbctopo;
  box Box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
  const auto shell1 = build_component_image_shifts(Box, 1);
  const auto shell2 = build_component_image_shifts(Box, 2);
  const auto shell3 = build_component_image_shifts(Box, 3);
  expect(shell1.size() == 27, "radius-1 shell is not the legacy 27 domain");
  expect(shell2.size() == 125, "radius-2 shell is not the complete 125 domain");
  expect(shell3.size() == 343, "radius-3 shell is not the complete 343 domain");
  const auto indices_match_storage = [](const auto &domain) {
    for (std::size_t index = 0; index < domain.size(); ++index) {
      if (domain[index].shift_index != index)
        return false;
    }
    return true;
  };
  expect(indices_match_storage(shell1) && indices_match_storage(shell2) &&
             indices_match_storage(shell3),
         "shift_index does not match the post-sort storage position");
  expect(shell1.front().ix == 0 && shell1.front().iy == 0 &&
             shell1.front().iz == 0,
         "identity shift is not ranked first");
  bool analytic_relative = false;
  const auto relative1 =
      build_component_relative_shift_domain(shell1, &analytic_relative);
  expect(analytic_relative, "complete shell did not use analytic S-S build");
  expect(relative1.size() == 125,
         "radius-1 relative domain is not the complete 125 domain");
  const auto zero_relative = std::find_if(
      relative1.begin(), relative1.end(), [](const auto &entry) {
        return entry.shift.ix == 0 && entry.shift.iy == 0 &&
               entry.shift.iz == 0;
      });
  expect(zero_relative != relative1.end() &&
             zero_relative->multiplicity == shell1.size(),
         "zero relative shift multiplicity is incorrect");

  std::set<std::tuple<int, int, int>> unique;
  for (const auto &shift : shell2) {
    unique.emplace(shift.ix, shift.iy, shift.iz);
    expect(std::max({std::abs(shift.ix), std::abs(shift.iy),
                     std::abs(shift.iz)}) <= 2,
           "shift lies outside requested shell");
  }
  expect(unique.size() == shell2.size(), "adaptive shell contains duplicates");
  const auto corner = unique.find({2, -2, 2});
  expect(corner != unique.end(), "radius-2 corner shift is missing");

  PbctopoShift minimum_shift;
  minimum_shift.ix = std::numeric_limits<int>::min();
  PbctopoShift maximum_shift;
  maximum_shift.ix = std::numeric_limits<int>::max();
  expect(build_component_relative_shift_domain(
             {minimum_shift, maximum_shift})
             .empty(),
         "unrepresentable relative shift domain did not fail closed");
  Int3 overflow_result{};
  expect(!checked_int3_subtract(
             Int3{std::numeric_limits<std::int64_t>::min(), 0, 0},
             Int3{std::numeric_limits<std::int64_t>::max(), 0, 0},
             overflow_result),
         "checked Int3 subtraction accepted signed overflow");
  std::cout << "PASS pbctopo_adaptive_image_shell_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\component_shift_domain.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO adaptive image shell smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO adaptive image shell smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
