param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_complete_uncertainty_band_policy_smoke'),
    [string]$Compiler = 'g++'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$source = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'
@'
#include "pbctopo/uncertainty_band.h"

int main() {
  using namespace titan_pbctopo;
  constexpr std::size_t feasible = 27;
  static_assert(pbctopo_should_certify_complete_uncertainty_band(true, true));
  static_assert(!pbctopo_should_certify_complete_uncertainty_band(false, true));
  static_assert(pbctopo_complete_uncertainty_band_exported(
      true, 0, feasible, feasible));
  static_assert(!pbctopo_complete_uncertainty_band_exported(
      true, 0, 1, feasible));
  return 0;
}
'@ | Set-Content -LiteralPath $source -Encoding ASCII

$compiler = (Get-Command $Compiler -ErrorAction Stop).Source
& $compiler -std=c++23 -Wall -Wextra -Werror -Isrc $source -o $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile uncertainty-band policy smoke: exit $LASTEXITCODE"
}
& $exe
if ($LASTEXITCODE -ne 0) {
    throw "Uncertainty-band policy smoke failed: exit $LASTEXITCODE"
}

Write-Host 'PASS pbctopo_complete_uncertainty_band_policy_smoke (27/27 complete; 1/27 incomplete)'
