param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_reference_strain_smoke')
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
#include "pbctopo/reference_strain.h"

#include <cmath>
#include <iostream>
#include <limits>

int main() {
  using titan_pbctopo::pbctopo_symmetric_reference_strain;
  if (pbctopo_symmetric_reference_strain(4.0, 4.0) != 0.0)
    return 1;
  const double stretched = pbctopo_symmetric_reference_strain(4.0, 1.0);
  const double compressed = pbctopo_symmetric_reference_strain(1.0, 4.0);
  if (std::fabs(stretched - compressed) > 1.0e-15)
    return 2;
  if (!(stretched > 0.0 && stretched < 1.0))
    return 3;
  if (!(pbctopo_symmetric_reference_strain(1.21, 1.0) < stretched))
    return 4;
  if (pbctopo_symmetric_reference_strain(
          std::numeric_limits<double>::infinity(), 1.0) != 0.0)
    return 5;
  std::cout << "PASS pbctopo_reference_strain_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO reference strain smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO reference strain smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
