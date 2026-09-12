param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\titan_simd_dispatch_override_smoke')
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
#include "titan_simd_dispatch.h"

int main(int argc, char **) {
  const bool expected = argc > 1;
  return titan_simd_dispatch::runtime_simd_forced_scalar() == expected ? 0 : 1;
}
'@ | Set-Content -LiteralPath $source -Encoding ASCII

& $Compiler '-std=c++23' '-O2' ('-I' + (Join-Path $repo 'src')) $source '-o' $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile SIMD override smoke"
}

Remove-Item Env:TITAN_FORCE_SCALAR -ErrorAction SilentlyContinue
& $exe
if ($LASTEXITCODE -ne 0) {
    throw 'SIMD override was active without an environment request'
}

$env:TITAN_FORCE_SCALAR = '1'
try {
    & $exe forced
    if ($LASTEXITCODE -ne 0) {
        throw 'TITAN_FORCE_SCALAR did not disable runtime SIMD dispatch'
    }
}
finally {
    Remove-Item Env:TITAN_FORCE_SCALAR -ErrorAction SilentlyContinue
}

Write-Host 'PASS titan_simd_dispatch_override_smoke'
