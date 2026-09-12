param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\gbsa_sasa_convergence_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = (Get-Command clang++.exe -ErrorAction SilentlyContinue).Source
if (-not $Compiler) {
    $Compiler = (Get-Command g++.exe -ErrorAction SilentlyContinue).Source
}
if (-not $Compiler) {
    throw 'No C++ compiler found for GBSA SASA convergence smoke.'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ExePath = Join-Path $WorkRoot 'titan_gbsa_sasa_convergence_smoke.exe'
$compileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_sasa.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_gbsa_sasa_convergence_smoke.cpp'),
    '-o',
    $ExePath
)

& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to build titan_gbsa_sasa_convergence_smoke.'
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw 'titan_gbsa_sasa_convergence_smoke failed.'
}
