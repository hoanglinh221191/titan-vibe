param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\energy_validation_unit_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = (Get-Command clang++.exe -ErrorAction SilentlyContinue).Source
if (-not $Compiler) {
    $Compiler = (Get-Command g++.exe -ErrorAction SilentlyContinue).Source
}
if (-not $Compiler) {
    throw 'No C++ compiler found for energy validation smoke.'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ExePath = Join-Path $WorkRoot 'titan_energy_validation_smoke.exe'
$compileArgs = @(
    '-std=c++23',
    '-O2',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\titan_energy_validation.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_energy_validation_smoke.cpp'),
    '-o',
    $ExePath
)

& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to build titan_energy_validation_smoke.'
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw 'titan_energy_validation_smoke failed.'
}
