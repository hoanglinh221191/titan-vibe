param(
  [string]$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Source = Join-Path $RepoRoot 'test\titan_pbe_grid_solver_api_smoke.cpp'
$Exe = Join-Path $RepoRoot 'test\tmp\pbe_grid_solver_api_smoke.exe'

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

New-Item -ItemType Directory -Path (Split-Path -Parent $Exe) -Force | Out-Null

& $Compiler -std=c++23 -O2 "-I$(Join-Path $RepoRoot 'src')" $Source -o $Exe
if ($LASTEXITCODE -ne 0) {
  throw "GridSolver API smoke compile failed with code $LASTEXITCODE"
}

& $Exe
if ($LASTEXITCODE -ne 0) {
  throw "GridSolver API smoke executable failed with code $LASTEXITCODE"
}

Write-Host 'PBE GridSolver API smoke passed.'
