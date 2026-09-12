param(
  [string]$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Source = Join-Path $RepoRoot 'test\titan_pbe_grid_ses_mask_smoke.cpp'
$Kernel = Join-Path $RepoRoot 'src\pbe\grid\ses_mask.cpp'
$Executor = Join-Path $RepoRoot 'src\titan_executor.cpp'
$Exe = Join-Path $RepoRoot 'test\tmp\pbe_grid_ses_mask_smoke.exe'

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

New-Item -ItemType Directory -Path (Split-Path -Parent $Exe) -Force | Out-Null

& $Compiler -std=c++23 -O2 -mavx2 "-I$(Join-Path $RepoRoot 'src')" `
  $Source $Kernel $Executor -o $Exe
if ($LASTEXITCODE -ne 0) {
  throw "Grid SES mask smoke compile failed with code $LASTEXITCODE"
}

& $Exe
if ($LASTEXITCODE -ne 0) {
  throw "Grid SES mask smoke executable failed with code $LASTEXITCODE"
}

Write-Host 'PBE Grid SES mask smoke passed.'
