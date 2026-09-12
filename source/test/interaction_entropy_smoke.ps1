param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\interaction_entropy_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = (Get-Command clang++.exe -ErrorAction SilentlyContinue).Source
if (-not $Compiler) {
    $Compiler = (Get-Command g++.exe -ErrorAction SilentlyContinue).Source
}
if (-not $Compiler) {
    throw 'No C++ compiler found for interaction entropy smoke.'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ExePath = Join-Path $WorkRoot 'titan_interaction_entropy_smoke.exe'
& $Compiler -std=c++23 -O2 ('-I' + (Join-Path $RepoRoot 'src')) `
    (Join-Path $RepoRoot 'src\titan_interaction_entropy.cpp') `
    (Join-Path $RepoRoot 'test\titan_interaction_entropy_smoke.cpp') `
    -o $ExePath
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to build titan_interaction_entropy_smoke.'
}
& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw 'titan_interaction_entropy_smoke failed.'
}
