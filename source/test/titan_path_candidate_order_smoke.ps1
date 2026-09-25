$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
$OutDir = Join-Path $PSScriptRoot 'tmp\titan_path_candidate_order_smoke'
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$Exe = Join-Path $OutDir 'titan_path_candidate_order_smoke.exe'
Push-Location $Root
try {
    & clang++ -std=c++20 -O2 -Isrc test\titan_path_candidate_order_smoke.cpp src\pathway_de_variants.cpp src\pull_de_variants.cpp -o $Exe
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $Exe
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
