$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
$OutDir = Join-Path $PSScriptRoot 'tmp\titan_path_navigation_smoke'
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$Exe = Join-Path $OutDir 'titan_path_navigation_smoke.exe'
& clang++ -std=c++20 -O2 -Isrc `
    test\titan_path_navigation_smoke.cpp `
    src\titan_path\navigation.cpp src\titan_path\geometry.cpp -o $Exe
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $Exe
exit $LASTEXITCODE
