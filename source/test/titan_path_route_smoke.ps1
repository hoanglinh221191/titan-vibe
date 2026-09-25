param(
    [string]$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\titan_path_route_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$resolvedTestRoot = [System.IO.Path]::GetFullPath((Join-Path $repo 'test'))
$resolvedWorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
if (-not $resolvedWorkRoot.StartsWith($resolvedTestRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "WorkRoot must remain under $resolvedTestRoot"
}
if (Test-Path -LiteralPath $resolvedWorkRoot) {
    Remove-Item -LiteralPath $resolvedWorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $resolvedWorkRoot | Out-Null
$exe = Join-Path $resolvedWorkRoot 'titan_path_route_smoke.exe'
$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source

Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc `
        test\titan_path_route_smoke.cpp src\titan_path\route.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile TITAN-Path route smoke test' }
    & $exe $resolvedWorkRoot
    if ($LASTEXITCODE -ne 0) { throw "TITAN-Path route smoke failed: $LASTEXITCODE" }
}
finally {
    Pop-Location
}
