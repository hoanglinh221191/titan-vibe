[CmdletBinding()]
param(
    [string]$Compiler = 'clang++',
    [string]$WorkRoot = 'test\tmp\p2b2b_retry_foundation',
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$source = Join-Path $RepoRoot 'test\titan_tabi_cartesian_retry_foundation_smoke.cpp'
$exe = Join-Path $WorkRoot 'titan_tabi_cartesian_retry_foundation_smoke.exe'
$needsBuild = $Force -or -not (Test-Path -LiteralPath $exe)
if (-not $needsBuild) {
    $needsBuild = (Get-Item -LiteralPath $source).LastWriteTimeUtc -gt
        (Get-Item -LiteralPath $exe).LastWriteTimeUtc
}
if ($needsBuild) {
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "Retry-foundation compiler was not found: $Compiler"
    }
    & $Compiler '-std=c++23' '-O2' "-I$RepoRoot\src" $source '-o' $exe
    if ($LASTEXITCODE -ne 0) { throw 'Retry-foundation unit build failed.' }
}

& $exe
if ($LASTEXITCODE -ne 0) { throw 'Retry-foundation unit smoke failed.' }
Write-Host 'PASS pbe_tabi_cartesian_retry_foundation_smoke'
