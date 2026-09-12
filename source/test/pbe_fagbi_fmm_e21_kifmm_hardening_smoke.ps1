param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e21_kifmm_hardening'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot

function Resolve-Cxx([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        if (-not (Test-Path -LiteralPath $Requested)) {
            throw "Requested C++ compiler not found: $Requested"
        }
        return $Requested
    }
    foreach ($candidate in @(
            'C:\msys64\ucrt64\bin\clang++.exe',
            'C:\msys64\mingw64\bin\clang++.exe',
            'C:\msys64\ucrt64\bin\g++.exe',
            'C:\msys64\mingw64\bin\g++.exe')) {
        if (Test-Path -LiteralPath $candidate) { return $candidate }
    }
    $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    $cmd = Get-Command g++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    throw 'No C++ compiler found for KIFMM hardening smoke.'
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$Compiler = Resolve-Cxx $Cxx
$ExePath = Join-Path $WorkRoot 'titan_fmm_kifmm_hardening_smoke.exe'
$Args = @(
    '-std=c++23',
    '-O2',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_hardening_smoke.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    '-o',
    $ExePath
)

Push-Location $RepoRoot
try {
    & $Compiler @Args
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build titan_fmm_kifmm_hardening_smoke'
    }
    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'titan_fmm_kifmm_hardening_smoke failed'
    }
} finally {
    Pop-Location
}
