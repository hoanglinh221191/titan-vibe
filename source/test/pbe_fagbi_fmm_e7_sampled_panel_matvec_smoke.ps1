param(
    [string]$Pqr1sh1 = '',
    [string]$Pqr1a63 = '',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e7_sampled_panel_matvec'),
    [string]$Cxx = '',
    [int]$SampleCount = 0,
    [switch]$Full
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
    $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
    $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}

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
    throw 'No C++ compiler found for KIFMM-lite E7 smoke.'
}

foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
    if (-not (Test-Path -LiteralPath $pqr)) {
        throw "Protein PQR fixture not found: $pqr"
    }
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$Compiler = Resolve-Cxx $Cxx
$ExePath = Join-Path $WorkRoot 'titan_fmm_kifmm_sampled_panel_matvec_smoke.exe'
$Args = @(
    '-std=c++23',
    '-O2',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
)

Push-Location $RepoRoot
try {
    & $Compiler @Args
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build titan_fmm_kifmm_sampled_panel_matvec_smoke'
    }
    if ($SampleCount -le 0) {
        $SampleCount = $(if ($Full) { 512 } else { 32 })
    }
    $RunArgs = @($Pqr1sh1, $Pqr1a63, "--sample-count=$SampleCount")
    if ($Full) { $RunArgs += '--full' }
    & $ExePath @RunArgs
    if ($LASTEXITCODE -ne 0) {
        throw 'titan_fmm_kifmm_sampled_panel_matvec_smoke failed'
    }
} finally {
    Pop-Location
}
