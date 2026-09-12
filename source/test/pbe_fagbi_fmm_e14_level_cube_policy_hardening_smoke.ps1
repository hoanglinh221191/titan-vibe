param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e14_level_cube_policy_smoke'),
  [string]$Cxx = ''
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
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) { throw "missing PQR fixture: $pqr" }
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) { throw "requested C++ compiler not found: $Requested" }
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
  throw 'No C++ compiler found for KIFMM-lite E14 smoke.'
}

$compiler = Resolve-Cxx $Cxx
$exe = Join-Path $WorkRoot 'titan_fmm_kifmm_sampled_panel_matvec_smoke.exe'
$args = @(
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
  $exe
)
Push-Location $RepoRoot
try {
  & $compiler @args
  if ($LASTEXITCODE -ne 0) { throw 'failed to build E14 smoke driver' }
  $output = & $exe $Pqr1sh1 $Pqr1a63 '--e14-far-floor-smoke' 2>&1
  $code = $LASTEXITCODE
} finally {
  Pop-Location
}
$log = Join-Path $WorkRoot 'e14_far_floor_smoke.log'
$output | Set-Content -LiteralPath $log -Encoding ascii
if ($code -ne 0) {
  throw "E14 far-fraction floor smoke failed; see $log"
}
if (-not (($output | ForEach-Object { [string]$_ }) -match 'fallback_dominated')) {
  throw 'E14 smoke did not exercise fallback_dominated status'
}
Write-Host 'FMM E14 level-cube policy hardening smoke passed'
