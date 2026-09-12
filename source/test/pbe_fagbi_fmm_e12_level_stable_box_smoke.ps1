param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e12_level_stable_box_smoke'),
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
  throw 'No C++ compiler found for KIFMM-lite E12 smoke.'
}

function Build-Driver {
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E12 smoke driver' }
  } finally {
    Pop-Location
  }
  return $exe
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @('status','source_box_policy','aggregate_rms','k1','k2','k3','k4','nonproduction_fallback_pairs','unique_source_node_half_widths','level_cube_geometry','worklist_hash','policy_hash','panel_geometry_hash','sampled_target_hash')) {
    if (-not $fields.Contains($required)) { throw "malformed E12 smoke row, missing $required`: $Line" }
  }
  [pscustomobject]$fields
}

function As-Double($Value) {
  [double]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}
function As-Long($Value) {
  [int64]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}
function Max-Channel($Row) {
  [Math]::Max([Math]::Max((As-Double $Row.k1), (As-Double $Row.k2)), [Math]::Max((As-Double $Row.k3), (As-Double $Row.k4)))
}

function Run-Profile([string]$Exe, [string]$Policy, [string]$Label) {
  $log = Join-Path $WorkRoot "$Label.log"
  $args = @(
    $Pqr1sh1,
    $Pqr1a63,
    '--e10-factor-profile',
    '--fixture=1SH1',
    '--sample=uniform_hash',
    '--sample-count=32',
    '--order=6',
    '--theta=0.6',
    '--condition-policy=qr-rank-only',
    "--source-box-policy=$Policy"
  )
  Push-Location $RepoRoot
  try {
    $output = & $Exe @args 2>&1
    $code = $LASTEXITCODE
  } finally {
    Pop-Location
  }
  $output | Set-Content -LiteralPath $log -Encoding ascii
  if ($code -ne 0) { throw "E12 smoke row failed: $Label`n$($output -join [Environment]::NewLine)" }
  $line = $output | Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } | Select-Object -Last 1
  if ($null -eq $line) { throw "missing E12 smoke result line: $Label" }
  Parse-CaseLine ([string]$line)
}

$exe = Build-Driver
$tight = Run-Profile $exe 'tight-bbox' 'tight_bbox'
$levelA = Run-Profile $exe 'level-cube' 'level_cube_a'
$levelB = Run-Profile $exe 'level-cube' 'level_cube_b'

if ($levelA.source_box_policy -ne 'level-cube') { throw 'level-cube row did not report source_box_policy=level-cube' }
if ($levelA.level_cube_geometry -ne 'true') { throw 'level-cube row did not report level_cube_geometry=true' }
if ((As-Long $levelA.unique_source_node_half_widths) -ge (As-Long $tight.unique_source_node_half_widths)) {
  throw "level-cube did not reduce source half-width diversity: tight=$($tight.unique_source_node_half_widths), level=$($levelA.unique_source_node_half_widths)"
}
if ((As-Double $levelA.aggregate_rms) -gt 1.0e-6) { throw "level-cube aggregate RMS gate failed: $($levelA.aggregate_rms)" }
if ((Max-Channel $levelA) -gt 5.0e-6) { throw "level-cube channel RMS gate failed" }
if ((As-Long $levelA.nonproduction_fallback_pairs) -ne 0) { throw 'level-cube used nonproduction fallback' }
foreach ($field in @('policy_hash','worklist_hash','panel_geometry_hash','sampled_target_hash','unique_source_node_half_widths')) {
  if ([string]$levelA.$field -ne [string]$levelB.$field) { throw "level-cube repeat changed $field" }
}

Write-Host 'FMM E12 level-stable box smoke passed'
