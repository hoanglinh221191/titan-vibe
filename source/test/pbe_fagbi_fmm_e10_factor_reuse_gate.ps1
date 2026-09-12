param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e10_factor_reuse'),
  [string]$BaselineTimingCsv = (Join-Path $PSScriptRoot 'tmp\fmm_e9_sampled_profile\fmm_e9_timing.csv'),
  [string]$Cxx = '',
  [switch]$Quick,
  [switch]$Full,
  [string[]]$RerunRows = @(),
  [switch]$RerunFailed,
  [switch]$Force,
  [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E10_FACTOR_REUSE_REPORT.txt'
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) { throw "missing PQR fixture: $pqr" }
}
if (-not $Quick -and -not $Full) { $Quick = $true }

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$BuildRoot = Join-Path $WorkRoot 'build'
$RunRoot = Join-Path $WorkRoot 'runs'
New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null

$MatrixPath = Join-Path $WorkRoot 'fmm_e10_matrix.csv'
$TimingPath = Join-Path $WorkRoot 'fmm_e10_timing.csv'
$CachePath = Join-Path $WorkRoot 'fmm_e10_cache.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e10_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e10_summary.dat'

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
  throw 'No C++ compiler found for KIFMM-lite E10 profile.'
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [double]$Theta, [string]$Role) {
  $sampleLabel = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  $thetaLabel = 'theta' + (($Theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture)) -replace '\.', 'p')
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_label = $sampleLabel
    sample_count = $SampleCount
    order = 6
    theta = $Theta
    role = $Role
    row_id = ('{0}_{1}_{2}_q6_cached_{3}' -f $Fixture, $sampleLabel, $SampleCount, $thetaLabel)
  }
}

# Literal row tokens for static contract:
# 1SH1_uniform_128_q6_cached, 1SH1_uniform_512_q6_cached.
$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 0.6 'quick_uniform'),
  (New-Row '1SH1' 'uniform_hash' 32 0.6 'sample_scaling')
)
$fullRows = @(
  (New-Row '1SH1' 'uniform_hash' 512 0.6 'uniform_512'),
  (New-Row '1SH1' 'inward_normal_groove' 512 0.6 'adversarial'),
  (New-Row '1SH1' 'closest_opposing_surface' 512 0.6 'adversarial'),
  (New-Row '1SH1' 'high_aspect' 512 0.6 'adversarial'),
  (New-Row '1A63' 'uniform_hash' 512 0.6 'uniform_512'),
  (New-Row '1A63' 'high_aspect' 512 0.6 'adversarial'),
  (New-Row '1SH1' 'uniform_hash' 128 0.4 'theta'),
  (New-Row '1SH1' 'uniform_hash' 128 0.6 'theta'),
  (New-Row '1SH1' 'uniform_hash' 128 0.7 'theta'),
  (New-Row '1SH1' 'uniform_hash' 32 0.6 'sample_scaling'),
  (New-Row '1SH1' 'uniform_hash' 128 0.6 'sample_scaling'),
  (New-Row '1SH1' 'uniform_hash' 512 0.6 'sample_scaling')
)
$rows = if ($Full) { $fullRows } else { $quickRows }
$seen = @{}
$rows = @($rows | Where-Object {
    if ($seen.ContainsKey($_.row_id)) { $false } else { $seen[$_.row_id] = $true; $true }
  })
if ($RerunRows.Count -gt 0) {
  $wanted = @{}
  foreach ($id in $RerunRows) { $wanted[$id] = $true }
  $rows = @($rows | Where-Object { $wanted.ContainsKey($_.row_id) })
  if ($rows.Count -eq 0) { throw 'RerunRows did not match any E10 row id' }
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @('fixture','sample','status','sampled_rows','order','aggregate_rms','equivalent_factor_cache_hits','equivalent_factor_cache_misses')) {
    if (-not $fields.Contains($required)) { throw "malformed E10 profile row, missing $required`: $Line" }
  }
  foreach ($optional in @(
      'equivalent_build_solve_seconds',
      'equivalent_matrix_build_seconds',
      'equivalent_condition_seconds',
      'equivalent_factor_seconds',
      'equivalent_rhs_build_seconds',
      'equivalent_solve_seconds')) {
    if (-not $fields.Contains($optional)) { $fields[$optional] = 'NaN' }
  }
  [pscustomobject]$fields
}

function Parse-ResultFromLog([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path)) { return $null }
  $line = Get-Content -LiteralPath $Path | Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } | Select-Object -Last 1
  if ($null -eq $line) { return $null }
  Parse-CaseLine ([string]$line)
}

function As-Double($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return [double]::NaN }
  [double]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}
function As-Long($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return [int64]0 }
  [int64]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}
function Is-AcceptableRow($Parsed) {
  if ($null -eq $Parsed) { return $false }
  return $Parsed.status -eq 'pass'
}

function Get-BaselineSeconds([string]$Fixture, [string]$Sample, [int]$SampleCount) {
  if (-not (Test-Path -LiteralPath $BaselineTimingCsv)) { return [double]::NaN }
  $label = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  $rowId = ('{0}_{1}_{2}_q6_theta0p6' -f $Fixture, $label, $SampleCount)
  $row = Import-Csv -LiteralPath $BaselineTimingCsv | Where-Object { $_.row_id -eq $rowId } | Select-Object -First 1
  if ($null -eq $row) { return [double]::NaN }
  return As-Double $row.equivalent_build_solve_seconds
}

$Compiler = Resolve-Cxx $Cxx
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_sampled_panel_matvec_smoke.exe'
$CompileArgs = @(
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
  & $Compiler @CompileArgs
  if ($LASTEXITCODE -ne 0) { throw 'failed to build E10 factor-reuse driver' }
} finally {
  Pop-Location
}

$parsedRows = @()
foreach ($row in $rows) {
  $rowDir = Join-Path $RunRoot $row.row_id
  New-Item -ItemType Directory -Path $rowDir -Force | Out-Null
  $rawLog = Join-Path $rowDir 'raw.log'
  $parsed = Parse-ResultFromLog $rawLog
  $shouldRun = $Force -or $null -eq $parsed
  if ($RerunFailed -and -not (Is-AcceptableRow $parsed)) { $shouldRun = $true }
  if ($RerunRows.Count -gt 0) { $shouldRun = $true }
  if (-not $shouldRun) {
    Write-Host "resume: $($row.row_id)"
  } else {
    Write-Host "run: $($row.row_id)"
    $runArgs = @(
      $Pqr1sh1,
      $Pqr1a63,
      '--e10-factor-profile',
      "--fixture=$($row.fixture)",
      "--sample=$($row.sample)",
      "--sample-count=$($row.sample_count)",
      "--order=$($row.order)",
      ("--theta={0}" -f $row.theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture))
    )
    Push-Location $RepoRoot
    try {
      $output = & $ExePath @runArgs 2>&1
      $code = $LASTEXITCODE
    } finally {
      Pop-Location
    }
    $output | Set-Content -LiteralPath $rawLog -Encoding ascii
    $parsed = Parse-ResultFromLog $rawLog
    if ($code -ne 0 -and $null -eq $parsed) {
      $parsed = [pscustomobject]@{
        fixture = $row.fixture
        sample = $row.sample
        status = 'runtime_failed'
        sampled_rows = $row.sample_count
        order = $row.order
        theta = $row.theta
        aggregate_rms = 'NaN'
        k1 = 'NaN'
        k2 = 'NaN'
        k3 = 'NaN'
        k4 = 'NaN'
        max_equiv_residual = 'NaN'
        max_condition = 'NaN'
        equivalent_build_solve_seconds = 'NaN'
        equivalent_matrix_build_seconds = 'NaN'
        equivalent_condition_seconds = 'NaN'
        equivalent_factor_seconds = 'NaN'
        equivalent_rhs_build_seconds = 'NaN'
        equivalent_solve_seconds = 'NaN'
        equivalent_factor_cache_hits = '0'
        equivalent_factor_cache_misses = '0'
      }
    }
  }
  if ($null -eq $parsed) { throw "no parseable result for $($row.row_id)" }
  $baseline = Get-BaselineSeconds $row.fixture $row.sample $row.sample_count
  $current = As-Double $parsed.equivalent_build_solve_seconds
  $speedup = if ([double]::IsNaN($baseline) -or $current -le 0.0) { [double]::NaN } else { $baseline / $current }
  $parsedRows += [pscustomobject]@{
    row_id = $row.row_id
    fixture = $row.fixture
    sample = $row.sample
    sample_count = $row.sample_count
    role = $row.role
    order = $row.order
    theta = $row.theta
    status = [string]$parsed.status
    aggregate_rms = As-Double $parsed.aggregate_rms
    max_channel_rms = [Math]::Max([Math]::Max((As-Double $parsed.k1), (As-Double $parsed.k2)), [Math]::Max((As-Double $parsed.k3), (As-Double $parsed.k4)))
    max_equiv_residual = As-Double $parsed.max_equiv_residual
    max_condition = As-Double $parsed.max_condition
    far_panel_pairs = As-Long $parsed.far_panel_pairs
    fallback_panel_pairs = As-Long $parsed.fallback_panel_pairs
    equivalent_build_solve_seconds = $current
    equivalent_matrix_build_seconds = As-Double $parsed.equivalent_matrix_build_seconds
    equivalent_condition_seconds = As-Double $parsed.equivalent_condition_seconds
    equivalent_factor_seconds = As-Double $parsed.equivalent_factor_seconds
    equivalent_rhs_build_seconds = As-Double $parsed.equivalent_rhs_build_seconds
    equivalent_solve_seconds = As-Double $parsed.equivalent_solve_seconds
    equivalent_factor_cache_hits = As-Long $parsed.equivalent_factor_cache_hits
    equivalent_factor_cache_misses = As-Long $parsed.equivalent_factor_cache_misses
    e9_baseline_equivalent_seconds = $baseline
    e10_vs_e9_equivalent_speedup = $speedup
    worklist_hash = [string]$parsed.worklist_hash
    policy_hash = [string]$parsed.policy_hash
  }
}

$parsedRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$parsedRows | Select-Object row_id,fixture,sample,sample_count,equivalent_build_solve_seconds,equivalent_matrix_build_seconds,equivalent_condition_seconds,equivalent_factor_seconds,equivalent_rhs_build_seconds,equivalent_solve_seconds,e9_baseline_equivalent_seconds,e10_vs_e9_equivalent_speedup | Export-Csv -LiteralPath $TimingPath -NoTypeInformation
$parsedRows | Select-Object row_id,equivalent_factor_cache_hits,equivalent_factor_cache_misses,max_condition,max_equiv_residual | Export-Csv -LiteralPath $CachePath -NoTypeInformation

$candidateRows = @($parsedRows | Where-Object { $_.order -eq 6 })
$accuracyPass = ($candidateRows.Count -gt 0) -and (@($candidateRows | Where-Object {
      $_.status -ne 'pass' -or $_.aggregate_rms -gt 1.0e-6 -or $_.max_channel_rms -gt 5.0e-6 -or $_.max_equiv_residual -gt 1.0e-10 -or $_.max_condition -gt 1.0e12
    }).Count -eq 0)
$coveragePass = (@($candidateRows | Where-Object { $_.far_panel_pairs -le 0 -or $_.fallback_panel_pairs -le 0 }).Count -eq 0)
$haveBaseline = @($candidateRows | Where-Object { -not [double]::IsNaN($_.e10_vs_e9_equivalent_speedup) }).Count -gt 0
$speedPass = $haveBaseline -and (@($candidateRows | Where-Object { -not [double]::IsNaN($_.e10_vs_e9_equivalent_speedup) -and $_.e10_vs_e9_equivalent_speedup -lt 3.0 }).Count -eq 0)
$decision = if (-not $coveragePass) {
  'fmm_e10_fallback_coverage_gap'
} elseif (-not $accuracyPass) {
  'fmm_e10_accuracy_gap'
} elseif ($Full -and $speedPass) {
  'fmm_e10_factor_reuse_pass'
} elseif ($accuracyPass -and $haveBaseline) {
  'fmm_e10_speed_gap'
} elseif ($accuracyPass) {
  'fmm_e10_profile_ready_without_baseline'
} else {
  'fmm_e10_profile_ready_with_accuracy_gaps'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Full) { 'full' } else { 'quick' }
  rows = $parsedRows.Count
  candidate_rows = $candidateRows.Count
  baseline_csv = $BaselineTimingCsv
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation

@(
  "decision $decision",
  "mode $(if ($Full) { 'full' } else { 'quick' })",
  "rows $($parsedRows.Count)",
  "candidate_rows $($candidateRows.Count)",
  "baseline_csv $BaselineTimingCsv"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

@"
TITAN FMM-E10 KIFMM-lite Equivalent Factor Reuse Gate

Decision: $decision
Mode: $(if ($Full) { 'full' } else { 'quick' })

Artifacts:
- Matrix: $MatrixPath
- Timing: $TimingPath
- Cache: $CachePath
- Decision: $DecisionPath
- Summary: $SummaryPath

Scope:
E10 is diagnostic-only. It optimizes and profiles q6 equivalent-system reuse
for the sampled KIFMM-lite path and does not wire KIFMM-lite into production
TABI matvec, GMRES, preconditioners, parser defaults, or manuscript claims.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E10 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Timing: $TimingPath"
Write-Host "Cache: $CachePath"
if (-not $AllowFailedDecision -and $decision -notin @('fmm_e10_factor_reuse_pass','fmm_e10_speed_gap','fmm_e10_profile_ready_without_baseline')) {
  throw "FMM E10 gate did not pass: $decision"
}
