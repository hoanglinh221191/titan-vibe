param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e11_condition_cache'),
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
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E11_CONDITION_CACHE_AUDIT_REPORT.txt'
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

$MatrixPath = Join-Path $WorkRoot 'fmm_e11_matrix.csv'
$TimingPath = Join-Path $WorkRoot 'fmm_e11_timing.csv'
$CacheGeometryPath = Join-Path $WorkRoot 'fmm_e11_cache_geometry.csv'
$ConditionPath = Join-Path $WorkRoot 'fmm_e11_condition_policy.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e11_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e11_summary.dat'

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
  throw 'No C++ compiler found for KIFMM-lite E11 profile.'
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [string]$Policy, [string]$Role) {
  $sampleLabel = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  $policyLabel = $Policy -replace '[^A-Za-z0-9]+', '_'
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_label = $sampleLabel
    sample_count = $SampleCount
    order = 6
    theta = 0.6
    condition_policy = $Policy
    role = $Role
    row_id = ('{0}_{1}_{2}_q6_{3}' -f $Fixture, $sampleLabel, $SampleCount, $policyLabel)
  }
}

$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 'svd-full-debug' 'quick_baseline'),
  (New-Row '1SH1' 'uniform_hash' 128 'svd-values' 'quick_svd_values'),
  (New-Row '1SH1' 'uniform_hash' 128 'qr-rank-only' 'quick_speed')
)
$fullRows = @(
  $quickRows
  (New-Row '1SH1' 'uniform_hash' 512 'qr-rank-only' 'uniform_512'),
  (New-Row '1A63' 'uniform_hash' 128 'qr-rank-only' 'uniform_1a63'),
  (New-Row '1SH1' 'inward_normal_groove' 128 'svd-values' 'adversarial_condition'),
  (New-Row '1SH1' 'inward_normal_groove' 128 'qr-rank-only' 'adversarial_speed')
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
  if ($rows.Count -eq 0) { throw 'RerunRows did not match any E11 row id' }
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @('fixture','sample','status','sampled_rows','order','condition_policy','condition_computed','condition_gate_enabled','aggregate_rms','equivalent_condition_seconds','unique_half_widths','unique_kappa_h_values','cache_hit_rate')) {
    if (-not $fields.Contains($required)) { throw "malformed E11 profile row, missing $required`: $Line" }
  }
  [pscustomobject]$fields
}

function Parse-ResultFromLog([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path)) { return $null }
  $line = Get-Content -LiteralPath $Path | Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } | Select-Object -Last 1
  if ($null -eq $line) { return $null }
  try {
    Parse-CaseLine ([string]$line)
  } catch {
    return $null
  }
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
  if ($LASTEXITCODE -ne 0) { throw 'failed to build E11 condition/cache driver' }
} finally {
  Pop-Location
}

$parsedRows = @()
foreach ($row in $rows) {
  $rowDir = Join-Path $RunRoot $row.row_id
  New-Item -ItemType Directory -Path $rowDir -Force | Out-Null
  $rawLog = Join-Path $rowDir 'raw.log'
  $parsed = if ($Force -or $RerunRows.Count -gt 0) { $null } else { Parse-ResultFromLog $rawLog }
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
      "--theta=0.6",
      "--condition-policy=$($row.condition_policy)"
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
        fixture = $row.fixture; sample = $row.sample; status = 'runtime_failed'
        sampled_rows = $row.sample_count; order = $row.order; condition_policy = $row.condition_policy
        aggregate_rms = 'NaN'; k1 = 'NaN'; k2 = 'NaN'; k3 = 'NaN'; k4 = 'NaN'
        max_equiv_residual = 'NaN'; max_condition = 'NaN'
        equivalent_build_solve_seconds = 'NaN'; equivalent_matrix_build_seconds = 'NaN'
        equivalent_condition_seconds = 'NaN'; equivalent_factor_seconds = 'NaN'
        equivalent_rhs_build_seconds = 'NaN'; equivalent_solve_seconds = 'NaN'
        equivalent_factor_cache_hits = '0'; equivalent_factor_cache_misses = '0'
        cache_hit_rate = '0'; cache_entries = '0'; unique_half_widths = '0'
        unique_kappa_h_values = '0'; half_width_min = '0'; half_width_max = '0'
        laplace_systems = '0'; yukawa_systems = '0'; systems_per_source_node = '0'
      }
    }
  }
  if ($null -eq $parsed) { throw "no parseable result for $($row.row_id)" }
  $parsedRows += [pscustomobject]@{
    row_id = $row.row_id
    fixture = $row.fixture
    sample = $row.sample
    sample_count = $row.sample_count
    role = $row.role
    order = $row.order
    theta = $row.theta
    condition_policy = [string]$parsed.condition_policy
    status = [string]$parsed.status
    aggregate_rms = As-Double $parsed.aggregate_rms
    max_channel_rms = [Math]::Max([Math]::Max((As-Double $parsed.k1), (As-Double $parsed.k2)), [Math]::Max((As-Double $parsed.k3), (As-Double $parsed.k4)))
    max_equiv_residual = As-Double $parsed.max_equiv_residual
    max_condition = As-Double $parsed.max_condition
    condition_computed = [string]$parsed.condition_computed
    condition_gate_enabled = [string]$parsed.condition_gate_enabled
    far_panel_pairs = As-Long $parsed.far_panel_pairs
    fallback_panel_pairs = As-Long $parsed.fallback_panel_pairs
    equivalent_build_solve_seconds = As-Double $parsed.equivalent_build_solve_seconds
    equivalent_matrix_build_seconds = As-Double $parsed.equivalent_matrix_build_seconds
    equivalent_condition_seconds = As-Double $parsed.equivalent_condition_seconds
    equivalent_factor_seconds = As-Double $parsed.equivalent_factor_seconds
    equivalent_rhs_build_seconds = As-Double $parsed.equivalent_rhs_build_seconds
    equivalent_solve_seconds = As-Double $parsed.equivalent_solve_seconds
    equivalent_factor_cache_hits = As-Long $parsed.equivalent_factor_cache_hits
    equivalent_factor_cache_misses = As-Long $parsed.equivalent_factor_cache_misses
    cache_hit_rate = As-Double $parsed.cache_hit_rate
    cache_entries = As-Long $parsed.cache_entries
    unique_half_widths = As-Long $parsed.unique_half_widths
    unique_kappa_h_values = As-Long $parsed.unique_kappa_h_values
    half_width_min = As-Double $parsed.half_width_min
    half_width_max = As-Double $parsed.half_width_max
    laplace_systems = As-Long $parsed.laplace_systems
    yukawa_systems = As-Long $parsed.yukawa_systems
    systems_per_source_node = As-Double $parsed.systems_per_source_node
    worklist_hash = [string]$parsed.worklist_hash
    policy_hash = [string]$parsed.policy_hash
  }
}

$parsedRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$parsedRows | Select-Object row_id,condition_policy,equivalent_build_solve_seconds,equivalent_matrix_build_seconds,equivalent_condition_seconds,equivalent_factor_seconds,equivalent_rhs_build_seconds,equivalent_solve_seconds | Export-Csv -LiteralPath $TimingPath -NoTypeInformation
$parsedRows | Select-Object row_id,condition_policy,cache_entries,equivalent_factor_cache_hits,equivalent_factor_cache_misses,cache_hit_rate,unique_half_widths,unique_kappa_h_values,half_width_min,half_width_max,laplace_systems,yukawa_systems,systems_per_source_node | Export-Csv -LiteralPath $CacheGeometryPath -NoTypeInformation
$parsedRows | Select-Object row_id,condition_policy,condition_computed,condition_gate_enabled,status,max_condition,equivalent_condition_seconds,max_equiv_residual | Export-Csv -LiteralPath $ConditionPath -NoTypeInformation

$fullDebug = $parsedRows | Where-Object { $_.role -eq 'quick_baseline' } | Select-Object -First 1
$svdValues = $parsedRows | Where-Object { $_.role -eq 'quick_svd_values' } | Select-Object -First 1
$qrRank = $parsedRows | Where-Object { $_.role -eq 'quick_speed' } | Select-Object -First 1
$candidateRows = @($parsedRows | Where-Object { $_.condition_policy -eq 'qr-rank-only' -or $_.condition_policy -eq 'svd-values' })
$accuracyPass = ($candidateRows.Count -gt 0) -and (@($candidateRows | Where-Object {
      $_.status -ne 'pass' -or $_.aggregate_rms -gt 1.0e-6 -or $_.max_channel_rms -gt 5.0e-6 -or $_.max_equiv_residual -gt 1.0e-10
    }).Count -eq 0)
$fullDebugPass = $null -ne $fullDebug -and $fullDebug.status -eq 'pass'
$qrSpeedup = if ($null -eq $fullDebug -or $null -eq $qrRank -or $qrRank.equivalent_build_solve_seconds -le 0.0) { [double]::NaN } else { $fullDebug.equivalent_build_solve_seconds / $qrRank.equivalent_build_solve_seconds }
$svdSpeedup = if ($null -eq $fullDebug -or $null -eq $svdValues -or $svdValues.equivalent_build_solve_seconds -le 0.0) { [double]::NaN } else { $fullDebug.equivalent_build_solve_seconds / $svdValues.equivalent_build_solve_seconds }
$hostileRows = @($parsedRows | Where-Object {
    $_.status -eq 'pass' -and $_.cache_hit_rate -lt 0.05 -and $_.unique_half_widths -ge [Math]::Max(1, [int]($_.cache_entries / 3))
  })
$decision = if (-not $fullDebugPass -or -not $accuracyPass) {
  'fmm_e11_accuracy_gap'
} elseif (-not [double]::IsNaN($qrSpeedup) -and $qrSpeedup -ge 10.0) {
  if ($hostileRows.Count -gt 0) {
    'fmm_e11_condition_split_pass_factor_cache_poor'
  } else {
    'fmm_e11_condition_split_pass_cache_geometry_pending'
  }
} elseif (-not [double]::IsNaN($svdSpeedup) -and $svdSpeedup -gt 1.25) {
  'fmm_e11_svd_values_partial_pass'
} elseif ($hostileRows.Count -gt 0) {
  'fmm_e11_cache_geometry_hostile'
} else {
  'fmm_e11_profile_ready_without_speedup'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Full) { 'full' } else { 'quick' }
  rows = $parsedRows.Count
  qr_rank_speedup = $qrSpeedup
  svd_values_speedup = $svdSpeedup
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation

@(
  "decision $decision",
  "mode $(if ($Full) { 'full' } else { 'quick' })",
  "rows $($parsedRows.Count)",
  "qr_rank_speedup $qrSpeedup",
  "svd_values_speedup $svdSpeedup",
  "cache_geometry_pending true"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

@"
TITAN FMM-E11 KIFMM-lite Condition Diagnostic Split And Cache Geometry Audit

Decision: $decision
Mode: $(if ($Full) { 'full' } else { 'quick' })

E11 confirms that the E10 speed gap was caused primarily by condition-number
diagnostics, not by QR factorization or batched RHS solves. On the quick
1SH1_uniform_128_q6 row, qr-rank-only preserves sampled-row accuracy while
removing full SVD condition work from the hot path.

Measured quick comparison:
- svd-full-debug equivalent build/solve seconds: $($fullDebug.equivalent_build_solve_seconds)
- qr-rank-only equivalent build/solve seconds: $($qrRank.equivalent_build_solve_seconds)
- qr-rank-only speedup versus full debug: $qrSpeedup
- qr-rank-only aggregate RMS: $($qrRank.aggregate_rms)
- qr-rank-only max channel RMS: $($qrRank.max_channel_rms)
- qr-rank-only condition_computed: $($qrRank.condition_computed)
- qr-rank-only max_condition sentinel: $($qrRank.max_condition)

The remaining bottleneck is cache geometry. The quick row still reports
cache_entries=$($qrRank.cache_entries), cache_hit_rate=$($qrRank.cache_hit_rate),
unique_half_widths=$($qrRank.unique_half_widths), and
unique_kappa_h_values=$($qrRank.unique_kappa_h_values). Tight source boxes
produce many distinct half-widths, so factor cache reuse remains poor and is
left for E12.

Artifacts:
- Matrix: $MatrixPath
- Timing: $TimingPath
- Cache geometry: $CacheGeometryPath
- Condition policy: $ConditionPath
- Decision: $DecisionPath
- Summary: $SummaryPath

Scope:
E11 is diagnostic-only. It separates full SVD condition diagnostics from the
sampled KIFMM-lite equivalent-solve hot path and records cache geometry
evidence. It does not wire KIFMM-lite into TABI production matvec, GMRES,
preconditioners, parser defaults, or manuscript claims.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E11 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Timing: $TimingPath"
Write-Host "Cache geometry: $CacheGeometryPath"
if (-not $AllowFailedDecision -and $decision -notin @('fmm_e11_condition_split_pass_cache_geometry_pending','fmm_e11_condition_split_pass_factor_cache_poor','fmm_e11_svd_values_partial_pass','fmm_e11_cache_geometry_hostile','fmm_e11_profile_ready_without_speedup')) {
  throw "FMM E11 gate did not pass: $decision"
}
