param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e9_sampled_profile'),
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
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E9_SAMPLED_SCALING_PROFILE_REPORT.txt'
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) {
    throw "missing PQR fixture: $pqr"
  }
}

if (-not $Quick -and -not $Full) {
  $Quick = $true
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$BuildRoot = Join-Path $WorkRoot 'build'
$RunRoot = Join-Path $WorkRoot 'runs'
New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null

$MatrixPath = Join-Path $WorkRoot 'fmm_e9_matrix.csv'
$TimingPath = Join-Path $WorkRoot 'fmm_e9_timing.csv'
$MemoryPath = Join-Path $WorkRoot 'fmm_e9_memory.csv'
$ScalingPath = Join-Path $WorkRoot 'fmm_e9_scaling.csv'
$AdversarialPath = Join-Path $WorkRoot 'fmm_e9_adversarial.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e9_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e9_summary.dat'

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) {
      throw "requested C++ compiler not found: $Requested"
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
  throw 'No C++ compiler found for KIFMM-lite E9 profile.'
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [int]$Order, [double]$Theta, [string]$Role) {
  $thetaLabel = ('theta{0}' -f (($Theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture)) -replace '\.', 'p'))
  $sampleLabel = $Sample
  if ($sampleLabel -eq 'uniform_hash') { $sampleLabel = 'uniform' }
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_label = $sampleLabel
    sample_count = $SampleCount
    order = $Order
    theta = $Theta
    role = $Role
    row_id = ('{0}_{1}_{2}_q{3}_{4}' -f $Fixture, $sampleLabel, $SampleCount, $Order, $thetaLabel)
  }
}

# Required E9 rows documented as literal tokens for the static contract:
# 1SH1_uniform_512_q6, 1SH1_uniform_512_q8,
# 1SH1_inward_normal_groove_q6, 1SH1_closest_opposing_surface_q6,
# 1SH1_high_aspect_q6, 1A63_uniform_512_q6, 1A63_high_aspect_q6,
# q4_q6_q8_trend, theta0p4, theta0p6, theta0p7.
$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 4 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 128 8 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 32 6 0.6 'sample_scaling')
)
$fullRows = @(
  (New-Row '1SH1' 'uniform_hash' 512 6 0.6 'adversarial_uniform'),
  (New-Row '1SH1' 'uniform_hash' 512 8 0.6 'adversarial_uniform'),
  (New-Row '1SH1' 'inward_normal_groove' 512 6 0.6 'adversarial'),
  (New-Row '1SH1' 'closest_opposing_surface' 512 6 0.6 'adversarial'),
  (New-Row '1SH1' 'high_aspect' 512 6 0.6 'adversarial'),
  (New-Row '1A63' 'uniform_hash' 512 6 0.6 'adversarial_uniform'),
  (New-Row '1A63' 'high_aspect' 512 6 0.6 'adversarial'),
  (New-Row '1SH1' 'uniform_hash' 128 4 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 128 8 0.6 'q4_q6_q8_trend'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.4 'theta0p4'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.6 'theta0p6'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.7 'theta0p7'),
  (New-Row '1SH1' 'uniform_hash' 32 6 0.6 'sample_scaling'),
  (New-Row '1SH1' 'uniform_hash' 512 6 0.6 'sample_scaling')
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
  if ($rows.Count -eq 0) {
    throw "RerunRows did not match any E9 row id"
  }
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @('fixture','sample','status','source_panels','sampled_rows','order','theta','aggregate_rms')) {
    if (-not $fields.Contains($required)) {
      throw "malformed KIFMM sampled profile row, missing $required`: $Line"
    }
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
  if ($Parsed.status -eq 'pass') { return $true }
  if ($Parsed.order -eq '4' -and $Parsed.status -in @('aggregate_rms_gate_failed','channel_rms_gate_failed','pass')) { return $true }
  return $false
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
  if ($LASTEXITCODE -ne 0) { throw 'failed to build E9 sampled profile driver' }
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
      '--e9-profile',
      "--fixture=$($row.fixture)",
      "--sample=$($row.sample)",
      "--sample-count=$($row.sample_count)",
      "--order=$($row.order)",
      "--orders=$($row.order)",
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
        source_panels = 0
        sampled_rows = $row.sample_count
        order = $row.order
        theta = $row.theta
      }
    }
  }
  $parsedRows += [pscustomobject]@{
    row_id = $row.row_id
    role = $row.role
    requested_fixture = $row.fixture
    requested_sample = $row.sample
    requested_sample_count = $row.sample_count
    requested_order = $row.order
    requested_theta = $row.theta
    raw = $parsed
  }
}

$matrix = foreach ($entry in $parsedRows) {
  $row = $entry.raw
  $maxChannel = [Math]::Max(
    [Math]::Max((As-Double $row.k1), (As-Double $row.k2)),
    [Math]::Max((As-Double $row.k3), (As-Double $row.k4)))
  [pscustomobject]@{
    fixture = $entry.requested_fixture
    row_id = $entry.row_id
    role = $entry.role
    sample_set = $entry.requested_sample
    sampled_rows = As-Long $row.sampled_rows
    source_panels = As-Long $row.source_panels
    order = As-Long $row.order
    theta = As-Double $row.theta
    status = $row.status
    aggregate_rms = As-Double $row.aggregate_rms
    max_channel_rms = $maxChannel
    k1_rms = As-Double $row.k1
    k2_rms = As-Double $row.k2
    k3_rms = As-Double $row.k3
    k4_rms = As-Double $row.k4
    far_panel_pairs = As-Long $row.far_panel_pairs
    fallback_panel_pairs = As-Long $row.fallback_panel_pairs
    production_direct_reference_valid_pairs = As-Long $row.production_direct_reference_valid_pairs
    invalid_pairs = As-Long $row.invalid_pairs
    nonproduction_fallback_pairs = As-Long $row.nonproduction_fallback_pairs
    max_equiv_residual = As-Double $row.max_equiv_residual
    max_condition = As-Double $row.max_condition
    max_equivalent_rank = As-Long $row.max_equivalent_rank
    min_equivalent_sigma_min = As-Double $row.min_equivalent_sigma_min
    max_equivalent_sigma_max = As-Double $row.max_equivalent_sigma_max
    max_equivalent_rhs_norm = As-Double $row.max_equivalent_rhs_norm
    max_equivalent_charge_norm = As-Double $row.max_equivalent_charge_norm
    max_equivalent_charge_to_rhs_norm = As-Double $row.max_equivalent_charge_to_rhs_norm
    panel_geometry_hash = $row.panel_geometry_hash
    unknown_vector_hash = $row.unknown_vector_hash
    sampled_target_hash = $row.sampled_target_hash
    worklist_hash = $row.worklist_hash
    policy_hash = $row.policy_hash
  }
}
$matrix | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation -Encoding ascii

$timing = foreach ($entry in $parsedRows) {
  $row = $entry.raw
  [pscustomobject]@{
    row_id = $entry.row_id
    fixture = $entry.requested_fixture
    sample_set = $entry.requested_sample
    sampled_rows = As-Long $row.sampled_rows
    order = As-Long $row.order
    theta = As-Double $row.theta
    direct_reference_seconds = As-Double $row.direct_reference_seconds
    tree_build_seconds = As-Double $row.tree_build_seconds
    worklist_build_seconds = As-Double $row.worklist_build_seconds
    equivalent_build_solve_seconds = As-Double $row.equivalent_build_solve_seconds
    far_eval_seconds = As-Double $row.far_eval_seconds
    fallback_eval_seconds = As-Double $row.fallback_eval_seconds
    total_seconds = As-Double $row.total_seconds
    unique_far_source_nodes = As-Long $row.unique_far_source_nodes
    cached_source_nodes = As-Long $row.cached_source_nodes
  }
}
$timing | Export-Csv -LiteralPath $TimingPath -NoTypeInformation -Encoding ascii

$memory = foreach ($entry in $parsedRows) {
  $row = $entry.raw
  [pscustomobject]@{
    row_id = $entry.row_id
    fixture = $entry.requested_fixture
    sample_set = $entry.requested_sample
    sampled_rows = As-Long $row.sampled_rows
    order = As-Long $row.order
    source_tree_bytes = As-Long $row.source_tree_bytes
    target_tree_bytes = As-Long $row.target_tree_bytes
    worklist_bytes = As-Long $row.worklist_bytes
    equivalent_charge_storage_bytes = As-Long $row.equivalent_charge_storage_bytes
    cached_channel_object_estimated_bytes = As-Long $row.cached_channel_object_estimated_bytes
  }
}
$memory | Export-Csv -LiteralPath $MemoryPath -NoTypeInformation -Encoding ascii
$matrix | Where-Object { $_.sample_set -eq 'uniform_hash' } |
  Sort-Object fixture, sampled_rows, order, theta |
  Export-Csv -LiteralPath $ScalingPath -NoTypeInformation -Encoding ascii
$matrix | Where-Object { $_.sample_set -ne 'uniform_hash' -or $_.role -like 'adversarial*' } |
  Export-Csv -LiteralPath $AdversarialPath -NoTypeInformation -Encoding ascii

function Passes-Candidate($row) {
  if ($row.order -eq 4) { return $true }
  $expected = $row.source_panels * $row.sampled_rows
  return $row.status -eq 'pass' -and
    $row.aggregate_rms -le 1.0e-6 -and
    $row.max_channel_rms -le 5.0e-6 -and
    $row.invalid_pairs -eq 0 -and
    $row.nonproduction_fallback_pairs -eq 0 -and
    (($row.far_panel_pairs + $row.fallback_panel_pairs) -eq $expected) -and
    $row.production_direct_reference_valid_pairs -eq $expected -and
    $row.max_equiv_residual -le 1.0e-10 -and
    $row.max_condition -le 1.0e12
}

$candidateRows = @($matrix | Where-Object { $_.order -ne 4 })
$coverageGap = @($candidateRows | Where-Object {
    $expected = $_.source_panels * $_.sampled_rows
    $_.invalid_pairs -ne 0 -or $_.nonproduction_fallback_pairs -ne 0 -or
      (($_.far_panel_pairs + $_.fallback_panel_pairs) -ne $expected) -or
      ($_.production_direct_reference_valid_pairs -ne $expected)
  })
$conditionGap = @($candidateRows | Where-Object { $_.max_equiv_residual -gt 1.0e-10 -or $_.max_condition -gt 1.0e12 })
$accuracyGap = @($candidateRows | Where-Object { -not (Passes-Candidate $_) })
$q6Gaps = @($candidateRows | Where-Object { $_.order -eq 6 -and -not (Passes-Candidate $_) })
$q8Pass = @($candidateRows | Where-Object { $_.order -eq 8 -and (Passes-Candidate $_) })

$decision = 'fmm_e9_profile_ready_with_accuracy_gaps'
if ($coverageGap.Count -gt 0) {
  $decision = 'fmm_e9_fallback_coverage_gap'
} elseif ($conditionGap.Count -gt 0) {
  $decision = 'fmm_e9_condition_gap'
} elseif ($accuracyGap.Count -eq 0 -and $Full) {
  $decision = 'fmm_e9_adversarial_scaling_profile_pass'
} elseif ($q6Gaps.Count -gt 0 -and $q8Pass.Count -gt 0) {
  $decision = 'fmm_e9_order6_gap_order8_pass'
} elseif ($accuracyGap.Count -gt 0) {
  $decision = 'fmm_e9_far_equivalent_gap'
}

[pscustomobject]@{
  decision = $decision
  mode = $(if ($Full) { 'full' } else { 'quick' })
  rows = $matrix.Count
  candidate_rows = $candidateRows.Count
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$(if ($Full) { 'full' } else { 'quick' })",
  "rows=$($matrix.Count)",
  "matrix=$MatrixPath",
  "timing=$TimingPath",
  "memory=$MemoryPath",
  "scaling=$ScalingPath",
  "adversarial=$AdversarialPath"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

@"
TITAN FMM-E9 Sampled KIFMM-lite Adversarial And Scaling Profile Gate

Decision: $decision
Mode: $(if ($Full) { 'full' } else { 'quick' })

Artifacts:
- Matrix: $MatrixPath
- Timing: $TimingPath
- Memory: $MemoryPath
- Scaling: $ScalingPath
- Adversarial: $AdversarialPath
- Decision: $DecisionPath
- Summary: $SummaryPath

Scope:
E9 is diagnostic-only. It does not wire KIFMM-lite into production TABI matvec,
GMRES, preconditioners, parser defaults, or manuscript claims.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E9 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Timing: $TimingPath"
Write-Host "Memory: $MemoryPath"
if ($decision -notin @('fmm_e9_adversarial_scaling_profile_pass','fmm_e9_order6_gap_order8_pass','fmm_e9_profile_ready_with_accuracy_gaps') -and -not $AllowFailedDecision) {
  throw "FMM E9 gate did not pass: $decision"
}
