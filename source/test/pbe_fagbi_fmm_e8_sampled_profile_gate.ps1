param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e8_sampled_profile'),
  [string]$Cxx = '',
  [int]$SampleCount = 32,
  [switch]$Full,
  [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$E7Smoke = Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e7_sampled_panel_matvec_smoke.ps1'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E8_SAMPLED_PROFILE_REPORT.txt'

if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
if (-not (Test-Path -LiteralPath $E7Smoke)) {
  throw "missing E7 sampled panel smoke: $E7Smoke"
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) {
    throw "missing PQR fixture: $pqr"
  }
}

if (Test-Path -LiteralPath $WorkRoot) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$BuildRoot = Join-Path $WorkRoot 'build'
$LogPath = Join-Path $WorkRoot 'fmm_e8_sampled_profile_raw.log'
$MatrixPath = Join-Path $WorkRoot 'fmm_e8_matrix.csv'
$TimingPath = Join-Path $WorkRoot 'fmm_e8_timing.csv'
$MemoryPath = Join-Path $WorkRoot 'fmm_e8_memory.csv'
$AdversarialPath = Join-Path $WorkRoot 'fmm_e8_adversarial.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e8_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e8_summary.dat'

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) {
      $fields[$parts[0]] = $parts[1]
    }
  }
  foreach ($required in @('fixture','sample','status','order','aggregate_rms','far_panel_pairs','fallback_panel_pairs')) {
    if (-not $fields.Contains($required)) {
      throw "malformed KIFMM sampled profile row, missing $required`: $Line"
    }
  }
  return [pscustomobject]$fields
}

function As-Double($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
    return [double]::NaN
  }
  return [double]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}

function As-Long($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
    return [int64]0
  }
  return [int64]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}

$runArgs = @(
  '-NoProfile',
  '-ExecutionPolicy', 'Bypass',
  '-File', $E7Smoke,
  '-Pqr1sh1', $Pqr1sh1,
  '-Pqr1a63', $Pqr1a63,
  '-WorkRoot', $BuildRoot,
  '-SampleCount', $SampleCount
)
if (-not [string]::IsNullOrWhiteSpace($Cxx)) {
  $runArgs += @('-Cxx', $Cxx)
}
if ($Full) {
  $runArgs += '-Full'
}

$output = & powershell @runArgs 2>&1
$exitCode = $LASTEXITCODE
$output | Set-Content -LiteralPath $LogPath -Encoding ascii
if ($exitCode -ne 0) {
  if (-not $AllowFailedDecision) {
    throw "E7 sampled profile driver failed; see $LogPath"
  }
}

$rows = @($output | ForEach-Object { [string]$_ } |
  Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } |
  ForEach-Object { Parse-CaseLine $_ })
if ($rows.Count -eq 0) {
  throw "no sampled profile rows parsed from $LogPath"
}

$matrix = @()
foreach ($row in $rows) {
  $maxChannel = [Math]::Max(
    [Math]::Max((As-Double $row.k1), (As-Double $row.k2)),
    [Math]::Max((As-Double $row.k3), (As-Double $row.k4)))
  $matrix += [pscustomobject]@{
    fixture = $row.fixture
    sample_set = $row.sample
    status = $row.status
    source_panels = As-Long $row.source_panels
    sampled_rows = As-Long $row.sampled_rows
    order = As-Long $row.order
    far_cell_pairs = As-Long $row.far_cell_pairs
    far_panel_pairs = As-Long $row.far_panel_pairs
    fallback_panel_pairs = As-Long $row.fallback_panel_pairs
    production_fallback_valid_pairs = As-Long $row.production_fallback_valid_pairs
    production_direct_reference_valid_pairs = As-Long $row.production_direct_reference_valid_pairs
    invalid_pairs = As-Long $row.invalid_pairs
    nonproduction_fallback_pairs = As-Long $row.nonproduction_fallback_pairs
    aggregate_rms = As-Double $row.aggregate_rms
    max_channel_rms = $maxChannel
    k1_rms = As-Double $row.k1
    k2_rms = As-Double $row.k2
    k3_rms = As-Double $row.k3
    k4_rms = As-Double $row.k4
    aggregate_abs_rms = As-Double $row.aggregate_abs_rms
    k1_abs_rms = As-Double $row.k1_abs_rms
    k2_abs_rms = As-Double $row.k2_abs_rms
    k3_abs_rms = As-Double $row.k3_abs_rms
    k4_abs_rms = As-Double $row.k4_abs_rms
    k1_reference_rms = As-Double $row.k1_reference_rms
    k2_reference_rms = As-Double $row.k2_reference_rms
    k3_reference_rms = As-Double $row.k3_reference_rms
    k4_reference_rms = As-Double $row.k4_reference_rms
    max_equiv_residual = As-Double $row.max_equiv_residual
    max_condition = As-Double $row.max_condition
    max_equivalent_rank = As-Long $row.max_equivalent_rank
    min_equivalent_sigma_min = As-Double $row.min_equivalent_sigma_min
    max_equivalent_sigma_max = As-Double $row.max_equivalent_sigma_max
    max_equivalent_rhs_norm = As-Double $row.max_equivalent_rhs_norm
    max_equivalent_charge_norm = As-Double $row.max_equivalent_charge_norm
    max_equivalent_charge_to_rhs_norm = As-Double $row.max_equivalent_charge_to_rhs_norm
    worklist_hash = $row.worklist_hash
    policy_hash = $row.policy_hash
  }
}
$matrix | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation -Encoding ascii

$timing = foreach ($row in $rows) {
  [pscustomobject]@{
    fixture = $row.fixture
    sample_set = $row.sample
    order = As-Long $row.order
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

$memory = foreach ($row in $rows) {
  [pscustomobject]@{
    fixture = $row.fixture
    sample_set = $row.sample
    order = As-Long $row.order
    source_tree_bytes = As-Long $row.source_tree_bytes
    target_tree_bytes = As-Long $row.target_tree_bytes
    worklist_bytes = As-Long $row.worklist_bytes
    equivalent_charge_storage_bytes = As-Long $row.equivalent_charge_storage_bytes
    cached_channel_object_estimated_bytes = As-Long $row.cached_channel_object_estimated_bytes
  }
}
$memory | Export-Csv -LiteralPath $MemoryPath -NoTypeInformation -Encoding ascii

$adversarial = $matrix | Where-Object {
  $_.sample_set -notlike 'uniform_hash*'
}
if (@($adversarial).Count -eq 0) {
  $adversarial = @([pscustomobject]@{
      fixture = 'not_run'
      sample_set = 'quick_mode_uniform_only'
      status = 'informational'
      order = 0
      aggregate_rms = ''
      max_channel_rms = ''
      far_panel_pairs = ''
      fallback_panel_pairs = ''
    })
}
$adversarial | Export-Csv -LiteralPath $AdversarialPath -NoTypeInformation -Encoding ascii

$badCoverage = @($matrix | Where-Object {
    $_.invalid_pairs -ne 0 -or $_.nonproduction_fallback_pairs -ne 0 -or
    $_.production_direct_reference_valid_pairs -ne ($_.source_panels * $_.sampled_rows)
  })
$conditionGap = @($matrix | Where-Object { $_.max_condition -gt 1.0e12 -or $_.max_equiv_residual -gt 1.0e-10 })
$candidateRows = @($matrix | Where-Object { $_.order -ne 4 -and $_.sample_set -notlike '*q4_diagnostic*' })
$farGap = @($candidateRows | Where-Object { $_.status -eq 'far_equivalent_failure' -or $_.aggregate_rms -gt 1.0e-6 -or $_.max_channel_rms -gt 5.0e-6 })
$q6Pass = @($matrix | Where-Object { $_.order -eq 6 -and $_.status -eq 'pass' -and $_.aggregate_rms -le 1.0e-6 -and $_.max_channel_rms -le 5.0e-6 })
$q8Pass = @($matrix | Where-Object { $_.order -eq 8 -and $_.status -eq 'pass' -and $_.aggregate_rms -le 1.0e-6 -and $_.max_channel_rms -le 5.0e-6 })

$decision = $(if ($Full) { 'fmm_e8_profile_ready_with_accuracy_gaps' } else { 'fmm_e8_uniform_quick_profile_pass_full_adversarial_pending' })
if ($badCoverage.Count -gt 0) {
  $decision = 'fmm_e8_fallback_coverage_gap'
} elseif ($conditionGap.Count -gt 0) {
  $decision = 'fmm_e8_condition_gap'
} elseif ($Full -and $q6Pass.Count -gt 0 -and $farGap.Count -eq 0) {
  $decision = 'fmm_e8_sampled_profile_pass'
} elseif ($q8Pass.Count -gt 0) {
  $decision = $(if ($Full) { 'fmm_e8_order6_gap_order8_pass' } else { 'fmm_e8_uniform_quick_profile_pass_full_adversarial_pending' })
} elseif ($farGap.Count -gt 0) {
  $decision = 'fmm_e8_far_equivalent_gap'
}
if ($exitCode -ne 0 -and $decision -eq 'fmm_e8_profile_ready_with_accuracy_gaps') {
  $decision = 'fmm_e8_runtime_limited'
}

[pscustomobject]@{
  decision = $decision
  mode = $(if ($Full) { 'full_adversarial' } else { 'quick_uniform' })
  rows = $matrix.Count
  sample_count = $SampleCount
  source = 'E7 sampled-panel driver'
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$(if ($Full) { 'full_adversarial' } else { 'quick_uniform' })",
  "rows=$($matrix.Count)",
  "matrix=$MatrixPath",
  "timing=$TimingPath",
  "memory=$MemoryPath",
  "adversarial=$AdversarialPath",
  "raw_log=$LogPath"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

@"
TITAN FMM-E8 Sampled Production-Surface KIFMM-lite Profile Gate

Decision: $decision
Mode: $(if ($Full) { 'full_adversarial' } else { 'quick_uniform' })
Sample count: $SampleCount

Artifacts:
- Matrix: $MatrixPath
- Timing: $TimingPath
- Memory: $MemoryPath
- Adversarial: $AdversarialPath
- Decision: $DecisionPath
- Summary: $SummaryPath
- Raw log: $LogPath

Scope:
- Diagnostic-only sampled production-surface KIFMM-lite profiling.
- No production TABI matvec, GMRES, preconditioner, parser-default, or policy change.
- Quick mode records uniform sampled rows only and is intentionally labelled
  full-adversarial-pending.
- Full mode delegates adversarial sampled sets to the E7 driver.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E8 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Timing: $TimingPath"
Write-Host "Memory: $MemoryPath"
if ($decision -notin @('fmm_e8_uniform_quick_profile_pass_full_adversarial_pending','fmm_e8_sampled_profile_pass','fmm_e8_order6_gap_order8_pass','fmm_e8_profile_ready_with_accuracy_gaps') -and -not $AllowFailedDecision) {
  throw "FMM E8 gate did not pass: $decision"
}
