param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e46_refactor_regression'),
  [int]$CaseTimeoutSeconds = 1800,
  [int]$Workers = 4,
  [switch]$AllowFailedDecision,
  [switch]$SkipRuntime
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$TitanExe = [IO.Path]::GetFullPath($TitanExe)
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$matrixCsv = Join-Path $WorkRoot 'fmm_e46_matrix.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e46_decision.csv'
$summaryDat = Join-Path $WorkRoot 'fmm_e46_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E46_REFACTOR_REPORT.txt'
$e45WorkRoot = Join-Path (Split-Path -Parent $PSScriptRoot) 'test\tmp\fmm_e45_staged_dual_worklist_gmres'

function Test-E46FileContains([string]$Path, [string[]]$Needles) {
  if (-not (Test-Path -LiteralPath $Path)) {
    throw "Missing required E46 file: $Path"
  }
  $text = Get-Content -Raw -LiteralPath $Path
  foreach ($needle in $Needles) {
    if (-not $text.Contains($needle)) {
      throw "Missing marker '$needle' in $Path"
    }
  }
}

function Assert-E46StaticContract {
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_real_panel_geometry.h') @(
    'KifmmRealPanelNodeEnvelope',
    'build_source_envelopes',
    'build_target_envelopes',
    'panel_radii',
    'valid_panel'
  )
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_real_panel_geometry.cpp') @(
    'KifmmRealPanelNodeEnvelope',
    'build_source_envelopes',
    'build_target_envelopes',
    'hash_panels'
  )
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_real_panel_worklist.h') @(
    'build_real_panel_guarded_worklist',
    'update_leaf_panel_telemetry',
    'update_l2l_containment_telemetry'
  )
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_real_panel_worklist.cpp') @(
    'strict-local-check',
    'build_real_panel_guarded_worklist',
    'negative_target_local_check_excess_pairs'
  )
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_channel_hybrid_operator.h') @(
    'DualWorklistChannelHybridPreparedOperator',
    'prepare_dual_worklist_channel_hybrid_operator',
    'apply_prepared_dual_worklist_channel_hybrid_operator',
    'assemble_dual_worklist_channel_hybrid_tabi_system_matvec'
  )
  Test-E46FileContains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_channel_hybrid_operator.cpp') @(
    'frozen_prepared_apply',
    'run_dual_worklist_channel_hybrid_with_geometry',
    'channel_hybrid_double_count_guard'
  )
  $solver = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp')
  if (-not $solver.Contains('tabi_fmm_kifmm_channel_hybrid_operator.h')) {
    throw 'tabi_solver.cpp must include the E46 channel-hybrid operator header.'
  }
}

function Get-E46Double([object]$Value, [double]$Default = [double]::NaN) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0.0
  if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Assert-E45Regression {
  $e45Script = Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e45_staged_dual_worklist_gmres_gate.ps1'
  $cmd = @(
    '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $e45Script,
    '-TitanExe', $TitanExe,
    '-Quick', '-Force', '-AllowFailedDecision',
    '-CaseTimeoutSeconds', $CaseTimeoutSeconds,
    '-Workers', $Workers
  )
  $childOutput = & powershell @cmd 2>&1
  foreach ($line in $childOutput) {
    Write-Host $line
  }
  if ($LASTEXITCODE -ne 0) {
    throw "E45 quick regression command failed with exit code $LASTEXITCODE"
  }
  $decisionPath = Join-Path $e45WorkRoot 'fmm_e45_decision.csv'
  $matrixPath = Join-Path $e45WorkRoot 'fmm_e45_matrix.csv'
  $invPath = Join-Path $e45WorkRoot 'fmm_e45_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $decisionPath)) { throw "Missing E45 decision: $decisionPath" }
  if (-not (Test-Path -LiteralPath $matrixPath)) { throw "Missing E45 matrix: $matrixPath" }
  if (-not (Test-Path -LiteralPath $invPath)) { throw "Missing E45 invariants: $invPath" }
  $decision = @(Import-Csv -LiteralPath $decisionPath)[0]
  if ($decision.decision -ne 'fmm_e45_short_gmres_pass_restart_pending') {
    throw "Unexpected E45 decision after refactor: $($decision.decision)"
  }
  $rows = @(Import-Csv -LiteralPath $matrixPath)
  foreach ($row in $rows) {
    if ($row.operator_status -ne 'pass' -or
        $row.prepared_freeze_status -ne 'pass' -or
        $row.frozen_worklist -ne 'true' -or
        $row.frozen_equivalent_system_cache -ne 'true' -or
        $row.prepared_rebuild_status -ne 'frozen_prepared_apply' -or
        [int64]$row.invalid_pairs -ne 0 -or
        [int64]$row.nonproduction_fallback_pairs -ne 0) {
      throw "E45 lifecycle regression on $($row.stage)/$($row.fixture)"
    }
  }
  $inv = @(Import-Csv -LiteralPath $invPath)
  $first = $inv[0]
  $drift = @($inv | Where-Object {
      $_.worklist_hash -ne $first.worklist_hash -or
      $_.policy_hash -ne $first.policy_hash -or
      $_.sampled_target_hash -ne $first.sampled_target_hash -or
      $_.cache_entries -ne $first.cache_entries
    })
  if ($drift.Count -ne 0) {
    throw 'E45 invariant hash/cache drift after E46 refactor.'
  }
  return [pscustomobject]@{
    decision=$decision.decision
    rows=$rows.Count
    invariant_rows=$inv.Count
    mean_apply_seconds=($rows | Measure-Object -Property mean_apply_seconds -Average).Average
    max_sampled_audit=($rows | ForEach-Object { Get-E46Double $_.sampled_audit_rms } | Measure-Object -Maximum).Maximum
  }
}

$staticStatus = 'not_run'
$runtimeStatus = 'not_run'
$runtimeSummary = $null
try {
  Assert-E46StaticContract
  $staticStatus = 'pass'
  if ($SkipRuntime) {
    $runtimeStatus = 'skipped'
  } else {
    $runtimeSummary = Assert-E45Regression
    $runtimeStatus = 'pass'
  }
  $decision = if ($runtimeStatus -eq 'pass' -or $runtimeStatus -eq 'skipped') {
    'fmm_e46_refactor_static_pass'
  } else {
    'fmm_e46_refactor_regression_gap'
  }
} catch {
  $decision = 'fmm_e46_refactor_regression_gap'
  $failure = $_.Exception.Message
}

$matrix = [pscustomobject]@{
  static_status=$staticStatus
  runtime_status=$runtimeStatus
  e45_decision=$(if ($runtimeSummary) { $runtimeSummary.decision } else { '' })
  e45_rows=$(if ($runtimeSummary) { $runtimeSummary.rows } else { 0 })
  e45_invariant_rows=$(if ($runtimeSummary) { $runtimeSummary.invariant_rows } else { 0 })
  e45_mean_apply_seconds=$(if ($runtimeSummary) { $runtimeSummary.mean_apply_seconds } else { [double]::NaN })
  e45_max_sampled_audit=$(if ($runtimeSummary) { $runtimeSummary.max_sampled_audit } else { [double]::NaN })
  failure=$(if ($decision -eq 'fmm_e46_refactor_regression_gap') { $failure } else { '' })
}
$matrix | Export-Csv -NoTypeInformation -LiteralPath $matrixCsv
@([pscustomobject]@{decision=$decision; static_status=$staticStatus; runtime_status=$runtimeStatus}) |
  Export-Csv -NoTypeInformation -LiteralPath $decisionCsv
@(
  "decision=$decision",
  "static_status=$staticStatus",
  "runtime_status=$runtimeStatus",
  "matrix=$matrixCsv"
) | Set-Content -Encoding ASCII -LiteralPath $summaryDat
@(
  'TITAN FMM-E46 channel-hybrid prepared-operator refactor gate',
  '================================================================',
  '',
  "Decision: $decision",
  "Static status: $staticStatus",
  "Runtime status: $runtimeStatus",
  '',
  'E46 is a surgical zero-behavior-change refactor. It only passes when the new geometry/worklist/channel-hybrid modules exist and E45 quick lifecycle invariants remain frozen and stable.',
  "Matrix: $matrixCsv"
) | Set-Content -Encoding ASCII -LiteralPath $reportPath

Write-Host "FMM-E46 decision: $decision"
Write-Host "Matrix: $matrixCsv"
if ($decision -ne 'fmm_e46_refactor_static_pass' -and -not $AllowFailedDecision) {
  throw "FMM-E46 gate did not pass: $decision"
}
