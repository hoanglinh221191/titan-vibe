param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e47_runtime_memory_attribution'),
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

$matrixCsv = Join-Path $WorkRoot 'fmm_e47_matrix.csv'
$phaseCsv = Join-Path $WorkRoot 'fmm_e47_phase_attribution.csv'
$memoryCsv = Join-Path $WorkRoot 'fmm_e47_memory_attribution.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e47_decision.csv'
$summaryDat = Join-Path $WorkRoot 'fmm_e47_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E47_RUNTIME_MEMORY_ATTRIBUTION_REPORT.txt'
$e45WorkRoot = Join-Path (Split-Path -Parent $PSScriptRoot) 'test\tmp\fmm_e45_staged_dual_worklist_gmres'

function Assert-Contains([string]$Path, [string[]]$Needles) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Missing file: $Path" }
  $text = Get-Content -Raw -LiteralPath $Path
  foreach ($needle in $Needles) {
    if (-not $text.Contains($needle)) {
      throw "Missing marker '$needle' in $Path"
    }
  }
}

function Get-Double([object]$Value, [double]$Default = [double]::NaN) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0.0
  if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Get-Int64([object]$Value, [int64]$Default = 0) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0L
  if ([int64]::TryParse([string]$Value, [ref]$parsed)) { return $parsed }
  return $Default
}

function Assert-E47StaticContract {
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.h') @(
    'apply_source_build_seconds',
    'apply_guard_audit_seconds',
    'apply_finalize_seconds',
    'apply_unattributed_seconds',
    'apply_rss_after_cartesian_k4_bytes',
    'apply_peak_working_set_after_finalize_bytes'
  )
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.cpp') @(
    'record_apply_memory_after_cartesian_k4',
    'compute_apply_unattributed_seconds',
    'current_process_memory_sample'
  )
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp') @(
    'apply_source_build_seconds',
    'apply_rss_after_cartesian_k4_bytes',
    'apply_peak_working_set_after_finalize_bytes'
  )
}

function Invoke-E45Quick {
  if ($SkipRuntime) { return }
  $e45Script = Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e45_staged_dual_worklist_gmres_gate.ps1'
  $cmd = @(
    '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $e45Script,
    '-TitanExe', $TitanExe,
    '-Quick', '-Force', '-AllowFailedDecision',
    '-CaseTimeoutSeconds', $CaseTimeoutSeconds,
    '-Workers', $Workers
  )
  & powershell @cmd
  if ($LASTEXITCODE -ne 0) {
    throw "E45 quick command failed with exit code $LASTEXITCODE"
  }
}

function Read-E47Rows {
  $path = Join-Path $e45WorkRoot 'fmm_e45_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    throw "Missing E45 invariant aggregate: $path"
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) { throw 'E45 invariant aggregate is empty.' }
  $required = @(
    'apply_source_build_seconds',
    'apply_hierarchy_seconds',
    'apply_guard_audit_seconds',
    'apply_production_fallback_seconds',
    'apply_cartesian_k4_seconds',
    'apply_dense_oracle_seconds',
    'apply_finalize_seconds',
    'apply_unattributed_seconds',
    'apply_rss_start_bytes',
    'apply_rss_after_cartesian_k4_bytes',
    'apply_peak_working_set_after_hierarchy_bytes',
    'apply_peak_working_set_after_cartesian_k4_bytes',
    'apply_peak_working_set_after_finalize_bytes',
    'apply_peak_working_set_delta_bytes'
  )
  foreach ($name in $required) {
    if (-not ($rows[0].PSObject.Properties.Name -contains $name)) {
      throw "Missing E47 invariant column: $name"
    }
  }
  return $rows
}

$decision = 'fmm_e47_runtime_memory_attribution_gap'
$failure = ''
$staticStatus = 'not_run'
$runtimeStatus = 'not_run'
try {
  Assert-E47StaticContract
  $staticStatus = 'pass'
  Invoke-E45Quick
  $runtimeStatus = if ($SkipRuntime) { 'skipped' } else { 'pass' }
  if ($SkipRuntime) {
    $decision = 'fmm_e47_static_attribution_contract_pass'
    $matrix = [pscustomobject]@{
      decision=$decision
      static_status=$staticStatus
      runtime_status=$runtimeStatus
      rows=0
      mean_apply_seconds=[double]::NaN
      mean_hierarchy_seconds=[double]::NaN
      mean_fallback_seconds=[double]::NaN
      mean_cartesian_k4_seconds=[double]::NaN
      max_peak_working_set_bytes=0
      max_peak_delta_bytes=0
      failure=''
    }
    @() | Export-Csv -NoTypeInformation -LiteralPath $phaseCsv
    @() | Export-Csv -NoTypeInformation -LiteralPath $memoryCsv
    throw '__E47_STATIC_PASS__'
  }
  $rows = @(Read-E47Rows)
  $bad = @($rows | Where-Object {
      (Get-Double $_.apply_seconds) -lt 0.0 -or
      (Get-Double $_.apply_source_build_seconds) -lt 0.0 -or
      (Get-Double $_.apply_hierarchy_seconds) -lt 0.0 -or
      (Get-Double $_.apply_guard_audit_seconds) -lt 0.0 -or
      (Get-Double $_.apply_production_fallback_seconds) -lt 0.0 -or
      (Get-Double $_.apply_cartesian_k4_seconds) -lt 0.0 -or
      (Get-Double $_.apply_finalize_seconds) -lt 0.0 -or
      (Get-Int64 $_.apply_peak_working_set_after_finalize_bytes) -le 0
    })
  if ($bad.Count -ne 0) {
    throw "Non-finite or invalid E47 attribution rows: $($bad.Count)"
  }
  $decision = if ($SkipRuntime) {
    'fmm_e47_static_attribution_contract_pass'
  } else {
    'fmm_e47_runtime_memory_attribution_pass'
  }

  $phaseRows = foreach ($row in $rows) {
    [pscustomobject]@{
      stage=$row.stage
      fixture=$row.fixture
      matvec_call=$row.matvec_call
      apply_seconds=$row.apply_seconds
      source_build_seconds=$row.apply_source_build_seconds
      hierarchy_seconds=$row.apply_hierarchy_seconds
      guard_audit_seconds=$row.apply_guard_audit_seconds
      production_fallback_seconds=$row.apply_production_fallback_seconds
      cartesian_k4_seconds=$row.apply_cartesian_k4_seconds
      dense_oracle_seconds=$row.apply_dense_oracle_seconds
      finalize_seconds=$row.apply_finalize_seconds
      unattributed_seconds=$row.apply_unattributed_seconds
    }
  }
  $memoryRows = foreach ($row in $rows) {
    [pscustomobject]@{
      stage=$row.stage
      fixture=$row.fixture
      matvec_call=$row.matvec_call
      rss_start_bytes=$row.apply_rss_start_bytes
      rss_after_source_build_bytes=$row.apply_rss_after_source_build_bytes
      rss_after_hierarchy_bytes=$row.apply_rss_after_hierarchy_bytes
      rss_after_guard_audit_bytes=$row.apply_rss_after_guard_audit_bytes
      rss_after_fallback_bytes=$row.apply_rss_after_fallback_bytes
      rss_after_cartesian_k4_bytes=$row.apply_rss_after_cartesian_k4_bytes
      rss_after_dense_oracle_bytes=$row.apply_rss_after_dense_oracle_bytes
      rss_after_finalize_bytes=$row.apply_rss_after_finalize_bytes
      peak_start_bytes=$row.apply_peak_working_set_start_bytes
      peak_after_source_build_bytes=$row.apply_peak_working_set_after_source_build_bytes
      peak_after_hierarchy_bytes=$row.apply_peak_working_set_after_hierarchy_bytes
      peak_after_guard_audit_bytes=$row.apply_peak_working_set_after_guard_audit_bytes
      peak_after_fallback_bytes=$row.apply_peak_working_set_after_fallback_bytes
      peak_after_cartesian_k4_bytes=$row.apply_peak_working_set_after_cartesian_k4_bytes
      peak_after_dense_oracle_bytes=$row.apply_peak_working_set_after_dense_oracle_bytes
      peak_after_finalize_bytes=$row.apply_peak_working_set_after_finalize_bytes
      peak_delta_bytes=$row.apply_peak_working_set_delta_bytes
    }
  }
  $phaseRows | Export-Csv -NoTypeInformation -LiteralPath $phaseCsv
  $memoryRows | Export-Csv -NoTypeInformation -LiteralPath $memoryCsv
  $matrix = [pscustomobject]@{
    decision=$decision
    static_status=$staticStatus
    runtime_status=$runtimeStatus
    rows=$rows.Count
    mean_apply_seconds=($rows | ForEach-Object { Get-Double $_.apply_seconds } | Measure-Object -Average).Average
    mean_hierarchy_seconds=($rows | ForEach-Object { Get-Double $_.apply_hierarchy_seconds } | Measure-Object -Average).Average
    mean_fallback_seconds=($rows | ForEach-Object { Get-Double $_.apply_production_fallback_seconds } | Measure-Object -Average).Average
    mean_cartesian_k4_seconds=($rows | ForEach-Object { Get-Double $_.apply_cartesian_k4_seconds } | Measure-Object -Average).Average
    max_peak_working_set_bytes=($rows | ForEach-Object { Get-Int64 $_.apply_peak_working_set_after_finalize_bytes } | Measure-Object -Maximum).Maximum
    max_peak_delta_bytes=($rows | ForEach-Object { Get-Int64 $_.apply_peak_working_set_delta_bytes } | Measure-Object -Maximum).Maximum
    failure=''
  }
} catch {
  if ($_.Exception.Message -eq '__E47_STATIC_PASS__') {
    $failure = ''
  } else {
    $failure = $_.Exception.Message
    $matrix = [pscustomobject]@{
      decision=$decision
      static_status=$staticStatus
      runtime_status=$runtimeStatus
      rows=0
      mean_apply_seconds=[double]::NaN
      mean_hierarchy_seconds=[double]::NaN
      mean_fallback_seconds=[double]::NaN
      mean_cartesian_k4_seconds=[double]::NaN
      max_peak_working_set_bytes=0
      max_peak_delta_bytes=0
      failure=$failure
    }
  }
}

$matrix | Export-Csv -NoTypeInformation -LiteralPath $matrixCsv
@([pscustomobject]@{decision=$decision; static_status=$staticStatus; runtime_status=$runtimeStatus; failure=$failure}) |
  Export-Csv -NoTypeInformation -LiteralPath $decisionCsv
@(
  "decision=$decision",
  "static_status=$staticStatus",
  "runtime_status=$runtimeStatus",
  "matrix=$matrixCsv",
  "phase_attribution=$phaseCsv",
  "memory_attribution=$memoryCsv"
) | Set-Content -Encoding ASCII -LiteralPath $summaryDat
@(
  'TITAN FMM-E47 runtime and memory attribution gate',
  '==================================================',
  '',
  "Decision: $decision",
  "Static status: $staticStatus",
  "Runtime status: $runtimeStatus",
  '',
  'E47 is instrumentation-only. It preserves the E45 channel-hybrid GMRES path and adds per-matvec phase timings plus RSS/peak-working-set snapshots so the next optimization can target the actual bottleneck.',
  "Matrix: $matrixCsv",
  "Phase attribution: $phaseCsv",
  "Memory attribution: $memoryCsv"
) | Set-Content -Encoding ASCII -LiteralPath $reportPath

Write-Host "FMM-E47 decision: $decision"
Write-Host "Matrix: $matrixCsv"
if ($decision -notmatch '_pass$' -and -not $AllowFailedDecision) {
  throw "FMM-E47 gate did not pass: $decision :: $failure"
}
