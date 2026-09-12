param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e48_k123_hierarchy_attribution'),
  [int]$CaseTimeoutSeconds = 1800,
  [int]$Workers = 4,
  [switch]$AllowFailedDecision,
  [switch]$SkipRuntime,
  [ValidateSet('canonical-offset-cache','static-translation-table-cache','target-node-streaming-cache','bounded-lru-cache')]
  [string]$MatrixCachePolicy = 'canonical-offset-cache',
  [double]$MatrixCacheLimitMb = 0.0
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$TitanExe = [IO.Path]::GetFullPath($TitanExe)
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$matrixCsv = Join-Path $WorkRoot 'fmm_e48_matrix.csv'
$stageCsv = Join-Path $WorkRoot 'fmm_e48_k123_stage_runtime.csv'
$memoryCsv = Join-Path $WorkRoot 'fmm_e48_k123_stage_memory.csv'
$retainedCsv = Join-Path $WorkRoot 'fmm_e48_k123_retained_bytes.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e48_decision.csv'
$summaryDat = Join-Path $WorkRoot 'fmm_e48_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E48_K123_HIERARCHY_ATTRIBUTION_REPORT.txt'
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

function Assert-E48StaticContract {
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.h') @(
    'hierarchy_runtime_memory_attribution_status',
    'hierarchy_peak_after_p2m_bytes',
    'hierarchy_peak_after_m2l_bytes',
    'k123_upward_equivalent_estimated_bytes',
    'k123_total_estimated_bytes'
  )
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp') @(
    'record_hierarchy_memory_after_p2m',
    'record_hierarchy_memory_after_m2l',
    'update_k123_retained_memory_estimates',
    'current_process_memory_sample'
  )
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp') @(
    'hierarchy_p2m_seconds',
    'hierarchy_peak_after_m2l_bytes',
    'k123_total_estimated_bytes',
    'hierarchy_runtime_memory_attribution_status'
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
    '-Workers', $Workers,
    '-MatrixCachePolicy', $MatrixCachePolicy,
    '-MatrixCacheLimitMb', $MatrixCacheLimitMb
  )
  & powershell @cmd
  if ($LASTEXITCODE -ne 0) {
    throw "E45 quick command failed with exit code $LASTEXITCODE"
  }
}

function Read-E48Rows {
  $path = Join-Path $e45WorkRoot 'fmm_e45_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    throw "Missing E45 invariant aggregate: $path"
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) { throw 'E45 invariant aggregate is empty.' }
  $required = @(
    'hierarchy_p2m_seconds',
    'hierarchy_m2m_seconds',
    'hierarchy_m2l_seconds',
    'hierarchy_local_solve_seconds',
    'hierarchy_l2l_seconds',
    'hierarchy_l2p_seconds',
    'hierarchy_peak_after_p2m_bytes',
    'hierarchy_peak_after_m2l_bytes',
    'hierarchy_peak_after_finalize_bytes',
    'hierarchy_peak_delta_bytes',
    'k123_upward_equivalent_estimated_bytes',
    'k123_local_rhs_estimated_bytes',
    'k123_equivalent_system_cache_estimated_bytes',
    'k123_m2l_matrix_cache_limit_bytes',
    'k123_m2l_matrix_cache_evictions',
    'k123_m2l_matrix_cache_estimated_bytes',
    'k123_m2l_matrix_cache_peak_estimated_bytes',
    'k123_m2l_matrix_cache_total_built_bytes',
    'k123_total_estimated_bytes',
    'hierarchy_runtime_memory_attribution_status'
  )
  foreach ($name in $required) {
    if (-not ($rows[0].PSObject.Properties.Name -contains $name)) {
      throw "Missing E48 invariant column: $name"
    }
  }
  return $rows
}

$decision = 'fmm_e48_k123_hierarchy_attribution_gap'
$failure = ''
$staticStatus = 'not_run'
$runtimeStatus = 'not_run'
try {
  Assert-E48StaticContract
  $staticStatus = 'pass'
  Invoke-E45Quick
  $runtimeStatus = if ($SkipRuntime) { 'skipped' } else { 'pass' }
  if ($SkipRuntime) {
    $decision = 'fmm_e48_static_attribution_contract_pass'
    $matrix = [pscustomobject]@{
      decision=$decision
      static_status=$staticStatus
      runtime_status=$runtimeStatus
      rows=0
      mean_hierarchy_seconds=[double]::NaN
      mean_k123_m2l_seconds=[double]::NaN
      mean_k123_l2p_seconds=[double]::NaN
      max_hierarchy_peak_bytes=0
      max_k123_total_estimated_bytes=0
      failure=''
    }
    @() | Export-Csv -NoTypeInformation -LiteralPath $stageCsv
    @() | Export-Csv -NoTypeInformation -LiteralPath $memoryCsv
    @() | Export-Csv -NoTypeInformation -LiteralPath $retainedCsv
    throw '__E48_STATIC_PASS__'
  }

  $rows = @(Read-E48Rows)
  $bad = @($rows | Where-Object {
      (Get-Double $_.hierarchy_p2m_seconds) -lt 0.0 -or
      (Get-Double $_.hierarchy_m2l_seconds) -lt 0.0 -or
      (Get-Double $_.hierarchy_l2p_seconds) -lt 0.0 -or
      (Get-Int64 $_.hierarchy_peak_after_finalize_bytes) -le 0 -or
      ([string]$_.hierarchy_runtime_memory_attribution_status) -ne 'pass'
    })
  if ($bad.Count -ne 0) {
    throw "Invalid E48 K123 attribution rows: $($bad.Count)"
  }

  $stageRows = foreach ($row in $rows) {
    [pscustomobject]@{
      stage=$row.stage
      fixture=$row.fixture
      matvec_call=$row.matvec_call
      hierarchy_total_seconds=$row.apply_hierarchy_seconds
      p2m_seconds=$row.hierarchy_p2m_seconds
      p2m_rhs_seconds=$row.hierarchy_p2m_rhs_seconds
      p2m_solve_seconds=$row.hierarchy_p2m_solve_seconds
      m2m_seconds=$row.hierarchy_m2m_seconds
      m2l_seconds=$row.hierarchy_m2l_seconds
      local_solve_seconds=$row.hierarchy_local_solve_seconds
      l2l_seconds=$row.hierarchy_l2l_seconds
      l2p_seconds=$row.hierarchy_l2p_seconds
      dense_oracle_seconds=$row.hierarchy_dense_oracle_seconds
    }
  }
  $memoryRows = foreach ($row in $rows) {
    [pscustomobject]@{
      stage=$row.stage
      fixture=$row.fixture
      matvec_call=$row.matvec_call
      rss_start_bytes=$row.hierarchy_rss_start_bytes
      rss_after_p2m_bytes=$row.hierarchy_rss_after_p2m_bytes
      rss_after_m2m_bytes=$row.hierarchy_rss_after_m2m_bytes
      rss_after_m2l_bytes=$row.hierarchy_rss_after_m2l_bytes
      rss_after_l2l_bytes=$row.hierarchy_rss_after_l2l_bytes
      rss_after_l2p_bytes=$row.hierarchy_rss_after_l2p_bytes
      rss_after_finalize_bytes=$row.hierarchy_rss_after_finalize_bytes
      peak_start_bytes=$row.hierarchy_peak_start_bytes
      peak_after_p2m_bytes=$row.hierarchy_peak_after_p2m_bytes
      peak_after_m2m_bytes=$row.hierarchy_peak_after_m2m_bytes
      peak_after_m2l_bytes=$row.hierarchy_peak_after_m2l_bytes
      peak_after_l2l_bytes=$row.hierarchy_peak_after_l2l_bytes
      peak_after_l2p_bytes=$row.hierarchy_peak_after_l2p_bytes
      peak_after_finalize_bytes=$row.hierarchy_peak_after_finalize_bytes
      peak_delta_bytes=$row.hierarchy_peak_delta_bytes
    }
  }
  $retainedRows = foreach ($row in $rows) {
    [pscustomobject]@{
      stage=$row.stage
      fixture=$row.fixture
      matvec_call=$row.matvec_call
      upward_equivalent_bytes=$row.k123_upward_equivalent_estimated_bytes
      local_rhs_bytes=$row.k123_local_rhs_estimated_bytes
      dense_local_rhs_bytes=$row.k123_dense_local_rhs_estimated_bytes
      local_equivalent_bytes=$row.k123_local_equivalent_estimated_bytes
      row_storage_bytes=$row.k123_rows_estimated_bytes
      equivalent_system_cache_bytes=$row.k123_equivalent_system_cache_estimated_bytes
      m2l_matrix_cache_limit_bytes=$row.k123_m2l_matrix_cache_limit_bytes
      m2l_matrix_cache_evictions=$row.k123_m2l_matrix_cache_evictions
      m2l_matrix_cache_bytes=$row.k123_m2l_matrix_cache_estimated_bytes
      m2l_matrix_cache_peak_bytes=$row.k123_m2l_matrix_cache_peak_estimated_bytes
      m2l_matrix_cache_total_built_bytes=$row.k123_m2l_matrix_cache_total_built_bytes
      total_estimated_bytes=$row.k123_total_estimated_bytes
    }
  }
  $stageRows | Export-Csv -NoTypeInformation -LiteralPath $stageCsv
  $memoryRows | Export-Csv -NoTypeInformation -LiteralPath $memoryCsv
  $retainedRows | Export-Csv -NoTypeInformation -LiteralPath $retainedCsv

  $decision = 'fmm_e48_k123_hierarchy_attribution_pass'
  $matrix = [pscustomobject]@{
    decision=$decision
    static_status=$staticStatus
    runtime_status=$runtimeStatus
    rows=$rows.Count
    mean_hierarchy_seconds=($rows | ForEach-Object { Get-Double $_.apply_hierarchy_seconds } | Measure-Object -Average).Average
    mean_k123_m2l_seconds=($rows | ForEach-Object { Get-Double $_.hierarchy_m2l_seconds } | Measure-Object -Average).Average
    mean_k123_l2p_seconds=($rows | ForEach-Object { Get-Double $_.hierarchy_l2p_seconds } | Measure-Object -Average).Average
    max_hierarchy_peak_bytes=($rows | ForEach-Object { Get-Int64 $_.hierarchy_peak_after_finalize_bytes } | Measure-Object -Maximum).Maximum
    max_k123_total_estimated_bytes=($rows | ForEach-Object { Get-Int64 $_.k123_total_estimated_bytes } | Measure-Object -Maximum).Maximum
    failure=''
  }
} catch {
  if ($_.Exception.Message -eq '__E48_STATIC_PASS__') {
    $failure = ''
  } else {
    $failure = $_.Exception.Message
    $matrix = [pscustomobject]@{
      decision=$decision
      static_status=$staticStatus
      runtime_status=$runtimeStatus
      rows=0
      mean_hierarchy_seconds=[double]::NaN
      mean_k123_m2l_seconds=[double]::NaN
      mean_k123_l2p_seconds=[double]::NaN
      max_hierarchy_peak_bytes=0
      max_k123_total_estimated_bytes=0
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
  "stage_runtime=$stageCsv",
  "stage_memory=$memoryCsv",
  "retained_bytes=$retainedCsv"
) | Set-Content -LiteralPath $summaryDat -Encoding ascii

@(
  'TITAN FMM-E48 K123 hierarchy internal runtime/memory attribution gate',
  '======================================================================',
  '',
  "Decision: $decision",
  "Static status: $staticStatus",
  "Runtime status: $runtimeStatus",
  '',
  'E48 is instrumentation-only. It preserves the E45/E47 channel-hybrid GMRES path and attributes the K123 hierarchy phase across P2M, M2M, M2L, local solve/L2L, L2P, retained bytes, and process peak-working-set snapshots.',
  "Matrix: $matrixCsv",
  "Stage runtime: $stageCsv",
  "Stage memory: $memoryCsv",
  "Retained bytes: $retainedCsv"
) | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "FMM-E48 decision: $decision"
Write-Host "Matrix: $matrixCsv"
if (-not $AllowFailedDecision -and $decision -notin @(
    'fmm_e48_static_attribution_contract_pass',
    'fmm_e48_k123_hierarchy_attribution_pass')) {
  exit 1
}
