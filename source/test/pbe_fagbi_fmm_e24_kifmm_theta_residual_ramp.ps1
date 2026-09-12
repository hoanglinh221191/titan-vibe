param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e24_theta_ramp'),
  [int]$CaseTimeoutSeconds = 0,
  [switch]$Quick,
  [switch]$Full,
  [switch]$Force,
  [switch]$RerunFailed,
  [string[]]$RerunStages = @(),
  [switch]$AllowFailedDecision,
  [int]$Workers = 4
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (-not $Quick -and -not $Full) { $Quick = $true }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$matrixCsv = Join-Path $WorkRoot 'fmm_e24_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fmm_e24_convergence_trace.csv'
$invariantCsv = Join-Path $WorkRoot 'fmm_e24_invariants.csv'
$thetaSummaryCsv = Join-Path $WorkRoot 'fmm_e24_theta_summary.csv'
$iluGraphCsv = Join-Path $WorkRoot 'fmm_e24_ilu_graph_summary.csv'
$runtimeCsv = Join-Path $WorkRoot 'fmm_e24_runtime.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e24_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e24_summary.dat'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E24_THETA_RAMP_REPORT.txt'

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }

$RuntimeRepoRoot = $RepoRoot
$ffProbe = 'lib_ff\ffAMBER99SBnonbonded.dat'
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot $ffProbe))) {
  $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
  if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
    $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
  }
  $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot $ffProbe))) {
  $siblingMain = 'C:\tmp\elec_vdw-main-sync'
  if (Test-Path -LiteralPath (Join-Path $siblingMain $ffProbe)) {
    $RuntimeRepoRoot = $siblingMain
  }
}
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot $ffProbe))) {
  throw "Missing AMBER99SB force-field file under runtime repo root: $RuntimeRepoRoot"
}

function Get-E24Double([object]$Value, [double]$Default = [double]::NaN) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0.0
  if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Get-E24Int64([object]$Value, [int64]$Default = 0) {
  if ($null -eq $Value) { return $Default }
  $parsed = [int64]0
  if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Test-E24Finite([double]$Value) {
  return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Get-E24ProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
  if ($null -eq $Profile) { return $Default }
  $prop = $Profile.PSObject.Properties[$Name]
  if ($null -eq $prop) { return $Default }
  return $prop.Value
}

function Assert-E24StaticContract {
  $args = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\arguments.cpp')
  foreach ($needle in @('kifmm_theta_allowed','0.75','0.8','pbe-bie-kifmm-ilu-graph')) {
    if (-not $args.Contains($needle)) { throw "Missing E24 arguments marker: $needle" }
  }
  $parse = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\cli\parse_arguments.cpp')
  foreach ($needle in @('pbe-bie-kifmm-ilu-graph','kifmm-fallback-capped')) {
    if (-not $parse.Contains($needle)) { throw "Missing E24 parser marker: $needle" }
  }
  $solver = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp')
  foreach ($needle in @(
      'pbe_tabi_kifmm_gmres_convergence_trace.csv',
      'reported_true_residual',
      'build_kifmm_lite_fallback_panel_graph',
      'fmm_nearfield_ilu0_kifmm_fallback_capped')) {
    if (-not $solver.Contains($needle)) { throw "Missing E24 solver marker: $needle" }
  }
}

function Get-E24Fixture([string]$Name) {
  $manifest = @(Import-Csv -LiteralPath $manifestPath)
  $row = @($manifest | Where-Object { $_.fixture -eq $Name })
  if ($row.Count -ne 1) { throw "Expected exactly one fixture '$Name', found $($row.Count)." }
  $pqr = Join-Path $FixtureRoot $row[0].pqr
  if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing PQR for $Name`: $pqr" }
  $row[0] | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqr -Force
  return $row[0]
}

function New-E24Stage([string]$Stage, [string]$Fixture, [double]$Theta,
    [string]$IluGraph, [int]$MaxIter, [int]$Restart, [int]$SampleCount,
    [int]$Timeout) {
  return [pscustomobject]@{
    stage=$Stage; fixture=$Fixture; theta=$Theta; ilu_graph=$IluGraph
    max_iter=$MaxIter; restart=$Restart; sample_count=$SampleCount
    timeout=$Timeout
  }
}

function Write-E24ArgFile([object]$Fixture, [object]$Stage, [string]$RunDir) {
  $arg = Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
    -PqrPath $Fixture.pqr_path -Backend tabi -Mode treecode `
    -FmmMode single-matvec-analytic -SurfacePoints 8 `
    -MaxIter ([int]$Stage.max_iter) -Restart ([int]$Stage.restart) `
    -Tolerance 1e-6 -Threads $Workers -SaltM 0.15 `
    -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
    -OpeningTheta ([double]$Stage.theta) -AuditSampleCount ([int]$Stage.sample_count) `
    -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
    -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
    -NearIntegration analytic_singular -TaylorOrder 4 `
    -FmmCartesianAudit off -FmmLocal off -FmmLocalRefine off `
    -FmmLocalRefineTolerance 0.002 -FmmLocalCache off `
    -FmmLocalCacheBudgetMb 256 -FmmLocalFallback scalar `
    -FmmLocalPlanBuilder legacy-rescan -FmmLocalRefineProbeCount 4 `
    -FmmLocalRefineProbeMode rademacher -FmmLocalGeometryGuard off `
    -FmmLocalGeometryTheta 0.0

  $lines = @(Get-Content -LiteralPath $arg)
  $libPath = Convert-FmmCliToTitanPath(
    (Resolve-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
  $lines = $lines | ForEach-Object {
    if ($_ -like '-path-a=*') {
      "-path-a=$libPath"
    } elseif ($_ -like '-path-n=*') {
      "-path-n=$libPath"
    } elseif ($_ -eq '-pbe-tabi-tree-leaf-size=64') {
      '-pbe-tabi-tree-leaf-size=32'
    } elseif ($_ -like '-pbe-tabi-true-residual-sample-count=*') {
      "-pbe-tabi-true-residual-sample-count=$($Stage.sample_count)"
    } else {
      $_
    }
  }
  $lines += @(
    '-pbe-tabi-fmm-engine=kifmm-lite-experimental',
    '-pbe-tabi-kifmm-order=6',
    '-pbe-tabi-kifmm-equiv-scale=1.05',
    '-pbe-tabi-kifmm-check-scale=2.40',
    '-pbe-tabi-kifmm-solver=qr',
    '-pbe-tabi-kifmm-profile=summary',
    '-pbe-tabi-kifmm-full-apply=target-node-parallel',
    "-pbe-tabi-kifmm-ilu-graph=$($Stage.ilu_graph)"
  )
  $bad = @($lines | Where-Object { -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*=' })
  if ($bad.Count -ne 0) {
    throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
  }
  Set-Content -LiteralPath $arg -Value $lines -Encoding ASCII
  return $arg
}

function Read-E24InvariantSummary([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_kifmm_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0; nonproduction_fallback_pairs=0; first_apply_seconds=[double]::NaN; mean_apply_seconds=[double]::NaN; far_pairs=0; fallback_pairs=0}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0; nonproduction_fallback_pairs=0; first_apply_seconds=[double]::NaN; mean_apply_seconds=[double]::NaN; far_pairs=0; fallback_pairs=0}
  }
  $first = $rows[0]
  $hashStable = @($rows | Where-Object {
      $_.worklist_hash -ne $first.worklist_hash -or
      $_.policy_hash -ne $first.policy_hash -or
      $_.panel_geometry_hash -ne $first.panel_geometry_hash -or
      $_.sampled_target_hash -ne $first.sampled_target_hash
    }).Count -eq 0
  $cacheStable = @($rows | Where-Object { $_.cache_entries -ne $first.cache_entries }).Count -eq 0
  $invalid = 0L; $nonprod = 0L; $applyTotal = 0.0
  foreach ($row in $rows) {
    $invalid += Get-E24Int64 $row.sampled_invalid_pairs
    $nonprod += Get-E24Int64 $row.nonproduction_fallback_pairs
    $applyTotal += Get-E24Double $row.apply_seconds 0.0
  }
  return [pscustomobject]@{
    rows=$rows.Count
    hash_stable=[string]$hashStable.ToString().ToLowerInvariant()
    cache_stable=[string]$cacheStable.ToString().ToLowerInvariant()
    invalid_pairs=$invalid
    nonproduction_fallback_pairs=$nonprod
    first_apply_seconds=Get-E24Double $first.apply_seconds
    mean_apply_seconds=$applyTotal / [double]$rows.Count
    far_pairs=Get-E24Int64 $first.far_panel_pairs
    fallback_pairs=Get-E24Int64 $first.fallback_panel_pairs
  }
}

function Read-E24Progress([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_kifmm_gmres_progress.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'}
  }
  return [pscustomobject]@{rows=$rows.Count; last_phase=[string]$rows[$rows.Count - 1].phase}
}

function Read-E24TraceSummary([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_kifmm_gmres_convergence_trace.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; last_preconditioned_residual=[double]::NaN; last_reported_true_residual=[double]::NaN}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; last_preconditioned_residual=[double]::NaN; last_reported_true_residual=[double]::NaN}
  }
  $last = $rows[$rows.Count - 1]
  return [pscustomobject]@{
    rows=$rows.Count
    last_preconditioned_residual=Get-E24Double $last.preconditioned_residual
    last_reported_true_residual=Get-E24Double $last.reported_true_residual
  }
}

function Get-E24OperatorStatus([object]$Inv) {
  if ([int64]$Inv.rows -gt 0 -and
      $Inv.hash_stable -eq 'true' -and
      $Inv.cache_stable -eq 'true' -and
      [int64]$Inv.invalid_pairs -eq 0 -and
      [int64]$Inv.nonproduction_fallback_pairs -eq 0 -and
      (Test-E24Finite ([double]$Inv.first_apply_seconds)) -and
      (Test-E24Finite ([double]$Inv.mean_apply_seconds))) {
    return 'pass'
  }
  return 'fail'
}

function Get-E24GmresStatus([int]$ExitCode, [string]$ProfileStatus, [string]$ConvergenceReason) {
  if ($ExitCode -ne 0) { return "exit_$ExitCode" }
  if ($ProfileStatus -eq 'ok') { return 'converged' }
  if ($ConvergenceReason -eq 'max_iter') { return 'max_iter_nonconverged' }
  if ([string]::IsNullOrWhiteSpace($ConvergenceReason)) { return 'not_observed' }
  return $ConvergenceReason
}

function Get-E24AuditStatus([string]$GmresStatus, [double]$SampledAuditRms, [double]$SampledDirectTrueResidual) {
  if ($GmresStatus -ne 'converged') { return 'not_interpreted_nonconverged' }
  if ((Test-E24Finite $SampledAuditRms) -and
      (Test-E24Finite $SampledDirectTrueResidual) -and
      $SampledAuditRms -le 8.0e-3 -and
      $SampledDirectTrueResidual -le 8.0e-3) {
    return 'pass'
  }
  return 'fail'
}

function New-E24FailureRow([object]$Fixture, [object]$Stage, [string]$Status, [object]$Run, [string]$RunDir, [string]$FailureLog = '') {
  $inv = Read-E24InvariantSummary $RunDir
  $progress = Read-E24Progress $RunDir
  $trace = Read-E24TraceSummary $RunDir
  $executionStatus = if ($Status -eq 'timeout') { 'timeout' } else { 'failed' }
  $gmresStatus = if ($Status -eq 'timeout') { 'not_observed_due_timeout' } else { 'not_observed_due_failure' }
  return [pscustomobject]@{
    stage=$Stage.stage; fixture=$Fixture.fixture; theta=$Stage.theta; ilu_graph=$Stage.ilu_graph
    status=$Status; completed='false'; execution_status=$executionStatus; gmres_status=$gmresStatus
    operator_status='not_evaluated'; audit_status='not_evaluated'; audit_interpreted='false'
    timeout_phase=$progress.last_phase; progress_rows=$progress.rows; trace_rows=$trace.rows
    max_iter=$Stage.max_iter; restart=$Stage.restart; sample_count=$Stage.sample_count
    iterations=0; convergence_reason=''; final_residual=[double]::NaN
    internal_residual=[double]::NaN; last_preconditioned_residual=$trace.last_preconditioned_residual
    last_reported_true_residual=$trace.last_reported_true_residual
    sampled_audit_rms=[double]::NaN; sampled_direct_true_residual=[double]::NaN
    preconditioner_effective=''; preconditioner_graph_source=$(if ($Stage.ilu_graph -eq 'kifmm-fallback-capped') { 'kifmm_fallback_capped_panel_graph' } else { 'cartesian_fmm_m2l_near_panel_graph' })
    preconditioner_nnz_blocks=0; preconditioner_avg_degree=0; preconditioner_max_degree=0
    plan_build_count=0; cache_build_count=0; kifmm_apply_calls=$inv.rows
    invariant_rows=$inv.rows; invariant_hash_stable=$inv.hash_stable; invariant_cache_stable=$inv.cache_stable
    invariant_invalid_pairs=$inv.invalid_pairs; invariant_nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
    far_pairs=$inv.far_pairs; fallback_pairs=$inv.fallback_pairs
    first_apply_seconds=$inv.first_apply_seconds; mean_apply_seconds=$inv.mean_apply_seconds
    wall_seconds=$Run.wall_seconds; peak_rss_mb=$Run.peak_working_set_mb
    failure_log=$FailureLog; run_dir=$RunDir
  }
}

function Invoke-E24Stage([object]$Stage) {
  $fixture = Get-E24Fixture $Stage.fixture
  $runName = "$($Stage.stage)_$($fixture.fixture)"
  $runDir = Join-Path $WorkRoot $runName
  $done = Join-Path $runDir 'fmm_e24_row.csv'
  if ((Test-Path -LiteralPath $done) -and -not $Force) {
    $existing = @(Import-Csv -LiteralPath $done)[0]
    $failed = $existing.status -eq 'timeout' -or $existing.status -eq 'failed' -or $existing.status -like 'exit_*'
    $stageRequested = $RerunStages -contains $Stage.stage
    if (-not (($RerunFailed -and $failed) -or $stageRequested)) {
      Write-Host "resume: $runName"
      return $existing
    }
  }
  if (Test-Path -LiteralPath $runDir) {
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
  } else {
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
  }
  $arg = Write-E24ArgFile -Fixture $fixture -Stage $Stage -RunDir $runDir
  $timeout = if ($CaseTimeoutSeconds -gt 0) { $CaseTimeoutSeconds } else { [int]$Stage.timeout }
  Write-Host "run: $runName (timeout ${timeout}s)"
  $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg -CaseTimeoutSeconds $timeout
  if ($run.timed_out) {
    $row = New-E24FailureRow $fixture $Stage 'timeout' $run $runDir
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
  }
  $energyPath = Join-Path $runDir 'pbe_energy.csv'
  $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
  if (-not (Test-Path -LiteralPath $profilePath)) {
    $stdout = Join-Path $runDir 'stdout.log'
    $stderr = Join-Path $runDir 'stderr.log'
    $log = ''
    if (Test-Path -LiteralPath $stdout) { $log += Get-Content -Raw -LiteralPath $stdout }
    if (Test-Path -LiteralPath $stderr) { $log += Get-Content -Raw -LiteralPath $stderr }
    $row = New-E24FailureRow $fixture $Stage 'failed' $run $runDir $log
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
  }
  $profile = @(Import-Csv -LiteralPath $profilePath)[0]
  $energy = if (Test-Path -LiteralPath $energyPath) { @(Import-Csv -LiteralPath $energyPath)[0] } else { $null }
  $inv = Read-E24InvariantSummary $runDir
  $progress = Read-E24Progress $runDir
  $trace = Read-E24TraceSummary $runDir
  $profileStatus = [string](Get-E24ProfileValue $profile 'status' '')
  $energyStatus = if ($null -ne $energy) { [string]$energy.status } else { 'missing' }
  $convergenceReason = [string](Get-E24ProfileValue $profile 'convergence_reason' '')
  $gmresStatus = Get-E24GmresStatus $run.exit_code $profileStatus $convergenceReason
  $operatorStatus = Get-E24OperatorStatus $inv
  $sampledAuditRms = Get-E24Double (Get-E24ProfileValue $profile 'tree_vs_direct_matvec_rel' ([double]::NaN))
  $sampledDirectTrueResidual = Get-E24Double (Get-E24ProfileValue $profile 'fmm_local_sampled_direct_true_residual' ([double]::NaN))
  $auditStatus = Get-E24AuditStatus $gmresStatus $sampledAuditRms $sampledDirectTrueResidual
  $rowStatus = if ($run.exit_code -eq 0 -and $energyStatus -eq 'ok' -and $profileStatus -eq 'ok') {
    'ok'
  } elseif ($run.exit_code -eq 0 -and $gmresStatus -eq 'max_iter_nonconverged') {
    'completed_nonconverged'
  } elseif ($run.exit_code -eq 0) {
    "completed_nonconverged:energy=$energyStatus;profile=$profileStatus"
  } else {
    "exit_$($run.exit_code):energy=$energyStatus;profile=$profileStatus"
  }
  $precondEffective = [string](Get-E24ProfileValue $profile 'preconditioner_effective' '')
  $graphSource = if ($precondEffective -like '*kifmm_fallback*' -or $Stage.ilu_graph -eq 'kifmm-fallback-capped') {
    'kifmm_fallback_capped_panel_graph'
  } else {
    'cartesian_fmm_m2l_near_panel_graph'
  }
  $row = [pscustomobject]@{
    stage=$Stage.stage; fixture=$fixture.fixture; theta=$Stage.theta; ilu_graph=$Stage.ilu_graph
    status=$rowStatus; completed=$(if ($run.exit_code -eq 0 -or $profileStatus -ne '') { 'true' } else { 'false' })
    execution_status=$(if ($run.exit_code -eq 0) { 'completed' } else { "exit_$($run.exit_code)" })
    gmres_status=$gmresStatus; operator_status=$operatorStatus; audit_status=$auditStatus
    audit_interpreted=$(if ($gmresStatus -eq 'converged') { 'true' } else { 'false' })
    timeout_phase=$progress.last_phase; progress_rows=$progress.rows; trace_rows=$trace.rows
    max_iter=$Stage.max_iter; restart=$Stage.restart; sample_count=$Stage.sample_count
    iterations=Get-E24Int64 (Get-E24ProfileValue $profile 'iterations' 0)
    convergence_reason=$convergenceReason
    final_residual=Get-E24Double (Get-E24ProfileValue $profile 'final_residual' ([double]::NaN))
    internal_residual=Get-E24Double (Get-E24ProfileValue $profile 'true_residual_rel' ([double]::NaN))
    last_preconditioned_residual=$trace.last_preconditioned_residual
    last_reported_true_residual=$trace.last_reported_true_residual
    sampled_audit_rms=$sampledAuditRms
    sampled_direct_true_residual=$sampledDirectTrueResidual
    preconditioner_effective=$precondEffective
    preconditioner_graph_source=$graphSource
    preconditioner_nnz_blocks=Get-E24Int64 (Get-E24ProfileValue $profile 'preconditioner_sparse_nnz_blocks' 0)
    preconditioner_avg_degree=Get-E24Double (Get-E24ProfileValue $profile 'preconditioner_sparse_avg_degree' 0.0)
    preconditioner_max_degree=Get-E24Int64 (Get-E24ProfileValue $profile 'preconditioner_sparse_max_degree' 0)
    plan_build_count=Get-E24Int64 (Get-E24ProfileValue $profile 'fmm_local_plan_build_count' 0)
    cache_build_count=Get-E24Int64 (Get-E24ProfileValue $profile 'fmm_local_cache_build_count' 0)
    kifmm_apply_calls=Get-E24Int64 (Get-E24ProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    invariant_rows=$inv.rows; invariant_hash_stable=$inv.hash_stable; invariant_cache_stable=$inv.cache_stable
    invariant_invalid_pairs=$inv.invalid_pairs; invariant_nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
    far_pairs=$inv.far_pairs; fallback_pairs=$inv.fallback_pairs
    first_apply_seconds=$inv.first_apply_seconds; mean_apply_seconds=$inv.mean_apply_seconds
    wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
    failure_log=''; run_dir=$runDir
  }
  $row | Export-Csv -LiteralPath $done -NoTypeInformation
  return $row
}

function Test-E24LifecyclePass([object]$Row) {
  return $Row.completed -eq 'true' -and
    $Row.execution_status -eq 'completed' -and
    $Row.operator_status -eq 'pass' -and
    [int64]$Row.plan_build_count -eq 1 -and
    [int64]$Row.cache_build_count -eq 1 -and
    [int64]$Row.invariant_rows -gt 0 -and
    $Row.invariant_hash_stable -eq 'true' -and
    $Row.invariant_cache_stable -eq 'true' -and
    [int64]$Row.invariant_invalid_pairs -eq 0 -and
    [int64]$Row.invariant_nonproduction_fallback_pairs -eq 0 -and
    (Test-E24Finite ([double]$Row.mean_apply_seconds))
}

function Test-E24CandidatePass([object]$Row, [double]$BaselineApply) {
  return (Test-E24LifecyclePass $Row) -and
    $Row.gmres_status -eq 'converged' -and
    $Row.audit_status -eq 'pass' -and
    $Row.ilu_graph -eq 'cartesian-near' -and
    (Test-E24Finite $BaselineApply) -and
    ([double]$Row.mean_apply_seconds -le 0.75 * $BaselineApply)
}

function Export-E24Aggregates([object[]]$Rows) {
  $Rows | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
  $Rows | Select-Object stage,fixture,theta,ilu_graph,status,gmres_status,iterations,convergence_reason,final_residual,internal_residual,last_preconditioned_residual,last_reported_true_residual,sampled_audit_rms,sampled_direct_true_residual |
    Export-Csv -LiteralPath $thetaSummaryCsv -NoTypeInformation
  $Rows | Select-Object stage,fixture,theta,ilu_graph,preconditioner_effective,preconditioner_graph_source,preconditioner_nnz_blocks,preconditioner_avg_degree,preconditioner_max_degree,peak_rss_mb |
    Export-Csv -LiteralPath $iluGraphCsv -NoTypeInformation
  $Rows | Select-Object stage,fixture,theta,ilu_graph,wall_seconds,first_apply_seconds,mean_apply_seconds,kifmm_apply_calls,far_pairs,fallback_pairs |
    Export-Csv -LiteralPath $runtimeCsv -NoTypeInformation
  $allInv = New-Object System.Collections.Generic.List[object]
  $allTrace = New-Object System.Collections.Generic.List[object]
  foreach ($row in $Rows) {
    $invPath = Join-Path $row.run_dir 'pbe_tabi_kifmm_iteration_invariants.csv'
    if (Test-Path -LiteralPath $invPath) {
      foreach ($inv in @(Import-Csv -LiteralPath $invPath)) {
        $inv | Add-Member -NotePropertyName stage -NotePropertyValue $row.stage -Force
        $inv | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
        $inv | Add-Member -NotePropertyName theta -NotePropertyValue $row.theta -Force
        $allInv.Add($inv) | Out-Null
      }
    }
    $tracePath = Join-Path $row.run_dir 'pbe_tabi_kifmm_gmres_convergence_trace.csv'
    if (Test-Path -LiteralPath $tracePath) {
      foreach ($trace in @(Import-Csv -LiteralPath $tracePath)) {
        $trace | Add-Member -NotePropertyName row_stage -NotePropertyValue $row.stage -Force
        $trace | Add-Member -NotePropertyName row_fixture -NotePropertyValue $row.fixture -Force
        $trace | Add-Member -NotePropertyName theta -NotePropertyValue $row.theta -Force
        $allTrace.Add($trace) | Out-Null
      }
    }
  }
  if ($allInv.Count -gt 0) { $allInv.ToArray() | Export-Csv -LiteralPath $invariantCsv -NoTypeInformation }
  if ($allTrace.Count -gt 0) { $allTrace.ToArray() | Export-Csv -LiteralPath $convergenceCsv -NoTypeInformation }
}

Assert-E24StaticContract

$quickStages = @(
  New-E24Stage 'theta0p6_1sh1' '1sh1' 0.6 'cartesian-near' 32 32 256 1800
  New-E24Stage 'theta0p7_1sh1' '1sh1' 0.7 'cartesian-near' 32 32 256 1800
  New-E24Stage 'theta0p75_1sh1' '1sh1' 0.75 'cartesian-near' 32 32 256 1800
  New-E24Stage 'theta0p8_1sh1' '1sh1' 0.8 'cartesian-near' 32 32 256 1800
)

$rows = New-Object System.Collections.Generic.List[object]
foreach ($stage in $quickStages) {
  $rows.Add((Invoke-E24Stage $stage)) | Out-Null
}

$baseline = @($rows | Where-Object { $_.fixture -eq '1sh1' -and [double]$_.theta -eq 0.6 -and $_.ilu_graph -eq 'cartesian-near' } | Select-Object -First 1)
$baselineApply = if ($baseline.Count -eq 1) { [double]$baseline[0].mean_apply_seconds } else { [double]::NaN }
$candidates = @($rows | Where-Object { [double]$_.theta -gt 0.6 -and (Test-E24CandidatePass $_ $baselineApply) } | Sort-Object @{Expression={[double]$_.mean_apply_seconds}; Ascending=$true})
$bestTheta = if ($candidates.Count -gt 0) { [double]$candidates[0].theta } else { 0.6 }
$observedThetaRows = @($rows | Where-Object {
    $_.fixture -eq '1sh1' -and $_.ilu_graph -eq 'cartesian-near' -and
    (Test-E24Finite ([double]$_.mean_apply_seconds))
  } | Sort-Object @{Expression={[double]$_.mean_apply_seconds}; Ascending=$true})
$bestObservedTheta = if ($observedThetaRows.Count -gt 0) { [double]$observedThetaRows[0].theta } else { [double]::NaN }
$bestObservedApply = if ($observedThetaRows.Count -gt 0) { [double]$observedThetaRows[0].mean_apply_seconds } else { [double]::NaN }
$bestObservedImprovementPct = if ((Test-E24Finite $baselineApply) -and (Test-E24Finite $bestObservedApply) -and $baselineApply -gt 0.0) {
  100.0 * ($baselineApply - $bestObservedApply) / $baselineApply
} else {
  [double]::NaN
}

if ($Full -and -not $Quick) {
  $fullStages = New-Object System.Collections.Generic.List[object]
  $fullStages.Add((New-E24Stage "theta$($bestTheta.ToString('0.##').Replace('.','p'))_1sh1_512" '1sh1' $bestTheta 'cartesian-near' 32 32 512 2400)) | Out-Null
  foreach ($theta in @(0.6,0.7,0.75,0.8)) {
    $fullStages.Add((New-E24Stage "theta$($theta.ToString('0.##').Replace('.','p'))_1a63_probe16" '1a63' $theta 'cartesian-near' 16 16 128 3600)) | Out-Null
  }
  if ($candidates.Count -gt 0 -or $AllowFailedDecision) {
    $graphTheta = $bestTheta
    $fullStages.Add((New-E24Stage "kifmm_graph_theta$($graphTheta.ToString('0.##').Replace('.','p'))_1sh1" '1sh1' $graphTheta 'kifmm-fallback-capped' 32 32 256 2400)) | Out-Null
    $fullStages.Add((New-E24Stage "kifmm_graph_theta$($graphTheta.ToString('0.##').Replace('.','p'))_1a63" '1a63' $graphTheta 'kifmm-fallback-capped' 16 16 128 3600)) | Out-Null
  }
  foreach ($stage in $fullStages) {
    $rows.Add((Invoke-E24Stage $stage)) | Out-Null
  }
}

$rowArray = $rows.ToArray()
Export-E24Aggregates $rowArray

$timeouts = @($rowArray | Where-Object { $_.status -eq 'timeout' })
$lifecycleGaps = @($rowArray | Where-Object { -not (Test-E24LifecyclePass $_) -and $_.status -ne 'timeout' })
$accuracyGaps = @($rowArray | Where-Object { $_.gmres_status -eq 'converged' -and $_.audit_status -eq 'fail' })
$quickCandidate = @($rowArray | Where-Object { [double]$_.theta -gt 0.6 -and (Test-E24CandidatePass $_ $baselineApply) } | Sort-Object @{Expression={[double]$_.mean_apply_seconds}; Ascending=$true})
$graphRows = @($rowArray | Where-Object { $_.ilu_graph -eq 'kifmm-fallback-capped' })
$graphCandidates = @($graphRows | Where-Object {
    (Test-E24LifecyclePass $_) -and
    $_.preconditioner_effective -eq 'fmm_nearfield_ilu0_kifmm_fallback_capped'
  })
$oneA63Rows = @($rowArray | Where-Object { $_.fixture -eq '1a63' -and $_.ilu_graph -eq 'cartesian-near' })
$oneA63Converged = @($oneA63Rows | Where-Object { $_.gmres_status -eq 'converged' -and $_.audit_status -eq 'pass' })

if ($timeouts.Count -gt 0) {
  $decision = 'fmm_e24_runtime_limited'
} elseif ($accuracyGaps.Count -gt 0) {
  $decision = 'fmm_e24_theta_accuracy_gap'
} elseif ($quickCandidate.Count -gt 0 -and $graphCandidates.Count -gt 0) {
  $decision = 'fmm_e24_kifmm_ilu_graph_candidate'
} elseif ($quickCandidate.Count -gt 0 -and $Full -and $oneA63Rows.Count -gt 0 -and $oneA63Converged.Count -eq 0) {
  $decision = 'fmm_e24_1a63_residual_ramp_pending'
} elseif ($quickCandidate.Count -gt 0) {
  $decision = 'fmm_e24_theta_candidate_pass'
} elseif ($graphRows.Count -gt 0 -and $graphCandidates.Count -eq 0) {
  $decision = 'fmm_e24_kifmm_ilu_graph_rejected'
} elseif ($lifecycleGaps.Count -gt 0) {
  $decision = 'fmm_e24_no_theta_speed_gain'
} else {
  $decision = 'fmm_e24_no_theta_speed_gain'
}

@([pscustomobject]@{
    decision=$decision
    rows=$rowArray.Count
    baseline_apply_seconds=$baselineApply
    best_observed_theta=$bestObservedTheta
    best_observed_apply_seconds=$bestObservedApply
    best_observed_apply_improvement_pct=$bestObservedImprovementPct
    recommended_theta=0.6
    theta_candidate_count=$quickCandidate.Count
    graph_rows=$graphRows.Count
    graph_candidate_count=$graphCandidates.Count
    timeout_rows=$timeouts.Count
  }) | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "rows=$($rowArray.Count)",
  "baseline_apply_seconds=$baselineApply",
  "best_observed_theta=$bestObservedTheta",
  "best_observed_apply_seconds=$bestObservedApply",
  "best_observed_apply_improvement_pct=$bestObservedImprovementPct",
  "recommended_theta=0.6",
  "candidate_theta_for_full_rows=$bestTheta",
  "matrix=$matrixCsv",
  "convergence=$convergenceCsv",
  "invariants=$invariantCsv"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

$report = @(
  'TITAN FMM-E24 KIFMM theta ramp and capped ILU-graph diagnostic',
  '================================================================',
  '',
  "Decision: $decision",
  '',
  'E24 is diagnostic-only. It relaxes only the explicit kifmm-lite-experimental theta whitelist and compares the existing Cartesian-near ILU graph against a capped KIFMM fallback-derived graph.',
  'The theta ramp is interpreted before ILU graph rows so runtime and preconditioner effects are not confounded.',
  'Closeout interpretation: fmm_e24_no_theta_speed_gain is a negative tuning result, not a numerical failure. It means flat KIFMM-lite has little remaining room from MAC-theta tuning because far work is still evaluated down to target panels.',
  ('Best observed theta was {0} with mean apply {1}s, an improvement of {2}% versus theta 0.6; this is below the 25% candidate gate. The recommended KIFMM-lite theta remains 0.6 for subsequent diagnostic baselines.' -f $bestObservedTheta,$bestObservedApply,$bestObservedImprovementPct),
  'Do not use additional ILU/preconditioner tuning to explain the measured bottleneck. The next useful work is an architecture change toward hierarchical box-level KIFMM, not another flat-worklist parameter sweep.',
  '',
  'Rows:'
)
foreach ($row in $rowArray) {
  $report += ('  {0}/{1}: theta={2}, graph={3}, gmres={4}, audit={5}, apply_mean={6}s, wall={7}s, precond={8}, nnz={9}, residual={10}, trace_rows={11}' -f
    $row.stage,$row.fixture,$row.theta,$row.ilu_graph,$row.gmres_status,
    $row.audit_status,$row.mean_apply_seconds,$row.wall_seconds,
    $row.preconditioner_effective,$row.preconditioner_nnz_blocks,
    $row.final_residual,$row.trace_rows)
}
$report += ''
$report += "Matrix: $matrixCsv"
$report += "Convergence trace: $convergenceCsv"
$report += "ILU graph summary: $iluGraphCsv"
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

Write-Host "FMM-E24 decision: $decision"
Write-Host "Matrix: $matrixCsv"
Write-Host "Convergence: $convergenceCsv"
if ($decision -notin @('fmm_e24_theta_candidate_pass','fmm_e24_1a63_residual_ramp_pending','fmm_e24_kifmm_ilu_graph_candidate','fmm_e24_kifmm_ilu_graph_rejected','fmm_e24_no_theta_speed_gain') -and -not $AllowFailedDecision) {
  throw "FMM-E24 gate did not pass: $decision"
}
