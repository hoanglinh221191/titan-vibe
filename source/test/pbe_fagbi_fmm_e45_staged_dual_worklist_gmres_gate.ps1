param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e45_staged_dual_worklist_gmres'),
  [int]$CaseTimeoutSeconds = 0,
  [switch]$Quick,
  [switch]$Full,
  [switch]$Force,
  [switch]$RerunFailed,
  [switch]$AllowFailedDecision,
  [switch]$Include1A63,
  [int]$Workers = 4,
  [ValidateSet('canonical-offset-cache','static-translation-table-cache','target-node-streaming-cache','bounded-lru-cache')]
  [string]$MatrixCachePolicy = 'canonical-offset-cache',
  [double]$MatrixCacheLimitMb = 0.0
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$TitanExe = [IO.Path]::GetFullPath($TitanExe)
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$matrixCsv = Join-Path $WorkRoot 'fmm_e45_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fmm_e45_convergence.csv'
$invariantCsv = Join-Path $WorkRoot 'fmm_e45_iteration_invariants.csv'
$bridgeCsv = Join-Path $WorkRoot 'fmm_e45_operator_bridge.csv'
$progressCsv = Join-Path $WorkRoot 'fmm_e45_progress.csv'
$runtimeCsv = Join-Path $WorkRoot 'fmm_e45_runtime.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e45_decision.csv'
$summaryDat = Join-Path $WorkRoot 'fmm_e45_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E45_STAGED_DUAL_WORKLIST_GMRES_REPORT.txt'
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

function Assert-E45StaticContract {
  $arguments = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\arguments.cpp')
  $parser = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\cli\parse_arguments.cpp')
  $solver = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp')
  $realPanelHeader = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.h')
  foreach ($needle in @(
      'kifmm-channel-hybrid-experimental',
      'pbe-bie-kifmm-hierarchy-m2l-eval-policy',
      'pbe-bie-kifmm-hierarchy-m2l-parallel-policy',
      'pbe-bie-cartesian-k4-geometry-theta')) {
    if (-not ($arguments.Contains($needle) -and $parser.Contains($needle))) {
      throw "Missing E45 parser contract marker: $needle"
    }
  }
  foreach ($needle in @(
      'prepare_dual_worklist_channel_hybrid_operator',
      'apply_prepared_dual_worklist_channel_hybrid_operator',
      'assemble_dual_worklist_channel_hybrid_tabi_system_matvec')) {
    if (-not ($solver.Contains($needle) -or $realPanelHeader.Contains($needle))) {
      throw "Missing E45 prepared operator marker: $needle"
    }
  }
  foreach ($needle in @(
      'pbe_tabi_channel_hybrid_iteration_invariants.csv',
      'pbe_tabi_channel_hybrid_gmres_progress.csv',
      'cartesian_fmm_m2l_near_panel_graph')) {
    if (-not $solver.Contains($needle)) {
      throw "Missing E45 solver contract marker: $needle"
    }
  }
}

function Get-E45Double([object]$Value, [double]$Default = [double]::NaN) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0.0
  if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Get-E45Int64([object]$Value, [int64]$Default = 0) {
  if ($null -eq $Value) { return $Default }
  $parsed = [int64]0
  if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Test-E45Finite([double]$Value) {
  return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Get-E45ProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
  if ($null -eq $Profile) { return $Default }
  $prop = $Profile.PSObject.Properties[$Name]
  if ($null -eq $prop) { return $Default }
  return $prop.Value
}

function Get-E45RuntimeRepoRoot {
  $root = $RepoRoot
  $ffProbe = 'lib_ff\ffAMBER99SBnonbonded.dat'
  if (Test-Path -LiteralPath (Join-Path $root $ffProbe)) { return $root }
  $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
  if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
    $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
  }
  $root = Split-Path -Parent $commonGitDir
  if (Test-Path -LiteralPath (Join-Path $root $ffProbe)) { return $root }
  $siblingMain = 'C:\tmp\elec_vdw-main-sync'
  if (Test-Path -LiteralPath (Join-Path $siblingMain $ffProbe)) { return $siblingMain }
  throw "Missing AMBER99SB force-field file under runtime repo root: $RepoRoot"
}

function Get-E45Fixture([string]$Name) {
  $manifest = @(Import-Csv -LiteralPath $manifestPath)
  $row = @($manifest | Where-Object { $_.fixture -eq $Name.ToLowerInvariant() })
  if ($row.Count -ne 1) { throw "Expected exactly one fixture '$Name', found $($row.Count)." }
  $pqr = Join-Path $FixtureRoot $row[0].pqr
  if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing PQR for $Name`: $pqr" }
  $row[0] | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqr -Force
  return $row[0]
}

function Write-E45ArgFile([object]$Fixture, [object]$Stage, [string]$RunDir, [string]$RuntimeRepoRoot) {
  $arg = Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
    -PqrPath $Fixture.pqr_path -Backend tabi -Mode treecode `
    -FmmMode single-matvec-analytic -SurfacePoints 8 `
    -MaxIter ([int]$Stage.max_iter) -Restart ([int]$Stage.restart) `
    -Tolerance 1e-6 -Threads $Workers -SaltM 0.15 `
    -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
    -OpeningTheta 0.6 -AuditSampleCount ([int]$Stage.sample_count) `
    -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
    -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
    -NearIntegration analytic_singular -TaylorOrder 4 `
    -FmmCartesianAudit off -FmmLocal off -FmmLocalRefine off `
    -FmmLocalRefineTolerance 0.002 -FmmLocalCache off `
    -FmmLocalCacheBudgetMb 256 -FmmLocalFallback scalar `
    -FmmLocalPlanBuilder legacy-rescan -FmmLocalRefineProbeCount 4 `
    -FmmLocalRefineProbeMode rademacher -FmmLocalGeometryGuard off `
    -FmmLocalGeometryTheta 0.0

  $libPath = Convert-FmmCliToTitanPath(
    (Resolve-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
  $lines = @(Get-Content -LiteralPath $arg) | ForEach-Object {
    if ($_ -like '-path-a=*') { "-path-a=$libPath" }
    elseif ($_ -like '-path-n=*') { "-path-n=$libPath" }
    elseif ($_ -eq '-pbe-tabi-tree-leaf-size=64') { '-pbe-tabi-tree-leaf-size=32' }
    elseif ($_ -like '-pbe-tabi-true-residual-sample-count=*') { "-pbe-tabi-true-residual-sample-count=$($Stage.sample_count)" }
    else { $_ }
  }
  $lines += @(
    '-pbe-tabi-fmm-engine=kifmm-channel-hybrid-experimental',
    '-pbe-tabi-kifmm-hierarchy-order=8',
    '-pbe-tabi-kifmm-hierarchy-theta=0.6',
    '-pbe-tabi-kifmm-hierarchy-real-panel-guard=strict-local-check',
    '-pbe-tabi-kifmm-hierarchy-guard-margin-factor=0',
    '-pbe-tabi-kifmm-hierarchy-tree-stop-policy=panel-radius-occupancy',
    '-pbe-tabi-kifmm-hierarchy-tree-panel-radius-factor=1.5',
    '-pbe-tabi-kifmm-hierarchy-tree-occupancy-stop=32',
    '-pbe-tabi-kifmm-hierarchy-m2l-eval-policy=target-check-batched',
    "-pbe-tabi-kifmm-hierarchy-matrix-cache-policy=$MatrixCachePolicy",
    "-pbe-tabi-kifmm-hierarchy-matrix-cache-limit-mb=$MatrixCacheLimitMb",
    '-pbe-tabi-kifmm-hierarchy-gemm-batch-policy=per-pair',
    '-pbe-tabi-kifmm-hierarchy-m2l-parallel-policy=target-node-key-group',
    "-pbe-tabi-kifmm-hierarchy-m2l-parallel-workers=$Workers",
    '-pbe-tabi-cartesian-k4-order=5',
    '-pbe-tabi-cartesian-k4-theta=0.6',
    '-pbe-tabi-cartesian-k4-geometry-theta=0.4',
    '-pbe-tabi-cartesian-k4-tail-tol=0.002'
  )
  $bad = @($lines | Where-Object { -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*=' })
  if ($bad.Count -ne 0) {
    throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
  }
  Set-Content -LiteralPath $arg -Value $lines -Encoding ASCII
  return $arg
}

function Read-E45InvariantSummary([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_channel_hybrid_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{
      rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0
      nonproduction_fallback_pairs=0; mean_apply_seconds=[double]::NaN
      preconditioner_graph_source='not_observed'; graph_source_stable='false'
      frozen_worklist='not_observed'; frozen_equivalent_system_cache='not_observed'
      prepared_rebuild_status='not_observed'; prepared_freeze_status='not_observed'
    }
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{
      rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0
      nonproduction_fallback_pairs=0; mean_apply_seconds=[double]::NaN
      preconditioner_graph_source='not_observed'; graph_source_stable='false'
      frozen_worklist='not_observed'; frozen_equivalent_system_cache='not_observed'
      prepared_rebuild_status='not_observed'; prepared_freeze_status='not_observed'
    }
  }
  $first = $rows[0]
  $hashStable = @($rows | Where-Object {
      $_.worklist_hash -ne $first.worklist_hash -or
      $_.policy_hash -ne $first.policy_hash -or
      $_.panel_geometry_hash -ne $first.panel_geometry_hash -or
      $_.sampled_target_hash -ne $first.sampled_target_hash
    }).Count -eq 0
  $cacheStable = @($rows | Where-Object { $_.cache_entries -ne $first.cache_entries }).Count -eq 0
  $firstGraphSource = [string](Get-E45ProfileValue $first 'preconditioner_graph_source' 'not_observed')
  $graphSourceStable = @($rows | Where-Object {
      [string](Get-E45ProfileValue $_ 'preconditioner_graph_source' 'not_observed') -ne $firstGraphSource
    }).Count -eq 0
  $firstFrozenWorklist = [string](Get-E45ProfileValue $first 'frozen_worklist' 'not_observed')
  $firstFrozenCache = [string](Get-E45ProfileValue $first 'frozen_equivalent_system_cache' 'not_observed')
  $firstRebuildStatus = [string](Get-E45ProfileValue $first 'prepared_rebuild_status' 'not_observed')
  $preparedFreezeStable = @($rows | Where-Object {
      [string](Get-E45ProfileValue $_ 'frozen_worklist' 'not_observed') -ne $firstFrozenWorklist -or
      [string](Get-E45ProfileValue $_ 'frozen_equivalent_system_cache' 'not_observed') -ne $firstFrozenCache -or
      [string](Get-E45ProfileValue $_ 'prepared_rebuild_status' 'not_observed') -ne $firstRebuildStatus
    }).Count -eq 0
  $preparedFreezeStatus = if ($preparedFreezeStable -and
      $firstFrozenWorklist -eq 'true' -and
      $firstFrozenCache -eq 'true' -and
      $firstRebuildStatus -eq 'frozen_prepared_apply') {
    'pass'
  } elseif ($firstRebuildStatus -eq 'apply_rebuilds_diagnostic_plan') {
    'apply_rebuilds_diagnostic_plan'
  } else {
    'fail'
  }
  $invalid = 0L
  $nonprod = 0L
  $applyTotal = 0.0
  foreach ($row in $rows) {
    $invalid += Get-E45Int64 $row.sampled_invalid_pairs
    $nonprod += Get-E45Int64 $row.nonproduction_fallback_pairs
    $applyTotal += Get-E45Double $row.apply_seconds 0.0
  }
  return [pscustomobject]@{
    rows=$rows.Count
    hash_stable=[string]$hashStable.ToString().ToLowerInvariant()
    cache_stable=[string]$cacheStable.ToString().ToLowerInvariant()
    invalid_pairs=$invalid
    nonproduction_fallback_pairs=$nonprod
    mean_apply_seconds=$applyTotal / [double]$rows.Count
    preconditioner_graph_source=$firstGraphSource
    graph_source_stable=[string]$graphSourceStable.ToString().ToLowerInvariant()
    frozen_worklist=$firstFrozenWorklist
    frozen_equivalent_system_cache=$firstFrozenCache
    prepared_rebuild_status=$firstRebuildStatus
    prepared_freeze_status=$preparedFreezeStatus
  }
}

function Read-E45Progress([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_channel_hybrid_gmres_progress.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'}
  }
  return [pscustomobject]@{rows=$rows.Count; last_phase=[string]$rows[-1].phase}
}

function Get-E45OperatorStatus([object]$Inv) {
  if ([int64]$Inv.rows -gt 0 -and $Inv.hash_stable -eq 'true' -and
      $Inv.cache_stable -eq 'true' -and [int64]$Inv.invalid_pairs -eq 0 -and
      [int64]$Inv.nonproduction_fallback_pairs -eq 0 -and
      $Inv.graph_source_stable -eq 'true' -and
      $Inv.preconditioner_graph_source -eq 'cartesian_fmm_m2l_near_panel_graph' -and
      $Inv.prepared_freeze_status -eq 'pass' -and
      (Test-E45Finite ([double]$Inv.mean_apply_seconds))) {
    return 'pass'
  }
  if ($Inv.prepared_freeze_status -eq 'apply_rebuilds_diagnostic_plan') {
    return 'freeze_pending'
  }
  return 'fail'
}

function Get-E45GmresStatus([int]$ExitCode, [string]$ProfileStatus, [string]$ConvergenceReason) {
  if ($ExitCode -ne 0) { return "exit_$ExitCode" }
  if ($ProfileStatus -eq 'ok') { return 'converged' }
  if ($ConvergenceReason -eq 'max_iter') { return 'max_iter_nonconverged' }
  if ([string]::IsNullOrWhiteSpace($ConvergenceReason)) { return 'not_observed' }
  return $ConvergenceReason
}

function Invoke-E45Stage([object]$Stage, [string]$RuntimeRepoRoot) {
  $fixture = Get-E45Fixture $Stage.fixture
  $runDir = Join-Path $WorkRoot "$($Stage.stage)_$($fixture.fixture)"
  $done = Join-Path $runDir 'fmm_e45_row.csv'
  if ((Test-Path -LiteralPath $done) -and -not $Force) {
    $existing = @(Import-Csv -LiteralPath $done)[0]
    $failed = $existing.status -eq 'timeout' -or $existing.status -like 'exit_*' -or $existing.operator_status -ne 'pass'
    if (-not ($RerunFailed -and $failed)) {
      Write-Host "resume: $($Stage.stage)_$($fixture.fixture)"
      return $existing
    }
  }
  if (Test-Path -LiteralPath $runDir) {
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
  } else {
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
  }
  $arg = Write-E45ArgFile -Fixture $fixture -Stage $Stage -RunDir $runDir -RuntimeRepoRoot $RuntimeRepoRoot
  $timeout = if ($CaseTimeoutSeconds -gt 0) { $CaseTimeoutSeconds } else { [int]$Stage.timeout }
  Write-Host "run: $($Stage.stage)_$($fixture.fixture) (timeout ${timeout}s)"
  $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg -CaseTimeoutSeconds $timeout
  $inv = Read-E45InvariantSummary $runDir
  $progress = Read-E45Progress $runDir
  if ($run.timed_out) {
    $row = [pscustomobject]@{
      stage=$Stage.stage; fixture=$fixture.fixture; status='timeout'; run_completed='false'
      execution_status='timeout'; gmres_status='not_observed_due_timeout'
      timeout_phase=$progress.last_phase; progress_rows=$progress.rows
      max_iter=$Stage.max_iter; restart=$Stage.restart; sampled_rows=$Stage.sample_count
      iterations='NA'; iteration_status='not_observed_due_timeout'
      internal_residual=[double]::NaN; sampled_audit_rms=[double]::NaN
      sampled_direct_true_residual=[double]::NaN
      preconditioner_effective='not_observed'
      preconditioner_graph_source=$inv.preconditioner_graph_source
      graph_source_stable=$inv.graph_source_stable
      frozen_worklist=$inv.frozen_worklist
      frozen_equivalent_system_cache=$inv.frozen_equivalent_system_cache
      prepared_rebuild_status=$inv.prepared_rebuild_status
      prepared_freeze_status=$inv.prepared_freeze_status
      plan_build_count=0; cache_build_count=0; kifmm_apply_calls=$inv.rows
      worklist_hash_stable=$inv.hash_stable; policy_hash_stable=$inv.hash_stable
      cache_entries_stable=$inv.cache_stable
      invalid_pairs=$inv.invalid_pairs; nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
      mean_apply_seconds=$inv.mean_apply_seconds
      wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
      operator_status='not_evaluated'; run_dir=$runDir
    }
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
  }
  $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
  if (-not (Test-Path -LiteralPath $profilePath)) {
    $row = [pscustomobject]@{
      stage=$Stage.stage; fixture=$fixture.fixture; status="exit_$($run.exit_code)"; run_completed='false'
      execution_status="exit_$($run.exit_code)"; gmres_status='not_observed_due_failure'
      timeout_phase=$progress.last_phase; progress_rows=$progress.rows
      max_iter=$Stage.max_iter; restart=$Stage.restart; sampled_rows=$Stage.sample_count
      iterations='NA'; iteration_status='not_observed_due_failure'
      internal_residual=[double]::NaN; sampled_audit_rms=[double]::NaN
      sampled_direct_true_residual=[double]::NaN
      preconditioner_effective='not_observed'
      preconditioner_graph_source=$inv.preconditioner_graph_source
      graph_source_stable=$inv.graph_source_stable
      frozen_worklist=$inv.frozen_worklist
      frozen_equivalent_system_cache=$inv.frozen_equivalent_system_cache
      prepared_rebuild_status=$inv.prepared_rebuild_status
      prepared_freeze_status=$inv.prepared_freeze_status
      plan_build_count=0; cache_build_count=0; kifmm_apply_calls=$inv.rows
      worklist_hash_stable=$inv.hash_stable; policy_hash_stable=$inv.hash_stable
      cache_entries_stable=$inv.cache_stable
      invalid_pairs=$inv.invalid_pairs; nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
      mean_apply_seconds=$inv.mean_apply_seconds
      wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
      operator_status='not_evaluated'; run_dir=$runDir
    }
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
  }
  $profile = @(Import-Csv -LiteralPath $profilePath)[0]
  $profileStatus = [string](Get-E45ProfileValue $profile 'status' '')
  $convergenceReason = [string](Get-E45ProfileValue $profile 'convergence_reason' '')
  $gmresStatus = Get-E45GmresStatus $run.exit_code $profileStatus $convergenceReason
  $operatorStatus = Get-E45OperatorStatus $inv
  $sampledAuditRms = Get-E45Double (Get-E45ProfileValue $profile 'tree_vs_direct_matvec_rel' ([double]::NaN))
  $sampledDirectTrueResidual = Get-E45Double (Get-E45ProfileValue $profile 'fmm_local_sampled_direct_true_residual' ([double]::NaN))
  $row = [pscustomobject]@{
    stage=$Stage.stage; fixture=$fixture.fixture
    status=$(if ($run.exit_code -eq 0 -and $profileStatus -eq 'ok') { 'ok' } elseif ($run.exit_code -eq 0) { 'completed_nonconverged' } else { "exit_$($run.exit_code)" })
    run_completed=$(if ($run.exit_code -eq 0 -or $profileStatus -ne '') { 'true' } else { 'false' })
    execution_status=$(if ($run.exit_code -eq 0) { 'completed' } else { "exit_$($run.exit_code)" })
    gmres_status=$gmresStatus; timeout_phase=$progress.last_phase; progress_rows=$progress.rows
    max_iter=$Stage.max_iter; restart=$Stage.restart; sampled_rows=$Stage.sample_count
    iterations=Get-E45Int64 (Get-E45ProfileValue $profile 'iterations' 0)
    iteration_status=$(if ($gmresStatus -eq 'not_observed_due_timeout') { 'not_observed_due_timeout' } else { $gmresStatus })
    internal_residual=Get-E45Double (Get-E45ProfileValue $profile 'true_residual_rel' ([double]::NaN))
    sampled_audit_rms=$sampledAuditRms
    sampled_direct_true_residual=$sampledDirectTrueResidual
    preconditioner_effective=[string](Get-E45ProfileValue $profile 'preconditioner_effective' '')
    preconditioner_graph_source=$inv.preconditioner_graph_source
    graph_source_stable=$inv.graph_source_stable
    frozen_worklist=$inv.frozen_worklist
    frozen_equivalent_system_cache=$inv.frozen_equivalent_system_cache
    prepared_rebuild_status=$inv.prepared_rebuild_status
    prepared_freeze_status=$inv.prepared_freeze_status
    plan_build_count=Get-E45Int64 (Get-E45ProfileValue $profile 'fmm_local_plan_build_count' 0)
    cache_build_count=Get-E45Int64 (Get-E45ProfileValue $profile 'fmm_local_cache_build_count' 0)
    kifmm_apply_calls=Get-E45Int64 (Get-E45ProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    worklist_hash_stable=$inv.hash_stable; policy_hash_stable=$inv.hash_stable
    cache_entries_stable=$inv.cache_stable
    invalid_pairs=$inv.invalid_pairs; nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
    mean_apply_seconds=$inv.mean_apply_seconds
    wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
    operator_status=$operatorStatus; run_dir=$runDir
  }
  $row | Export-Csv -LiteralPath $done -NoTypeInformation
  return $row
}

function Test-E45LifecyclePass([object]$Row) {
  return $Row.run_completed -eq 'true' -and
    $Row.execution_status -eq 'completed' -and
    $Row.operator_status -eq 'pass' -and
    $Row.preconditioner_graph_source -eq 'cartesian_fmm_m2l_near_panel_graph' -and
    $Row.graph_source_stable -eq 'true' -and
    $Row.prepared_freeze_status -eq 'pass' -and
    [int64]$Row.kifmm_apply_calls -gt 0 -and
    $Row.worklist_hash_stable -eq 'true' -and
    $Row.policy_hash_stable -eq 'true' -and
    $Row.cache_entries_stable -eq 'true' -and
    [int64]$Row.invalid_pairs -eq 0 -and
    [int64]$Row.nonproduction_fallback_pairs -eq 0
}

Assert-E45StaticContract
$RuntimeRepoRoot = Get-E45RuntimeRepoRoot

$stages = @(
  [pscustomobject]@{stage='matvec1'; fixture='1sh1'; max_iter=1; restart=8; sample_count=128; timeout=900},
  [pscustomobject]@{stage='short8'; fixture='1sh1'; max_iter=8; restart=8; sample_count=128; timeout=1800},
  [pscustomobject]@{stage='restart32'; fixture='1sh1'; max_iter=32; restart=32; sample_count=256; timeout=3600}
)
if ($Quick -and -not $Full) {
  $stages = @($stages | Where-Object { $_.stage -in @('matvec1', 'short8') })
}
if ($Include1A63) {
  $stages += [pscustomobject]@{stage='probe1a63'; fixture='1a63'; max_iter=1; restart=8; sample_count=128; timeout=3600}
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($stage in $stages) {
  $rows.Add((Invoke-E45Stage $stage $RuntimeRepoRoot)) | Out-Null
}

$rows.ToArray() | Export-Csv -NoTypeInformation -LiteralPath $matrixCsv
$rows.ToArray() | Select-Object stage,fixture,execution_status,gmres_status,iterations,iteration_status,internal_residual |
  Export-Csv -NoTypeInformation -LiteralPath $convergenceCsv
$rows.ToArray() | Select-Object stage,fixture,wall_seconds,peak_rss_mb,mean_apply_seconds,kifmm_apply_calls |
  Export-Csv -NoTypeInformation -LiteralPath $runtimeCsv

$allInv = New-Object System.Collections.Generic.List[object]
$allProgress = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
  $invPath = Join-Path $row.run_dir 'pbe_tabi_channel_hybrid_iteration_invariants.csv'
  if (Test-Path -LiteralPath $invPath) {
    foreach ($inv in @(Import-Csv -LiteralPath $invPath)) {
      $inv | Add-Member -NotePropertyName stage -NotePropertyValue $row.stage -Force
      $inv | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
      $allInv.Add($inv) | Out-Null
    }
  }
  $progressPath = Join-Path $row.run_dir 'pbe_tabi_channel_hybrid_gmres_progress.csv'
  if (Test-Path -LiteralPath $progressPath) {
    foreach ($pr in @(Import-Csv -LiteralPath $progressPath)) {
      $pr | Add-Member -NotePropertyName stage -NotePropertyValue $row.stage -Force
      $pr | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
      $allProgress.Add($pr) | Out-Null
    }
  }
}
if ($allInv.Count -gt 0) { $allInv.ToArray() | Export-Csv -NoTypeInformation -LiteralPath $invariantCsv }
if ($allProgress.Count -gt 0) { $allProgress.ToArray() | Export-Csv -NoTypeInformation -LiteralPath $progressCsv }

$bridgeRows = foreach ($row in $rows) {
  $status = if ($row.operator_status -eq 'pass') {
    'pass'
  } elseif ($row.prepared_freeze_status -eq 'apply_rebuilds_diagnostic_plan') {
    'lifecycle_gap'
  } else {
    'fail'
  }
  [pscustomobject]@{
    stage=$row.stage
    fixture=$row.fixture
    status=$status
    gmres_path='assemble_dual_worklist_channel_hybrid_tabi_system_matvec'
    full_panel_path='run_kifmm_hierarchy_real_panel_full_apply_sampled_audit'
    prepared_rebuild_status=$row.prepared_rebuild_status
    frozen_worklist=$row.frozen_worklist
    frozen_equivalent_system_cache=$row.frozen_equivalent_system_cache
    worklist_hash_stable=$row.worklist_hash_stable
    policy_hash_stable=$row.policy_hash_stable
    note=$(if ($status -eq 'pass') {
        'GMRES path used frozen prepared channel-hybrid apply.'
      } elseif ($status -eq 'lifecycle_gap') {
        'GMRES path delegates to deterministic full-panel apply, but the dual-worklist plan/cache are rebuilt per apply; frozen bridge parity remains pending.'
      } else {
        'Bridge evidence failed or was not observed.'
      })
  }
}
$bridgeRows | Export-Csv -NoTypeInformation -LiteralPath $bridgeCsv

$timeouts = @($rows | Where-Object { $_.status -eq 'timeout' })
$precondGaps = @($rows | Where-Object {
    $_.preconditioner_graph_source -ne 'cartesian_fmm_m2l_near_panel_graph' -or
    $_.graph_source_stable -ne 'true'
  })
$lifecycleGaps = @($rows | Where-Object { -not (Test-E45LifecyclePass $_) -and $_.status -ne 'timeout' })
$accuracyGaps = @($rows | Where-Object {
    $_.gmres_status -eq 'converged' -and
    ((Get-E45Double $_.sampled_audit_rms) -gt 8.0e-3 -or
     (Get-E45Double $_.sampled_direct_true_residual) -gt 8.0e-3)
  })
$quickPass = @($rows | Where-Object { $_.stage -in @('matvec1', 'short8') -and (Test-E45LifecyclePass $_) })
$restartPass = @($rows | Where-Object { $_.stage -eq 'restart32' -and (Test-E45LifecyclePass $_) })

if ($timeouts.Count -gt 0) {
  $decision = 'fmm_e45_runtime_limited'
} elseif ($precondGaps.Count -gt 0) {
  $decision = 'fmm_e45_preconditioner_graph_gap'
} elseif ($lifecycleGaps.Count -gt 0) {
  $decision = 'fmm_e45_lifecycle_gap'
} elseif ($accuracyGaps.Count -gt 0) {
  $decision = 'fmm_e45_operator_accuracy_gap'
} elseif ($Quick -and $quickPass.Count -eq 2) {
  $decision = 'fmm_e45_short_gmres_pass_restart_pending'
} elseif ($restartPass.Count -gt 0) {
  $decision = 'fmm_e45_staged_gmres_wiring_pass'
} else {
  $decision = 'fmm_e45_gmres_convergence_gap'
}

@([pscustomobject]@{
    decision=$decision
    rows=$rows.Count
    timeout_rows=$timeouts.Count
    lifecycle_gaps=$lifecycleGaps.Count
    preconditioner_gaps=$precondGaps.Count
  }) | Export-Csv -NoTypeInformation -LiteralPath $decisionCsv

@(
  "decision=$decision",
  "rows=$($rows.Count)",
  "matrix=$matrixCsv",
  "invariants=$invariantCsv",
  "progress=$progressCsv"
) | Set-Content -LiteralPath $summaryDat -Encoding ASCII

$report = @(
  'TITAN FMM-E45 staged dual-worklist channel-hybrid GMRES gate',
  '==============================================================',
  '',
  "Decision: $decision",
  '',
  'E45 is diagnostic-only. It wires the dual-worklist channel-hybrid operator into staged TABI GMRES without changing defaults or manuscript-facing policy.',
  'K1/K2/K3 use the q8 hierarchy channel path; K4 uses Cartesian p5 geometry0p4. ILU(0) remains sourced from cartesian_fmm_m2l_near_panel_graph.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}/{1}: execution={2}, gmres={3}, operator={4}, freeze={5}, graph={6}, iter={7}, apply_mean={8}s, phase={9}' -f
    $row.stage, $row.fixture, $row.execution_status, $row.gmres_status,
    $row.operator_status, $row.prepared_freeze_status,
    $row.preconditioner_graph_source, $row.iterations, $row.mean_apply_seconds,
    $row.timeout_phase)
}
$report += ''
$report += 'Bridge: GMRES-path channel-hybrid matvec uses a frozen prepared dual-worklist plan when freeze=pass. Rows marked lifecycle_gap indicate the frozen prepared-apply lifecycle was not observed.'
$report += "Bridge CSV: $bridgeCsv"
$report += "Matrix: $matrixCsv"
$report += "Invariants: $invariantCsv"
$report | Set-Content -LiteralPath $reportPath -Encoding ASCII

Write-Host "FMM-E45 decision: $decision"
Write-Host "Matrix: $matrixCsv"
Write-Host "Resume: rerun the same command to continue; use -RerunFailed for failed/timeout rows or -Force for all rows."
if ($decision -notin @('fmm_e45_short_gmres_pass_restart_pending','fmm_e45_staged_gmres_wiring_pass') -and -not $AllowFailedDecision) {
  throw "FMM-E45 gate did not pass: $decision"
}
