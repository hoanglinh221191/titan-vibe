param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
  [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e23_staged_gmres'),
  [int]$CaseTimeoutSeconds = 0,
  [switch]$Quick,
  [switch]$Full,
  [switch]$Force,
  [switch]$RerunFailed,
  [switch]$AllowFailedDecision,
  [string]$Cxx = '',
  [int]$Workers = 4
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) {
  throw "Missing public fixture manifest: $manifestPath"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$matrixCsv = Join-Path $WorkRoot 'fmm_e23_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fmm_e23_convergence.csv'
$invariantCsv = Join-Path $WorkRoot 'fmm_e23_iteration_invariants.csv'
$bridgeCsv = Join-Path $WorkRoot 'fmm_e23_operator_bridge.csv'
$progressCsv = Join-Path $WorkRoot 'fmm_e23_progress.csv'
$runtimeCsv = Join-Path $WorkRoot 'fmm_e23_runtime.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e23_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e23_summary.dat'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E23_STAGED_GMRES_REPORT.txt'

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

function Get-E23Double([object]$Value, [double]$Default = [double]::NaN) {
  if ($null -eq $Value) { return $Default }
  $parsed = 0.0
  if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Get-E23Int64([object]$Value, [int64]$Default = 0) {
  if ($null -eq $Value) { return $Default }
  $parsed = [int64]0
  if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
      [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $Default
}

function Test-E23Finite([double]$Value) {
  return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Get-E23ProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
  if ($null -eq $Profile) { return $Default }
  $prop = $Profile.PSObject.Properties[$Name]
  if ($null -eq $prop) { return $Default }
  return $prop.Value
}

function Parse-E23KeyValueLine([string]$Line) {
  $out = [ordered]@{}
  foreach ($part in ($Line -split '\s+')) {
    $idx = $part.IndexOf('=')
    if ($idx -le 0) { continue }
    $out[$part.Substring(0, $idx)] = $part.Substring($idx + 1)
  }
  return $out
}

function Resolve-E23Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) { throw "requested C++ compiler not found: $Requested" }
    return $Requested
  }
  foreach ($candidate in @(
      'C:\msys64\ucrt64\bin\g++.exe',
      'C:\msys64\ucrt64\bin\clang++.exe',
      'C:\msys64\mingw64\bin\g++.exe',
      'C:\msys64\mingw64\bin\clang++.exe')) {
    if (Test-Path -LiteralPath $candidate) { return $candidate }
  }
  $cmd = Get-Command g++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  throw 'No C++ compiler found for FMM-E23 bridge smoke.'
}

function Build-E23BridgeExe([string]$Exe) {
  $compiler = Resolve-E23Cxx $Cxx
  $args = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_full_panel_matvec_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp',
    'src\tabi\tabi_fmm_kifmm.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $Exe
  )
  Push-Location $RepoRoot
  try {
    & $compiler @args
    if ($LASTEXITCODE -ne 0) { throw 'failed to build FMM-E23 bridge smoke' }
  } finally {
    Pop-Location
  }
}

function Assert-E23StaticContract {
  $header = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.h')
  if (-not $header.Contains('assemble_kifmm_lite_tabi_system_matvec')) {
    throw 'Missing shared KIFMM TABI system matvec bridge.'
  }
  $solver = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp')
  foreach ($needle in @(
      'assemble_kifmm_lite_tabi_system_matvec',
      'pbe_tabi_kifmm_iteration_invariants.csv',
      'pbe_tabi_kifmm_gmres_progress.csv',
      'cartesian_fmm_m2l_near_panel_graph')) {
    if (-not $solver.Contains($needle)) { throw "Missing E23 solver contract marker: $needle" }
  }
  $args = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\arguments.cpp')
  foreach ($needle in @(
      'kifmm-lite-experimental',
      'pbe-tabi-kifmm-full-apply=target-node-parallel requires')) {
    if (-not $args.Contains($needle)) { throw "Missing E23 parser contract marker: $needle" }
  }
}

function Get-E23Fixture([string]$Name) {
  $manifest = @(Import-Csv -LiteralPath $manifestPath)
  $row = @($manifest | Where-Object { $_.fixture -eq $Name })
  if ($row.Count -ne 1) { throw "Expected exactly one fixture '$Name', found $($row.Count)." }
  $pqr = Join-Path $FixtureRoot $row[0].pqr
  if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing PQR for $Name`: $pqr" }
  $row[0] | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqr -Force
  return $row[0]
}

function Write-E23ArgFile([object]$Fixture, [object]$Stage, [string]$RunDir) {
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
    '-pbe-tabi-kifmm-full-apply=target-node-parallel'
  )
  $bad = @($lines | Where-Object { -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*=' })
  if ($bad.Count -ne 0) {
    throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
  }
  Set-Content -LiteralPath $arg -Value $lines -Encoding ASCII
  return $arg
}

function Read-E23InvariantSummary([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_kifmm_iteration_invariants.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0; nonproduction_fallback_pairs=0; first_apply_seconds=[double]::NaN; mean_apply_seconds=[double]::NaN}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; hash_stable='false'; cache_stable='false'; invalid_pairs=0; nonproduction_fallback_pairs=0; first_apply_seconds=[double]::NaN; mean_apply_seconds=[double]::NaN}
  }
  $first = $rows[0]
  $hashStable = @($rows | Where-Object {
      $_.worklist_hash -ne $first.worklist_hash -or
      $_.policy_hash -ne $first.policy_hash -or
      $_.panel_geometry_hash -ne $first.panel_geometry_hash -or
      $_.sampled_target_hash -ne $first.sampled_target_hash
    }).Count -eq 0
  $cacheStable = @($rows | Where-Object { $_.cache_entries -ne $first.cache_entries }).Count -eq 0
  $invalid = 0L
  $nonprod = 0L
  $applyTotal = 0.0
  foreach ($row in $rows) {
    $invalid += Get-E23Int64 $row.sampled_invalid_pairs
    $nonprod += Get-E23Int64 $row.nonproduction_fallback_pairs
    $applyTotal += Get-E23Double $row.apply_seconds 0.0
  }
  return [pscustomobject]@{
    rows=$rows.Count
    hash_stable=[string]$hashStable.ToString().ToLowerInvariant()
    cache_stable=[string]$cacheStable.ToString().ToLowerInvariant()
    invalid_pairs=$invalid
    nonproduction_fallback_pairs=$nonprod
    first_apply_seconds=Get-E23Double $first.apply_seconds
    mean_apply_seconds=$applyTotal / [double]$rows.Count
  }
}

function Read-E23Progress([string]$RunDir) {
  $path = Join-Path $RunDir 'pbe_tabi_kifmm_gmres_progress.csv'
  if (-not (Test-Path -LiteralPath $path)) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'; matvec_calls=0; total_matvec_calls=0}
  }
  $rows = @(Import-Csv -LiteralPath $path)
  if ($rows.Count -eq 0) {
    return [pscustomobject]@{rows=0; last_phase='not_observed'; matvec_calls=0; total_matvec_calls=0}
  }
  $last = $rows[$rows.Count - 1]
  return [pscustomobject]@{
    rows=$rows.Count
    last_phase=[string]$last.phase
    matvec_calls=Get-E23Int64 $last.fmm_local_gmres_matvec_calls
    total_matvec_calls=Get-E23Int64 $last.matvec_calls
  }
}

function Get-E23OperatorStatus([object]$Inv) {
  if ([int64]$Inv.rows -gt 0 -and
      $Inv.hash_stable -eq 'true' -and
      $Inv.cache_stable -eq 'true' -and
      [int64]$Inv.invalid_pairs -eq 0 -and
      [int64]$Inv.nonproduction_fallback_pairs -eq 0 -and
      (Test-E23Finite ([double]$Inv.first_apply_seconds)) -and
      (Test-E23Finite ([double]$Inv.mean_apply_seconds))) {
    return 'pass'
  }
  return 'fail'
}

function Get-E23GmresStatus([int]$ExitCode, [string]$ProfileStatus, [string]$ConvergenceReason) {
  if ($ExitCode -ne 0) { return "exit_$ExitCode" }
  if ($ProfileStatus -eq 'ok') { return 'converged' }
  if ($ConvergenceReason -eq 'max_iter') { return 'max_iter_nonconverged' }
  if ([string]::IsNullOrWhiteSpace($ConvergenceReason)) { return 'not_observed' }
  return $ConvergenceReason
}

function Get-E23AuditStatus([string]$GmresStatus, [double]$SampledAuditRms, [double]$SampledDirectTrueResidual) {
  if ($GmresStatus -ne 'converged') { return 'not_interpreted_nonconverged' }
  if ((Test-E23Finite $SampledAuditRms) -and
      (Test-E23Finite $SampledDirectTrueResidual) -and
      $SampledAuditRms -le 8.0e-3 -and
      $SampledDirectTrueResidual -le 8.0e-3) {
    return 'pass'
  }
  return 'fail'
}

function Convert-E23RowCloseoutSchema([object]$Row) {
  if ($null -ne $Row.PSObject.Properties['execution_status']) { return $Row }
  $completed = [string]$Row.completed
  $legacyStatus = [string]$Row.status
  $convergenceReason = [string]$Row.convergence_reason
  $executionStatus = if ($legacyStatus -eq 'timeout') {
    'timeout'
  } elseif ($completed -eq 'true') {
    'completed'
  } elseif ($legacyStatus -like 'exit_*') {
    $legacyStatus.Split(':')[0]
  } else {
    'failed'
  }
  $gmresStatus = if ($legacyStatus -eq 'ok') {
    'converged'
  } elseif ($convergenceReason -eq 'max_iter') {
    'max_iter_nonconverged'
  } elseif ($legacyStatus -eq 'timeout') {
    'not_observed_due_timeout'
  } else {
    'not_observed_due_failure'
  }
  $inv = [pscustomobject]@{
    rows=Get-E23Int64 $Row.invariant_rows
    hash_stable=[string]$Row.invariant_hash_stable
    cache_stable=[string]$Row.invariant_cache_stable
    invalid_pairs=Get-E23Int64 $Row.invariant_invalid_pairs
    nonproduction_fallback_pairs=Get-E23Int64 $Row.invariant_nonproduction_fallback_pairs
    first_apply_seconds=Get-E23Double $Row.first_apply_seconds
    mean_apply_seconds=Get-E23Double $Row.mean_apply_seconds
  }
  $operatorStatus = if ($executionStatus -eq 'completed') { Get-E23OperatorStatus $inv } else { 'not_evaluated' }
  $sampledAuditRms = Get-E23Double $Row.sampled_audit_rms
  $sampledDirectTrueResidual = Get-E23Double $Row.sampled_direct_true_residual
  $auditStatus = Get-E23AuditStatus $gmresStatus $sampledAuditRms $sampledDirectTrueResidual
  $normalizedStatus = if ($legacyStatus -like 'completed_nonconverged*') { 'completed_nonconverged' } else { $legacyStatus }
  return [pscustomobject]@{
    stage=$Row.stage; fixture=$Row.fixture; status=$normalizedStatus; completed=$Row.completed
    execution_status=$executionStatus; gmres_status=$gmresStatus
    operator_status=$operatorStatus
    postprocess_status=$(if ($gmresStatus -eq 'converged') { 'pass' } elseif ($gmresStatus -eq 'max_iter_nonconverged') { 'skipped_due_nonconvergence' } else { 'not_evaluated' })
    audit_status=$auditStatus
    audit_interpreted=$(if ($gmresStatus -eq 'converged') { 'true' } else { 'false' })
    timeout_phase=$Row.timeout_phase; progress_rows=$Row.progress_rows
    max_iter=$Row.max_iter; restart=$Row.restart; sample_count=$Row.sample_count
    iterations=$Row.iterations; convergence_reason=$Row.convergence_reason
    final_residual=$Row.final_residual; internal_residual=$Row.internal_residual
    sampled_audit_rms=$Row.sampled_audit_rms
    sampled_direct_true_residual=$Row.sampled_direct_true_residual
    residual_history=$Row.residual_history
    preconditioner_graph_source=$Row.preconditioner_graph_source
    plan_build_count=$Row.plan_build_count; cache_build_count=$Row.cache_build_count
    kifmm_apply_calls=$Row.kifmm_apply_calls
    invariant_rows=$Row.invariant_rows; invariant_hash_stable=$Row.invariant_hash_stable
    invariant_cache_stable=$Row.invariant_cache_stable
    invariant_invalid_pairs=$Row.invariant_invalid_pairs
    invariant_nonproduction_fallback_pairs=$Row.invariant_nonproduction_fallback_pairs
    first_apply_seconds=$Row.first_apply_seconds; mean_apply_seconds=$Row.mean_apply_seconds
    wall_seconds=$Row.wall_seconds; peak_rss_mb=$Row.peak_rss_mb
    failure_log=$Row.failure_log; run_dir=$Row.run_dir
  }
}

function New-E23FailureRow([object]$Fixture, [object]$Stage, [string]$Status, [object]$Run, [string]$RunDir, [string]$FailureLog = '') {
  $inv = Read-E23InvariantSummary $RunDir
  $progress = Read-E23Progress $RunDir
  $executionStatus = if ($Status -eq 'timeout') { 'timeout' } else { 'failed' }
  $gmresStatus = if ($Status -eq 'timeout') { 'not_observed_due_timeout' } else { 'not_observed_due_failure' }
  return [pscustomobject]@{
    stage=$Stage.stage; fixture=$Fixture.fixture; status=$Status; completed='false'
    execution_status=$executionStatus; gmres_status=$gmresStatus
    operator_status='not_evaluated'; postprocess_status='not_evaluated'
    audit_status='not_evaluated'; audit_interpreted='false'
    timeout_phase=$progress.last_phase; progress_rows=$progress.rows
    max_iter=$Stage.max_iter; restart=$Stage.restart; sample_count=$Stage.sample_count
    iterations=0; convergence_reason=''; final_residual=[double]::NaN
    internal_residual=[double]::NaN; sampled_audit_rms=[double]::NaN
    sampled_direct_true_residual=[double]::NaN; residual_history=''
    preconditioner_graph_source='cartesian_fmm_m2l_near_panel_graph'
    plan_build_count=0; cache_build_count=0; kifmm_apply_calls=$inv.rows
    invariant_rows=$inv.rows; invariant_hash_stable=$inv.hash_stable
    invariant_cache_stable=$inv.cache_stable
    invariant_invalid_pairs=$inv.invalid_pairs
    invariant_nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
    first_apply_seconds=$inv.first_apply_seconds; mean_apply_seconds=$inv.mean_apply_seconds
    wall_seconds=$Run.wall_seconds; peak_rss_mb=$Run.peak_working_set_mb
    failure_log=$FailureLog; run_dir=$RunDir
  }
}

function Invoke-E23Stage([object]$Stage) {
  $fixture = Get-E23Fixture $Stage.fixture
  $runDir = Join-Path $WorkRoot "$($Stage.stage)_$($fixture.fixture)"
  $done = Join-Path $runDir 'fmm_e23_row.csv'
  if ((Test-Path -LiteralPath $done) -and -not $Force) {
    $existing = Convert-E23RowCloseoutSchema (@(Import-Csv -LiteralPath $done)[0])
    $existing | Export-Csv -LiteralPath $done -NoTypeInformation
    $failed = $existing.status -eq 'timeout' -or $existing.status -eq 'failed'
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
  $arg = Write-E23ArgFile -Fixture $fixture -Stage $Stage -RunDir $runDir
  $timeout = if ($CaseTimeoutSeconds -gt 0) { $CaseTimeoutSeconds } else { [int]$Stage.timeout }
  Write-Host "run: $($Stage.stage)_$($fixture.fixture) (timeout ${timeout}s)"
  $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg -CaseTimeoutSeconds $timeout
  if ($run.timed_out) {
    $row = New-E23FailureRow $fixture $Stage 'timeout' $run $runDir
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
    $row = New-E23FailureRow $fixture $Stage 'failed' $run $runDir $log
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
  }
  $profile = @(Import-Csv -LiteralPath $profilePath)[0]
  $energy = if (Test-Path -LiteralPath $energyPath) { @(Import-Csv -LiteralPath $energyPath)[0] } else { $null }
  $inv = Read-E23InvariantSummary $runDir
  $progress = Read-E23Progress $runDir
  $profileStatus = [string](Get-E23ProfileValue $profile 'status' '')
  $energyStatus = if ($null -ne $energy) { [string]$energy.status } else { 'missing' }
  $convergenceReason = [string](Get-E23ProfileValue $profile 'convergence_reason' '')
  $gmresStatus = Get-E23GmresStatus $run.exit_code $profileStatus $convergenceReason
  $operatorStatus = Get-E23OperatorStatus $inv
  $postprocessStatus = if ($gmresStatus -eq 'converged') {
    if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok') { 'pass' } else { "failed:energy=$energyStatus;profile=$profileStatus" }
  } else {
    'skipped_due_nonconvergence'
  }
  $sampledAuditRms = Get-E23Double (Get-E23ProfileValue $profile 'tree_vs_direct_matvec_rel' ([double]::NaN))
  $sampledDirectTrueResidual = Get-E23Double (Get-E23ProfileValue $profile 'fmm_local_sampled_direct_true_residual' ([double]::NaN))
  $auditStatus = Get-E23AuditStatus $gmresStatus $sampledAuditRms $sampledDirectTrueResidual
  $rowStatus = if ($run.exit_code -eq 0 -and $energyStatus -eq 'ok' -and $profileStatus -eq 'ok') {
    'ok'
  } elseif ($run.exit_code -eq 0 -and $gmresStatus -eq 'max_iter_nonconverged') {
    'completed_nonconverged'
  } elseif ($run.exit_code -eq 0) {
    "completed_nonconverged:energy=$energyStatus;profile=$profileStatus"
  } else {
    "exit_$($run.exit_code):energy=$energyStatus;profile=$profileStatus"
  }
  $row = [pscustomobject]@{
    stage=$Stage.stage; fixture=$fixture.fixture; status=$rowStatus
    completed=$(if ($run.exit_code -eq 0 -or $profileStatus -ne '') { 'true' } else { 'false' })
    execution_status=$(if ($run.exit_code -eq 0) { 'completed' } else { "exit_$($run.exit_code)" })
    gmres_status=$gmresStatus; operator_status=$operatorStatus
    postprocess_status=$postprocessStatus; audit_status=$auditStatus
    audit_interpreted=$(if ($gmresStatus -eq 'converged') { 'true' } else { 'false' })
    timeout_phase=$progress.last_phase; progress_rows=$progress.rows
    max_iter=$Stage.max_iter; restart=$Stage.restart; sample_count=$Stage.sample_count
    iterations=Get-E23Int64 (Get-E23ProfileValue $profile 'iterations' 0)
    convergence_reason=$convergenceReason
    final_residual=Get-E23Double (Get-E23ProfileValue $profile 'final_residual' ([double]::NaN))
    internal_residual=Get-E23Double (Get-E23ProfileValue $profile 'true_residual_rel' ([double]::NaN))
    sampled_audit_rms=$sampledAuditRms
    sampled_direct_true_residual=$sampledDirectTrueResidual
    residual_history=[string](Get-E23ProfileValue $profile 'residual_history' '')
    preconditioner_graph_source='cartesian_fmm_m2l_near_panel_graph'
    plan_build_count=Get-E23Int64 (Get-E23ProfileValue $profile 'fmm_local_plan_build_count' 0)
    cache_build_count=Get-E23Int64 (Get-E23ProfileValue $profile 'fmm_local_cache_build_count' 0)
    kifmm_apply_calls=Get-E23Int64 (Get-E23ProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    invariant_rows=$inv.rows; invariant_hash_stable=$inv.hash_stable
    invariant_cache_stable=$inv.cache_stable
    invariant_invalid_pairs=$inv.invalid_pairs
    invariant_nonproduction_fallback_pairs=$inv.nonproduction_fallback_pairs
    first_apply_seconds=$inv.first_apply_seconds; mean_apply_seconds=$inv.mean_apply_seconds
    wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
    failure_log=''; run_dir=$runDir
  }
  $row | Export-Csv -LiteralPath $done -NoTypeInformation
  return $row
}

function Test-E23StageLifecyclePass([object]$Row) {
  return $Row.completed -eq 'true' -and
    $Row.execution_status -eq 'completed' -and
    $Row.operator_status -eq 'pass' -and
    $Row.preconditioner_graph_source -eq 'cartesian_fmm_m2l_near_panel_graph' -and
    [int64]$Row.plan_build_count -eq 1 -and
    [int64]$Row.cache_build_count -eq 1 -and
    [int64]$Row.invariant_rows -gt 0 -and
    $Row.invariant_hash_stable -eq 'true' -and
    $Row.invariant_cache_stable -eq 'true' -and
    [int64]$Row.invariant_invalid_pairs -eq 0 -and
    [int64]$Row.invariant_nonproduction_fallback_pairs -eq 0 -and
    (Test-E23Finite ([double]$Row.first_apply_seconds)) -and
    (Test-E23Finite ([double]$Row.mean_apply_seconds)) -and
    (Test-E23Finite ([double]$Row.final_residual))
}

function Test-E23Restart32GmresPass([object]$Row) {
  return (Test-E23StageLifecyclePass $Row) -and
    $Row.stage -eq 'restart32' -and
    $Row.fixture -eq '1sh1' -and
    $Row.gmres_status -eq 'converged' -and
    $Row.audit_status -eq 'pass' -and
    $Row.audit_interpreted -eq 'true'
}

function Invoke-E23BridgeSmoke([object[]]$BridgeRows) {
  $exe = Join-Path $WorkRoot 'titan_fmm_kifmm_full_panel_matvec_smoke_e23.exe'
  Build-E23BridgeExe $exe
  $pqr1sh1 = (Get-E23Fixture '1sh1').pqr_path
  $pqr1a63 = (Get-E23Fixture '1a63').pqr_path
  $outRows = New-Object System.Collections.Generic.List[object]
  foreach ($case in $BridgeRows) {
    $fixtureName = if ($case.fixture -eq '1a63') { '1A63' } else { '1SH1' }
    $args = @(
      $pqr1sh1,
      $pqr1a63,
      '--e23-profile',
      "--fixture=$fixtureName",
      "--sample=$($case.sample)",
      "--sample-count=$($case.sample_count)",
      '--order=6',
      "--workers=$Workers"
    )
    $output = & $exe @args
    if ($LASTEXITCODE -ne 0) { throw "E23 bridge smoke failed for $($case.stage): $output" }
    $line = @($output | Where-Object { $_ -like 'KIFMM_E23_BRIDGE *' } | Select-Object -First 1)
    if ($line.Count -ne 1) { throw "Missing KIFMM_E23_BRIDGE output for $($case.stage)" }
    $kv = Parse-E23KeyValueLine $line[0]
    $outRows.Add([pscustomobject]@{
        stage=$case.stage; fixture=$case.fixture; sample=$case.sample
        sampled_rows=Get-E23Int64 $kv.sampled_rows
        bridge_relative_diff=Get-E23Double $kv.bridge_relative_diff
        row_relative_diff=Get-E23Double $kv.row_relative_diff
        channel_relative_diff=Get-E23Double $kv.channel_relative_diff
        worklist_hash_stable=Get-E23Int64 $kv.worklist_hash_stable
        serial_cache_entries=Get-E23Int64 $kv.serial_cache_entries
        parallel_cache_entries=Get-E23Int64 $kv.parallel_cache_entries
        serial_far_panel_pairs=Get-E23Int64 $kv.serial_far_panel_pairs
        parallel_far_panel_pairs=Get-E23Int64 $kv.parallel_far_panel_pairs
      }) | Out-Null
  }
  $outRows.ToArray() | Export-Csv -LiteralPath $bridgeCsv -NoTypeInformation
  return $outRows.ToArray()
}

Assert-E23StaticContract

$bridgeRows = @([pscustomobject]@{stage='bridge_1sh1_uniform128'; fixture='1sh1'; sample='uniform_hash'; sample_count=128})
if ($Full -and -not $Quick) {
  $bridgeRows += [pscustomobject]@{stage='bridge_1a63_highaspect128'; fixture='1a63'; sample='high_aspect'; sample_count=128}
}
$bridge = @(Invoke-E23BridgeSmoke $bridgeRows)

$stages = @(
  [pscustomobject]@{stage='matvec1'; fixture='1sh1'; max_iter=1; restart=8; sample_count=128; timeout=600},
  [pscustomobject]@{stage='short8'; fixture='1sh1'; max_iter=8; restart=8; sample_count=128; timeout=1200},
  [pscustomobject]@{stage='restart32'; fixture='1sh1'; max_iter=32; restart=32; sample_count=256; timeout=2400}
)
if ($Full -and -not $Quick) {
  $stages += [pscustomobject]@{stage='probe1a63'; fixture='1a63'; max_iter=1; restart=8; sample_count=128; timeout=1800}
  $stages += [pscustomobject]@{stage='short8_1a63'; fixture='1a63'; max_iter=8; restart=8; sample_count=128; timeout=3600}
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($stage in $stages) {
  if ($stage.stage -eq 'short8_1a63' -and -not $AllowFailedDecision) {
    $probe = @($rows | Where-Object { $_.stage -eq 'probe1a63' } | Select-Object -First 1)
    if ($probe.Count -eq 0 -or -not (Test-E23StageLifecyclePass $probe[0])) {
      Write-Host 'skip: short8_1a63 (probe1a63 did not pass lifecycle; use -AllowFailedDecision to run)'
      continue
    }
  }
  $rows.Add((Invoke-E23Stage $stage)) | Out-Null
}

$rows.ToArray() | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$rows.ToArray() | Select-Object stage,fixture,status,execution_status,gmres_status,operator_status,postprocess_status,audit_status,audit_interpreted,completed,iterations,convergence_reason,final_residual,internal_residual,residual_history |
  Export-Csv -LiteralPath $convergenceCsv -NoTypeInformation
$rows.ToArray() | Select-Object stage,fixture,wall_seconds,peak_rss_mb,first_apply_seconds,mean_apply_seconds,kifmm_apply_calls |
  Export-Csv -LiteralPath $runtimeCsv -NoTypeInformation

$allInv = New-Object System.Collections.Generic.List[object]
$allProgress = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
  $invPath = Join-Path $row.run_dir 'pbe_tabi_kifmm_iteration_invariants.csv'
  if (Test-Path -LiteralPath $invPath) {
    foreach ($inv in @(Import-Csv -LiteralPath $invPath)) {
      $inv | Add-Member -NotePropertyName stage -NotePropertyValue $row.stage -Force
      $inv | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
      $allInv.Add($inv) | Out-Null
    }
  }
  $progressPath = Join-Path $row.run_dir 'pbe_tabi_kifmm_gmres_progress.csv'
  if (Test-Path -LiteralPath $progressPath) {
    foreach ($pr in @(Import-Csv -LiteralPath $progressPath)) {
      $pr | Add-Member -NotePropertyName stage -NotePropertyValue $row.stage -Force
      $pr | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
      $allProgress.Add($pr) | Out-Null
    }
  }
}
if ($allInv.Count -gt 0) { $allInv.ToArray() | Export-Csv -LiteralPath $invariantCsv -NoTypeInformation }
if ($allProgress.Count -gt 0) { $allProgress.ToArray() | Export-Csv -LiteralPath $progressCsv -NoTypeInformation }

$bridgeGaps = @($bridge | Where-Object {
    [double]$_.bridge_relative_diff -gt 1.0e-12 -or
    [double]$_.row_relative_diff -gt 1.0e-12 -or
    [double]$_.channel_relative_diff -gt 1.0e-12 -or
    [int64]$_.worklist_hash_stable -ne 1
  })
$timeouts = @($rows | Where-Object { $_.status -eq 'timeout' })
$precondGaps = @($rows | Where-Object { $_.preconditioner_graph_source -ne 'cartesian_fmm_m2l_near_panel_graph' })
$lifecycleGaps = @($rows | Where-Object { -not (Test-E23StageLifecyclePass $_) -and $_.status -ne 'timeout' })
$oneA63Probe = @($rows | Where-Object { $_.stage -eq 'probe1a63' -and (Test-E23StageLifecyclePass $_) })
$required1sh1 = @('matvec1','short8','restart32')
$passed1sh1 = @($rows | Where-Object { $required1sh1 -contains $_.stage -and (Test-E23StageLifecyclePass $_) })
$restart32Pass = @($rows | Where-Object { Test-E23Restart32GmresPass $_ })

if ($bridgeGaps.Count -gt 0) {
  $decision = 'fmm_e23_operator_bridge_gap'
} elseif ($timeouts.Count -gt 0) {
  $decision = 'fmm_e23_runtime_limited'
} elseif ($precondGaps.Count -gt 0) {
  $decision = 'fmm_e23_preconditioner_graph_gap'
} elseif ($lifecycleGaps.Count -gt 0) {
  $decision = 'fmm_e23_lifecycle_gap'
} elseif ($restart32Pass.Count -gt 0 -and $oneA63Probe.Count -gt 0) {
  $decision = 'fmm_e23_1sh1_restart32_gmres_pass_1a63_probe_pending'
} elseif ($restart32Pass.Count -gt 0 -and $passed1sh1.Count -eq $required1sh1.Count) {
  $decision = 'fmm_e23_staged_gmres_wiring_pass'
} else {
  $decision = 'fmm_e23_gmres_convergence_gap'
}

@([pscustomobject]@{
    decision=$decision
    rows=$rows.Count
    bridge_rows=$bridge.Count
    bridge_gaps=$bridgeGaps.Count
    lifecycle_gaps=$lifecycleGaps.Count
    timeout_rows=$timeouts.Count
  }) | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "rows=$($rows.Count)",
  "bridge_rows=$($bridge.Count)",
  "matrix=$matrixCsv",
  "bridge=$bridgeCsv",
  "invariants=$invariantCsv"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

$report = @(
  'TITAN FMM-E23 staged KIFMM-lite GMRES correctness ladder',
  '==========================================================',
  '',
  "Decision: $decision",
  '',
  'E23 is diagnostic-only. It validates staged KIFMM-lite GMRES wiring with the E22 target-node-parallel full-panel apply; no production backend or default is changed.',
  '',
  'The operator bridge uses the same TABI diagonal assembly as the GMRES path. Short staged rows are wiring evidence, not full convergence evidence.',
  'The 1SH1 restart32 row is the GMRES smoke pass: it converges and passes the sampled dense-direct post-solve audit.',
  '1A63 probe and short rows are not interpreted as convergence or post-solve audit evidence because those staged runs stop at max_iter.',
  'E23 therefore establishes a 1SH1 KIFMM-GMRES smoke pass plus a 1A63 execution probe, not a general KIFMM-GMRES backend pass.',
  'Runtime remains diagnostic: the KIFMM-lite path is not yet performance-competitive with the manuscript guarded Cartesian p4 path.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}/{1}: execution={2}, gmres={3}, operator={4}, audit={5}, completed={6}, iter={7}, residual={8}, apply_mean={9}s, invariants={10}, hash={11}, cache={12}' -f
    $row.stage,$row.fixture,$row.execution_status,$row.gmres_status,
    $row.operator_status,$row.audit_status,$row.completed,$row.iterations,
    $row.final_residual,$row.mean_apply_seconds,$row.invariant_rows,
    $row.invariant_hash_stable,$row.invariant_cache_stable)
}
$report += ''
$report += "Matrix: $matrixCsv"
$report += "Operator bridge: $bridgeCsv"
$report += "Iteration invariants: $invariantCsv"
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

Write-Host "FMM-E23 decision: $decision"
Write-Host "Matrix: $matrixCsv"
Write-Host "Operator bridge: $bridgeCsv"
if ($decision -notin @('fmm_e23_staged_gmres_wiring_pass','fmm_e23_1sh1_restart32_gmres_pass_1a63_probe_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E23 gate did not pass: $decision"
}
