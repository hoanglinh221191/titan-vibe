param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e21_kifmm_gmres'),
    [int]$CaseTimeoutSeconds = 1800,
    [int]$MaxIter = 720,
    [int]$Restart = 128,
    [switch]$Quick,
    [switch]$Full,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
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

$matrixCsv = Join-Path $WorkRoot 'fmm_e21_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fmm_e21_convergence.csv'
$invariantCsv = Join-Path $WorkRoot 'fmm_e21_iteration_invariants.csv'
$auditCsv = Join-Path $WorkRoot 'fmm_e21_operator_audit.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e21_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e21_summary.dat'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E21_KIFMM_GMRES_REPORT.txt'

function Get-E21Double([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-E21Int64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-E21ProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Test-E21Finite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Set-E21RowProperty([object]$Row, [string]$Name, [object]$Value) {
    if ($null -ne $Row.PSObject.Properties[$Name]) {
        $Row.$Name = $Value
    } else {
        $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $Value -Force
    }
}

function Get-E21Fixture([string]$Name) {
    $manifest = @(Import-Csv -LiteralPath $manifestPath)
    $row = @($manifest | Where-Object { $_.fixture -eq $Name })
    if ($row.Count -ne 1) { throw "Expected exactly one fixture '$Name', found $($row.Count)." }
    $pqr = Join-Path $FixtureRoot $row[0].pqr
    if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing PQR for $Name`: $pqr" }
    $row[0] | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqr -Force
    return $row[0]
}

function Write-E21ArgFile([object]$Fixture, [string]$RunDir) {
    $arg = Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $Fixture.pqr_path -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter $MaxIter `
        -Restart $Restart -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off -FmmLocal off `
        -FmmLocalRefine off -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache off -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback scalar -FmmLocalPlanBuilder legacy-rescan `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard off -FmmLocalGeometryTheta 0.0

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
        '-pbe-tabi-kifmm-profile=summary'
    )
    $bad = @($lines | Where-Object { -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*=' })
    if ($bad.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
    }
    Set-Content -LiteralPath $arg -Value $lines -Encoding ASCII
    return $arg
}

function Read-E21ProgressSummary([string]$RunDir) {
    $path = Join-Path $RunDir 'pbe_tabi_kifmm_gmres_progress.csv'
    if (-not (Test-Path -LiteralPath $path)) {
        return [pscustomobject]@{
            rows=0; timeout_phase='not_observed'; last_phase='not_observed'
            kifmm_matvec_calls=0; total_matvec_calls=0
        }
    }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -eq 0) {
        return [pscustomobject]@{
            rows=0; timeout_phase='not_observed'; last_phase='not_observed'
            kifmm_matvec_calls=0; total_matvec_calls=0
        }
    }
    $last = $rows[$rows.Count - 1]
    return [pscustomobject]@{
        rows=$rows.Count
        timeout_phase=[string]$last.phase
        last_phase=[string]$last.phase
        kifmm_matvec_calls=Get-E21Int64 $last.fmm_local_gmres_matvec_calls
        total_matvec_calls=Get-E21Int64 $last.matvec_calls
    }
}

function New-E21FailureRow([object]$Fixture, [string]$Status, [object]$Run, [string]$RunDir, [string]$FailureLog = '') {
    $invariants = Read-E21InvariantSummary $RunDir
    $progress = Read-E21ProgressSummary $RunDir
    return [pscustomobject]@{
        fixture=$Fixture.fixture; pdb_id=$Fixture.pdb_id; protein_name=$Fixture.protein_name
        status=$Status; convergence_reason=''; panels=0; iterations=0
        run_completed='false'; timeout_phase=$progress.timeout_phase
        preconditioner_graph_source='cartesian_fmm_m2l_near_panel_graph'
        progress_rows=$progress.rows; iterations_observed=0
        iteration_status=$(if ($Status -eq 'timeout') { 'not_observed_due_timeout' } else { 'not_observed_due_failure' })
        energy_kj_mol=[double]::NaN; internal_residual=[double]::NaN
        final_residual=[double]::NaN; sampled_audit_rms=[double]::NaN
        sampled_direct_true_residual=[double]::NaN; wall_seconds=$Run.wall_seconds
        peak_rss_mb=$Run.peak_working_set_mb; matvec_seconds_avg=[double]::NaN
        matvec_calls=0; gmres_seconds=[double]::NaN
        matvec_mode=''; fmm_m2l_mode=''; fmm_m2l_setup_seconds=[double]::NaN
        fmm_m2l_apply_seconds=[double]::NaN; fmm_m2l_near_pairs=0; fmm_m2l_far_pairs=0
        fmm_local_gmres_active='false'; kifmm_prepare_seconds=0.0
        kifmm_apply_calls=0; kifmm_apply_seconds_total=0.0
        kifmm_apply_seconds_mean=0.0; kifmm_m2l_seconds_total=0.0
        kifmm_fallback_seconds_total=0.0; plan_build_count=0; cache_build_count=0
        operator_finite=$(if ($invariants.rows -gt 0 -and $invariants.invalid_pairs -eq 0 -and $invariants.nonproduction_fallback_pairs -eq 0) { 'true' } else { 'false' })
        invariant_rows=$invariants.rows
        invariant_hash_stable=$invariants.hash_stable
        invariant_cache_stable=$invariants.cache_stable
        invariant_invalid_pairs=$invariants.invalid_pairs
        invariant_nonproduction_fallback_pairs=$invariants.nonproduction_fallback_pairs
        failure_log=$FailureLog; run_dir=$RunDir
    }
}

function Read-E21InvariantSummary([string]$RunDir) {
    $path = Join-Path $RunDir 'pbe_tabi_kifmm_iteration_invariants.csv'
    if (-not (Test-Path -LiteralPath $path)) {
        return [pscustomobject]@{
            rows=0; hash_stable='false'; cache_stable='false'
            invalid_pairs=0; nonproduction_fallback_pairs=0
        }
    }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -eq 0) {
        return [pscustomobject]@{
            rows=0; hash_stable='false'; cache_stable='false'
            invalid_pairs=0; nonproduction_fallback_pairs=0
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
    $invalid = 0L
    $nonprod = 0L
    foreach ($row in $rows) {
        $invalid += Get-E21Int64 $row.sampled_invalid_pairs
        $nonprod += Get-E21Int64 $row.nonproduction_fallback_pairs
    }
    return [pscustomobject]@{
        rows=$rows.Count; hash_stable=[string]$hashStable.ToString().ToLowerInvariant()
        cache_stable=[string]$cacheStable.ToString().ToLowerInvariant()
        invalid_pairs=$invalid; nonproduction_fallback_pairs=$nonprod
    }
}

function Update-E21ExistingRowWithInvariants([object]$Row) {
    $invariants = Read-E21InvariantSummary $Row.run_dir
    $progress = Read-E21ProgressSummary $Row.run_dir
    Set-E21RowProperty $Row 'invariant_rows' $invariants.rows
    Set-E21RowProperty $Row 'invariant_hash_stable' $invariants.hash_stable
    Set-E21RowProperty $Row 'invariant_cache_stable' $invariants.cache_stable
    Set-E21RowProperty $Row 'invariant_invalid_pairs' $invariants.invalid_pairs
    Set-E21RowProperty $Row 'invariant_nonproduction_fallback_pairs' $invariants.nonproduction_fallback_pairs
    Set-E21RowProperty $Row 'timeout_phase' $progress.timeout_phase
    Set-E21RowProperty $Row 'preconditioner_graph_source' 'cartesian_fmm_m2l_near_panel_graph'
    Set-E21RowProperty $Row 'progress_rows' $progress.rows
    $observedIterations = Get-E21Int64 $Row.iterations 0
    if ($Row.status -eq 'ok') {
        Set-E21RowProperty $Row 'run_completed' 'true'
        Set-E21RowProperty $Row 'iterations_observed' $observedIterations
        Set-E21RowProperty $Row 'iteration_status' 'observed'
    } elseif ($observedIterations -gt 0) {
        Set-E21RowProperty $Row 'run_completed' 'false'
        Set-E21RowProperty $Row 'iterations_observed' $observedIterations
        Set-E21RowProperty $Row 'iteration_status' 'observed_nonconverged'
    } else {
        Set-E21RowProperty $Row 'run_completed' 'false'
        Set-E21RowProperty $Row 'iterations_observed' 0
        Set-E21RowProperty $Row 'iteration_status' $(if ($Row.status -eq 'timeout') { 'not_observed_due_timeout' } else { 'not_observed_due_failure' })
    }
    if ($invariants.rows -gt 0 -and $invariants.invalid_pairs -eq 0 -and
        $invariants.nonproduction_fallback_pairs -eq 0) {
        Set-E21RowProperty $Row 'operator_finite' 'true'
        Set-E21RowProperty $Row 'fmm_local_gmres_active' 'true'
        Set-E21RowProperty $Row 'kifmm_apply_calls' $invariants.rows
    }
    return $Row
}

function Invoke-E21Case([object]$Fixture) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_kifmm_gmres"
    $done = Join-Path $runDir 'fmm_e21_row.csv'
    if ((Test-Path -LiteralPath $done) -and -not $Force) {
        $existing = @(Import-Csv -LiteralPath $done)[0]
        $failed = $existing.status -ne 'ok'
        if (-not ($RerunFailed -and $failed)) {
            Write-Host "resume: $($Fixture.fixture)_kifmm_gmres"
            $updated = Update-E21ExistingRowWithInvariants $existing
            $updated | Export-Csv -LiteralPath $done -NoTypeInformation
            return $updated
        }
    }

    if (Test-Path -LiteralPath $runDir) {
        $resolvedRunDir = [IO.Path]::GetFullPath($runDir)
        $resolvedWorkRoot = [IO.Path]::GetFullPath($WorkRoot)
        if (-not $resolvedRunDir.StartsWith($resolvedWorkRoot, [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove run directory outside WorkRoot: $resolvedRunDir"
        }
        Remove-Item -LiteralPath $resolvedRunDir -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg = Write-E21ArgFile -Fixture $Fixture -RunDir $runDir
    Write-Host "run: $($Fixture.fixture)_kifmm_gmres (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        $row = New-E21FailureRow $Fixture 'timeout' $run $runDir
        $row | Export-Csv -LiteralPath $done -NoTypeInformation
        return $row
    }
    if ($run.exit_code -ne 0) {
        $stdoutPath = Join-Path $runDir 'stdout.log'
        $stderrPath = Join-Path $runDir 'stderr.log'
        $stdoutText = ''
        $stderrText = ''
        if (Test-Path -LiteralPath $stdoutPath) {
            $stdoutText = Get-Content -Raw -LiteralPath $stdoutPath
        }
        if (Test-Path -LiteralPath $stderrPath) {
            $stderrText = Get-Content -Raw -LiteralPath $stderrPath
        }
        $log = $stdoutText + $stderrText
        $row = New-E21FailureRow $Fixture 'failed' $run $runDir $log
        $row | Export-Csv -LiteralPath $done -NoTypeInformation
        return $row
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        $row = New-E21FailureRow $Fixture 'missing_energy_output' $run $runDir
        $row | Export-Csv -LiteralPath $done -NoTypeInformation
        return $row
    }
    if (-not (Test-Path -LiteralPath $profilePath)) {
        $row = New-E21FailureRow $Fixture 'missing_tabi_profile' $run $runDir
        $row | Export-Csv -LiteralPath $done -NoTypeInformation
        return $row
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    $invariants = Read-E21InvariantSummary $runDir
    $progress = Read-E21ProgressSummary $runDir
    $energyStatus = [string]$energy.status
    $profileStatus = [string](Get-E21ProfileValue $profile 'status' '')
    $rowStatus = if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok') {
        'ok'
    } else {
        "energy=$energyStatus;profile=$profileStatus"
    }
    $observedIterations = Get-E21Int64 (Get-E21ProfileValue $profile 'iterations' 0)
    $row = [pscustomobject]@{
        fixture=$Fixture.fixture; pdb_id=$Fixture.pdb_id; protein_name=$Fixture.protein_name
        status=$rowStatus; convergence_reason=[string](Get-E21ProfileValue $profile 'convergence_reason' '')
        panels=Get-E21Int64 (Get-E21ProfileValue $profile 'panels' 0)
        iterations=$observedIterations
        run_completed=$(if ($rowStatus -eq 'ok') { 'true' } else { 'false' })
        timeout_phase=$progress.timeout_phase
        preconditioner_graph_source='cartesian_fmm_m2l_near_panel_graph'
        progress_rows=$progress.rows
        iterations_observed=$observedIterations
        iteration_status=$(if ($rowStatus -eq 'ok') { 'observed' } elseif ($observedIterations -gt 0) { 'observed_nonconverged' } else { 'not_observed_due_failure' })
        energy_kj_mol=Get-E21Double $energy.E_rxn_kJmol
        internal_residual=Get-E21Double (Get-E21ProfileValue $profile 'true_residual_rel' ([double]::NaN))
        final_residual=Get-E21Double (Get-E21ProfileValue $profile 'final_residual' ([double]::NaN))
        sampled_audit_rms=Get-E21Double (Get-E21ProfileValue $profile 'tree_vs_direct_matvec_rel' ([double]::NaN))
        sampled_direct_true_residual=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_sampled_direct_true_residual' 0.0) 0.0
        wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
        matvec_seconds_avg=Get-E21Double (Get-E21ProfileValue $profile 'matvec_seconds_avg' ([double]::NaN))
        matvec_calls=Get-E21Int64 (Get-E21ProfileValue $profile 'matvec_calls' 0)
        gmres_seconds=Get-E21Double (Get-E21ProfileValue $profile 'gmres_seconds' 0.0) 0.0
        matvec_mode=[string](Get-E21ProfileValue $profile 'matvec_mode' '')
        fmm_m2l_mode=[string](Get-E21ProfileValue $profile 'fmm_m2l_mode' '')
        fmm_m2l_setup_seconds=Get-E21Double (Get-E21ProfileValue $profile 'fmm_m2l_setup_seconds' 0.0) 0.0
        fmm_m2l_apply_seconds=Get-E21Double (Get-E21ProfileValue $profile 'fmm_m2l_apply_seconds' 0.0) 0.0
        fmm_m2l_near_pairs=Get-E21Int64 (Get-E21ProfileValue $profile 'fmm_m2l_near_pairs' 0)
        fmm_m2l_far_pairs=Get-E21Int64 (Get-E21ProfileValue $profile 'fmm_m2l_far_pairs' 0)
        fmm_local_gmres_active=[string](Get-E21ProfileValue $profile 'fmm_local_gmres_active' 'false')
        kifmm_prepare_seconds=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_prepare_seconds' 0.0) 0.0
        kifmm_apply_calls=Get-E21Int64 (Get-E21ProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
        kifmm_apply_seconds_total=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_gmres_apply_seconds_total' 0.0) 0.0
        kifmm_apply_seconds_mean=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_gmres_apply_seconds_mean' 0.0) 0.0
        kifmm_m2l_seconds_total=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_m2l_seconds_total' 0.0) 0.0
        kifmm_fallback_seconds_total=Get-E21Double (Get-E21ProfileValue $profile 'fmm_local_fallback_seconds_total' 0.0) 0.0
        plan_build_count=Get-E21Int64 (Get-E21ProfileValue $profile 'fmm_local_plan_build_count' 0)
        cache_build_count=Get-E21Int64 (Get-E21ProfileValue $profile 'fmm_local_cache_build_count' 0)
        operator_finite=[string](Get-E21ProfileValue $profile 'fmm_local_operator_finite' 'false')
        invariant_rows=$invariants.rows
        invariant_hash_stable=$invariants.hash_stable
        invariant_cache_stable=$invariants.cache_stable
        invariant_invalid_pairs=$invariants.invalid_pairs
        invariant_nonproduction_fallback_pairs=$invariants.nonproduction_fallback_pairs
        failure_log=''; run_dir=$runDir
    }
    $row | Export-Csv -LiteralPath $done -NoTypeInformation
    return $row
}

function Test-E21Pass([object]$Row) {
    return $Row.status -eq 'ok' -and
        (Test-E21Finite ([double]$Row.internal_residual)) -and
        (Test-E21Finite ([double]$Row.sampled_audit_rms)) -and
        (Test-E21Finite ([double]$Row.sampled_direct_true_residual)) -and
        (Test-E21Finite ([double]$Row.energy_kj_mol)) -and
        [double]$Row.internal_residual -le 1.0e-6 -and
        [double]$Row.sampled_audit_rms -le 8.0e-3 -and
        [double]$Row.sampled_direct_true_residual -le 8.0e-3 -and
        $Row.fmm_local_gmres_active -eq 'true' -and
        [int64]$Row.plan_build_count -eq 1 -and
        [int64]$Row.cache_build_count -eq 1 -and
        [int64]$Row.kifmm_apply_calls -gt 0 -and
        $Row.operator_finite -eq 'true' -and
        [int64]$Row.invariant_rows -gt 0 -and
        $Row.invariant_hash_stable -eq 'true' -and
        $Row.invariant_cache_stable -eq 'true' -and
        [int64]$Row.invariant_invalid_pairs -eq 0 -and
        [int64]$Row.invariant_nonproduction_fallback_pairs -eq 0
}

$fixtures = @('1sh1')
if ($Full -and -not $Quick) { $fixtures = @('1sh1','1a63') }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($name in $fixtures) {
    if ($name -eq '1a63' -and -not $AllowFailedDecision) {
        $first = @($rows | Where-Object { $_.fixture -eq '1sh1' })[0]
        if (-not (Test-E21Pass $first)) {
            Write-Host 'skip: 1a63 (1SH1 did not pass; use -AllowFailedDecision to run confirm row)'
            continue
        }
    }
    $rows.Add((Invoke-E21Case (Get-E21Fixture $name))) | Out-Null
}

$rows.ToArray() | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$rows.ToArray() | Select-Object fixture,status,run_completed,timeout_phase,preconditioner_graph_source,progress_rows,iteration_status,iterations_observed,convergence_reason,iterations,internal_residual,final_residual,residual_history |
    Export-Csv -LiteralPath $convergenceCsv -NoTypeInformation
$rows.ToArray() | Select-Object fixture,sampled_audit_rms,sampled_direct_true_residual,matvec_mode,fmm_m2l_mode,matvec_calls,kifmm_apply_calls |
    Export-Csv -LiteralPath $auditCsv -NoTypeInformation

$allInvariantRows = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
    $path = Join-Path $row.run_dir 'pbe_tabi_kifmm_iteration_invariants.csv'
    if (Test-Path -LiteralPath $path) {
        foreach ($inv in @(Import-Csv -LiteralPath $path)) {
            $inv | Add-Member -NotePropertyName fixture -NotePropertyValue $row.fixture -Force
            $allInvariantRows.Add($inv) | Out-Null
        }
    }
}
if ($allInvariantRows.Count -gt 0) {
    $allInvariantRows.ToArray() | Export-Csv -LiteralPath $invariantCsv -NoTypeInformation
} else {
    @([pscustomobject]@{fixture='none'; matvec_call=0; status='not_available'}) |
        Export-Csv -LiteralPath $invariantCsv -NoTypeInformation
}

$passRows = @($rows | Where-Object { Test-E21Pass $_ })
$decision = 'fmm_e21_gmres_convergence_gap'
if (@($rows | Where-Object { $_.status -eq 'timeout' }).Count -gt 0) {
    $decision = 'fmm_e21_runtime_limited'
} elseif (@($rows | Where-Object {
        $_.status -eq 'ok' -and
        (([double]$_.sampled_audit_rms -gt 8.0e-3) -or
         ([double]$_.sampled_direct_true_residual -gt 8.0e-3))
    }).Count -gt 0) {
    $decision = 'fmm_e21_operator_accuracy_gap'
} elseif (@($rows | Where-Object {
        $_.status -eq 'ok' -and
        ($_.invariant_hash_stable -ne 'true' -or $_.invariant_cache_stable -ne 'true' -or
         [int64]$_.invariant_invalid_pairs -ne 0 -or [int64]$_.invariant_nonproduction_fallback_pairs -ne 0)
    }).Count -gt 0) {
    $decision = 'fmm_e21_lifecycle_gap'
} elseif ($fixtures.Count -eq 1 -and $passRows.Count -eq 1) {
    $decision = 'fmm_e21_kifmm_gmres_1sh1_pass'
} elseif ($fixtures.Count -eq 2 -and $passRows.Count -eq 2) {
    $decision = 'fmm_e21_kifmm_gmres_pass'
}

@([pscustomobject]@{
    decision=$decision; rows=$rows.Count; pass_rows=$passRows.Count
    quick=[string]($Quick -or -not $Full).ToString().ToLowerInvariant()
    max_iter=$MaxIter; restart=$Restart
}) | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation

@(
    "decision=$decision",
    "rows=$($rows.Count)",
    "pass_rows=$($passRows.Count)",
    "max_iter=$MaxIter",
    "restart=$Restart",
    "matrix=$matrixCsv",
    "invariants=$invariantCsv"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

$report = @(
    'TITAN FMM-E21 KIFMM-lite GMRES smoke',
    '=====================================',
    '',
    "Decision: $decision",
    '',
    'E21 is diagnostic-only. It enables the explicit hidden kifmm-lite-experimental engine for TABI GMRES under the frozen E20.1 prepared-apply protocol; no default or production backend is changed.',
    '',
    'Timeout rows are not convergence evidence. The first opt-in GMRES attempt was runtime-limited before a complete convergence/audit record was produced.',
    '',
    $(if ($decision -eq 'fmm_e21_runtime_limited') {
        'This run is runtime-limited; no convergence conclusion is drawn from rows without a completed post-solve audit.'
      } else {
        'Completed rows are interpreted only when convergence, lifecycle invariants, and post-solve sampled dense-direct audit records are all present.'
      }),
    '',
    'For this E21.1 smoke, fmm-nearfield-ilu0 still uses the existing Cartesian FMM/M2L near-panel graph, not the KIFMM fallback worklist; preconditioner phases are recorded separately so timeout evidence can distinguish setup from Krylov/apply time.',
    '',
    'Rows:'
)
foreach ($row in $rows) {
    $report += ('  {0}: status={1}, completed={2}, phase={3}, precond_graph={4}, observed_iter={5}, residual={6}, audit={7}, direct_residual={8}, apply_mean={9}s, invariant_hash={10}, cache={11}' -f
        $row.fixture,$row.status,$row.run_completed,$row.timeout_phase,
        $row.preconditioner_graph_source,$row.iterations_observed,$row.internal_residual,
        $row.sampled_audit_rms,$row.sampled_direct_true_residual,
        $row.kifmm_apply_seconds_mean,$row.invariant_hash_stable,
        $row.invariant_cache_stable)
}
$report += ''
$report += "Matrix: $matrixCsv"
$report += "Iteration invariants: $invariantCsv"
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

Write-Host "FMM-E21 decision: $decision"
Write-Host "Matrix: $matrixCsv"
Write-Host "Iteration invariants: $invariantCsv"
if ($decision -notin @('fmm_e21_kifmm_gmres_1sh1_pass','fmm_e21_kifmm_gmres_pass') -and -not $AllowFailedDecision) {
    throw "FMM-E21 failed decision: $decision"
}
