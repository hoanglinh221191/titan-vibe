[CmdletBinding()]
param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string[]]$Fixtures = @('binding-complex','1sh1'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5d1r',
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [int]$MaxIter = 0,
    [string[]]$Stages = @(),
    [switch]$Quick,
    [switch]$RerunFailed,
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$Fixtures = @($Fixtures | ForEach-Object {
    @([string]$_ -split ',') | ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})
if ($Quick) { $Fixtures = @('1sh1') }
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$MatrixPath = Join-Path $WorkRoot 'calderon_c5d1r_regularized_ilu_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5d1r_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5d1r_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_CALDERON_C5D1R_REGULARIZED_ILU_REPORT.txt'

$stageDefs = @(
    [pscustomobject]@{ name='ilu0_off'; mode='off'; max=0.01 },
    [pscustomobject]@{ name='adaptive_1e6'; mode='adaptive'; max=1.0e-6 },
    [pscustomobject]@{ name='adaptive_1e4'; mode='adaptive'; max=1.0e-4 },
    [pscustomobject]@{ name='adaptive_1e2'; mode='adaptive'; max=1.0e-2 }
)
if ($Quick) {
    $stageDefs = @($stageDefs | Where-Object {
        $_.name -in @('ilu0_off','adaptive_1e2')
    })
}
$Stages = @($Stages | ForEach-Object {
    @([string]$_ -split ',') | ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})
if ($Stages.Count -gt 0) {
    $stageDefs = @($stageDefs | Where-Object { $_.name -in $Stages })
    if ($stageDefs.Count -ne $Stages.Count) {
        throw "Unknown C5D.1R stage. Valid stages: ilu0_off, adaptive_1e6, adaptive_1e4, adaptive_1e2."
    }
}
if ($MaxIter -le 0) { $MaxIter = if ($Quick) { 80 } else { 240 } }

function Resolve-C5D1rFixture([string]$Name) {
    switch ($Name) {
        'binding-complex' {
            return Join-Path $RepoRoot `
                'test\pbe_inputs_binding_jpcb\1brs_complex_amber_ph7.pqr'
        }
        '1sh1' {
            return Join-Path $RepoRoot `
                'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
        }
        default { throw "Unsupported C5D.1R fixture '$Name'." }
    }
}

function Get-ProfileValue([object]$Profile, [string]$Name,
        [string]$Default = '') {
    if ($null -eq $Profile -or
        -not ($Profile.PSObject.Properties.Name -contains $Name)) {
        return $Default
    }
    return [string]$Profile.$Name
}

function Read-C5D1rRow([string]$Fixture, [object]$Stage,
        [object]$Invocation, [string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profile = if (Test-Path -LiteralPath $profilePath) {
        @(Import-Csv -LiteralPath $profilePath)[0]
    } else { $null }
    $energy = if (Test-Path -LiteralPath $energyPath) {
        @(Import-Csv -LiteralPath $energyPath)[0]
    } else { $null }
    $runStatus = if ($Invocation.timed_out) { 'timeout' }
        elseif ($Invocation.exit_code -ne 0) { 'failed' }
        elseif ($null -eq $profile) { 'missing_profile' }
        else { 'ok' }
    return [pscustomobject]@{
        fixture = $Fixture
        stage = $Stage.name
        regularization_mode = $Stage.mode
        regularization_max = $Stage.max
        run_status = $runStatus
        exit_code = $Invocation.exit_code
        wall_seconds = $Invocation.wall_seconds
        peak_rss_mb = $Invocation.peak_working_set_mb
        status = Get-ProfileValue $profile 'status' 'missing'
        panels = Get-ProfileValue $profile 'panels' '0'
        iterations = Get-ProfileValue $profile 'iterations' '0'
        true_residual = Get-ProfileValue $profile 'true_residual_rel' 'NaN'
        audit_rms = Get-ProfileValue $profile 'tree_vs_direct_matvec_rel' 'NaN'
        energy_kj_mol = if ($null -eq $energy) { 'NaN' } else {
            Get-ProfileValue $energy 'E_rxn_kJmol' 'NaN'
        }
        preconditioner_effective = Get-ProfileValue $profile `
            'preconditioner_effective' 'missing'
        quality_status = Get-ProfileValue $profile `
            'preconditioner_quality_status' 'missing'
        unshifted_failure_stage = Get-ProfileValue $profile `
            'preconditioner_ilu_unshifted_failure_stage' 'missing'
        selected_relative_shift = Get-ProfileValue $profile `
            'preconditioner_ilu_regularization_selected' 'NaN'
        selected_absolute_shift_max = Get-ProfileValue $profile `
            'preconditioner_ilu_regularization_absolute_max' 'NaN'
        regularization_attempts = Get-ProfileValue $profile `
            'preconditioner_ilu_regularization_attempts' '0'
        rejected_attempts = Get-ProfileValue $profile `
            'preconditioner_ilu_regularization_rejected_attempts' '0'
        shifted_rows = Get-ProfileValue $profile `
            'preconditioner_ilu_regularization_shifted_rows' '0'
        pivot_ratio_min = Get-ProfileValue $profile `
            'preconditioner_pivot_ratio_min' 'NaN'
        factor_growth_max = Get-ProfileValue $profile `
            'preconditioner_factor_growth_max' 'NaN'
        quality_probe_gain = Get-ProfileValue $profile `
            'preconditioner_quality_probe_max_gain' 'NaN'
        quality_probe_residual = Get-ProfileValue $profile `
            'preconditioner_quality_probe_rel_residual' 'NaN'
        setup_seconds = Get-ProfileValue $profile `
            'preconditioner_setup_seconds' 'NaN'
        apply_seconds = Get-ProfileValue $profile `
            'preconditioner_apply_seconds' 'NaN'
    }
}

$existing = @{}
if ((Test-Path -LiteralPath $MatrixPath) -and -not $Force) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $existing["$($row.fixture)|$($row.stage)"] = $row
    }
}
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in $Fixtures) {
    $pqr = [IO.Path]::GetFullPath((Resolve-C5D1rFixture $fixture))
    if (-not (Test-Path -LiteralPath $pqr)) {
        throw "C5D.1R fixture is missing: $pqr"
    }
    foreach ($stage in $stageDefs) {
        $key = "$fixture|$($stage.name)"
        $prior = $existing[$key]
        if ($null -ne $prior -and -not $Force -and
            (-not $RerunFailed -or [string]$prior.run_status -eq 'ok')) {
            Write-Host "resume: $($fixture)_$($stage.name)"
            $rows.Add($prior)
            continue
        }
        $runDir = Join-Path $WorkRoot "$($fixture)_$($stage.name)"
        New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
        $auditRows = if ($Quick) { 128 } else { 256 }
        $salt = if ($fixture -eq 'binding-complex') { 0.05 } else { 0.15 }
        $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
            -PqrPath $pqr -Backend tabi -Mode treecode `
            -FmmMode single-matvec-analytic -SurfacePoints 8 `
            -MaxIter $MaxIter -Restart 64 -Tolerance 1e-6 -Threads $Threads `
            -SaltM $salt -InternalDielectric 1 -ExternalDielectric 80 `
            -TemperatureK 298.15 -OpeningTheta 0.6 `
            -AuditSampleCount $auditRows -AuditMode sampled `
            -Preconditioner fmm-nearfield-ilu0 `
            -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
            -PreconditionerIluRegularization $stage.mode `
            -PreconditionerIluRegularizationMax $stage.max `
            -NearIntegration analytic_singular -TaylorOrder 4
        Write-Host "run: $($fixture)_$($stage.name) (timeout ${CaseTimeoutSeconds}s)"
        $invocation = Invoke-FmmCliTitanCase -TitanExe $TitanExe `
            -RunDir $runDir -ArgPath $arg `
            -CaseTimeoutSeconds $CaseTimeoutSeconds
        $rows.Add((Read-C5D1rRow -Fixture $fixture -Stage $stage `
            -Invocation $invocation -RunDir $runDir))
        $rows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
    }
}
$rows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$timedOut = @($rows | Where-Object { $_.run_status -eq 'timeout' }).Count -gt 0
$failed = @($rows | Where-Object {
    $_.run_status -notin @('ok','timeout')
}).Count -gt 0
$rescues = @($rows | Where-Object {
    $_.preconditioner_effective -eq 'fmm_nearfield_ilu0_regularized' -and
    $_.quality_status -eq 'pass' -and [double]$_.selected_relative_shift -gt 0
})
$baselinePasses = @($rows | Where-Object {
    $_.stage -eq 'ilu0_off' -and
    $_.preconditioner_effective -eq 'fmm_nearfield_ilu0' -and
    $_.quality_status -eq 'pass'
})
$rescueConverged = @($rescues | Where-Object {
    $_.status -eq 'ok' -and [double]$_.true_residual -le 1.0e-4
})
$decision = if ($timedOut) { 'c5d1r_runtime_limited' }
    elseif ($failed) { 'c5d1r_harness_gap' }
    elseif ($rescues.Count -gt 0 -and $rescueConverged.Count -eq $rescues.Count) {
        'c5d1r_regularized_ilu_rescue'
    } elseif ($rescues.Count -gt 0) {
        'c5d1r_factor_rescue_convergence_gap'
    } elseif ($baselinePasses.Count -eq $Fixtures.Count) {
        'c5d1r_regularization_not_needed'
    } else { 'c5d1r_regularization_quality_gap' }

[pscustomobject]@{
    decision = $decision
    fixtures = $Fixtures.Count
    rows = $rows.Count
    regularized_rescues = $rescues.Count
    regularized_converged = $rescueConverged.Count
    baseline_passes = $baselinePasses.Count
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    "fixtures=$($Fixtures.Count)",
    "rows=$($rows.Count)",
    "regularized_rescues=$($rescues.Count)",
    "regularized_converged=$($rescueConverged.Count)",
    "baseline_passes=$($baselinePasses.Count)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$table = @($rows | ForEach-Object {
    "$($_.fixture),$($_.stage),$($_.run_status),$($_.status)," +
    "$($_.preconditioner_effective),$($_.quality_status)," +
    "$($_.unshifted_failure_stage),$($_.selected_relative_shift)," +
    "$($_.regularization_attempts),$($_.true_residual)," +
    "$($_.setup_seconds),$($_.peak_rss_mb)"
})
@(
    'TITAN C5D.1R Adaptive Block-Regularized ILU(0) Gate',
    '',
    "Decision: $decision",
    '',
    'The adaptive path retries the same 2x2 block-CSR ILU(0) sparsity with',
    'bounded row-scaled diagonal shifts. Every accepted factor is quality-',
    'probed against the original unshifted near-field matrix.',
    '',
    'fixture,stage,run_status,status,effective,quality,unshifted_failure,' +
        'selected_shift,attempts,true_residual,setup_seconds,peak_rss_mb',
    $table,
    '',
    'Sparse ILU defaults to the automatic ILU(0) -> regularized ILU ->',
    'diagonal pipeline; the existing off and maximum-shift options remain.'
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5D.1R decision: $decision"
Write-Host "Matrix: $MatrixPath"
if (-not $AllowFailedDecision -and $decision -notin @(
        'c5d1r_regularized_ilu_rescue','c5d1r_regularization_not_needed')) {
    throw "C5D.1R gate did not pass: $decision"
}
