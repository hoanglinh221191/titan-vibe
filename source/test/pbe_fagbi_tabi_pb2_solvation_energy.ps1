[CmdletBinding()]
param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string[]]$Fixtures = @('1a63'),
    [string]$WorkRoot = 'test\tmp\fagbi_tabi_pb2_solvation_energy',
    [ValidateSet('sphere-patch','grid-ses','analytic-ses')]
    [string]$PrimarySurface = 'sphere-patch',
    [ValidateRange(0.25,2.0)]
    [double]$SesGridSpacing = 0.75,
    [ValidateSet('green-jet-compressed','green-jet-tiled')]
    [string]$FmmLocalCache = 'green-jet-compressed',
    [ValidateRange(64,8192)]
    [int]$FmmLocalCacheBudgetMb = 256,
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$RerunFailed,
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ($Fixtures.Count -eq 1 -and $Fixtures[0] -match ',') {
    $Fixtures = @($Fixtures[0] -split ',' | ForEach-Object { $_.Trim() } |
        Where-Object { $_ -ne '' })
}

$fixtureMap = @{
    'barnase' = 'test\pbe_inputs_public\barnase_full_amber_ph7.pqr'
    '1a63' = 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
    '1sh1' = 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
    '1ubq' = 'test\pbe_inputs_public_z4q\1ubq_full_amber_ph7.pqr'
    '2ci2' = 'test\pbe_inputs_public_z4q\2ci2_full_amber_ph7.pqr'
    '1vii' = 'test\pbe_inputs_public_z4q\1vii_full_amber_ph7.pqr'
}

function Read-Double([object]$Value) {
    $parsed = 0.0
    if (-not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float, $InvariantCulture,
            [ref]$parsed)) { return [double]::NaN }
    return $parsed
}

function Test-Finite([double]$Value) {
    return -not ([double]::IsNaN($Value) -or [double]::IsInfinity($Value))
}

function New-FailureRow([string]$Fixture, [string]$Failure,
        [object]$Run = $null) {
    return [pscustomobject][ordered]@{
        fixture = $Fixture; pdb_id = $Fixture.ToUpperInvariant()
        atoms = 0; panels = 0; status = 'failed'; failure_class = $Failure
        deltaG_solv_kcalmol = [double]::NaN
        deltaG_solv_kJmol = [double]::NaN
        reaction_energy_kcalmol = [double]::NaN
        iterations = 0; true_residual = [double]::NaN
        sampled_operator_audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        wall_seconds = if ($null -ne $Run) { $Run.wall_seconds } else { 0.0 }
        peak_rss_mb = if ($null -ne $Run) { $Run.peak_working_set_mb } else { 0.0 }
        preconditioner = ''; convergence_reason = ''
        dgsolv_status = 'failed'; salt_M = 0.15; temperature_K = 298.15
        plan_build_count = 0; cache_build_count = 0
        primary_surface = $PrimarySurface; ses_grid_spacing_A = $SesGridSpacing
        fmm_local_cache = $FmmLocalCache
        fmm_local_cache_budget_mb = $FmmLocalCacheBudgetMb
        run_dir = ''
    }
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureRaw in $Fixtures) {
    $fixture = $fixtureRaw.ToLowerInvariant()
    if (-not $fixtureMap.ContainsKey($fixture)) {
        throw "Unsupported fixture '$fixtureRaw'. Choose: $($fixtureMap.Keys -join ', ')."
    }
    $pqr = Join-Path $RepoRoot $fixtureMap[$fixture]
    if (-not (Test-Path -LiteralPath $pqr)) {
        throw "Missing full public PQR fixture for ${fixture}: $pqr"
    }
    $runDir = Join-Path $WorkRoot $fixture
    $rowPath = Join-Path $runDir 'solvation_energy_row.csv'
    if (-not $Force -and (Test-Path -LiteralPath $rowPath)) {
        $saved = @(Import-Csv -LiteralPath $rowPath)[0]
        $savedOk = [string]$saved.status -eq 'ok'
        $savedPrimary = if ($saved.PSObject.Properties.Name -contains 'primary_surface') {
            [string]$saved.primary_surface
        } else { 'sphere-patch' }
        $savedSpacing = if ($saved.PSObject.Properties.Name -contains 'ses_grid_spacing_A') {
            Read-Double $saved.ses_grid_spacing_A
        } else { 0.75 }
        $savedCache = if ($saved.PSObject.Properties.Name -contains 'fmm_local_cache') {
            [string]$saved.fmm_local_cache
        } else { 'green-jet-compressed' }
        $savedBudget = if ($saved.PSObject.Properties.Name -contains 'fmm_local_cache_budget_mb') {
            [int]$saved.fmm_local_cache_budget_mb
        } else { 256 }
        $protocolMatches = $savedPrimary -eq $PrimarySurface -and
            [Math]::Abs($savedSpacing - $SesGridSpacing) -le 1.0e-12 -and
            $savedCache -eq $FmmLocalCache -and
            $savedBudget -eq $FmmLocalCacheBudgetMb
        if ($protocolMatches -and ($savedOk -or -not $RerunFailed)) {
            Write-Host "resume: $fixture"
            $rows.Add($saved)
            continue
        }
    }

    Write-Host "run: $fixture guarded-p4 deltaG_solv at 0.15 M (timeout ${CaseTimeoutSeconds}s)"
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath $pqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 `
        -Restart 128 -Tolerance 1.0e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmLocal gmres-p4 -FmmLocalRefine operator-tail-static `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalPlanBuilder metric-cache `
        -FmmLocalCache $FmmLocalCache `
        -FmmLocalCacheBudgetMb $FmmLocalCacheBudgetMb `
        -FmmLocalFallback packed-avx2 `
        -FmmLocalGeometryGuard strict-geometry-mac `
        -FmmLocalGeometryTheta 0.5 -DgSolv $true `
        -PrimarySurface $PrimarySurface -SesGridSpacing $SesGridSpacing
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds

    if ($run.exit_code -ne 0) {
        $failure = if ($run.timed_out) { 'timeout' } else { "exit_$($run.exit_code)" }
        $row = New-FailureRow -Fixture $fixture -Failure $failure -Run $run
        $row.run_dir = $runDir
    } else {
        try {
            $energy = @(Import-Csv (Join-Path $runDir 'pbe_energy.csv'))[0]
            $profile = @(Import-Csv (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
            $dgKcal = Read-Double $energy.deltaG_solv_kcalmol
            $dgKj = Read-Double $energy.deltaG_solv_kJmol
            $rxnKcal = Read-Double $energy.E_rxn_kcalmol
            $residual = Read-Double $profile.true_residual_rel
            $audit = Read-Double $profile.tree_vs_direct_matvec_rel
            $directResidual = Read-Double $profile.fmm_local_sampled_direct_true_residual
            $statusOk = [string]$energy.status -eq 'ok' -and
                [string]$profile.status -eq 'ok' -and
                [string]$energy.dgsolv_status -eq 'ok'
            $accuracyOk = (Test-Finite $dgKcal) -and (Test-Finite $residual) -and
                $residual -le 1.0e-6 -and (Test-Finite $audit) -and
                $audit -le 8.0e-3 -and (Test-Finite $directResidual) -and
                $directResidual -le 8.0e-3
            $closureOk = [Math]::Abs($dgKcal - $rxnKcal) -le
                1.0e-10 * [Math]::Max(1.0, [Math]::Abs($rxnKcal))
            $row = [pscustomobject][ordered]@{
                fixture = $fixture; pdb_id = $fixture.ToUpperInvariant()
                atoms = [int64]$energy.atom_count; panels = [int64]$profile.panels
                status = if ($statusOk -and $accuracyOk -and $closureOk) { 'ok' } else { 'failed' }
                failure_class = if (-not $statusOk) { 'solver_status' } elseif (-not $closureOk) { 'energy_closure' } elseif (-not $accuracyOk) { 'accuracy_gate' } else { 'pass' }
                deltaG_solv_kcalmol = $dgKcal; deltaG_solv_kJmol = $dgKj
                reaction_energy_kcalmol = $rxnKcal
                iterations = [int]$profile.iterations; true_residual = $residual
                sampled_operator_audit_rms = $audit
                sampled_direct_true_residual = $directResidual
                wall_seconds = $run.wall_seconds; peak_rss_mb = $run.peak_working_set_mb
                preconditioner = [string]$profile.preconditioner_effective
                convergence_reason = [string]$profile.convergence_reason
                dgsolv_status = [string]$energy.dgsolv_status
                salt_M = 0.15; temperature_K = 298.15
                plan_build_count = [int]$profile.fmm_local_plan_build_count
                cache_build_count = [int]$profile.fmm_local_cache_build_count
                primary_surface = $PrimarySurface
                ses_grid_spacing_A = $SesGridSpacing
                fmm_local_cache = $FmmLocalCache
                fmm_local_cache_budget_mb = $FmmLocalCacheBudgetMb
                run_dir = $runDir
            }
        } catch {
            $row = New-FailureRow -Fixture $fixture -Failure "schema: $($_.Exception.Message)" -Run $run
            $row.run_dir = $runDir
        }
    }
    @($row) | Export-Csv -LiteralPath $rowPath -NoTypeInformation
    $rows.Add($row)
}

$matrixPath = Join-Path $WorkRoot 'fagbi_tabi_pb2_solvation_energy_matrix.csv'
$rowArray = @($rows | ForEach-Object { $_ })
$rowArray | Export-Csv -LiteralPath $matrixPath -NoTypeInformation
$passed = @($rowArray | Where-Object { [string]$_.status -eq 'ok' }).Count
$decision = if ($passed -eq $rows.Count) {
    'tabi_pb2_solvation_energy_pass'
} elseif ($passed -gt 0) {
    'tabi_pb2_solvation_energy_partial'
} else {
    'tabi_pb2_solvation_energy_failed'
}
$decisionPath = Join-Path $WorkRoot 'fagbi_tabi_pb2_solvation_energy_decision.csv'
@([pscustomobject]@{ decision = $decision; passed = $passed; required = $rows.Count }) |
    Export-Csv -LiteralPath $decisionPath -NoTypeInformation
$summaryPath = Join-Path $WorkRoot 'fagbi_tabi_pb2_solvation_energy_summary.dat'
@(
    "decision $decision"
    "salt_M 0.15"
    "temperature_K 298.15"
    "primary_surface $PrimarySurface"
    "ses_grid_spacing_A $SesGridSpacing"
    "fmm_local_cache $FmmLocalCache"
    "fmm_local_cache_budget_mb $FmmLocalCacheBudgetMb"
    "fixtures $($Fixtures -join ',')"
    "passed $passed"
    "required $($rows.Count)"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

Write-Host "Decision: $decision"
Write-Host "Matrix: $matrixPath"
Write-Host 'Resume: rerun unchanged; add -RerunFailed for failed/timeout rows or -Force for all rows.'
if (-not $AllowFailedDecision -and $decision -ne 'tabi_pb2_solvation_energy_pass') {
    throw "Solvation-energy gate did not pass: $decision"
}
