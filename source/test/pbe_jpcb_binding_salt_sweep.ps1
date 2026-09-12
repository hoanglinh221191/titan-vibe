param(
    [string]$TitanExe = '.\titan.exe',
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_binding_jpcb'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\jpcb_binding_salt'),
    [string]$Salts = '0.05,0.15,0.30,0.50',
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [int]$SurfacePoints = 8,
    [ValidateSet('guarded-p4','degree1-grid-ses')]
    [string]$SolverPolicy = 'guarded-p4',
    [ValidateRange(0.5,2.5)][double]$SesGridSpacing = 1.25,
    [ValidateSet('green-jet-compressed','green-jet-tiled')]
    # The binding protocol performs hundreds of p4 matvecs per component.
    # Rebuilding tiled jets in every accepted far pair is diagnostic-only and
    # dominates wall time; use the reusable lossless compressed cache here.
    [string]$CacheMode = 'green-jet-compressed',
    [ValidateRange(64,8192)][int]$CacheBudgetMb = 4096,
    [switch]$Quick,
    [switch]$ContractOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$culture = [Globalization.CultureInfo]::InvariantCulture
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$manifestPath = Join-Path $FixtureRoot 'binding_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "jpcb_binding_fixture_gap: run test\pbe_jpcb_prepare_barnase_barstar.ps1 -RefreshFixtures"
}
$fixtureRows = @(Import-Csv -LiteralPath $manifestPath)
if ($fixtureRows.Count -ne 1) { throw 'Expected exactly one JPCB binding fixture.' }
$fixture = $fixtureRows[0]
$protocolId = if ($SolverPolicy -eq 'degree1-grid-ses') {
    "jpcb-binding-v2-degree1-grid-ses-h$($SesGridSpacing.ToString('0.00',$culture).Replace('.','p'))-s$SurfacePoints-nt$Threads"
} else {
    "jpcb-binding-v1-p4-$CacheMode-b$CacheBudgetMb-s$SurfacePoints-nt$Threads"
}
$expectedCacheBuildCount = if ($SolverPolicy -eq 'degree1-grid-ses') {
    0
} elseif ($CacheMode -eq 'green-jet-tiled') {
    0
} else {
    1
}
$targetTolerance = if ($SolverPolicy -eq 'degree1-grid-ses') { 1.0e-4 } else { 1.0e-6 }
$complexPqr = Join-Path $FixtureRoot $fixture.complex_pqr
$receptorPqr = Join-Path $FixtureRoot $fixture.receptor_pqr
$ligandPqr = Join-Path $FixtureRoot $fixture.ligand_pqr
foreach ($path in @($complexPqr,$receptorPqr,$ligandPqr)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing cached PQR: $path" }
}
if ((Get-FileHash -Algorithm SHA256 -LiteralPath $complexPqr).Hash.ToLowerInvariant() -ne
        $fixture.complex_pqr_sha256 -or
    (Get-FileHash -Algorithm SHA256 -LiteralPath $receptorPqr).Hash.ToLowerInvariant() -ne
        $fixture.receptor_pqr_sha256 -or
    (Get-FileHash -Algorithm SHA256 -LiteralPath $ligandPqr).Hash.ToLowerInvariant() -ne
        $fixture.ligand_pqr_sha256) {
    throw 'jpcb_binding_fixture_gap: cached PQR checksum drift.'
}

function Convert-JpcbSaltList([string]$Text) {
    $values = New-Object Collections.Generic.List[double]
    foreach ($token in ($Text -split ',')) {
        $value = 0.0
        if (-not [double]::TryParse($token.Trim(), [Globalization.NumberStyles]::Float,
                $culture, [ref]$value) -or $value -le 0.0 -or $value -gt 1.0) {
            throw "Invalid ionic strength for guarded FMM/ILU0 binding runs: $token. Use a value in (0, 1]."
        }
        $values.Add($value) | Out-Null
    }
    return @($values.ToArray() | Sort-Object -Unique)
}

function Test-JpcbFinite([double]$Value) {
    return -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}

function Get-JpcbCaseId([double]$Salt) {
    return 'salt_' + $Salt.ToString('0.00', $culture).Replace('.','p')
}

function Write-JpcbCaseArgs([string]$RunDir, [double]$Salt) {
    if ($SolverPolicy -eq 'degree1-grid-ses') {
        return Write-FmmCliArgFile `
            -RunDir $RunDir -RepoRoot $RepoRoot -PqrPath $complexPqr `
            -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
            -SurfacePoints $SurfacePoints -MaxIter 720 -Restart 128 `
            -Tolerance $targetTolerance -Threads $Threads -SaltM $Salt `
            -InternalDielectric 1.0 -ExternalDielectric 80.0 `
            -TemperatureK 298.15 -OpeningTheta 0.45 `
            -AuditSampleCount 512 -AuditMode sampled `
            -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
            -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
            -TaylorOrder 4 -FmmLocal off -FmmLocalRefine off `
            -FmmLocalCache off -FmmLocalFallback scalar `
            -FmmLocalPlanBuilder legacy-rescan `
            -PrimarySurface grid-ses -SesGridSpacing $SesGridSpacing `
            -AtomPotentials csv -BindingCycle pqr-components `
            -BindingReceptorPqr $receptorPqr -BindingLigandPqr $ligandPqr
    }
    return Write-FmmCliArgFile `
        -RunDir $RunDir -RepoRoot $RepoRoot -PqrPath $complexPqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints $SurfacePoints -MaxIter 720 -Restart 128 `
        -Tolerance 1.0e-6 -Threads $Threads -SaltM $Salt `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -OpeningTheta 0.6 `
        -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmLocal gmres-p4 `
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $CacheMode -FmmLocalCacheBudgetMb $CacheBudgetMb `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard strict-geometry-mac -FmmLocalGeometryTheta 0.5 `
        -AtomPotentials csv -BindingCycle pqr-components `
        -BindingReceptorPqr $receptorPqr -BindingLigandPqr $ligandPqr
}

function Test-JpcbComponentGate([object]$Row) {
    $internal = Convert-FmmCliDouble $Row.internal_residual
    $trueResidual = Convert-FmmCliDouble $Row.true_residual
    $audit = Convert-FmmCliDouble $Row.sampled_direct_operator_rms
    $directResidual = Convert-FmmCliDouble $Row.sampled_direct_true_residual
    $lifecycleGate = if ($SolverPolicy -eq 'degree1-grid-ses') {
        $Row.fmm_local_gmres_active -eq 'false' -and
            [int64]$Row.plan_build_count -eq 0 -and
            [int64]$Row.cache_build_count -eq 0
    } else {
        $Row.fmm_local_gmres_active -eq 'true' -and
            [int64]$Row.plan_build_count -eq 1 -and
            [int64]$Row.cache_build_count -eq $expectedCacheBuildCount
    }
    return $Row.status -eq 'ok' -and
        (Test-JpcbFinite $internal) -and $internal -le $targetTolerance -and
        (Test-JpcbFinite $trueResidual) -and $trueResidual -le $targetTolerance -and
        (Test-JpcbFinite $audit) -and $audit -le 8.0e-3 -and
        (Test-JpcbFinite $directResidual) -and $directResidual -le 8.0e-3 -and
        $lifecycleGate
}

function Write-JpcbCaseFailure([string]$RunDir, [string]$CaseId,
                               [double]$Salt, [string]$Status,
                               [object]$Run) {
    [pscustomobject]@{
        protocol_id=$protocolId; case_id=$CaseId; salt_M=$Salt
        status=$Status; exit_code=$Run.exit_code; timed_out=$Run.timed_out
        wall_seconds=$Run.wall_seconds; process_peak_rss_mb=$Run.peak_working_set_mb
    } | Export-Csv -LiteralPath (Join-Path $RunDir 'case_result.csv') -NoTypeInformation
}

function Invoke-JpcbSaltCase([double]$Salt) {
    $caseId = Get-JpcbCaseId $Salt
    $runDir = Join-Path $WorkRoot $caseId
    $casePath = Join-Path $runDir 'case_result.csv'
    if (-not $Force -and (Test-Path -LiteralPath $casePath)) {
        $prior = @(Import-Csv -LiteralPath $casePath)[0]
        $rerun = $RerunFailed -and $prior.status -ne 'ok'
        if (-not $rerun -and $prior.protocol_id -eq $protocolId) {
            Write-Host "resume: $caseId"
            return
        }
    }
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $argPath = Write-JpcbCaseArgs -RunDir $runDir -Salt $Salt
    Write-Host "run: $caseId (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        $status = if ($run.timed_out) { 'timeout' } else { 'failed' }
        Write-JpcbCaseFailure -RunDir $runDir -CaseId $caseId -Salt $Salt `
            -Status $status -Run $run
        return
    }

    $componentPath = Join-Path $runDir 'pbe_component_electrostatics.csv'
    $bindingPath = Join-Path $runDir 'pbe_binding_electrostatics.csv'
    $atomPath = Join-Path $runDir 'pbe_atom_reaction_potential.csv'
    $residuePath = Join-Path $runDir 'pbe_residue_electrostatics.csv'
    foreach ($path in @($componentPath,$bindingPath,$atomPath,$residuePath)) {
        if (-not (Test-Path -LiteralPath $path)) {
            Write-JpcbCaseFailure -RunDir $runDir -CaseId $caseId -Salt $Salt `
                -Status 'missing_output' -Run $run
            return
        }
    }
    $components = @(Import-Csv -LiteralPath $componentPath)
    $bindings = @(Import-Csv -LiteralPath $bindingPath)
    if ($components.Count -ne 3 -or $bindings.Count -ne 1) {
        Write-JpcbCaseFailure -RunDir $runDir -CaseId $caseId -Salt $Salt `
            -Status 'schema_gap' -Run $run
        return
    }
    $binding = $bindings[0]
    $componentGate = @($components | Where-Object { Test-JpcbComponentGate $_ }).Count -eq 3
    $finiteBinding = @('Gpol_complex_kJmol','Gpol_receptor_kJmol',
        'Gpol_ligand_kJmol','Delta_Gpol_bind_kJmol',
        'E_coul_receptor_ligand_kJmol','Delta_Gelec_bind_kJmol') | ForEach-Object {
            Test-JpcbFinite (Convert-FmmCliDouble $binding.$_)
        }
    $bindingGate = $binding.status -eq 'ok' -and
        @($finiteBinding | Where-Object { -not $_ }).Count -eq 0 -and
        (Convert-FmmCliDouble $binding.closure_error_max_kJmol) -le 1.0e-8 -and
        (Convert-FmmCliDouble $binding.sampled_audit_max) -le 8.0e-3 -and
        (Convert-FmmCliDouble $binding.sampled_direct_true_residual_max) -le 8.0e-3
    $status = if ($componentGate -and $bindingGate) { 'ok' } else { 'gate_failed' }

    $componentRows = foreach ($row in $components) {
        [pscustomobject]@{
            case_id=$caseId; salt_M=$Salt; component=$row.component
            atom_count=[int64]$row.atom_count; pqr_sha256=$row.pqr_sha256
            energy_kJmol=Convert-FmmCliDouble $row.reaction_energy_kJmol
            iterations=[int64]$row.iterations
            internal_residual=Convert-FmmCliDouble $row.internal_residual
            true_residual=Convert-FmmCliDouble $row.true_residual
            sampled_audit_rms=Convert-FmmCliDouble $row.sampled_direct_operator_rms
            sampled_audit_rows=[int64]$row.sampled_operator_rows
            sampled_direct_true_residual=Convert-FmmCliDouble $row.sampled_direct_true_residual
            plan_build_count=[int64]$row.plan_build_count
            cache_build_count=[int64]$row.cache_build_count
            prepare_seconds=Convert-FmmCliDouble $row.fmm_prepare_seconds
            p4_matvec_calls=[int64]$row.fmm_local_matvec_calls
            p4_apply_seconds=Convert-FmmCliDouble $row.fmm_local_apply_seconds
            p4_moment_seconds=Convert-FmmCliDouble $row.fmm_local_moment_seconds
            p4_m2l_seconds=Convert-FmmCliDouble $row.fmm_local_m2l_seconds
            p4_fallback_seconds=Convert-FmmCliDouble $row.fmm_local_fallback_seconds
            p4_l2p_seconds=Convert-FmmCliDouble $row.fmm_local_l2p_seconds
            gmres_seconds=Convert-FmmCliDouble $row.gmres_seconds
            gmres_orthogonalization_seconds=Convert-FmmCliDouble $row.gmres_orthogonalization_seconds
            gmres_preconditioner_seconds=Convert-FmmCliDouble $row.gmres_preconditioner_seconds
            matvec_seconds_mean=Convert-FmmCliDouble $row.matvec_seconds_mean
            peak_memory_mb=(Convert-FmmCliDouble $row.peak_memory_bytes) / 1MB
            status=$row.status; gate_pass=(Test-JpcbComponentGate $row)
        }
    }
    $componentRows | Export-Csv -LiteralPath `
        (Join-Path $runDir 'case_component_result.csv') -NoTypeInformation

    [pscustomobject]@{
        case_id=$caseId; salt_M=$Salt
        delta_g_pol_bind_kJmol=Convert-FmmCliDouble $binding.Delta_Gpol_bind_kJmol
        cross_coulomb_kJmol=Convert-FmmCliDouble $binding.E_coul_receptor_ligand_kJmol
        delta_g_elec_bind_kJmol=Convert-FmmCliDouble $binding.Delta_Gelec_bind_kJmol
        sampled_audit_max=Convert-FmmCliDouble $binding.sampled_audit_max
        sampled_direct_true_residual_max=Convert-FmmCliDouble $binding.sampled_direct_true_residual_max
        closure_error_max_kJmol=Convert-FmmCliDouble $binding.closure_error_max_kJmol
        wall_seconds=$run.wall_seconds
        process_peak_rss_mb=$run.peak_working_set_mb
        status=$status
    } | Export-Csv -LiteralPath (Join-Path $runDir 'case_binding_result.csv') -NoTypeInformation

    $atomRows = @(Import-Csv -LiteralPath $atomPath)
    $residueRows = @(Import-Csv -LiteralPath $residuePath)
    $atomSummary = foreach ($component in @('complex','receptor','ligand')) {
        $selectedAtoms = @($atomRows | Where-Object component -eq $component)
        $selectedResidues = @($residueRows | Where-Object component -eq $component)
        $sum = ($selectedAtoms | ForEach-Object {
            Convert-FmmCliDouble $_.reaction_energy_contribution_kJmol
        } | Measure-Object -Sum).Sum
        [pscustomobject]@{
            case_id=$caseId; salt_M=$Salt; component=$component
            atom_count=$selectedAtoms.Count; residue_count=$selectedResidues.Count
            atomwise_energy_sum_kJmol=$sum
            pqr_sha256=@($selectedAtoms)[0].pqr_sha256
        }
    }
    $atomSummary | Export-Csv -LiteralPath `
        (Join-Path $runDir 'case_atom_summary.csv') -NoTypeInformation
    [pscustomobject]@{
        protocol_id=$protocolId; case_id=$caseId; salt_M=$Salt
        status=$status; exit_code=$run.exit_code; timed_out=$run.timed_out
        wall_seconds=$run.wall_seconds; process_peak_rss_mb=$run.peak_working_set_mb
    } | Export-Csv -LiteralPath $casePath -NoTypeInformation
}

function Write-JpcbPlot([object[]]$Rows, [string]$Path) {
    $valid = @($Rows | Where-Object status -eq 'ok' | Sort-Object { [double]$_.salt_M })
    if ($valid.Count -eq 0) { return }
    $width=820; $height=500; $left=90; $right=30; $top=50; $bottom=75
    $plotW=$width-$left-$right; $plotH=$height-$top-$bottom
    $allY = @($valid | ForEach-Object {
        [double]$_.delta_g_pol_bind_kJmol; [double]$_.delta_g_elec_bind_kJmol
    })
    $yMin=($allY | Measure-Object -Minimum).Minimum
    $yMax=($allY | Measure-Object -Maximum).Maximum
    if ($yMax -le $yMin) { $yMax=$yMin+1.0 }
    $xMax=[Math]::Max(0.5,($valid | ForEach-Object {[double]$_.salt_M} |
        Measure-Object -Maximum).Maximum)
    function X([double]$v){ return $left+$plotW*$v/$xMax }
    function Y([double]$v){ return $top+$plotH*(1.0-($v-$yMin)/($yMax-$yMin)) }
    $polPoints=($valid | ForEach-Object { "$(X ([double]$_.salt_M)),$(Y ([double]$_.delta_g_pol_bind_kJmol))" }) -join ' '
    $elecPoints=($valid | ForEach-Object { "$(X ([double]$_.salt_M)),$(Y ([double]$_.delta_g_elec_bind_kJmol))" }) -join ' '
    $svg=@(
        "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<line x1='$left' y1='$top' x2='$left' y2='$($top+$plotH)' stroke='black' stroke-width='2'/>",
        "<line x1='$left' y1='$($top+$plotH)' x2='$($left+$plotW)' y2='$($top+$plotH)' stroke='black' stroke-width='2'/>",
        "<polyline points='$polPoints' fill='none' stroke='#0072B2' stroke-width='3'/>",
        "<polyline points='$elecPoints' fill='none' stroke='#D55E00' stroke-width='3'/>",
        "<text x='$($width/2)' y='$($height-20)' text-anchor='middle' font-family='Arial' font-size='20'>Ionic strength (M)</text>",
        "<text x='25' y='$($height/2)' text-anchor='middle' font-family='Arial' font-size='20' transform='rotate(-90 25 $($height/2))'>Binding contribution (kJ mol-1)</text>",
        "<line x1='470' y1='28' x2='505' y2='28' stroke='#0072B2' stroke-width='3'/><text x='515' y='34' font-family='Arial' font-size='17'>Polar</text>",
        "<line x1='610' y1='28' x2='645' y2='28' stroke='#D55E00' stroke-width='3'/><text x='655' y='34' font-family='Arial' font-size='17'>Total electrostatic</text>",
        '</svg>'
    )
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$saltValues = @(Convert-JpcbSaltList $Salts)
if ($Quick) { $saltValues = @(0.15) }
if ($ContractOnly) {
    $runDir = Join-Path $WorkRoot 'contract_only'
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $argPath = Write-JpcbCaseArgs -RunDir $runDir -Salt 0.15
    $bad = @(Get-Content -LiteralPath $argPath | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($bad.Count -ne 0) { throw "Contract arg.in has bare options: $($bad -join '; ')" }
    Write-Host "PASS JPCB binding salt contract: $argPath"
    exit 0
}
foreach ($salt in $saltValues) { Invoke-JpcbSaltCase -Salt $salt }

$requestedCaseIds = @($saltValues | ForEach-Object { Get-JpcbCaseId $_ })
$caseDirs = @(Get-ChildItem -LiteralPath $WorkRoot -Directory |
    Where-Object { $_.Name -in $requestedCaseIds } | Sort-Object Name)
$caseRows = @($caseDirs | ForEach-Object {
    $path=Join-Path $_.FullName 'case_result.csv'
    if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$componentRows = @($caseDirs | ForEach-Object {
    $path=Join-Path $_.FullName 'case_component_result.csv'
    if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$bindingRows = @($caseDirs | ForEach-Object {
    $path=Join-Path $_.FullName 'case_binding_result.csv'
    if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$atomRows = @($caseDirs | ForEach-Object {
    $path=Join-Path $_.FullName 'case_atom_summary.csv'
    if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$matrixPath=Join-Path $WorkRoot 'jpcb_binding_matrix.csv'
$thermoPath=Join-Path $WorkRoot 'jpcb_binding_thermocycle.csv'
$atomSummaryPath=Join-Path $WorkRoot 'jpcb_binding_atomwise_summary.csv'
$decisionPath=Join-Path $WorkRoot 'jpcb_binding_decision.csv'
$summaryPath=Join-Path $WorkRoot 'jpcb_binding_summary.dat'
$plotPath=Join-Path $WorkRoot 'jpcb_binding_salt_plot.svg'
$componentRows | Export-Csv -LiteralPath $matrixPath -NoTypeInformation
$bindingRows | Export-Csv -LiteralPath $thermoPath -NoTypeInformation
$atomRows | Export-Csv -LiteralPath $atomSummaryPath -NoTypeInformation
Write-JpcbPlot -Rows $bindingRows -Path $plotPath

$requiredCount=$saltValues.Count
$okCount=@($caseRows | Where-Object status -eq 'ok').Count
$decision = if (@($caseRows | Where-Object status -eq 'timeout').Count -gt 0) {
    'jpcb_binding_salt_runtime_limited'
} elseif ($okCount -ne $requiredCount) {
    'jpcb_binding_salt_accuracy_gap'
} elseif ($Quick) {
    'jpcb_binding_salt_quick_pass'
} else {
    'jpcb_binding_salt_pass'
}
[pscustomobject]@{ decision=$decision; required_cases=$requiredCount; ok_cases=$okCount;
    fixture=$fixture.fixture; protocol_id=$protocolId; cache_mode=$CacheMode;
    cache_budget_mb=$CacheBudgetMb; generated_utc=[DateTime]::UtcNow.ToString('o') } |
    Export-Csv -LiteralPath $decisionPath -NoTypeInformation
Set-Content -LiteralPath $summaryPath -Encoding ASCII -Value @(
    "decision=$decision","fixture=$($fixture.fixture)","required_cases=$requiredCount",
    "ok_cases=$okCount","matrix=$matrixPath","thermocycle=$thermoPath",
    "atomwise_summary=$atomSummaryPath","plot=$plotPath")
$reportPath=Join-Path $RepoRoot 'TITAN_JPCB_BINDING_SALT_REPORT.txt'
Set-Content -LiteralPath $reportPath -Encoding ASCII -Value @(
    'TITAN JPCB Barnase-Barstar Salt-Binding Gate','',"Decision: $decision",
    "Fixture: $($fixture.biological_name) ($($fixture.pdb_id), assembly $($fixture.assembly_id))",
    "Protocol: guarded local p4, strict geometry 0.5, near-field ILU0 cap 64, cache layout $CacheMode ($CacheBudgetMb MB)",
    "Cache lifecycle: expected global cache builds per solve = $expectedCacheBuildCount (tiled mode rebuilds per-pair scratch jets and reports tile work separately)",
    "Cases: $okCount / $requiredCount passed","Matrix: $matrixPath",
    "Thermodynamic cycle: $thermoPath","Plot: $plotPath")
Write-Host "JPCB binding decision: $decision"
Write-Host "Thermodynamic cycle: $thermoPath"
if (-not $AllowFailedDecision -and $decision -notin @('jpcb_binding_salt_pass','jpcb_binding_salt_quick_pass')) {
    throw "JPCB binding salt gate did not pass: $decision"
}
