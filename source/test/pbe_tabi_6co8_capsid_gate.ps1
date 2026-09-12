[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$PqrPath = 'test\pbe_inputs_capsid\6co8_assembly1_charmm.pqr',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_6co8_capsid',
    [string[]]$Stages = @('surface-only','coarse-solve'),
    [ValidateScript({ $_ -ge 0.25 -and $_ -le 4.0 })]
    [double]$SesGridSpacing = 2.0,
    [ValidateRange(4000000,100000000)][int64]$SesMaxGridVertices = 32000000,
    [ValidateRange(0.05,2.0)][double]$OpeningTheta = 0.45,
    [ValidateRange(1,64)][int]$Threads = 4,
    [ValidateRange(1,20000)][int]$MaxIter = 720,
    [ValidateRange(1,1000)][int]$Restart = 10,
    [ValidateRange(1.0e-8,1.0e-2)][double]$Tolerance = 1.0e-4,
    [ValidateRange(32,1024)][int]$IluNeighborCap = 32,
    [ValidateRange(1,10000000)][int64]$MinimumAtomCount = 1000000,
    [ValidateRange(1,20000000)][int64]$MaximumCoarsePanels = 6000000,
    [ValidateRange(1024,65536)][int]$MemoryBudgetMb = 24576,
    [ValidateRange(512,16384)][int]$ProjectedSolverBytesPerPanel = 4096,
    [ValidateScript({ $_ -ge 1000000 })]
    [int64]$MaximumDirectRhsPairs = 10000000000,
    [ValidateRange(1,86400)][int]$CaseTimeoutSeconds = 7200,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function Split-Values([string[]]$Values) {
    return @($Values | ForEach-Object { ([string]$_).Split(',') } |
        ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
}
function F([double]$Value) { $Value.ToString('G15', $Inv) }
function D($Value, [double]$Fallback = [double]::NaN) {
    $parsed = 0.0
    if ($null -ne $Value -and [double]::TryParse(
            [string]$Value, [Globalization.NumberStyles]::Float,
            $Inv, [ref]$parsed)) { return $parsed }
    return $Fallback
}
function I($Value, [int64]$Fallback = 0) {
    $parsed = 0L
    if ($null -ne $Value -and [int64]::TryParse([string]$Value,
            [ref]$parsed)) { return $parsed }
    return $Fallback
}
function Count-PqrAtoms([string]$Path) {
    $count = 0L
    $reader = [IO.File]::OpenText($Path)
    try {
        while ($null -ne ($line = $reader.ReadLine())) {
            if ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM')) {
                ++$count
            }
        }
    } finally { $reader.Dispose() }
    return $count
}
function Read-LastCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    return @(Import-Csv -LiteralPath $Path) | Select-Object -Last 1
}
function Ensure-EqualsSyntax([string]$ArgPath) {
    $bad = @((Get-Content -LiteralPath $ArgPath) | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($bad.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
    }
}

$stageSet = @(Split-Values $Stages)
foreach ($stage in $stageSet) {
    if ($stage -notin @('surface-only','coarse-solve')) {
        throw "Unknown stage '$stage'; expected surface-only or coarse-solve."
    }
}
$TitanExe = [IO.Path]::GetFullPath((Join-Path $RepoRoot $TitanExe))
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Titan executable not found: $TitanExe"
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

if ($Quick) {
    $MinimumAtomCount = 1
    $PqrPath = Join-Path $WorkRoot 'quick_capsid_contract.pqr'
    @(
        'ATOM      1  C1  MOL A   1      -2.000   0.000   0.000   0.5000 1.7000',
        'ATOM      2  C2  MOL A   1       2.000   0.000   0.000  -0.5000 1.7000',
        'END'
    ) | Set-Content -LiteralPath $PqrPath -Encoding ASCII
} else {
    if (-not [IO.Path]::IsPathRooted($PqrPath)) {
        $PqrPath = Join-Path $RepoRoot $PqrPath
    }
    if (-not (Test-Path -LiteralPath $PqrPath)) {
        throw "capsid_fixture_gap: missing full 6CO8 PQR: $PqrPath. Supply -PqrPath to a prepared biological-assembly PQR; the asymmetric unit is not accepted."
    }
    $PqrPath = (Resolve-Path -LiteralPath $PqrPath).Path
}

$atomCount = Count-PqrAtoms $PqrPath
if ($atomCount -lt $MinimumAtomCount) {
    throw "capsid_fixture_gap: PQR has $atomCount atoms, below required $MinimumAtomCount; use the expanded 6CO8 biological assembly, not the asymmetric unit."
}
$pqrSha = (Get-FileHash -LiteralPath $PqrPath -Algorithm SHA256).Hash.ToLowerInvariant()
$titanSha = (Get-FileHash -LiteralPath $TitanExe -Algorithm SHA256).Hash.ToLowerInvariant()
$matrixPath = Join-Path $WorkRoot 'tabi_6co8_capsid_matrix.csv'
$feasibilityPath = Join-Path $WorkRoot 'tabi_6co8_surface_feasibility.csv'
$decisionPath = Join-Path $WorkRoot 'tabi_6co8_capsid_decision.csv'
$summaryPath = Join-Path $WorkRoot 'tabi_6co8_capsid_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_TABI_6CO8_CAPSID_FEASIBILITY_REPORT.txt'

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
    foreach ($row in @(Import-Csv -LiteralPath $matrixPath)) {
        $existing[[string]$row.row_id] = $row
    }
}
$rows = [Collections.Generic.List[object]]::new()
$surfaceRow = $null
foreach ($stage in @('surface-only','coarse-solve')) {
    if ($stage -notin $stageSet) { continue }
    $rowId = "6co8_$($stage.Replace('-','_'))_h$(([string](F $SesGridSpacing)).Replace('.','p'))"
    $old = $existing[$rowId]
    $sameBinary = $null -ne $old -and
        $old.PSObject.Properties.Name -contains 'titan_sha256' -and
        [string]$old.titan_sha256 -eq $titanSha -and
        [string]$old.pqr_sha256 -eq $pqrSha
    $reuse = -not $Force -and $sameBinary -and
        (-not $RerunFailed -or [string]$old.status -eq 'ok')
    if ($reuse) {
        Write-Host "resume: $rowId"
        if ($old.PSObject.Properties.Name -notcontains 'estimated_rhs_pairs') {
            $old | Add-Member -NotePropertyName estimated_rhs_pairs `
                -NotePropertyValue ([int64](I $old.panels) * [int64]$atomCount)
        }
        $rows.Add($old)
        if ($stage -eq 'surface-only') { $surfaceRow = $old }
        continue
    }

    if ($stage -eq 'coarse-solve') {
        if ($null -eq $surfaceRow) {
            $surfaceRow = @($rows | Where-Object { $_.stage -eq 'surface-only' } |
                Select-Object -Last 1)
        }
        if ($null -eq $surfaceRow -or [string]$surfaceRow.status -ne 'ok') {
            throw 'capsid_surface_probe_required: run the surface-only stage successfully before coarse-solve.'
        }
        $surfacePanels = I $surfaceRow.panels
        $estimatedRhsPairs = [int64]$surfacePanels * [int64]$atomCount
        $surfacePeakMb = D $surfaceRow.peak_rss_mb 0.0
        $projectedMb = $surfacePeakMb +
            ($surfacePanels * [double]$ProjectedSolverBytesPerPanel / 1MB)
        if ($surfacePanels -gt $MaximumCoarsePanels -or
            $projectedMb -gt $MemoryBudgetMb) {
            Write-Host "skip: $rowId (memory gate: panels=$surfacePanels projected=$([Math]::Round($projectedMb,1)) MB)"
            $skip = [pscustomobject]@{
                row_id=$rowId; stage=$stage; atom_count=$atomCount
                panels=$surfacePanels; status='skipped_memory_gate'; exit_code=0
                timed_out=$false; wall_seconds=0; peak_rss_mb=$surfacePeakMb
                projected_solve_mb=$projectedMb; surface_build_seconds=D $surfaceRow.surface_build_seconds 0
                estimated_rhs_pairs=$estimatedRhsPairs
                iterations=0; final_residual='NA'; energy_kcal_mol='NA'
                matvec_calls=0; fmm_mode='off'; preconditioner_effective='not_built'
                convergence_reason='surface_memory_gate'; memory_budget_mb=$MemoryBudgetMb
                titan_sha256=$titanSha; pqr_sha256=$pqrSha
            }
            $rows.Add($skip)
            continue
        }
        if ($estimatedRhsPairs -gt $MaximumDirectRhsPairs) {
            Write-Host "skip: $rowId (direct RHS gate: pairs=$estimatedRhsPairs limit=$MaximumDirectRhsPairs)"
            $skip = [pscustomobject]@{
                row_id=$rowId; stage=$stage; atom_count=$atomCount
                panels=$surfacePanels; status='skipped_rhs_cost_gate'; exit_code=0
                timed_out=$false; wall_seconds=0; peak_rss_mb=$surfacePeakMb
                projected_solve_mb=$projectedMb; surface_build_seconds=D $surfaceRow.surface_build_seconds 0
                estimated_rhs_pairs=$estimatedRhsPairs
                iterations=0; final_residual='NA'; energy_kcal_mol='NA'
                matvec_calls=0; fmm_mode='off'; preconditioner_effective='not_built'
                convergence_reason='direct_rhs_acceleration_required'
                memory_budget_mb=$MemoryBudgetMb; titan_sha256=$titanSha; pqr_sha256=$pqrSha
            }
            $rows.Add($skip)
            continue
        }
    }

    $runDir = Join-Path $WorkRoot $rowId
    if (Test-Path -LiteralPath $runDir) {
        Remove-Item -LiteralPath $runDir -Recurse -Force
    }
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $argPath = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
        -PqrPath $PqrPath -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -Threads $Threads -MaxIter $MaxIter -Restart $Restart `
        -Tolerance $Tolerance -SaltM 0.15 -InternalDielectric 4.0 `
        -ExternalDielectric 80.0 -TemperatureK 298.15 `
        -ProbeRadius 1.4 -OpeningTheta $OpeningTheta `
        -AuditMode off -AuditSampleCount 0 `
        -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap $IluNeighborCap `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -PrimarySurface grid-ses `
        -SesGridSpacing $SesGridSpacing -DgSolv $true
    $argLines = @(Get-Content -LiteralPath $argPath)
    $argLines = @($argLines | ForEach-Object {
        if ($_ -match '^-pbe-tabi-true-residual-probe=') {
            '-pbe-tabi-true-residual-probe=off'
        } else { $_ }
    })
    if ($stage -eq 'surface-only') {
        $argLines += '-pbe-tabi-surface-only=true'
    }
    $argLines += "-pbe-tabi-ses-max-grid-vertices=$SesMaxGridVertices"
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    Ensure-EqualsSyntax $argPath

    Write-Host "run: $rowId (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) {
        'failed'
    } else { 'ok' }
    $profile = Read-LastCsv (Join-Path $runDir 'pbe_tabi_profile.csv')
    $memory = Read-LastCsv (Join-Path $runDir 'pbe_memory_profile.csv')
    $energy = Read-LastCsv (Join-Path $runDir 'pbe_energy.csv')
    if ($status -eq 'ok' -and $null -eq $profile) { $status = 'profile_missing' }
    if ($status -eq 'ok' -and [string]$profile.status -ne 'ok') {
        $status = [string]$profile.status
    }
    $panels = if ($null -ne $profile) { I $profile.panels } else { 0 }
    $profilePeakMb = if ($null -ne $memory) {
        (D $memory.peak_working_set_bytes 0) / 1MB
    } else { 0.0 }
    $peakMb = [Math]::Max([double]$run.peak_working_set_mb, $profilePeakMb)
    $projectedMb = $peakMb +
        ($panels * [double]$ProjectedSolverBytesPerPanel / 1MB)
    $row = [pscustomobject]@{
        row_id=$rowId; stage=$stage; atom_count=$atomCount; panels=$panels
        status=$status; exit_code=$run.exit_code; timed_out=$run.timed_out
        wall_seconds=$run.wall_seconds; peak_rss_mb=$peakMb
        projected_solve_mb=$projectedMb
        estimated_rhs_pairs=([int64]$panels * [int64]$atomCount)
        surface_build_seconds=$(if($null -ne $profile){D $profile.surface_build_seconds 0}else{0})
        iterations=$(if($null -ne $profile){I $profile.iterations}else{0})
        final_residual=$(if($null -ne $profile){D $profile.final_residual}else{'NA'})
        energy_kcal_mol=$(if($null -ne $energy){D $energy.E_rxn_kcalmol}else{'NA'})
        matvec_calls=$(if($null -ne $profile){I $profile.matvec_calls}else{0})
        fmm_mode=$(if($null -ne $profile){[string]$profile.fmm_m2l_mode}else{'off'})
        preconditioner_effective=$(if($null -ne $memory){[string]$memory.preconditioner_effective}else{'NA'})
        convergence_reason=$(if($null -ne $profile){[string]$profile.convergence_reason}else{$status})
        memory_budget_mb=$MemoryBudgetMb; titan_sha256=$titanSha; pqr_sha256=$pqrSha
    }
    $rows.Add($row)
    if ($stage -eq 'surface-only') { $surfaceRow = $row }
}

$rows.ToArray() | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding UTF8
$surfaceRows = @($rows | Where-Object { $_.stage -eq 'surface-only' })
$surfaceRows | Select-Object atom_count,panels,surface_build_seconds,wall_seconds,
    peak_rss_mb,projected_solve_mb,memory_budget_mb,status,pqr_sha256 |
    Export-Csv -LiteralPath $feasibilityPath -NoTypeInformation -Encoding UTF8
$surface = @($surfaceRows | Select-Object -Last 1)
$coarse = @($rows | Where-Object { $_.stage -eq 'coarse-solve' } |
    Select-Object -Last 1)
$decision = if ($surface.Count -eq 0 -or [string]$surface[0].status -ne 'ok') {
    'capsid_surface_probe_failed'
} elseif ($coarse.Count -eq 0) {
    'capsid_surface_probe_pass'
} elseif ([string]$coarse[0].status -eq 'skipped_memory_gate') {
    'capsid_surface_memory_limited'
} elseif ([string]$coarse[0].status -eq 'skipped_rhs_cost_gate') {
    'capsid_rhs_acceleration_required'
} elseif ([string]$coarse[0].status -eq 'timeout') {
    'capsid_runtime_limited'
} elseif ([string]$coarse[0].status -eq 'ok' -and
    (D $coarse[0].final_residual 1.0) -le $Tolerance -and
    [string]$coarse[0].fmm_mode -eq 'single-matvec-analytic') {
    'capsid_coarse_solve_pass'
} else { 'capsid_coarse_solve_gap' }

@([pscustomobject]@{ decision=$decision; fixture='6CO8 biological assembly 1'
    atom_count=$atomCount; ses_grid_spacing_A=$SesGridSpacing
    tolerance=$Tolerance; restart=$Restart; ilu_neighbor_cap=$IluNeighborCap
    note='Feasibility gate only; no direct-operator audit is claimed for the capsid row.' }) |
    Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding UTF8
@(
    "decision $decision",
    "atom_count $atomCount",
    "pqr_sha256 $pqrSha",
    "ses_grid_spacing_A $(F $SesGridSpacing)",
    "memory_budget_mb $MemoryBudgetMb",
    "maximum_direct_rhs_pairs $MaximumDirectRhsPairs",
    "matrix $matrixPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
@(
    'TITAN 6CO8 CAPSID FEASIBILITY GATE',
    '===================================',
    '',
    "Decision: $decision",
    "Input atoms: $atomCount",
    "PQR SHA256: $pqrSha",
    "Grid-SES spacing: $(F $SesGridSpacing) A",
    "Coarse policy: degree-one analytical FMM/M2L, theta $(F $OpeningTheta), GMRES restart $Restart, tolerance $(F $Tolerance), sparse ILU(0) neighbor cap $IluNeighborCap.",
    "Memory gate: maximum $MaximumCoarsePanels panels and projected process memory <= $MemoryBudgetMb MB.",
    "Direct RHS gate: at most $MaximumDirectRhsPairs panel-atom pairs before a charge-to-surface tree/FMM is required.",
    '',
    'This gate demonstrates capsid surface feasibility and preflights the coarse solve. A coarse solve is not claimed when the direct panel-atom RHS gate requires charge-to-surface acceleration. It does not reproduce TABI-PB 2.0 Table 7 or claim a sampled dense direct-operator audit for the million-atom row.',
    '',
    "Matrix: $matrixPath",
    "Surface feasibility: $feasibilityPath"
) | Set-Content -LiteralPath $reportPath -Encoding ASCII

Write-Host "6CO8 decision: $decision"
Write-Host "Matrix: $matrixPath"
if (-not $AllowFailedDecision -and
    $decision -notin @('capsid_surface_probe_pass','capsid_coarse_solve_pass')) {
    throw "6CO8 gate decision: $decision"
}
