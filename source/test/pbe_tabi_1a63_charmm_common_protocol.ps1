param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_1a63_charmm_common_protocol',
    [string]$SourcePdb = 'test\pbe_inputs_public\1a63_model1_chainA.pdb',
    [string]$FixtureRoot = 'test\pbe_inputs_public',
    [string]$PythonExe = 'python',
    [string]$Pdb2PqrTarget = '.codex_tmp\pdb2pqr-3.7.1',
    [double[]]$GridSpacings = @(1.50,1.25,1.00,0.80),
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$RefreshFixture,
    [switch]$PrepareOnly,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pdb2PqrVersion = '3.7.1'
$PreparationArguments = '--ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0'
$PaperFagbiEnergy = -2382.50
$PaperTabiPb2Energy = -2393.50

function Resolve-RepoPath([string]$Path) {
    if ([IO.Path]::IsPathRooted($Path)) { return [IO.Path]::GetFullPath($Path) }
    return [IO.Path]::GetFullPath((Join-Path $RepoRoot $Path))
}
function Parse-Double([string]$Value) {
    return [double]::Parse($Value, [Globalization.NumberStyles]::Float, $Inv)
}
function Format-Double([double]$Value) { return $Value.ToString('G15', $Inv) }
function Get-PqrStats([string]$Path) {
    $atoms = 0
    $charge = 0.0
    $minRadius = [double]::PositiveInfinity
    $maxRadius = 0.0
    $malformed = 0
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        ++$atoms
        $fields = @($line -split '\s+' | Where-Object { $_ })
        try {
            if ($fields.Count -lt 10) { throw 'short row' }
            $q = Parse-Double ($fields[-2])
            $r = Parse-Double ($fields[-1])
            if ([double]::IsNaN($q) -or [double]::IsInfinity($q) -or
                [double]::IsNaN($r) -or [double]::IsInfinity($r) -or $r -le 0.0) {
                throw 'non-finite or non-positive PQR field'
            }
            $charge += $q
            $minRadius = [math]::Min($minRadius, $r)
            $maxRadius = [math]::Max($maxRadius, $r)
        } catch {
            ++$malformed
        }
    }
    return [pscustomobject]@{
        atoms=$atoms; total_charge=$charge; minimum_radius_A=$minRadius
        maximum_radius_A=$maxRadius; malformed_rows=$malformed
    }
}
function Invoke-Pdb2Pqr([string]$InputPath, [string]$OutputPath) {
    $target = Resolve-RepoPath $Pdb2PqrTarget
    $marker = Join-Path $target 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        throw "Pinned PDB2PQR $Pdb2PqrVersion is missing at $target. Run test\pbe_fagbi_scope_z4g_prepare_public_fixtures.ps1 -RefreshFixtures once, or install pdb2pqr==$Pdb2PqrVersion into that target."
    }
    $oldPythonPath = $env:PYTHONPATH
    try {
        $env:PYTHONPATH = $target
        & $PythonExe -m pdb2pqr --ff=CHARMM --ffout=CHARMM --keep-chain `
            --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath
        if ($LASTEXITCODE -ne 0) { throw "PDB2PQR CHARMM preparation failed with exit code $LASTEXITCODE" }
    } finally {
        $env:PYTHONPATH = $oldPythonPath
    }
}

$SourcePdb = Resolve-RepoPath $SourcePdb
$FixtureRoot = Resolve-RepoPath $FixtureRoot
$WorkRoot = Resolve-RepoPath $WorkRoot
$TitanExe = Resolve-RepoPath $TitanExe
$PqrPath = Join-Path $FixtureRoot '1a63_full_charmm_ph7.pqr'
$FixtureManifestPath = Join-Path $FixtureRoot '1a63_charmm_common_protocol_manifest.csv'
New-Item -ItemType Directory -Force $FixtureRoot,$WorkRoot | Out-Null

if ($RefreshFixture -or -not (Test-Path -LiteralPath $PqrPath)) {
    Invoke-Pdb2Pqr -InputPath $SourcePdb -OutputPath $PqrPath
}
if (-not (Test-Path -LiteralPath $PqrPath)) { throw "Missing CHARMM PQR fixture: $PqrPath" }

$stats = Get-PqrStats $PqrPath
if ($stats.atoms -ne 2065 -or $stats.malformed_rows -ne 0 -or
    [math]::Abs($stats.total_charge - (-1.0)) -gt 1.0e-8) {
    throw "CHARMM PQR quality gate failed: atoms=$($stats.atoms), charge=$($stats.total_charge), malformed=$($stats.malformed_rows)."
}
$sourceSha = (Get-FileHash $SourcePdb -Algorithm SHA256).Hash.ToLowerInvariant()
$pqrSha = (Get-FileHash $PqrPath -Algorithm SHA256).Hash.ToLowerInvariant()
$fixtureRow = [pscustomobject][ordered]@{
    fixture='1A63'; pdb_id='1A63'; protein='RNA-binding domain of E. coli Rho factor'
    source_pdb=$SourcePdb; source_pdb_sha256=$sourceSha; pqr=$PqrPath
    pqr_sha256=$pqrSha; atom_count=$stats.atoms
    total_charge_e=Format-Double $stats.total_charge
    minimum_radius_A=Format-Double $stats.minimum_radius_A
    maximum_radius_A=Format-Double $stats.maximum_radius_A
    preparation_tool='PDB2PQR'; preparation_version=$Pdb2PqrVersion
    force_field='CHARMM'; preparation_arguments=$PreparationArguments
    protocol_scope='paper-aligned input; original paper did not report PDB2PQR version or protonation flags'
}
$fixtureRow | Export-Csv -LiteralPath $FixtureManifestPath -NoTypeInformation

if ($PrepareOnly) {
    Write-Host "Prepared CHARMM fixture: $PqrPath"
    Write-Host "Fixture manifest: $FixtureManifestPath"
    exit 0
}
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing TITAN executable: $TitanExe" }

$calibrationScript = Join-Path $PSScriptRoot 'pbe_tabi_1a63_policy_calibration.ps1'
$calibrationReport = Join-Path $WorkRoot 'TITAN_TABI_1A63_CHARMM_CALIBRATION_REPORT.txt'
$calibrationArgs = @{
    TitanExe=$TitanExe
    WorkRoot=$WorkRoot
    PqrPath=$PqrPath
    ReportPath=$calibrationReport
    Stages=@('surface')
    GridSpacings=$GridSpacings
    PolicyProbeRadius=1.40
    PolicyOpeningTheta=0.45
    PolicyIluNeighborCap=64
    Threads=$Threads
    CaseTimeoutSeconds=$CaseTimeoutSeconds
    AllowFailedDecision=$true
}
if ($Quick) { $calibrationArgs.Quick = $true }
if ($Force) { $calibrationArgs.Force = $true }
if ($RerunFailed) { $calibrationArgs.RerunFailed = $true }
& $calibrationScript @calibrationArgs

$matrixPath = Join-Path $WorkRoot 'tabi_1a63_policy_matrix.csv'
$comparisonPath = Join-Path $WorkRoot 'tabi_1a63_charmm_literature_context.csv'
$decisionPath = Join-Path $WorkRoot 'tabi_1a63_charmm_common_protocol_decision.csv'
$summaryPath = Join-Path $WorkRoot 'tabi_1a63_charmm_common_protocol_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_TABI_1A63_CHARMM_COMMON_PROTOCOL_REPORT.txt'
if (-not (Test-Path -LiteralPath $matrixPath)) { throw "Calibration matrix is missing: $matrixPath" }
$matrix = @(Import-Csv -LiteralPath $matrixPath)
$completed = @($matrix | Where-Object { $_.status -eq 'ok' -and $_.backend -eq 'fmm' })
$valid = @($completed | Where-Object {
    (Parse-Double $_.true_residual) -le 1.0e-4 -and
    (Parse-Double $_.audit_rms) -le 8.0e-3
})
$comparison = [Collections.Generic.List[object]]::new()
foreach ($row in $completed) {
    $energy = Parse-Double $row.deltaG_solv_kcal_mol
    $operatorGatePass = (Parse-Double $row.true_residual) -le 1.0e-4 -and
        (Parse-Double $row.audit_rms) -le 8.0e-3
    foreach ($reference in @(
        [pscustomobject]@{source='FAGBI Cartesian FMM';kind='finite mesh, N=132133';energy=$PaperFagbiEnergy},
        [pscustomobject]@{source='TABI-PB 2.0';kind='node-patch extrapolation';energy=$PaperTabiPb2Energy}
    )) {
        $delta = $energy - $reference.energy
        $comparison.Add([pscustomobject][ordered]@{
            fixture='1A63'; force_field='CHARMM'; grid_spacing_A=$row.grid_spacing_A
            titan_panels=$row.panels; titan_energy_kcal_mol=$energy
            literature_source=$reference.source; literature_reference_kind=$reference.kind
            literature_energy_kcal_mol=$reference.energy; signed_delta_kcal_mol=$delta
            absolute_delta_kcal_mol=[math]::Abs($delta)
            relative_delta_percent=100.0*[math]::Abs($delta)/[math]::Abs($reference.energy)
            operator_gate_pass=$operatorGatePass
            comparison_scope='context only; surface generator and discretization are not matched'
        })
    }
}
$comparison | Export-Csv -LiteralPath $comparisonPath -NoTypeInformation

$timedOut = @($matrix | Where-Object { $_.status -eq 'timeout' }).Count
$failed = @($matrix | Where-Object { $_.status -notin @('ok','timeout') }).Count
$sorted = @($valid | Sort-Object { Parse-Double $_.grid_spacing_A } -Descending)
$surfaceSequence = @($completed | Sort-Object { Parse-Double $_.grid_spacing_A } -Descending)
$adjacentDelta = [double]::NaN
$surfaceMonotonic = $true
if ($sorted.Count -ge 2) {
    $a = Parse-Double $sorted[-1].deltaG_solv_kcal_mol
    $b = Parse-Double $sorted[-2].deltaG_solv_kcal_mol
    $adjacentDelta = [math]::Abs($a-$b)/[math]::Max(1.0e-300,[math]::Abs($a))
}
if ($surfaceSequence.Count -ge 2) {
    $previousSign = 0
    for ($i=1; $i -lt $surfaceSequence.Count; ++$i) {
        $coarse = Parse-Double $surfaceSequence[$i-1].deltaG_solv_kcal_mol
        $fine = Parse-Double $surfaceSequence[$i].deltaG_solv_kcal_mol
        $delta = $fine - $coarse
        $sign = if ($delta -gt 0.0) { 1 } elseif ($delta -lt 0.0) { -1 } else { 0 }
        if ($sign -ne 0 -and $previousSign -ne 0 -and $sign -ne $previousSign) {
            $surfaceMonotonic = $false
        }
        if ($sign -ne 0) { $previousSign = $sign }
    }
}
$decision = if ($timedOut -gt 0) {
    'charmm_common_protocol_runtime_limited'
} elseif ($failed -gt 0 -or $valid.Count -eq 0) {
    if ($completed.Count -gt 0) { 'charmm_common_protocol_operator_accuracy_gap' }
    else { 'charmm_common_protocol_harness_gap' }
} elseif (-not $surfaceMonotonic) {
    'charmm_common_protocol_surface_nonmonotonic'
} elseif (-not [double]::IsNaN($adjacentDelta) -and $adjacentDelta -le 0.03) {
    'charmm_common_protocol_mesh_stable'
} else {
    'charmm_common_protocol_ready_with_mesh_gap'
}
[pscustomobject]@{
    decision=$decision; fixture_sha256=$pqrSha; completed_rows=$completed.Count
    valid_rows=$valid.Count
    timeout_rows=$timedOut; failed_rows=$failed
    finest_adjacent_relative_energy_delta=$adjacentDelta
    surface_energy_monotonic=$surfaceMonotonic
    literature_comparison_is_oracle=$false
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation
@(
    "decision=$decision"
    "fixture_pqr=$PqrPath"
    "fixture_sha256=$pqrSha"
    "atom_count=$($stats.atoms)"
    "total_charge_e=$(Format-Double $stats.total_charge)"
    "finest_adjacent_relative_energy_delta=$(Format-Double $adjacentDelta)"
    "surface_energy_monotonic=$surfaceMonotonic"
    "matrix=$matrixPath"
    "literature_context=$comparisonPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
@(
    'TITAN TABI 1A63 CHARMM Common-Protocol Report'
    "Decision: $decision"
    ''
    'Input protocol:'
    "- PDB2PQR $Pdb2PqrVersion with CHARMM, PROPKA at pH 7.0, chain A."
    '- Dielectric constants 1/80, ionic concentration 0.15 M, GMRES tolerance 1e-4.'
    '- TITAN watertight grid SES, degree-one analytical FMM/M2L, theta 0.45, sparse near-field ILU(0) cap 64.'
    ''
    'Interpretation boundary:'
    '- TABI-PB 2.0 reports PDB2PQR with CHARMM but does not report its PDB2PQR version or protonation flags.'
    '- Published energies are included as context, not as an oracle, because NanoShaper/node-patch geometry is not reproduced here.'
    '- A non-monotonic surface-density sequence is classified explicitly and is not presented as mesh convergence.'
    '- Direct TITAN lib_ff preparation is excluded from the scientific matrix until its unmatched-atom gate reaches zero.'
    ''
    'Outputs:'
    "- $FixtureManifestPath"
    "- $matrixPath"
    "- $comparisonPath"
    "- $decisionPath"
) | Set-Content -LiteralPath $reportPath -Encoding UTF8

Write-Host "CHARMM common-protocol decision: $decision"
Write-Host "Matrix: $matrixPath"
Write-Host "Literature context: $comparisonPath"
if ($decision -notin @('charmm_common_protocol_mesh_stable','charmm_common_protocol_ready_with_mesh_gap') -and
    -not $AllowFailedDecision) { exit 2 }
