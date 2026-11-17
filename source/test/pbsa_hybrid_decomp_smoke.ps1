param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_hybrid_decomp_smoke'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'lib_ff')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Read-SummaryMap([string]$Path) {
    $map = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ([string]::IsNullOrWhiteSpace($line) -or $line.StartsWith('#')) {
            continue
        }
        $parts = $line -split '\s+', 2
        if ($parts.Count -eq 2) {
            $map[$parts[0]] = $parts[1]
        }
    }
    return $map
}

function New-SmokeInput([string]$RunDir, [string]$Mol2Type = 'lig',
                        [string]$ReuseMode = 'validate',
                        [string[]]$ExtraArgLines = @()) {
    New-Item -ItemType Directory -Path $RunDir -Force | Out-Null
    # Receptor chain A, ligand chain B, full topology: see the fixture notes in
    # pbe_grid_fixture_helpers.ps1. -common_topo comes before the
    # -mol-N-type lines, which then keep the protein/ligand roles.
    $pdb = Join-Path $RunDir 'ala2.pdb'
    $top = Join-Path $RunDir 'ala2.top'
    Write-PbsaTwoChainPdb $pdb
    Write-PbsaTwoChainTopology $top
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
    $arg = Join-Path $RunDir 'arg.in'
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdb)",
        '-1pdb',
        '-ff=AMBER99SB',
        "-g=$PbsaTwoChainGroups",
        "-common_topo=$(Convert-ToTitanPath $top)",
        '-mol-1-type=protein',
        "-mol-2-type=$Mol2Type",
        "-path-a=$lib",
        "-path-n=$lib",
        '-pbsa=true',
        '-pbsa-components=complex,receptor,ligand',
        '-pbsa-residue-frames=all',
        '-pbsa-ie=true',
        '-pbsa-ie-min-samples=1',
        '-pbsa-ie-block=1',
        '-pbsa-perf=detail',
        "-pbsa-grid-reuse=$ReuseMode",
        '-pbe-scale=0.35',
        '-pbe-margin=2.0',
        '-pbe-max-iter=80',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-3',
        '-pbe-max-grid=200000',
        '-pbe-map=off',
        '-nt=1'
    )
    $argLines += $ExtraArgLines
    $argLines | Set-Content -LiteralPath $arg -Encoding ASCII
    return $arg
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (-not (Test-Path -LiteralPath $LibRoot)) {
    throw "Missing TITAN force-field library directory: $LibRoot"
}
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
$argPath = New-SmokeInput -RunDir $WorkRoot

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for PBSA hybrid smoke"
    }
} finally {
    Pop-Location
}

$summary = Read-SummaryMap (Join-Path $WorkRoot 'mm_pbsa_summary.dat')
if ($summary['pbsa_backend'] -ne 'grid') {
    throw "Expected PBSA production backend grid, observed $($summary['pbsa_backend'])"
}
if ($summary['pbsa_decomp_requested'] -ne 'auto' -or
    $summary['pbsa_decomp_effective'] -ne 'hybrid') {
    throw "Expected default decomposition auto -> hybrid for protein/ligand"
}
if ($summary['interaction_entropy_source'] -ne 'E_MM_elec_plus_vdw') {
    throw "Expected PBSA IE source E_MM_elec_plus_vdw, observed $($summary['interaction_entropy_source'])"
}
if ($summary['pbsa_grid_reuse_requested'] -ne 'validate' -or
    $summary['pbsa_grid_reuse_effective'] -ne 'validate') {
    throw "Expected PBSA grid reuse validation mode in summary"
}

$energyRows = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'mm_pbsa_energy.csv'))
if ($energyRows.Count -ne 1 -or $energyRows[0].backend -ne 'grid') {
    throw "Expected one PBSA energy row reporting backend=grid"
}

$unitSummary = Join-Path $WorkRoot 'mm_pbsa_unit_summary.csv'
$unitFrame = Join-Path $WorkRoot 'mm_pbsa_unit_frame.csv'
$decompPath = Join-Path $WorkRoot 'mm_pbsa_decomposition.dat'
$endpointTiming = Join-Path $WorkRoot 'mm_pbsa_endpoint_timing.csv'
$reuseValidation = Join-Path $WorkRoot 'mm_pbsa_grid_reuse_validation.csv'
$performance = Join-Path $WorkRoot 'mm_pbsa_performance.dat'
$performanceFrames = Join-Path $WorkRoot 'mm_pbsa_performance_frames.csv'
foreach ($required in @($unitSummary, $unitFrame, $decompPath,
                         $endpointTiming, $reuseValidation,
                         $performance, $performanceFrames)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing hybrid decomposition output: $required"
    }
}

$perf = Read-SummaryMap $performance
if ($perf['timing_granularity'] -ne 'natural_boundaries' -or
    -not $perf.ContainsKey('pb_solve_ms_mean') -or
    -not $perf.ContainsKey('mm_sasa_components_ms_mean') -or
    -not $perf.ContainsKey('entropy_output_ms') -or
    -not $perf.ContainsKey('total_ms_mean')) {
    throw "PBSA performance report must expose natural-boundary MM/SASA and entropy/output timing"
}
$performanceFrameRows = @(Import-Csv -LiteralPath $performanceFrames)
if ($performanceFrameRows.Count -ne 1 -or
    $null -eq $performanceFrameRows[0].PSObject.Properties['pb_solve_ms'] -or
    $null -eq $performanceFrameRows[0].PSObject.Properties['mm_sasa_components_ms'] -or
    $null -ne $performanceFrameRows[0].PSObject.Properties['entropy_output_ms']) {
    throw "PBSA detail timing must report per-frame component timing without assigning aggregate entropy time to frames"
}

$unitRows = @(Import-Csv -LiteralPath $unitSummary)
if (-not ($unitRows | Where-Object { $_.molecule_role -eq 'receptor' -and $_.unit_type -eq 'residue' })) {
    throw "Hybrid decomposition lacks receptor residue rows"
}
if (-not ($unitRows | Where-Object { $_.molecule_role -eq 'ligand' -and $_.unit_type -eq 'atom' -and $_.atom_name.Length -gt 0 })) {
    throw "Hybrid decomposition lacks ligand atom rows"
}
foreach ($column in @('E_elec_kJmol_mean', 'E_vdw_kJmol_mean',
                      'DeltaDeltaG_solv_PB_kJmol_mean',
                      'DeltaG_nonpolar_SASA_kJmol_mean',
                      'DeltaG_no_entropy_kJmol_mean',
                      'minus_TDeltaS_IE_kJmol',
                      'DeltaG_bind_IE_kJmol')) {
    if ($null -eq $unitRows[0].PSObject.Properties[$column]) {
        throw "Missing hybrid output column: $column"
    }
}

$decomp = Read-SummaryMap $decompPath
$closure = Convert-ToDouble $decomp['unit_closure_total_max_abs_kJmol']
if ([math]::Abs($closure) -gt 1e-5) {
    throw "Hybrid no-entropy closure exceeds tolerance: $closure"
}
if ($decomp['unit_entropy_closure'] -ne 'not_applicable_nonlinear') {
    throw "Hybrid output must state that unit entropy is non-additive"
}

$timingRows = @(Import-Csv -LiteralPath $endpointTiming)
if ($timingRows.Count -ne 12 -or
    -not ($timingRows | Where-Object { $_.endpoint -eq 'complex' -and $_.phase -eq 'solvent' -and [double]$_.surface_ms -ge 0.0 }) -or
    -not ($timingRows | Where-Object { $_.endpoint -eq 'ligand' -and $_.phase -eq 'vacuum' -and $_.timing_variant -eq 'baseline' -and [int]$_.reuse_surface_hits -eq 0 }) -or
    -not ($timingRows | Where-Object { $_.endpoint -eq 'ligand' -and $_.phase -eq 'vacuum' -and $_.timing_variant -eq 'reuse' -and [int]$_.reuse_surface_hits -eq 1 })) {
    throw "Expected baseline and reuse timing rows for each validated PBSA endpoint"
}
$reuseRows = @(Import-Csv -LiteralPath $reuseValidation)
if ($reuseRows.Count -ne 3 -or
    ($reuseRows | Where-Object { [math]::Abs((Convert-ToDouble $_.reaction_energy_delta_kJmol)) -gt 1e-4 })) {
    throw "PBSA reuse validation failed its numerical gate"
}

$nonLigDir = Join-Path $WorkRoot 'auto_non_ligand'
$nonLigArg = New-SmokeInput -RunDir $nonLigDir -Mol2Type 'protein' -ReuseMode 'off'
Push-Location $nonLigDir
try {
    & $TitanExe -a $nonLigArg | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for PBSA auto non-ligand smoke"
    }
} finally {
    Pop-Location
}
$nonLigSummary = Read-SummaryMap (Join-Path $nonLigDir 'mm_pbsa_summary.dat')
if ($nonLigSummary['pbsa_decomp_requested'] -ne 'auto' -or
    $nonLigSummary['pbsa_decomp_effective'] -ne 'off' -or
    (Test-Path -LiteralPath (Join-Path $nonLigDir 'mm_pbsa_unit_summary.csv'))) {
    throw "PBSA auto decomposition must remain off when no ligand role is declared"
}

$vcycleDir = Join-Path $WorkRoot 'vcycle'
$vcycleArg = New-SmokeInput -RunDir $vcycleDir -ReuseMode 'off' -ExtraArgLines @(
    '-pbe-multigrid=vcycle',
    '-pbe-mg-levels=2',
    '-pbe-mg-pre=2',
    '-pbe-mg-post=2',
    '-pbe-mg-coarse-iter=120',
    # one V-cycle costs 124 sweeps; 80 allowed one cycle, which cannot
    # converge (see pbsa_vcycle_smoke.ps1)
    '-pbe-max-iter=1240'
)
Push-Location $vcycleDir
try {
    & $TitanExe -a $vcycleArg | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for PBSA vcycle smoke"
    }
} finally {
    Pop-Location
}
$vcycleRows = @(Import-Csv -LiteralPath (Join-Path $vcycleDir 'mm_pbsa_energy.csv'))
if ($vcycleRows.Count -ne 1 -or $vcycleRows[0].backend -ne 'grid' -or
    $vcycleRows[0].solver_mode -ne 'vcycle') {
    throw "PBSA vcycle smoke did not preserve grid backend and vcycle solver"
}

Write-Host "PBSA hybrid decomposition smoke PASS"
