param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_multigrid_smoke')
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBE multigrid smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   34.000   34.000   34.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -3.200   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -1.800   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1      -1.160   0.500  -1.320  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -1.720   1.260  -2.100  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Write-PbsaPdb([string]$Path) {
    $lines = @(
        'CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
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

function Invoke-PbeCase([string]$Fixture, [string]$Mode) {
    $caseName = "${Fixture}_${Mode}"
    $runDir = Join-Path $WorkRoot $caseName
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath

    $pqrPath = Get-PbeGridFixturePath -FixtureRoot $FixtureRoot -FixtureName $Fixture

    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ForceFieldRoot).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        "-common_topo=$(Convert-ToTitanPath (Write-PbeDriverTopology $pdbPath))",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=grid',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-scale=0.45',
        '-pbe-margin=3.0',
        '-pbe-gridbox=cuboidal',
        '-pbe-surface=ses',
        '-pbe-probe=1.4',
        '-pbe-max-iter=1600',
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-5',
        '-pbe-max-grid=900000',
        '-pbe-map=off',
        '-nt=1'
    )
    if ($Mode -eq 'focus') {
        $argLines += @('-pbe-multigrid=focus', '-pbe-focus-levels=2')
    } elseif ($Mode -eq 'nl_default') {
        # GRID-04: nonlinear PB with the default solver (auto = Newton)
        $argLines += @('-pbe-multigrid=off', '-pbe-model=nonlinear', '-pbe-salt=0.15')
    } elseif ($Mode -eq 'nl_picard') {
        $argLines += @('-pbe-multigrid=off', '-pbe-model=nonlinear', '-pbe-salt=0.15',
                       '-pbe-nl-solver=picard')
    } elseif ($Mode -eq 'vcycle') {
        $argLines += @(
            '-pbe-multigrid=vcycle',
            '-pbe-mg-levels=2',
            '-pbe-mg-pre=2',
            '-pbe-mg-post=2',
            '-pbe-mg-coarse-iter=200'
        )
    } else {
        $argLines += '-pbe-multigrid=off'
    }

    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $caseName"
        }
    } finally {
        Pop-Location
    }

    $csvPath = Join-Path $runDir 'pbe_energy.csv'
    $summaryPath = Join-Path $runDir 'pbe_summary.dat'
    foreach ($required in @($csvPath, $summaryPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing PBE multigrid output for ${caseName}: $required"
        }
    }
    $rows = @(Import-Csv -LiteralPath $csvPath)
    if ($rows.Count -ne 1) {
        throw "Expected one PBE row for ${caseName}, observed $($rows.Count)"
    }
    $row = $rows[0]
    if ($Mode -eq 'nl_picard') {
        # GRID-04: Picard may stall; ok must then mean residual <= tolerance
        $residual = Convert-ToDouble $row.nl_residual
        $within = $residual -le 1.0e-5
        if (($row.status -eq 'ok') -ne $within -or ($row.nl_status -eq 'ok') -ne $within) {
            throw "Picard status disagrees with its residual for ${caseName}: status=$($row.status) nl_status=$($row.nl_status) residual=$residual"
        }
        return $row
    }
    if ($row.status -ne 'ok') {
        throw "PBE multigrid case did not converge for ${caseName}: $($row.status)"
    }
    if ($Mode -eq 'nl_default') {
        $residual = Convert-ToDouble $row.nl_residual
        if ($row.nl_solver -ne 'newton' -or $row.nl_status -ne 'ok' -or $residual -gt 1.0e-5) {
            throw "Default nonlinear PB for ${caseName}: solver=$($row.nl_solver) nl_status=$($row.nl_status) residual=$residual (expected newton, ok, <= 1e-5)"
        }
    }
    $energy = Convert-ToDouble $row.E_rxn_kJmol
    if ([double]::IsNaN($energy) -or [double]::IsInfinity($energy)) {
        throw "Non-finite energy for ${caseName}: $energy"
    }
    if ([int]$row.grid_nx -eq [int]$row.grid_ny -and
        [int]$row.grid_ny -eq [int]$row.grid_nz) {
        throw "Expected cuboidal grid for ${caseName}; observed cubic $($row.grid_nx)x$($row.grid_ny)x$($row.grid_nz)"
    }

    $summary = Read-SummaryMap $summaryPath
    if ($Mode -eq 'vcycle') {
        if ($row.solver_mode -ne 'vcycle') {
            throw "Expected solver_mode=vcycle for ${caseName}; observed $($row.solver_mode)"
        }
        if ($summary['multigrid_effective'] -ne 'vcycle') {
            throw "Expected multigrid_effective=vcycle for ${caseName}"
        }
        if ([int]$row.mg_levels -lt 2 -or [int]$row.mg_cycles -le 0) {
            throw "Expected active V-cycle diagnostics for ${caseName}"
        }
        if ([int]$row.iterations_solvent -gt 1600 -or
            [int]$row.iterations_vacuum -gt 1600) {
            throw "V-cycle exceeded -pbe-max-iter for ${caseName}: solvent=$($row.iterations_solvent) vacuum=$($row.iterations_vacuum)"
        }
    }
    return $row
}

function Invoke-PbsaVcycleCase() {
    $runDir = Join-Path $WorkRoot 'pbsa_vcycle'
    New-CleanDir $runDir
    # ENE-21: MM/PBSA refuses -g=x&x, so the case uses the two-chain fixture
    # of pbe_grid_fixture_helpers.ps1 (receptor chain A, ligand chain B).
    # One V-cycle costs (2 + 2) + 80 sweeps here and convergence is judged
    # between cycles, so the budget allows several cycles.
    $pdbPath = Join-Path $runDir 'ala2.pdb'
    $topPath = Join-Path $runDir 'ala2.top'
    Write-PbsaTwoChainPdb $pdbPath
    Write-PbsaTwoChainTopology $topPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ForceFieldRoot).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        "-g=$PbsaTwoChainGroups",
        "-common_topo=$(Convert-ToTitanPath $topPath)",
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbsa=true',
        '-pbe-scale=0.35',
        '-pbe-margin=2.0',
        '-pbe-max-iter=1240',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-3',
        '-pbe-max-grid=200000',
        '-pbe-map=off',
        '-pbe-multigrid=vcycle',
        '-pbe-mg-levels=2',
        '-pbe-mg-coarse-iter=80',
        '-nt=1'
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for PBSA V-cycle"
        }
    } finally {
        Pop-Location
    }

    foreach ($required in @('pbe_energy.csv', 'pbe_summary.dat', 'mm_pbsa_summary.dat')) {
        $path = Join-Path $runDir $required
        if (-not (Test-Path -LiteralPath $path)) {
            throw "Missing expected PBSA V-cycle output: $path"
        }
    }
    $rows = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))
    if ($rows.Count -ne 1 -or $rows[0].solver_mode -ne 'vcycle' -or
        $rows[0].dgsolv_status -ne 'ok') {
        throw "PBSA V-cycle did not report vcycle dgsolv success"
    }
    if ([int]$rows[0].iterations_solvent -gt 1240 -or
        [int]$rows[0].iterations_vacuum -gt 1240) {
        throw "PBSA V-cycle exceeded -pbe-max-iter."
    }
    $pbsa = Read-SummaryMap (Join-Path $runDir 'mm_pbsa_summary.dat')
    if ($pbsa['pbsa_enabled'] -ne 'true') {
        throw "PBSA V-cycle summary did not report pbsa_enabled=true"
    }
    Write-Host 'PBE multigrid smoke PASS: pbsa_vcycle'
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$FixtureRoot = Join-Path $WorkRoot 'fixtures'
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridFixtures -FixtureRoot $FixtureRoot -FixtureNames @(
    'mini-protein',
    'charged-peptide',
    'compact-mini-protein'
)
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot

foreach ($fixture in @('mini-protein', 'charged-peptide', 'compact-mini-protein')) {
    $baseline = Invoke-PbeCase -Fixture $fixture -Mode 'sor'
    [void](Invoke-PbeCase -Fixture $fixture -Mode 'focus')
    $vcycle = Invoke-PbeCase -Fixture $fixture -Mode 'vcycle'
    $baseEnergy = Convert-ToDouble $baseline.E_rxn_kJmol
    $vcycleEnergy = Convert-ToDouble $vcycle.E_rxn_kJmol
    $delta = [math]::Abs($baseEnergy - $vcycleEnergy)
    $limit = [math]::Max(0.01, 1.0e-5 * [math]::Max([math]::Abs($baseEnergy), 1.0))
    if ($delta -gt $limit) {
        throw "V-cycle parity mismatch for ${fixture}: baseline=$baseEnergy vcycle=$vcycleEnergy limit=$limit"
    }
    if ([int]$vcycle.mg_cycles -ge [int]$baseline.iterations_solvent) {
        throw "V-cycle cycle diagnostic did not improve over SOR iterations for ${fixture}: cycles=$($vcycle.mg_cycles) sor_iter=$($baseline.iterations_solvent)"
    }
    Write-Host ("PBE multigrid smoke PASS: {0} sor_iter={1} mg_cycles={2} E_kJ={3}" -f
        $fixture, $baseline.iterations_solvent, $vcycle.mg_cycles, $vcycle.E_rxn_kJmol)
    $nonlinear = Invoke-PbeCase -Fixture $fixture -Mode 'nl_default'
    $picard = Invoke-PbeCase -Fixture $fixture -Mode 'nl_picard'
    Write-Host ("PBE multigrid smoke PASS: {0} nonlinear default newton residual={1} E_kJ={2}; picard {3} residual={4}" -f
        $fixture, $nonlinear.nl_residual, $nonlinear.E_rxn_kJmol, $picard.status, $picard.nl_residual)
}

Invoke-PbsaVcycleCase
Write-Host "PBE multigrid smoke PASS: all cases"
