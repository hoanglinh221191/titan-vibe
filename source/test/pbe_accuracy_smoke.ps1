param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_accuracy_smoke')
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
        throw "Refusing to clean path outside PBE accuracy smoke root: $full"
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
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -3.200   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -1.800   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1      -1.160   0.500  -1.320  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -1.720   1.260  -2.100  1.00  0.00           O',
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

function Invoke-PbeAccuracyCase([string]$Name, [string]$FocusMode) {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath

    $pqrPath = Get-PbeGridFixturePath -FixtureRoot $FixtureRoot -FixtureName 'mini-protein'

    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ForceFieldRoot).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        "-common_topo=$(Convert-ToTitanPath (Write-PbeDriverTopology $pdbPath))",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-scale=0.45',
        '-pbe-margin=3.0',
        '-pbe-gridbox=cuboidal',
        '-pbe-surface=ses',
        '-pbe-probe=1.4',
        '-pbe-max-iter=500',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-3',
        '-pbe-max-grid=500000',
        '-pbe-map=off',
        '-nt=1'
    )
    if ($FocusMode -eq 'auto') {
        $argLines += @(
            '-pbe-multigrid=focus',
            '-pbe-focus-levels=2',
            '-pbe-focus-margin=4.0'
        )
    } else {
        $argLines += '-pbe-focus=off'
    }

    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $Name"
        }
    } finally {
        Pop-Location
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $summaryPath = Join-Path $runDir 'pbe_summary.dat'
    foreach ($required in @($energyPath, $summaryPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing PBE output for ${Name}: $required"
        }
    }

    $rows = @(Import-Csv -LiteralPath $energyPath)
    if ($rows.Count -ne 1) {
        throw "Expected one PBE row for ${Name}, observed $($rows.Count)"
    }
    $row = $rows[0]
    $summary = Read-SummaryMap $summaryPath

    if ($row.status -ne 'ok') {
        throw "PBE accuracy case did not converge for ${Name}: $($row.status)"
    }
    $energyKJ = Convert-ToDouble $row.E_rxn_kJmol
    if ([double]::IsNaN($energyKJ) -or [double]::IsInfinity($energyKJ)) {
        throw "Non-finite PBE energy for ${Name}: $energyKJ"
    }
    if ([int]$row.grid_nx -eq [int]$row.grid_ny -and
        [int]$row.grid_ny -eq [int]$row.grid_nz) {
        throw "Expected cuboidal PBE grid for ${Name}; observed cubic $($row.grid_nx)x$($row.grid_ny)x$($row.grid_nz)"
    }

    if ($FocusMode -eq 'auto') {
        if ($row.focus_enabled -ne 'true' -or [int]$row.focus_levels -lt 2) {
            throw "Focused PBE row did not report active focusing for $Name"
        }
        if ($summary['focus_mode'] -ne 'auto') {
            throw "Focused PBE summary did not report focus_mode=auto for $Name"
        }
        if ((Convert-ToDouble $row.focus_initial_scale) -ge 0.45) {
            throw "Focused PBE initial scale was not coarser for $Name"
        }
        if ((Convert-ToDouble $row.focus_initial_margin) -le 3.0) {
            throw "Focused PBE initial margin was not expanded for $Name"
        }
    } else {
        if ($row.focus_enabled -ne 'false') {
            throw "Baseline PBE row unexpectedly enabled focusing for $Name"
        }
    }

    Write-Host ("PBE accuracy smoke PASS: {0} focus={1} grid={2}x{3}x{4} E_kJ={5}" -f
        $Name, $FocusMode, $row.grid_nx, $row.grid_ny, $row.grid_nz, $row.E_rxn_kJmol)
    return $row
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$FixtureRoot = Join-Path $WorkRoot 'fixtures'
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridFixtures -FixtureRoot $FixtureRoot -FixtureNames @('mini-protein')
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot
[void](Invoke-PbeAccuracyCase -Name 'mini_protein_baseline' -FocusMode 'off')
[void](Invoke-PbeAccuracyCase -Name 'mini_protein_focus' -FocusMode 'auto')
Write-Host "PBE accuracy smoke PASS: all cases"
