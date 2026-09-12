param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_surface_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBE surface smoke root: $full"
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

function Invoke-SurfaceCase([string]$Fixture, [string]$Surface) {
    $caseName = "${Fixture}_${Surface}"
    $runDir = Join-Path $WorkRoot $caseName
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath

    $fixtureName = if ($Fixture -eq 'twoatoms') { 'two-atoms' } else { $Fixture }
    $pqrPath = Get-PbeGridFixturePath -FixtureRoot $FixtureRoot -FixtureName $fixtureName

    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ForceFieldRoot).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-scale=2.0',
        '-pbe-margin=3.0',
        "-pbe-surface=$Surface",
        '-pbe-probe=1.4',
        '-pbe-max-iter=80',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-3',
        '-pbe-max-grid=500000',
        '-pbe-map=off',
        '-nt=1'
    )
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
    if (-not (Test-Path -LiteralPath $csvPath)) {
        throw "Missing pbe_energy.csv for $caseName"
    }
    $rows = @(Import-Csv -LiteralPath $csvPath)
    if ($rows.Count -ne 1) {
        throw "Expected one PBE row for $caseName, observed $($rows.Count)"
    }
    if ($rows[0].surface -ne $Surface) {
        throw "Surface mismatch for ${caseName}: $($rows[0].surface)"
    }
    if ($rows[0].status -ne 'ok') {
        throw "PBE surface case did not converge for ${caseName}: $($rows[0].status)"
    }
    return [int]$rows[0].solute_points
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$FixtureRoot = Join-Path $WorkRoot 'fixtures'
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridFixtures -FixtureRoot $FixtureRoot -FixtureNames @('sphere', 'two-atoms')
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot

foreach ($fixture in @('sphere', 'twoatoms')) {
    $vdw = Invoke-SurfaceCase -Fixture $fixture -Surface 'vdw'
    $sas = Invoke-SurfaceCase -Fixture $fixture -Surface 'sas'
    $ses = Invoke-SurfaceCase -Fixture $fixture -Surface 'ses'
    if ($sas -lt $vdw) {
        throw "SAS solute points regressed below VDW for ${fixture}: vdw=$vdw sas=$sas"
    }
    if ($ses -lt $vdw -or $ses -ge $sas) {
        throw "SES closing contract failed for ${fixture}: vdw=$vdw ses=$ses sas=$sas"
    }
    if ($fixture -eq 'sphere' -and ($ses - $vdw) -ge ($sas - $ses)) {
        throw "Single-sphere SES is not closer to VDW than SAS: vdw=$vdw ses=$ses sas=$sas"
    }
    Write-Host "PBE surface smoke PASS: $fixture vdw=$vdw sas=$sas ses=$ses"
}

Write-Host "PBE surface smoke PASS: all cases"
