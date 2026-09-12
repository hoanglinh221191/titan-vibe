param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_grid_dgsolv_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Culture = [System.Globalization.CultureInfo]::InvariantCulture
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function Parse-Number([object]$Value) {
    return [double]::Parse([string]$Value, $Culture)
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$fixtureRoot = Join-Path $WorkRoot 'fixtures'
$ffRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridFixtures -FixtureRoot $fixtureRoot -FixtureNames @('mini-protein')
Initialize-PbeGridForceField -ForceFieldRoot $ffRoot
$pqr = Get-PbeGridFixturePath -FixtureRoot $fixtureRoot -FixtureName 'mini-protein'
$driver = Join-Path $WorkRoot 'driver.pdb'
Set-Content -LiteralPath $driver -Encoding ASCII -Value @(
    'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
    'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
    'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
    'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
    'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
    'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
    'END'
)

$lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ffRoot).Path.TrimEnd('\') + '\')
$arg = Join-Path $WorkRoot 'arg.in'
Set-Content -LiteralPath $arg -Encoding ASCII -Value @(
    "-f=$(Convert-ToTitanPath $driver)",
    '-1pdb=true',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$lib",
    "-path-n=$lib",
    '-pbe=grid',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqr)",
    '-pbe-dgsolv=true',
    '-pbe-surface=ses',
    '-pbe-probe=1.4',
    '-pbe-scale=1.0',
    '-pbe-margin=5.0',
    '-pbe-max-grid=500000',
    '-pbe-max-iter=500',
    '-pbe-rms=1e-5',
    '-pbe-dphi=1e-4',
    '-pbe-map=off',
    '-nt=2'
)

Push-Location $WorkRoot
try {
    $env:TITAN_NO_PAUSE = '1'
    & $TitanExe -a $arg | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "Grid dgsolv run failed with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$energyPath = Join-Path $WorkRoot 'pbe_energy.csv'
$summaryPath = Join-Path $WorkRoot 'pbe_summary.dat'
if (-not (Test-Path -LiteralPath $energyPath)) { throw 'Missing pbe_energy.csv' }
if (-not (Test-Path -LiteralPath $summaryPath)) { throw 'Missing pbe_summary.dat' }
$rows = @(Import-Csv -LiteralPath $energyPath)
if ($rows.Count -ne 1) { throw "Expected one energy row, observed $($rows.Count)" }
$row = $rows[0]
if ($row.dgsolv_status -ne 'ok') { throw "Grid dgsolv status is $($row.dgsolv_status)" }
$reaction = Parse-Number $row.E_rxn_kJmol
$dgsolv = Parse-Number $row.deltaG_solv_kJmol
if ([double]::IsNaN($reaction) -or [double]::IsInfinity($reaction) -or
    [double]::IsNaN($dgsolv) -or [double]::IsInfinity($dgsolv)) {
    throw "Grid dgsolv produced a non-finite energy: reaction=$reaction dgsolv=$dgsolv"
}
if ([math]::Abs($reaction - $dgsolv) -gt 1e-8 * [math]::Max(1.0, [math]::Abs($reaction))) {
    throw "Grid single-system dgsolv does not match reaction energy: $reaction vs $dgsolv"
}
$summary = Get-Content -LiteralPath $summaryPath -Raw
if ($summary -notmatch '(?m)^deltaG_solv_frames\s+1\s*$') {
    throw 'Grid dgsolv summary did not record exactly one frame.'
}

Write-Host "PBE Grid dgsolv smoke passed: $dgsolv kJ/mol"
