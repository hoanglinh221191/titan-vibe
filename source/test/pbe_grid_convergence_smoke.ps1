param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_grid_convergence_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

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
    "-f=$(Convert-ToTitanPath $driver)", '-1pdb=true', '-ff=AMBER99SB', '-g=1&1',
    "-path-a=$lib", "-path-n=$lib", '-pbe=grid',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqr)", '-pbe-surface=vdw',
    '-pbe-scale=1.0', '-pbe-margin=5.0', '-pbe-max-grid=500000',
    '-pbe-max-iter=1', '-pbe-rms=1e9', '-pbe-dphi=1e-20',
    '-pbe-map=off', '-nt=2'
)

Push-Location $WorkRoot
try {
    $env:TITAN_NO_PAUSE = '1'
    & $TitanExe -a $arg | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "Grid convergence run failed with code $LASTEXITCODE" }
} finally { Pop-Location }

$row = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbe_energy.csv'))[0]
if ($row.status -eq 'ok') {
    throw "Grid accepted convergence when only RMS passed: rms=$($row.rms_solvent) dphi=$($row.max_dphi_solvent)"
}
if ([int]$row.iterations_solvent -ne 1 -or [int]$row.iterations_vacuum -ne 1) {
    throw "Grid convergence fixture did not exhaust its one-iteration budget."
}

Write-Host 'PBE Grid convergence smoke passed: both update criteria are required.'
