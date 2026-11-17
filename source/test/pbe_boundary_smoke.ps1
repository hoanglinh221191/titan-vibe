param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_boundary_smoke')
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function Convert-ToDouble([object]$Value) { return [double]::Parse([string]$Value, $InvariantCulture) }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBE boundary smoke root: $full"
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
    @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'END'
    ) | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Invoke-PbeBoundaryCase([string]$Name, [string]$Boundary) {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Get-PbeGridFixturePath -FixtureRoot $FixtureRoot -FixtureName 'two-atoms'
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
        '-pbe-scale=0.6',
        '-pbe-margin=5.0',
        '-pbe-surface=vdw',
        "-pbe-boundary=$Boundary",
        '-pbe-salt=0.0',
        '-pbe-max-iter=1000',
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-5',
        '-pbe-max-grid=300000',
        '-pbe-map=off',
        '-nt=1'
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE for $Name" }
    } finally {
        Pop-Location
    }
    $row = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    if ($row.boundary_mode -ne $Boundary) {
        throw "Boundary mode mismatch for ${Name}: expected=$Boundary observed=$($row.boundary_mode)"
    }
    foreach ($field in @('E_rxn_kJmol', 'boundary_rms', 'boundary_max')) {
        $value = Convert-ToDouble $row.$field
        if ([double]::IsNaN($value) -or [double]::IsInfinity($value)) {
            throw "Non-finite $field for ${Name}: $value"
        }
    }
    return $row
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$FixtureRoot = Join-Path $WorkRoot 'fixtures'
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridFixtures -FixtureRoot $FixtureRoot -FixtureNames @('two-atoms')
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot

$coulombic = Invoke-PbeBoundaryCase -Name 'coulombic' -Boundary 'coulombic'
$coulombicRepeat = Invoke-PbeBoundaryCase -Name 'coulombic_repeat' -Boundary 'coulombic'
$dipolar = Invoke-PbeBoundaryCase -Name 'dipolar' -Boundary 'dipolar'
$zero = Invoke-PbeBoundaryCase -Name 'zero' -Boundary 'zero'

$repeatDelta = [math]::Abs((Convert-ToDouble $coulombic.E_rxn_kJmol) - (Convert-ToDouble $coulombicRepeat.E_rxn_kJmol))
if ($repeatDelta -gt 0.000000001) {
    throw "Coulombic boundary was not deterministic; repeat delta=$repeatDelta kJ/mol"
}
if ((Convert-ToDouble $coulombic.boundary_max) -le 0.0) {
    throw "Coulombic boundary did not report nonzero boundary potential"
}
if ((Convert-ToDouble $dipolar.boundary_max) -le 0.0) {
    throw "Dipolar boundary did not report nonzero boundary potential"
}
if ([math]::Abs((Convert-ToDouble $zero.boundary_rms)) -gt 0.000000000001 -or
    [math]::Abs((Convert-ToDouble $zero.boundary_max)) -gt 0.000000000001) {
    throw "Zero boundary should report zero boundary diagnostics"
}

Write-Host ("PBE boundary smoke PASS: repeat_delta={0} zero_E={1} dipolar_E={2}" -f
    $repeatDelta, $zero.E_rxn_kJmol, $dipolar.E_rxn_kJmol)
