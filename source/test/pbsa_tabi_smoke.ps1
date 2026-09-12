param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_tabi_smoke'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'example\1\lib_ff')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

# Windows Start-Process rejects an inherited environment containing both
# `Path` and `PATH`; normalize it before the negative subprocess check.
if ([Environment]::OSVersion.Platform -eq [PlatformID]::Win32NT) {
    $processPath = $env:Path
    Remove-Item Env:Path -ErrorAction SilentlyContinue
    $env:Path = $processPath
}

function Read-SummaryMap([string]$Path) {
    $map = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ([string]::IsNullOrWhiteSpace($line) -or $line.StartsWith('#')) { continue }
        $parts = $line -split '\s+', 2
        if ($parts.Count -eq 2) { $map[$parts[0]] = $parts[1] }
    }
    return $map
}

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing TITAN executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $LibRoot)) { throw "Missing force-field library: $LibRoot" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$LibRoot = (Resolve-Path -LiteralPath $LibRoot).Path
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)

if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$pdb = Join-Path $WorkRoot 'two_molecule_fixture.pdb'
@(
    'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
    'ATOM      1  N   ALA A   1      -4.000   0.000   0.000  1.00  0.00           N',
    'ATOM      2  CA  ALA A   1      -2.550   0.000   0.000  1.00  0.00           C',
    'ATOM      3  C   ALA A   1      -1.950   1.350   0.000  1.00  0.00           C',
    'ATOM      4  O   ALA A   1      -2.550   2.350   0.000  1.00  0.00           O',
    'ATOM      5  CB  ALA A   1      -2.050  -0.750   1.250  1.00  0.00           C',
    'TER',
    'ATOM      6  N   ALA B   2       2.000   0.000   0.000  1.00  0.00           N',
    'ATOM      7  CA  ALA B   2       3.450   0.000   0.000  1.00  0.00           C',
    'ATOM      8  C   ALA B   2       4.050   1.350   0.000  1.00  0.00           C',
    'ATOM      9  O   ALA B   2       3.450   2.350   0.000  1.00  0.00           O',
    'ATOM     10  CB  ALA B   2       3.950  -0.750   1.250  1.00  0.00           C',
    'END'
) | Set-Content -LiteralPath $pdb -Encoding ASCII

$ndx = Join-Path $WorkRoot 'two_molecule_fixture.ndx'
@(
    '[ System ]',
    '1 2 3 4 5 6 7 8 9 10',
    '[ Receptor ]',
    '1 2 3 4 5',
    '[ Ligand ]',
    '6 7 8 9 10'
) | Set-Content -LiteralPath $ndx -Encoding ASCII

$arg = Join-Path $WorkRoot 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $pdb)",
    '-1pdb',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $ndx)",
    '-g=1-2',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$(Convert-ToTitanPath $LibRoot)/",
    "-path-n=$(Convert-ToTitanPath $LibRoot)/",
    '-pbsa',
    '-pbsa-components=complex,receptor,ligand',
    '-pbsa-decomp=hybrid',
    '-pbsa-residue-frames=all',
    '-pbsa-ie=false',
    '-pbsa-perf=detail',
    '-pbsa-grid-reuse=off',
    '-edecomp',
    '-edecomp-source=pbsa',
    '-edecomp-level=both',
    '-pbe=tabi',
    '-pbe-model=linear',
    '-pbe-biomodel=pbe',
    '-pbe-dielectric=two-dielectric',
    '-pbe-int-diel=1',
    '-pbe-ext-diel=80',
    '-pbe-salt=0.15',
    '-pbe-temp=298.15',
    '-pbe-surface=vdw',
    '-pbe-map=off',
    '-pbe-focus=off',
    '-pbe-multigrid=off',
    '-pbe-bie-sum=direct',
    '-pbe-bie-fmm-m2l=off',
    '-pbe-bie-primary-surface=sphere-patch',
    '-pbe-bie-discretization=centroid',
    '-pbe-bie-near-integration=analytic_singular',
    '-pbe-bie-quad-buffer=off',
    '-pbe-bie-precond=diagonal',
    '-pbe-bie-surface-points=8',
    '-pbe-bie-gmres-restart=32',
    '-pbe-max-iter=120',
    '-pbe-rms=1e-4',
    '-nt=1'
) | Set-Content -LiteralPath $arg -Encoding ASCII

Push-Location $WorkRoot
try {
    & $TitanExe -a $arg | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$summary = Read-SummaryMap (Join-Path $WorkRoot 'mm_pbsa_summary.dat')
if ($summary['pbsa_backend'] -ne 'tabi') { throw "Expected pbsa_backend=tabi" }
if ($summary['endpoint_solves_per_frame'] -ne '3') { throw "Expected three TABI endpoint solves" }
if ($summary['pbsa_grid_reuse_effective'] -ne 'not_applicable') {
    throw "Expected TABI grid reuse to be not_applicable"
}
if ($summary['pbsa_decomp_effective'] -ne 'hybrid') { throw "Expected hybrid decomposition" }

$energy = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'mm_pbsa_energy.csv'))
if ($energy.Count -ne 1 -or $energy[0].backend -ne 'tabi' -or $energy[0].status -ne 'ok') {
    throw "Expected one successful TABI MM/PBSA energy row"
}

$units = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'mm_pbsa_unit_summary.csv'))
if (-not ($units | Where-Object { $_.molecule_role -eq 'receptor' -and $_.unit_type -eq 'residue' })) {
    throw 'Missing receptor residue decomposition'
}
if (-not ($units | Where-Object { $_.molecule_role -eq 'ligand' -and $_.unit_type -eq 'atom' })) {
    throw 'Missing ligand atom decomposition'
}

foreach ($name in @('edecomp_mmpbsa_residue_map.csv', 'edecomp_mmpbsa_atom_map.csv',
                     'mm_pbsa_endpoint_timing.csv', 'mm_pbsa_performance.dat')) {
    if (-not (Test-Path -LiteralPath (Join-Path $WorkRoot $name))) {
        throw "Missing TABI MM/PBSA output: $name"
    }
}

$decomp = Read-SummaryMap (Join-Path $WorkRoot 'mm_pbsa_decomposition.dat')
$closure = [double]::Parse($decomp['unit_closure_total_max_abs_kJmol'], $InvariantCulture)
if ([math]::Abs($closure) -gt 1e-5) { throw "Decomposition closure failed: $closure kJ/mol" }

$rejectRoot = Join-Path $WorkRoot 'grid_reuse_validate_reject'
New-Item -ItemType Directory -Path $rejectRoot -Force | Out-Null
$rejectArg = Join-Path $rejectRoot 'arg.in'
@((Get-Content -LiteralPath $arg) | ForEach-Object {
    if ($_ -eq '-pbsa-grid-reuse=off') { '-pbsa-grid-reuse=validate' } else { $_ }
}) | Set-Content -LiteralPath $rejectArg -Encoding ASCII

$rejectStdout = Join-Path $rejectRoot 'stdout.txt'
$rejectStderr = Join-Path $rejectRoot 'stderr.txt'
$rejectProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a', $rejectArg) `
    -WorkingDirectory $rejectRoot `
    -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $rejectStdout `
    -RedirectStandardError $rejectStderr
$rejectExit = $rejectProcess.ExitCode
if ($rejectExit -eq 0) {
    throw 'Expected TABI MM/PBSA to reject -pbsa-grid-reuse=validate'
}

$nonconvergedRoot = Join-Path $WorkRoot 'strict_nonconvergence'
New-Item -ItemType Directory -Path $nonconvergedRoot -Force | Out-Null
$nonconvergedArg = Join-Path $nonconvergedRoot 'arg.in'
@((Get-Content -LiteralPath $arg) | ForEach-Object {
    if ($_ -eq '-pbe-max-iter=120') { '-pbe-max-iter=1' }
    elseif ($_ -eq '-pbe-rms=1e-4') { '-pbe-rms=1e-14' }
    else { $_ }
}) | Set-Content -LiteralPath $nonconvergedArg -Encoding ASCII

$nonconvergedStdout = Join-Path $nonconvergedRoot 'stdout.txt'
$nonconvergedStderr = Join-Path $nonconvergedRoot 'stderr.txt'
$nonconvergedProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a', $nonconvergedArg) `
    -WorkingDirectory $nonconvergedRoot `
    -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $nonconvergedStdout `
    -RedirectStandardError $nonconvergedStderr
if ($nonconvergedProcess.ExitCode -eq 0) {
    throw 'Expected a strict MM/PBSA non-convergence failure'
}
$failedSummary = Read-SummaryMap (Join-Path $nonconvergedRoot 'mm_pbsa_summary.dat')
if ($failedSummary['frames_valid'] -ne '0' -or
    $failedSummary['frames_failed'] -ne '1' -or
    $failedSummary['statistics_status'] -ne 'no_valid_frames') {
    throw 'Failed MM/PBSA frame leaked into aggregate statistics'
}
$failedEnergy = @(Import-Csv -LiteralPath (Join-Path $nonconvergedRoot 'mm_pbsa_energy.csv'))
if ($failedEnergy.Count -ne 1 -or
    $failedEnergy[0].included_in_statistics -ne 'false' -or
    $failedEnergy[0].status -eq 'ok' -or
    $failedEnergy[0].DeltaG_bind_kJmol -ne 'nan') {
    throw 'Failed MM/PBSA frame did not retain the strict diagnostic-output contract'
}
$failureText = (Get-Content -LiteralPath $nonconvergedStdout -Raw) + "`n" +
               (Get-Content -LiteralPath $nonconvergedStderr -Raw)
if ($failureText -notmatch 'strict failure policy') {
    throw 'Strict MM/PBSA failure did not report its failure policy'
}

Write-Host 'PBSA TABI smoke PASS'
