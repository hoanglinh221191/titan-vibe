param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\gbsa_perf_entropy_smoke'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'lib_ff')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

# Windows Start-Process treats Path/PATH case-insensitively and rejects an
# inherited environment containing both spellings.
if ([Environment]::OSVersion.Platform -eq [PlatformID]::Win32NT) {
    $processPath = $env:Path
    Remove-Item Env:Path -ErrorAction SilentlyContinue
    $env:Path = $processPath
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
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

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $LibRoot)) {
    throw "Missing TITAN force-field library directory: $LibRoot"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$pdb = Join-Path $WorkRoot 'two_molecule_fixture.pdb'
@(
    'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
    'TITLE     GBSA smoke t= 0.00000',
    'MODEL        1',
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
    'ENDMDL',
    'TITLE     GBSA smoke t= 1.00000',
    'MODEL        2',
    'ATOM      1  N   ALA A   1      -4.000   0.000   0.000  1.00  0.00           N',
    'ATOM      2  CA  ALA A   1      -2.550   0.000   0.000  1.00  0.00           C',
    'ATOM      3  C   ALA A   1      -1.950   1.350   0.000  1.00  0.00           C',
    'ATOM      4  O   ALA A   1      -2.550   2.350   0.000  1.00  0.00           O',
    'ATOM      5  CB  ALA A   1      -2.050  -0.750   1.250  1.00  0.00           C',
    'TER',
    'ATOM      6  N   ALA B   2       2.100   0.000   0.000  1.00  0.00           N',
    'ATOM      7  CA  ALA B   2       3.550   0.000   0.000  1.00  0.00           C',
    'ATOM      8  C   ALA B   2       4.150   1.350   0.000  1.00  0.00           C',
    'ATOM      9  O   ALA B   2       3.550   2.350   0.000  1.00  0.00           O',
    'ATOM     10  CB  ALA B   2       4.050  -0.750   1.250  1.00  0.00           C',
    'ENDMDL',
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
$lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
$argPath = Join-Path $WorkRoot 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $pdb)",
    '-1pdb',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $ndx)",
    '-g=1-2',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$lib",
    "-path-n=$lib",
    '-gbsa',
    '-gbsa-perf=detail',
    '-nt=1'
) | Set-Content -LiteralPath $argPath -Encoding ASCII

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for GBSA performance smoke"
    }
} finally {
    Pop-Location
}

$summaryPath = Join-Path $WorkRoot 'gbsa_summary.dat'
$performancePath = Join-Path $WorkRoot 'gbsa_performance.dat'
$framesPath = Join-Path $WorkRoot 'gbsa_performance_frames.csv'
$snapshotsPath = Join-Path $WorkRoot 'traj_0_gbsa_snapshots.dat'
foreach ($required in @($summaryPath, $performancePath, $framesPath, $snapshotsPath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing GBSA output: $required"
    }
}

$summaryText = Get-Content -LiteralPath $summaryPath -Raw
if ($summaryText -notmatch 'interaction_entropy_source E_MM_elec_plus_vdw') {
    throw "GBSA summary must state that IE is computed from E_MM_elec_plus_vdw"
}

$perf = Read-SummaryMap $performancePath
if ($perf['timing_granularity'] -ne 'natural_boundaries_fused_energy_kernel' -or
    -not $perf.ContainsKey('gbsa_energy_kernel_ms_total') -or
    -not $perf.ContainsKey('entropy_output_ms') -or
    -not $perf.ContainsKey('total_ms_total')) {
    throw "GBSA performance report lacks fused-kernel or entropy/output timing contract"
}

$rows = @(Import-Csv -LiteralPath $framesPath)
if ($rows.Count -ne 2 -or
    $null -eq $rows[0].PSObject.Properties['gbsa_energy_kernel_ms'] -or
    $null -ne $rows[0].PSObject.Properties['entropy_output_ms']) {
    throw "GBSA must analyze both the reference snapshot and the subsequent frame while keeping timing scopes separate"
}
$snapshotRows = @(Get-Content -LiteralPath $snapshotsPath | Select-Object -Skip 1)
if ($snapshotRows.Count -ne 2 -or
    $snapshotRows[0] -notmatch '^0\.0000\s+\|' -or
    $snapshotRows[1] -notmatch '^1\.0000\s+\|') {
    throw "GBSA snapshot output must include frame 0 with its original time"
}
if (($snapshotRows -join "`n") -notmatch '-53\.3133' -or
    ($snapshotRows -join "`n") -notmatch '-52\.8034') {
    throw "ALPB-off GBSA polar energies changed from the established fixture"
}
if ($summaryText -match '(?m)^alpb_' -or
    (Test-Path -LiteralPath (Join-Path $WorkRoot 'gbsa_alpb_sizes.csv'))) {
    throw "ALPB-off output must retain its existing summary and artifact contract"
}

$alpbDir = Join-Path $WorkRoot 'alpb'
New-Item -ItemType Directory -Path $alpbDir -Force | Out-Null
$alpbPdb = Join-Path $alpbDir 'two_molecule_fixture.pdb'
$alpbNdx = Join-Path $alpbDir 'two_molecule_fixture.ndx'
Copy-Item -LiteralPath $pdb -Destination $alpbPdb
Copy-Item -LiteralPath $ndx -Destination $alpbNdx
$alpbArg = Join-Path $alpbDir 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $alpbPdb)",
    '-1pdb',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $alpbNdx)",
    '-g=1-2',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$lib",
    "-path-n=$lib",
    '-alpb',
    '-gbsa-perf=detail',
    '-nt=1'
) | Set-Content -LiteralPath $alpbArg -Encoding ASCII
Push-Location $alpbDir
try {
    & $TitanExe -a $alpbArg | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for ALPB integration smoke"
    }
} finally {
    Pop-Location
}

$alpbSummary = Join-Path $alpbDir 'gbsa_summary.dat'
$alpbSizes = Join-Path $alpbDir 'gbsa_alpb_sizes.csv'
$alpbSnapshots = Join-Path $alpbDir 'traj_0_gbsa_snapshots.dat'
$alpbLog = Join-Path $alpbDir 'log.txt'
foreach ($required in @($alpbSummary, $alpbSizes, $alpbSnapshots, $alpbLog)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing ALPB provenance output: $required"
    }
}
$alpbSummaryText = Get-Content -LiteralPath $alpbSummary -Raw
foreach ($field in @(
    'alpb_enabled true',
    'alpb_alpha_dimensionless 0.571412',
    'alpb_formula Sigalov_Eq25_salt_screened',
    'alpb_size_method ADet_Sigalov_Eq23',
    'alpb_sizes_export_status ok',
    'alpb_sizes_file gbsa_alpb_sizes.csv',
    'alpb_A_receptor_nm_mean',
    'alpb_A_ligand_nm_mean',
    'alpb_A_complex_nm_mean'
)) {
    if ($alpbSummaryText -notmatch [regex]::Escape($field)) {
        throw "ALPB summary is missing provenance field: $field"
    }
}
$alpbSizeRows = @(Import-Csv -LiteralPath $alpbSizes)
if ($alpbSizeRows.Count -ne 2) {
    throw "ALPB size artifact must contain one row per analyzed frame"
}
$invariant = [Globalization.CultureInfo]::InvariantCulture
foreach ($row in $alpbSizeRows) {
    foreach ($column in @('A_receptor_nm', 'A_ligand_nm', 'A_complex_nm')) {
        $value = [double]::Parse($row.$column, $invariant)
        if ([double]::IsNaN($value) -or [double]::IsInfinity($value) -or
            $value -le 0.0) {
            throw "ALPB size artifact contains invalid $column=$($row.$column)"
        }
    }
}
$alpbSnapshotText = Get-Content -LiteralPath $alpbSnapshots -Raw
if ($alpbSnapshotText -notmatch '-53\.1605' -or
    $alpbSnapshotText -notmatch '-52\.6429') {
    throw "ALPB integration polar-energy oracle changed"
}
$alpbLogText = Get-Content -LiteralPath $alpbLog -Raw
if ($alpbLogText -notmatch 'MM/GBSA ALPB PROVENANCE' -or
    $alpbLogText -notmatch 'Sigalov Eq\.25 salt-screened' -or
    $alpbLogText -notmatch 'ADet, Sigalov Eq\.23') {
    throw "ALPB runtime log is missing scientific provenance"
}

$alpbNt8Dir = Join-Path $WorkRoot 'alpb_nt8'
New-Item -ItemType Directory -Path $alpbNt8Dir -Force | Out-Null
$alpbNt8Pdb = Join-Path $alpbNt8Dir 'two_molecule_fixture.pdb'
$alpbNt8Ndx = Join-Path $alpbNt8Dir 'two_molecule_fixture.ndx'
Copy-Item -LiteralPath $pdb -Destination $alpbNt8Pdb
Copy-Item -LiteralPath $ndx -Destination $alpbNt8Ndx
$alpbNt8Arg = Join-Path $alpbNt8Dir 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $alpbNt8Pdb)",
    '-1pdb',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $alpbNt8Ndx)",
    '-g=1-2',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$lib",
    "-path-n=$lib",
    '-alpb',
    '-gbsa-perf=detail',
    '-nt=8'
) | Set-Content -LiteralPath $alpbNt8Arg -Encoding ASCII
Push-Location $alpbNt8Dir
try {
    & $TitanExe -a $alpbNt8Arg | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for ALPB nt=8 smoke"
    }
} finally {
    Pop-Location
}
foreach ($name in @(
    'traj_0_gbsa_snapshots.dat',
    'gbsa_summary.dat',
    'gbsa_alpb_sizes.csv'
)) {
    $nt1Hash = (Get-FileHash -LiteralPath (Join-Path $alpbDir $name) `
        -Algorithm SHA256).Hash
    $nt8Hash = (Get-FileHash -LiteralPath (Join-Path $alpbNt8Dir $name) `
        -Algorithm SHA256).Hash
    if ($nt1Hash -cne $nt8Hash) {
        throw "ALPB scientific output is not nt=1/8 deterministic: $name"
    }
}

$offDir = Join-Path $WorkRoot 'perf_off'
New-Item -ItemType Directory -Path $offDir -Force | Out-Null
$offPdb = Join-Path $offDir 'two_molecule_fixture.pdb'
Copy-Item -LiteralPath $pdb -Destination $offPdb
$offNdx = Join-Path $offDir 'two_molecule_fixture.ndx'
Copy-Item -LiteralPath $ndx -Destination $offNdx
$offArgPath = Join-Path $offDir 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $offPdb)",
    '-1pdb',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $offNdx)",
    '-g=1-2',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$lib",
    "-path-n=$lib",
    '-gbsa',
    '-gbsa-perf=off',
    '-nt=1'
) | Set-Content -LiteralPath $offArgPath -Encoding ASCII
Push-Location $offDir
try {
    & $TitanExe -a $offArgPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for GBSA perf=off smoke"
    }
} finally {
    Pop-Location
}
if ((Test-Path -LiteralPath (Join-Path $offDir 'gbsa_performance.dat')) -or
    (Test-Path -LiteralPath (Join-Path $offDir 'gbsa_performance_frames.csv'))) {
    throw "GBSA perf=off must not create performance output files"
}

$identicalDir = Join-Path $WorkRoot 'identical_selection_reject'
New-Item -ItemType Directory -Path $identicalDir -Force | Out-Null
$identicalPdb = Join-Path $identicalDir 'two_molecule_fixture.pdb'
$identicalNdx = Join-Path $identicalDir 'two_molecule_fixture.ndx'
Copy-Item -LiteralPath $pdb -Destination $identicalPdb
Copy-Item -LiteralPath $ndx -Destination $identicalNdx
$identicalArg = Join-Path $identicalDir 'arg.in'
@(
    "-f=$(Convert-ToTitanPath $identicalPdb)",
    '-1pdb=true',
    '-ff=AMBER99SB',
    "-n=$(Convert-ToTitanPath $identicalNdx)",
    '-g=1-1',
    '-mol-1-type=protein',
    '-mol-2-type=lig',
    "-path-a=$lib",
    "-path-n=$lib",
    '-gbsa=true',
    '-nt=1'
) | Set-Content -LiteralPath $identicalArg -Encoding ASCII
$identicalStdout = Join-Path $identicalDir 'stdout.txt'
$identicalStderr = Join-Path $identicalDir 'stderr.txt'
$identicalProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a', $identicalArg) `
    -WorkingDirectory $identicalDir `
    -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $identicalStdout `
    -RedirectStandardError $identicalStderr
$identicalExit = $identicalProcess.ExitCode
$identicalOutput = (Get-Content -LiteralPath $identicalStdout -Raw) + "`n" +
                   (Get-Content -LiteralPath $identicalStderr -Raw)
if ($identicalExit -eq 0 -or
    $identicalOutput -notmatch 'invalid_gbsa_identical_selections') {
    throw 'GBSA did not clearly reject identical receptor/ligand selections'
}

Write-Host "GBSA entropy and performance smoke PASS"
