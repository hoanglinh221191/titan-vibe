param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_llps_local_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$ffPath = Join-Path $repo 'example\2\lib_ff'
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
    $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

function Write-CommonIndex {
    param([string]$Path)
    @'
[ MOL_A ]
1 2 3

[ MOL_B ]
4 5 6
'@ | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Get-FirstModelAtomCoordinates {
    param([string]$Path)
    $coords = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line.StartsWith('ENDMDL')) { break }
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
            continue
        }
        $serial = [int]$line.Substring(6, 5).Trim()
        $coords[$serial] = [pscustomobject]@{
            X = [double]$line.Substring(30, 8).Trim()
            Y = [double]$line.Substring(38, 8).Trim()
            Z = [double]$line.Substring(46, 8).Trim()
        }
    }
    return $coords
}

function Assert-Near {
    param(
        [double]$Actual,
        [double]$Expected,
        [double]$Tolerance,
        [string]$Message
    )
    if ([Math]::Abs($Actual - $Expected) -gt $Tolerance) {
        throw "$Message (actual=$Actual expected=$Expected tolerance=$Tolerance)"
    }
}

function Distance {
    param($A, $B)
    $dx = $A.X - $B.X
    $dy = $A.Y - $B.Y
    $dz = $A.Z - $B.Z
    return [Math]::Sqrt($dx * $dx + $dy * $dy + $dz * $dz)
}

function Invoke-TopologyRangeCase {
    $caseRoot = Join-Path $WorkRoot 'topology_range'
    New-Item -ItemType Directory -Force -Path $caseRoot | Out-Null
    $pdb = Join-Path $caseRoot 'llps_topology_range_fixture.pdb'
    $ndx = Join-Path $caseRoot 'index.ndx'
    $top = Join-Path $caseRoot 'topol.top'
    $arg = Join-Path $caseRoot 'arg.in'
    $outName = 'llps_topology_range_out.pdb'

@'
REMARK    PBCTOPO LLPS-local topology range smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1      19.200  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       0.800  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3       2.400  10.000  10.000  1.00  0.00           C
ATOM      4  CA  ALA A   1       2.500  11.000  10.000  1.00  0.00           C
ATOM      5  CA  GLY A   2       4.100  11.000  10.000  1.00  0.00           C
ATOM      6  CA  SER A   3       5.700  11.000  10.000  1.00  0.00           C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA A   1      19.200  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       0.800  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3       2.400  10.000  10.000  1.00  0.00           C
ATOM      4  CA  ALA A   1       2.500  11.000  10.000  1.00  0.00           C
ATOM      5  CA  GLY A   2       4.100  11.000  10.000  1.00  0.00           C
ATOM      6  CA  SER A   3       5.700  11.000  10.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

    Write-CommonIndex -Path $ndx

@'
[ moleculetype ]
PEP_A 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
PEP_B 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ system ]
pbctopo llps-local topology range smoke

[ molecules ]
PEP_A 1
PEP_B 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

    $argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo-llps-local=true
-pbctopo-bond-cutoff=2.20
-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-path-a=$ffPath
-path-n=$ffPath
-opdb=$outName
"@
    [System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

    Push-Location $caseRoot
    try {
        & $titanPath -a $arg
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN topology range case exited with code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }

    $debugLog = Join-Path $caseRoot 'pbctopo_debug.log'
    $report = Join-Path $caseRoot 'pbctopo_frame_report.csv'
    $outPdb = Join-Path $caseRoot $outName
    foreach ($path in @($debugLog, $report, $outPdb)) {
        if (-not (Test-Path -LiteralPath $path)) {
            throw "Missing expected output: $path"
        }
    }

    $logText = Get-Content -Raw -LiteralPath $debugLog
    foreach ($pattern in @(
            'pbctopo_mode=llps-local',
            'pbctopo_topology_source=gromacs_topology',
            'pbctopo_llps_local_frames=2',
            'pbctopo_llps_local_units=4',
            'pbctopo_llps_local_unit_source=topology_range',
            'pbctopo_llps_local_edges_skipped_cross_unit=\d+'
        )) {
        if ($logText -notmatch $pattern) {
            throw "Missing expected LLPS-local topology log pattern: $pattern"
        }
    }

    $rows = @(Import-Csv -LiteralPath $report)
    if ($rows.Count -ne 2 -or @($rows | Where-Object { $_.mode -ne 'llps-local' -or $_.status -ne 'local_unwrap' -or $_.hard_certified -ne '0' }).Count -ne 0) {
        throw "Unexpected LLPS-local frame report row: $($rows | ConvertTo-Json -Compress)"
    }

    $coords = Get-FirstModelAtomCoordinates -Path $outPdb
    Assert-Near $coords[1].X 19.2 0.25 'Topology case should keep the first molecule anchor image'
    Assert-Near (Distance $coords[1] $coords[2]) 1.6 0.25 'Topology case should unwrap bond 1-2 internally'
    Assert-Near (Distance $coords[2] $coords[3]) 1.6 0.25 'Topology case should unwrap bond 2-3 internally'
    if ((Distance $coords[3] $coords[4]) -lt 10.0) {
        throw 'Topology case compacted two topology molecule ranges together'
    }
}

function Invoke-ChainFallbackCase {
    $caseRoot = Join-Path $WorkRoot 'chain_fallback'
    New-Item -ItemType Directory -Force -Path $caseRoot | Out-Null
    $sourcePdb = Join-Path $repo 'example\2\md1_short.pdb'
    $sourceIndex = Join-Path $repo 'example\2\index.ndx'
    $pdb = Join-Path $caseRoot 'llps_chain_fallback_two_frames.pdb'
    $arg = Join-Path $caseRoot 'arg.in'
    $outName = 'llps_chain_fallback_out.pdb'

    $frameOne = [System.IO.File]::ReadAllText($sourcePdb)
    $frameTwo = $frameOne.Replace('MODEL        1', 'MODEL        2').Replace(
        't= 100000.00000', 't= 100001.00000')
    [System.IO.File]::WriteAllText(
        $pdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)

    $argsText = @"
-ff=AMBER14SB
-f=$pdb
-pbctopo-mode=llps-local
-pbctopo-bond-cutoff=2.20
-g=3-4
-nt=1
-1pdb=true
-nocal=true
-n=$sourceIndex
-no-dssp=true
-noprint-side_chain_contact=true
-path-a=$ffPath
-path-n=$ffPath
-opdb=$outName
"@
    [System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

    Push-Location $caseRoot
    try {
        & $titanPath -a $arg
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN chain fallback case exited with code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }

    $debugLog = Join-Path $caseRoot 'pbctopo_debug.log'
    $outPdb = Join-Path $caseRoot $outName
    if (-not (Test-Path -LiteralPath $debugLog)) {
        throw "Missing chain fallback debug log"
    }
    if (-not (Test-Path -LiteralPath $outPdb)) {
        throw "Missing chain fallback output PDB"
    }

    $logText = Get-Content -Raw -LiteralPath $debugLog
    foreach ($pattern in @(
            'pbctopo_mode=llps-local',
            'pbctopo_topology_source=metadata_chain_residue',
            'pbctopo_llps_local_frames=2',
            'pbctopo_llps_local_units=([1-9]\d*)',
            'pbctopo_llps_local_unit_source=chain'
        )) {
        if ($logText -notmatch $pattern) {
            throw "Missing expected LLPS-local chain fallback log pattern: $pattern"
        }
    }

    if ((Select-String -Path $outPdb -Pattern '^ATOM').Count -le 0) {
        throw 'Chain fallback did not write atom coordinates'
    }
}

Invoke-TopologyRangeCase
Invoke-ChainFallbackCase

Write-Host 'PASS pbctopo_llps_local_smoke'
