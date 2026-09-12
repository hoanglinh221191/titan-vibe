param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_mixed_local_graph_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$ffRoot = (Resolve-Path -LiteralPath (Join-Path $repo 'example\2\lib_ff')).Path
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
    $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$pdb = Join-Path $WorkRoot 'mixed_graph_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO mixed topology and metadata graph smoke
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1 A1   RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 A2   RES A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      3 A3   RES A   1       6.000   0.000   0.000  1.00  0.00           C
ATOM      4 CA   ALA B   2      15.000   0.000   0.000  1.00  0.00           C
ATOM      5 C    ALA B   2      18.000   0.000   0.000  1.00  0.00           C
ENDMDL
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1 A1   RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 A2   RES A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      3 A3   RES A   1       6.000   0.000   0.000  1.00  0.00           C
ATOM      4 CA   ALA B   2      15.000   0.000   0.000  1.00  0.00           C
ATOM      5 C    ALA B   2      18.000   0.000   0.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ TEST ]
1 2 3

[ OTHER ]
4 5
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
TEST 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0
3 CT 1 RES A3 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
OTHER 3

[ atoms ]
1 CT 1 OTH B1 1 0.0 12.0
2 CT 1 OTH B2 1 0.0 12.0

[ bonds ]
1 2 1

[ system ]
pbctopo mixed graph smoke

[ molecules ]
TEST 1
OTHER 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-path-a=$ffRoot
-path-n=$ffRoot
-f=$pdb
-pbctopo=true
-pbctopo-bond-cutoff=2.15
-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=mixed_graph_out.pdb
"@
[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a $arg
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

$debugLog = Join-Path $WorkRoot 'pbctopo_debug.log'
$reportPath = Join-Path $WorkRoot 'pbctopo_frame_report.csv'
if (-not (Test-Path -LiteralPath $debugLog) -or
    -not (Test-Path -LiteralPath $reportPath)) {
    throw 'Missing mixed local-graph output'
}
$logText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'pbctopo_topology_source=mixed',
        'pbctopo_topology_static_graph_frames=2',
        'pbctopo_metadata_graph_frames=2',
        'pbctopo_geometry_edges_added=0',
        'pbctopo_hard_connectivity_edges=6',
        'pbctopo_steric_exclusion_pairs=2',
        'pbctopo_topology_component_source=topology\+metadata'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing mixed local-graph contract: $pattern"
    }
}

$rows = @(Import-Csv -LiteralPath $reportPath)
if ($rows.Count -ne 2) {
    throw "Expected two report rows, found $($rows.Count)"
}
foreach ($row in $rows) {
    if ($row.status -ne 'weak_observation' -or
        $row.component_count -ne '2' -or
        $row.hard_certified -ne '1' -or
        $row.identified -ne '0' -or
        $row.fallback_reason -ne 'evidence_disconnected') {
        throw "Mixed graph component contract failed at frame $($row.frame)"
    }
}

Write-Host 'PASS pbctopo_mixed_local_graph_smoke'
