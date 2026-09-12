param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_no_topology_geometry_graph_smoke')
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

$pdb = Join-Path $WorkRoot 'no_topology_geometry_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$arg = Join-Path $WorkRoot 'arg.in'
$outName = 'no_topology_geometry_out.pdb'

@'
REMARK    PBCTOPO no-topology geometry graph smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      25.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3      40.000  10.000  10.000  1.00  0.00           C
ATOM      4  CA  DUM B   1      60.000  20.000  20.000  1.00  0.00           C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      25.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3      40.000  10.000  10.000  1.00  0.00           C
ATOM      4  CA  DUM B   1      60.000  20.000  20.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ CHAIN_A ]
1 2 3

[ DUMMY ]
4
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

$argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
-pbctopo-bond-cutoff=1.10
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=64
-pbctopo-anchor-backend=auto
-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-no-dssp=true
-noprint-side_chain_contact=true
-path-a=$ffPath
-path-n=$ffPath
-opdb=$outName
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
$report = Join-Path $WorkRoot 'pbctopo_frame_report.csv'
foreach ($path in @($debugLog, $report)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing expected output: $path"
    }
}

$logText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'pbctopo_topology_source=metadata_chain_residue',
        'pbctopo_topology_component_source=metadata',
        'pbctopo_geometry_fallback_frames=0',
        'component_graph_failed=0'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected no-topology geometry graph log pattern: $pattern"
    }
}
if ($logText -match 'fell back to pbci=0') {
    throw 'No-topology geometry graph fixture unexpectedly fell back to pbci=0'
}

$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -ne 2 -or
    @($rows | Where-Object { $_.status -ne 'weak_observation' }).Count -ne 0 -or
    @($rows | Where-Object {
            [int]$_.component_count -ne 4 -or
            $_.certificate_scope -ne 'image_intercomponent_steric' -or
            $_.certificate_graph_source -ne 'none' -or
            $_.search_evidence_source -ne 'none' -or
            $_.identified -ne '0' -or
            $_.hard_certified -ne '1' -or
            $_.fallback_reason -ne 'evidence_disconnected' -or
            [int]$_.hard_connectivity_edges_available -ne 0 -or
            [int]$_.metadata_edges_available -ne 0
        }).Count -ne 0) {
    throw "Unexpected no-topology geometry report row: $($rows | ConvertTo-Json -Compress)"
}

Write-Host 'PASS pbctopo_no_topology_geometry_graph_smoke'
