param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_glycan_continuity_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
    $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$pdb = Join-Path $WorkRoot 'glycan_continuity_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'
$outName = 'glycan_continuity_out.pdb'

@'
REMARK    PBCTOPO glycan continuity smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      12.000  10.000  10.000  1.00  0.00           C
HETATM    3  C1  UYB A   1      50.000  50.000  50.000  1.00  0.00           C
HETATM    4  C2  0SA A   2      52.000  50.000  50.000  1.00  0.00           C
ATOM      5  CA  DUM B   1      20.000  20.000  20.000  1.00  0.00           C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      12.000  10.000  10.000  1.00  0.00           C
HETATM    3  C1  UYB A   1      50.000  50.000  50.000  1.00  0.00           C
HETATM    4  C2  0SA A   2      52.000  50.000  50.000  1.00  0.00           C
ATOM      5  CA  DUM B   1      20.000  20.000  20.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ GLYCOPROTEIN ]
1 2 3 4

[ DUMMY ]
5
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
GLYCOPROTEIN 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 1 UYB C1 1 0.0 12.0
4 CT 2 0SA C2 1 0.0 12.0

[ bonds ]
1 2 1
3 4 1

[ moleculetype ]
DUMMY 3

[ atoms ]
1 CT 1 DUM CA 1 0.0 12.0

[ system ]
pbctopo glycan continuity smoke

[ molecules ]
GLYCOPROTEIN 1
DUMMY 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

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
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
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
if (-not (Test-Path -LiteralPath $debugLog)) {
    throw "Missing pbctopo debug log: $debugLog"
}
$logText = Get-Content -Raw -LiteralPath $debugLog

foreach ($pattern in @(
        'pbctopo_topology_source=gromacs_topology',
        'pbctopo_topology_component_source=topology',
        'pbctopo_topology_static_graph_frames=2',
        'pbctopo_geometry_edges_added=0'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected glycan-continuity pbctopo log pattern: $pattern"
    }
}

if ($logText -match 'component_search_prefix_continuity_hard_rejects=[1-9]') {
    throw 'Same-chain glycan residues created a hard sequence-continuity reject'
}

if ($logText -match 'fell back to pbci=0') {
    throw 'Glycan continuity fixture unexpectedly fell back to pbci=0'
}

Write-Host 'PASS pbctopo_glycan_continuity_smoke'
