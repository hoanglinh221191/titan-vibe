param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_topology_single_selection_smoke')
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

$pdb = Join-Path $WorkRoot 'single_selection_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'
$runLog = Join-Path $WorkRoot 'run.log'
$stdoutLog = Join-Path $WorkRoot 'stdout.log'
$stderrLog = Join-Path $WorkRoot 'stderr.log'

@'
REMARK    PBCTOPO single-selection topology smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1 A1   RES A   1       1.000   1.000   1.000  1.00  0.00           C
ATOM      2 A2   RES A   1       2.500   1.000   1.000  1.00  0.00           C
ATOM      3 B1   OTH B   2      15.000  15.000  15.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ SELECTED ]
1 2

[ UNUSED ]
3
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
SELECTED 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
UNUSED 3

[ atoms ]
1 CT 1 OTH B1 1 0.0 12.0

[ system ]
pbctopo single-selection topology smoke

[ molecules ]
SELECTED 1
UNUSED 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
-g=0-0
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=single_selection_out.pdb
"@
[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a $arg 1> $stdoutLog 2> $stderrLog
    $exitCode = $LASTEXITCODE
    $stdoutText = Get-Content -Raw -LiteralPath $stdoutLog
    $stderrText = Get-Content -Raw -LiteralPath $stderrLog
    [System.IO.File]::WriteAllText(
        $runLog, $stdoutText + [Environment]::NewLine + $stderrText,
        [System.Text.Encoding]::UTF8)
    if ($exitCode -ne 0) {
        throw "TITAN exited with code $exitCode"
    }
}
finally {
    Pop-Location
}

$runText = Get-Content -Raw -LiteralPath $runLog
if ($runText -notmatch 'selection_atoms=2, projected_bonds=1') {
    throw 'Selected source atom IDs were not projected into the topology graph'
}
if ($runText -match 'MolTopology graph unavailable') {
    throw 'Single selected group emitted a topology-unavailable warning'
}
if (($runText | Select-String -Pattern 'Parsing topo' -AllMatches).Matches.Count -ne 1) {
    throw 'Single selected group should parse the common topology exactly once'
}
if ($runText -match 'Loading Common Topology 2') {
    throw 'Single selected group incorrectly attempted to load ligand topology'
}

$debugLog = Join-Path $WorkRoot 'pbctopo_debug.log'
$debugText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'local graph prepared: source=gromacs_topology, selected_atoms=2, topology_bonds=1, projected_topology_bonds=1',
        'pbctopo_topology_bonds_available=1',
        'pbctopo_topology_bond_edges_added=1',
        'pbctopo_topology_source=gromacs_topology'
    )) {
    if ($debugText -notmatch $pattern) {
        throw "Missing expected single-selection topology pattern: $pattern"
    }
}

Write-Host 'PASS pbctopo_topology_single_selection_smoke'
