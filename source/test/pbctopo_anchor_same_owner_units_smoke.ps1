param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_anchor_same_owner_units_smoke')
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

$pdb = Join-Path $WorkRoot 'same_owner_units_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO same-owner molecular-unit anchor smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1       8.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      10.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       8.000  13.000  10.000  1.00  0.00           C
ATOM      5  CA  GLY B   2      10.000  13.000  10.000  1.00  0.00           C
ATOM      6  CA  SER B   3      12.000  13.000  10.000  1.00  0.00           C
ATOM      7  CA  ASN C   1      30.000  30.000  30.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ OWNER_AB ]
1 2 3 4 5 6

[ CONTEXT ]
7
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
CHAIN_A 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
CHAIN_B 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
CONTEXT 3

[ atoms ]
1 CT 1 ASN CA 1 0.0 12.0

[ system ]
pbctopo same-owner molecular-unit anchor smoke

[ molecules ]
CHAIN_A 1
CHAIN_B 1
CONTEXT 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
-pbctopo-bond-cutoff=2.50
-pbctopo-contact-cutoff=4.50
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=64
-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=same_owner_units_out.pdb
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
$logText = Get-Content -Raw -LiteralPath $debugLog
if ($logText -notmatch 'pbctopo_anchor_intra_owner_contact_pairs=([1-9]\d*)') {
    throw 'Expected contact anchors between topology molecules in the same owner'
}
if ($logText -notmatch 'pbctopo_topology_source=gromacs_topology') {
    throw 'Expected GROMACS topology-backed molecular units'
}

$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -ne 1 -or
    $rows[0].status -ne 'weak_observation' -or
    [int]$rows[0].identified -ne 0 -or
    [int]$rows[0].hard_certified -ne 1 -or
    [int]$rows[0].hard_clashes -ne 0 -or
    [int]$rows[0].anchor_contact_map_lost -ne 0 -or
    $rows[0].fallback_reason -ne 'evidence_disconnected') {
    throw "Unexpected same-owner molecular-unit report: $($rows | ConvertTo-Json -Compress)"
}

Write-Host 'PASS pbctopo_anchor_same_owner_units_smoke'
