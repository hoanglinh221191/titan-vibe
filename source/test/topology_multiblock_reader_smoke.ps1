param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\topology_multiblock_reader_smoke')
)
# The common-topology reader (Multi_MOL::topo_protein) must expand the system in
# [ molecules ] order whatever the file layout:
#   - two moleculetypes in one itp, plus one that [ molecules ] does not use;
#   - include order different from [ molecules ] order;
#   - a molecule name listed twice in [ molecules ];
#   - an [atoms] header without spaces, a moleculetype without the
#     "; name nrexcl" comment, and [ atoms ] rows without a mass column
#     (mass from [ atomtypes ], as in GROMACS).
# With the positions right, the ENE-12 name check finds no mismatch.

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

$pdb = Join-Path $WorkRoot 'multiblock_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$itpAb = Join-Path $WorkRoot 'ab.itp'
$itpC = Join-Path $WorkRoot 'c.itp'

@'
REMARK    topology multi-block reader smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1 C1   ALA A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2 C2   ALA A   1      11.500  10.000  10.000  1.00  0.00           C
ATOM      3 O1   SER B   2      13.500  10.000  10.000  1.00  0.00           O
ATOM      4 N1   GLY C   3      15.500  10.000  10.000  1.00  0.00           N
ATOM      5 N2   GLY C   3      17.000  10.000  10.000  1.00  0.00           N
ATOM      6 C1   ALA D   4      19.000  10.000  10.000  1.00  0.00           C
ATOM      7 C2   ALA D   4      20.500  10.000  10.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ ALL ]
1 2 3 4 5 6 7
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
; moleculetype UNUSED is defined but not in [ molecules ] (as in martini.itp)
[ moleculetype ]
; name	nrexcl
UNUSED 1

[ atoms ]
1 CT2 1 XXX X1 1 0.0 12.011

[ moleculetype ]
; name	nrexcl
AAA 3

[ atoms ]
1 CT2 1 ALA C1 1  0.10 12.011
2 CT2 1 ALA C2 2 -0.10 12.011

[ bonds ]
1 2 1

[ moleculetype ]
; molname  nrexcl
BBB 3

[ atoms ]
; no mass column: the mass of NH1 comes from [ atomtypes ]
1 NH1 1 GLY N1 1  0.20
2 NH1 1 GLY N2 2 -0.20

[ bonds ]
1 2 1
'@ | Set-Content -LiteralPath $itpAb -Encoding ASCII

@'
[moleculetype]
CCC 3

[atoms]
1 OH1 1 SER O1 1 0.0 15.999
'@ | Set-Content -LiteralPath $itpC -Encoding ASCII

@'
[ atomtypes ]
; name  at.num  mass     charge  ptype  sigma      epsilon
CT2     6       12.011   0.0     A      0.3875     0.2343
NH1     7       14.007   0.0     A      0.3296     0.8368
OH1     8       15.999   0.0     A      0.3154     0.6364

#include "c.itp"
#include "ab.itp"

[ system ]
topology multi-block reader smoke

[ molecules ]
AAA 1
CCC 1
BBB 1
AAA 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

# Two runs: with the index (source atom IDs, direct subset of the expanded
# system) and without it (the whole structure against the expanded system).
function Invoke-Case([string]$label, [bool]$withIndex) {
    $caseDir = Join-Path $WorkRoot $label
    New-Item -ItemType Directory -Force -Path $caseDir | Out-Null
    $indexLine = if ($withIndex) { "-n=$ndx`n" } else { '' }
    $argsText = "-ff=charmm36`n-f=$pdb`n-pbctopo=true`n-g=0-0`n-nt=1`n-1pdb=true`n" +
        "-nocal=true`n$indexLine-common_topo=$top`n-no-dssp=true`n" +
        "-noprint-side_chain_contact=true`n-opdb=multiblock_out.pdb`n"
    $caseArg = Join-Path $caseDir 'arg.in'
    $caseOut = Join-Path $caseDir 'stdout.log'
    $caseErr = Join-Path $caseDir 'stderr.log'
    [System.IO.File]::WriteAllText($caseArg, $argsText, [System.Text.Encoding]::ASCII)
    Push-Location $caseDir
    try {
        & $titanPath -a $caseArg 1> $caseOut 2> $caseErr
        $exitCode = $LASTEXITCODE
    }
    finally {
        Pop-Location
    }
    $runText = (Get-Content -Raw -LiteralPath $caseOut) + "`n" + (Get-Content -Raw -LiteralPath $caseErr)
    if ($exitCode -ne 0) {
        throw "${label}: TITAN exited with code $exitCode"
    }
    foreach ($bad in @('names differ', 'Direct topology subset parse missed',
                       'no included \[ moleculetype \]', 'could not be read')) {
        if ($runText -match $bad) {
            throw "${label}: unexpected '$bad' in the run log"
        }
    }
    if ($runText -notmatch 'selection_atoms=7, projected_bonds=3') {
        throw "${label}: MolTopology graph did not project the 7 atoms and 3 bonds"
    }
    $reportPath = Join-Path $caseDir 'pbctopo_frame_report.csv'
    if (-not (Test-Path -LiteralPath $reportPath)) {
        throw "${label}: the VIBE frame report is missing"
    }
    $report = Import-Csv -LiteralPath $reportPath
    if ($report.component_count -ne '4') {
        throw "${label}: expected 4 components, found $($report.component_count)"
    }
}

Invoke-Case 'with_index' $true
Invoke-Case 'without_index' $false
Write-Host 'PASS topology_multiblock_reader_smoke'
