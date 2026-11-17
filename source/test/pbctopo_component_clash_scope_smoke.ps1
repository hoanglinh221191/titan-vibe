param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_component_clash_scope_smoke')
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

$pdb = Join-Path $WorkRoot 'component_clash_scope_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO component clash scope smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  C1  AAA A   1       1.000   1.000   1.000  1.00  0.00           C
ATOM      2  O1  AAA A   1       2.300   1.000   1.000  1.00  0.00           O
ATOM      3  N1  AAA A   1       2.800   1.000   1.000  1.00  0.00           N
ATOM      4  C2  BBB B   2       6.000   1.000   1.000  1.00  0.00           C
ATOM      5  H2  BBB B   2       1.200   1.000   1.000  1.00  0.00           H
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ SYSTEM ]
1 2 3 4 5
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
AAA 3

[ atoms ]
1 CT 1 AAA C1 1 0.0 12.0
2 CT 1 AAA O1 1 0.0 16.0
3 CT 1 AAA N1 1 0.0 14.0

[ bonds ]
1 2 1
1 3 1

[ moleculetype ]
BBB 3

[ atoms ]
1 CT 1 BBB C2 1 0.0 12.0
2 CT 1 BBB H2 1 0.0 1.0

[ bonds ]
1 2 1

[ system ]
pbctopo component clash scope smoke

[ molecules ]
AAA 1
BBB 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
-pbctopo-steric=filter
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=64
-g=0-0
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=component_clash_scope_out.pdb
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

$debugText = Get-Content -Raw -LiteralPath (Join-Path $WorkRoot 'pbctopo_debug.log')
foreach ($pattern in @(
        'component_graph_used=0',
        'component_graph_ultraclash=1',
        'component_search_frames=1',
        'pbctopo_output_certification_band_recoveries=1',
        # The search winner sits one cell from its partner, so its offsets
        # come back rebased; certifying them against unrebased points once
        # flagged every atom.
        'pbctopo_certification_assignment_point_residuals=0,',
        'pbci_fallback_frames=0'
    )) {
    if ($debugText -notmatch $pattern) {
        throw "Missing expected component clash-scope pattern: $pattern"
    }
}
if ($debugText -notmatch 'component graph fell through to search') {
    throw 'Hydrogen-only graph winner did not trigger exact-search recovery'
}

$report = Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbctopo_frame_report.csv')
if ($report.component_count -ne '2' -or $report.hard_certified -ne '1' -or
    $report.hard_clashes -ne '0') {
    throw 'Certification band did not recover an all-atom-safe image assignment'
}
# The only contact-supported relative image puts H2 on C1, so the recovered
# image has no contact support. Since R2B3 (2654075a) that is weak_ambiguous
# evidence rather than a certified frame.
if ($report.status -ne 'weak_observation' -or
    $report.evidence_state -ne 'weak_ambiguous' -or
    $report.fallback_reason -ne 'evidence_underdetermined') {
    throw "Unsupported recovered image was not reported as weak evidence: status=$($report.status), evidence_state=$($report.evidence_state), fallback_reason=$($report.fallback_reason)"
}

Write-Host 'PASS pbctopo_component_clash_scope_smoke'
