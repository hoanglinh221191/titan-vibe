param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_complete_uncertainty_band_smoke')
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

function Invoke-Case([string]$Name) {
    $root = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $root | Out-Null
    $pdb = Join-Path $root 'two_component.pdb'
    $ndx = Join-Path $root 'index.ndx'
    $top = Join-Path $root 'topol.top'
    $arg = Join-Path $root 'arg.in'

    @'
TITLE     PBCTOPO complete uncertainty band smoke
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY A   2      11.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA A   1      15.000  10.000  10.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY A   2      16.000  10.000  10.000  1.00  0.00      MOLB C
ENDMDL
TITLE     PBCTOPO discontact frame t= 1.00000 step= 1
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY A   2      11.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA A   1      30.000  30.000  30.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY A   2      31.000  30.000  30.000  1.00  0.00      MOLB C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

    @'
[ SYSTEM ]
1 2 3 4
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

    @'
[ moleculetype ]
MOLA 3
[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
[ bonds ]
1 2 1

[ moleculetype ]
MOLB 3
[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
[ bonds ]
1 2 1

[ system ]
PBCTOPO complete uncertainty band smoke
[ molecules ]
MOLA 1
MOLB 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

    $argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
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
-opdb=complete_band_out.pdb
"@
    [System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

    Push-Location $root
    try {
        $savedPreference = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        & $titanPath -a $arg 1> stdout.log 2> stderr.log
        $exitCode = $LASTEXITCODE
        $ErrorActionPreference = $savedPreference
        if ($exitCode -ne 0) {
            throw "TITAN exited with code $exitCode"
        }
    }
    finally {
        Pop-Location
    }
    return $root
}

$first = Invoke-Case 'run1'
$second = Invoke-Case 'run2'

$frameRows = @(Import-Csv -LiteralPath (Join-Path $first 'pbctopo_frame_report.csv'))
$candidateRows = @(Import-Csv -LiteralPath (Join-Path $first 'pbctopo_assignment_candidates.csv'))
if ($frameRows.Count -ne 2) {
    throw "Expected two frame reports, found $($frameRows.Count)"
}
$frame = $frameRows[0]
if ($frame.schema_version -ne '8' -or $frame.component_count -ne '2' -or
    $frame.component_feasible_assignment_set_enumerated -ne '1' -or
    $frame.component_feasible_assignments_within_domain -ne '125' -or
    $frame.feasible_assignment_set_enumerated -ne '1' -or
    $frame.feasible_assignments_within_domain -ne '125') {
    throw "Two-component full-domain enumeration was not reported: $($frame | ConvertTo-Json -Compress)"
}
$frame0Candidates = @($candidateRows | Where-Object { $_.frame -eq '0' })
$frame1Candidates = @($candidateRows | Where-Object { $_.frame -eq '1' })
if ($frame0Candidates.Count -ne 125 -or $frame1Candidates.Count -ne 125) {
    throw "Expected 125 certified assignments per frame, found $($frame0Candidates.Count)/$($frame1Candidates.Count)"
}
if (@($candidateRows | Where-Object { $_.credible_alternative_set_complete -ne '1' }).Count -ne 0) {
    throw 'Complete uncertainty band was not marked complete'
}
$target = @($frame0Candidates | Where-Object {
    $_.component_assignment -match '(^|\|)[^|]+:2:2:2($|\|)'
})
if ($target.Count -ne 1) {
    throw "Expected one separated H0 assignment at relative shift (2,2,2), found $($target.Count)"
}
if ([int]$target[0].candidate_rank -le 48) {
    throw "Separated H0 assignment was not below the historical cap: rank=$($target[0].candidate_rank)"
}
$selected = @($frame0Candidates | Where-Object { $_.selected_framewise -eq '1' })
if ($selected.Count -ne 1 -or $selected[0].assignment_signature -eq $target[0].assignment_signature) {
    throw 'Framewise winner was forced to the separated H0 assignment'
}

$discontactFrame = $frameRows[1]
$knownSeparated = @($frame1Candidates | Where-Object {
    $_.component_assignment -match '(^|\|)[^|]+:1:0:0($|\|)'
})
if ($discontactFrame.status -ne 'weak_observation' -or
    $discontactFrame.evidence_state -notin @('weak_disconnected', 'weak_ambiguous') -or
    $discontactFrame.mode -eq 'pbci' -or
    $knownSeparated.Count -ne 1) {
    throw 'No-contact frame was not retained as a weak complete-band observation'
}
$h0Alternatives = @($frame0Candidates | Where-Object {
    $_.evidence_no_support_relation_count -eq '1'
})
if ($h0Alternatives.Count -eq 0 -or
    @($h0Alternatives | Where-Object {
        $_.evidence_compatible_hypothesis_count -ne '0' -or
        $_.evidence_compatible_contact_count -ne '0' -or
        $_.evidence_compatible_hypothesis_hash -ne '00000000000000000000000000000000' -or
        $_.evidence_compatible_pair_hash -ne '00000000000000000000000000000000'
    }).Count -ne 0) {
    throw 'Single-cluster H0 alternatives were not retained with zero support identity'
}
$debugLog = Get-Content -LiteralPath (Join-Path $first 'pbctopo_debug.log') -Raw
if ($debugLog -notmatch 'pbctopo_warm_state_commits=1(?:\D|$)' -or
    $debugLog -notmatch 'pbctopo_weak_observation_warm_preserved=[1-9][0-9]*') {
    throw 'Weak no-contact frame contaminated or bypassed warm-state accounting'
}

$firstPdb = Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $first 'complete_band_out.pdb')
$secondPdb = Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $second 'complete_band_out.pdb')
$secondSelected = @(Import-Csv -LiteralPath (Join-Path $second 'pbctopo_assignment_candidates.csv') |
    Where-Object { $_.frame -eq '0' -and $_.selected_framewise -eq '1' })
if ($firstPdb.Hash -ne $secondPdb.Hash -or $secondSelected.Count -ne 1 -or
    $secondSelected[0].assignment_signature -ne $selected[0].assignment_signature) {
    throw 'Complete-band certification changed deterministic framewise output'
}

Write-Host "PASS pbctopo_complete_uncertainty_band_smoke (2x125; target rank=$($target[0].candidate_rank); weak frame keeps known separated alternative; pdb_sha256=$($firstPdb.Hash))"
