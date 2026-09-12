param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_interface_scope_integration_smoke')
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

$pdb = Join-Path $WorkRoot 'interface_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$intactContacts = Join-Path $WorkRoot 'intact.contacts'
$impossibleContacts = Join-Path $WorkRoot 'impossible.contacts'

@'
REMARK    PBCTOPO interface certificate scope smoke
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1 A1   RES A   1       9.000  15.000  15.000  1.00  0.00           C
ATOM      2 A2   RES A   1      10.200  15.000  15.000  1.00  0.00           C
ATOM      3 B1   OTH B   2      13.000  15.000  15.000  1.00  0.00           C
ATOM      4 B2   OTH B   2      14.200  15.000  15.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ RECEPTOR ]
1 2

[ LIGAND ]
3 4
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
RECEPTOR 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
LIGAND 3

[ atoms ]
1 CT 1 OTH B1 1 0.0 12.0
2 CT 1 OTH B2 1 0.0 12.0

[ bonds ]
1 2 1

[ system ]
pbctopo interface scope smoke

[ molecules ]
RECEPTOR 1
LIGAND 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

'1 3 4.0 2.5' | Set-Content -LiteralPath $intactContacts -Encoding ASCII
'1 3 0.5 2.5' | Set-Content -LiteralPath $impossibleContacts -Encoding ASCII

function Invoke-Case {
    param(
        [Parameter(Mandatory = $true)][string]$Name,
        [string]$ContactPath
    )
    $runDir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg = Join-Path $runDir 'arg.in'
    $extra = if ([string]::IsNullOrWhiteSpace($ContactPath)) {
        ''
    } else {
        "-pbctopo-required-contacts=$ContactPath`r`n"
    }
    $argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
${extra}-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=$Name.pdb
"@
    [System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)
    $console = Join-Path $runDir 'titan_console.txt'
    $stderr = Join-Path $runDir 'titan_stderr.txt'
    Push-Location $runDir
    try {
        & $titanPath -a $arg 1> $console 2> $stderr
        $exitCode = $LASTEXITCODE
    }
    finally {
        Pop-Location
    }
    if ($exitCode -ne 0) {
        throw "TITAN failed for $Name with code $exitCode"
    }
    return [pscustomobject]@{
        Row = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbctopo_frame_report.csv'))[0]
        Log = Get-Content -Raw -LiteralPath (Join-Path $runDir 'pbctopo_debug.log')
    }
}

$generic = Invoke-Case -Name 'generic'
if ($generic.Row.status -ne 'certified' -or
    $generic.Row.certificate_scope -ne 'image_intercomponent_steric' -or
    $generic.Row.hard_explicit_interface_contacts_checked -ne '0' -or
    $generic.Row.hard_feasible -ne '1' -or
    $generic.Row.identified -ne '1' -or
    $generic.Row.hard_graph_source -ne 'explicit_topology' -or
    $generic.Row.search_evidence_source -ne 'scored_contacts') {
    throw 'Generic scored-contact certificate scope is incorrect'
}
$expectedVirtualExhaustive = if (
    $generic.Row.virtual_seam_enumeration_exhaustive -eq '1' -and
    $generic.Row.virtual_projection_scoring_exhaustive -eq '1' -and
    $generic.Row.virtual_full_atom_verification_exhaustive -eq '1') {
    '1'
} else {
    '0'
}
if ($generic.Row.virtual_search_exhaustive -ne $expectedVirtualExhaustive -or
    ($generic.Row.output_certificate_search_exhaustive -eq '1' -and
     $generic.Row.component_certificate_search_exhaustive -ne '1' -and
     $generic.Row.virtual_full_atom_verification_exhaustive -ne '1')) {
    throw 'Virtual and component certificate-proof layers were conflated'
}

$intact = Invoke-Case -Name 'interface_intact' -ContactPath $intactContacts
if ($intact.Row.status -ne 'certified' -or
    $intact.Row.certificate_scope -ne 'image_intercomponent_steric_interface' -or
    $intact.Row.hard_explicit_interface_contacts_checked -ne '1' -or
    $intact.Row.hard_explicit_interface_contacts_lost -ne '0') {
    throw 'Intact explicit interface was not certified with interface scope'
}
if ($intact.Log -notmatch 'pbctopo_required_interface_contacts_loaded=1') {
    throw 'Explicit interface load telemetry is missing'
}

$impossible = Invoke-Case -Name 'interface_impossible' -ContactPath $impossibleContacts
if ($impossible.Row.status -ne 'fallback' -or
    $impossible.Row.certificate_scope -ne 'none' -or
    $impossible.Row.evidence_state -ne 'not_evaluated' -or
    $impossible.Row.evidence_consistency_evaluated -ne '0' -or
    $impossible.Row.evidence_consistent -ne '0' -or
    $impossible.Row.assignment_signature -ne '0' -or
    $impossible.Row.assignment_component_count -ne '0' -or
    $impossible.Row.committed_candidate_count -ne '0' -or
    $impossible.Row.soft_observed_contact_pair_hash -eq ('0' * 32) -or
    $impossible.Row.soft_observed_lost_pair_hash -ne ('0' * 32) -or
    $impossible.Row.soft_observed_selected_hypothesis_hash -ne ('0' * 32) -or
    $impossible.Row.soft_observed_selected_compatible_pair_hash -ne ('0' * 32) -or
    [int]$impossible.Row.soft_observed_contacts_checked -le 0 -or
    $impossible.Row.soft_observed_contacts_lost -ne '0' -or
    [double]$impossible.Row.soft_observed_contact_cutoff_A -le 0.0 -or
    $impossible.Row.soft_observed_contact_construction_attempted -ne '1' -or
    $impossible.Row.soft_observed_contact_construction_complete -ne '1' -or
    $impossible.Row.soft_observed_contact_cell_setup_failures -ne '0' -or
    $impossible.Row.soft_observed_contact_atom_index_failures -ne '0' -or
    $impossible.Row.soft_observed_contact_mic_query_failures -ne '0') {
    throw 'Impossible explicit interface did not fail closed'
}
if ($impossible.Log -notmatch
    'pbctopo_certification_explicit_interface_contacts_lost=[1-9][0-9]*') {
    throw 'Lost explicit-interface telemetry is missing'
}

Write-Host 'PASS pbctopo_interface_scope_integration_smoke'
