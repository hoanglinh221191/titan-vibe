param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$SourceRoot = 'D:\data_MD\vibe\Npro\tetramer_WT',
    [string]$WorkRoot = '',
    [string]$BaselineSha256 = '9ECF815F4AB92DC335BB2349F085FF8A74CB35BCA576C0139CAB17CBF2EF9751'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$requiredInputs = @('md1_snap0.pdb', 'md1_1.xtc', 'index.ndx', 'topol.top')
foreach ($inputName in $requiredInputs) {
    $inputPath = Join-Path $SourceRoot $inputName
    if (-not (Test-Path -LiteralPath $inputPath)) {
        throw "Missing Npro pbctopo smoke input: $inputPath"
    }
}

if ([string]::IsNullOrWhiteSpace($WorkRoot)) {
    $WorkRoot = Join-Path $SourceRoot 'test_pbctopo_component_graph_precompute_lazy_smoke'
}
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$argPath = Join-Path $WorkRoot 'arg.in'
$outputName = 'tetramer_pbctopo_component_graph_precompute_lazy_smoke.pdb'
$outputPath = Join-Path $WorkRoot $outputName

$argsText = @"
-ff=charmm36
-s=$SourceRoot\md1_snap0.pdb
-traj=$SourceRoot\md1_1.xtc
-pbctopo=true
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=128
-pbctopo-anchor-backend=auto
-g=0-0
-nt=8
-eu=kCal/mol
-1pdb=true
-nocal=true
-dt=400
-b=0
-e=0
-n=$SourceRoot\index.ndx
-common_topo=$SourceRoot\topol.top
-noprint-side_chain_contact=true
-no-dssp=true
-opdb=$outputName
"@
[System.IO.File]::WriteAllText($argPath, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

if (-not (Test-Path -LiteralPath $outputPath)) {
    throw "Missing pbctopo output PDB: $outputPath"
}
$sha = (Get-FileHash -LiteralPath $outputPath -Algorithm SHA256).Hash.ToUpperInvariant()
if ($sha -ne $BaselineSha256.ToUpperInvariant()) {
    throw "PDB SHA256 changed: got $sha expected $BaselineSha256"
}

$debugLog = Join-Path $WorkRoot 'pbctopo_debug.log'
if (-not (Test-Path -LiteralPath $debugLog)) {
    throw "Missing pbctopo debug log: $debugLog"
}
$logText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'pbci_fallback_frames=0',
        'component_graph_shift_descriptors=\d+',
        'component_graph_shift_materialized=\d+',
        'component_graph_bbox_only_scores=\d+',
        'component_graph_spatial_shortlist_builds=\d+',
        'component_graph_spatial_shortlist_fallbacks=\d+',
        'component_graph_spatial_shortlist_pairs=\d+',
        'component_graph_spatial_shortlist_shift_scores_avoided=\d+',
        'component_graph_lazy_materialize=([0-9.]+) s',
        'component_graph_refine_trial_build_incremental=\d+',
        'component_graph_refine_trial_build_full=\d+',
        'component_graph_refine_trial_build_global_rebase=\d+',
        'component_graph_refine_stats_cached=\d+',
        'component_graph_refine_stats_full=\d+',
        'component_graph_refine_stats=([0-9.]+) s',
        'component_cycle_sync_attempts=\d+',
        'component_cycle_sync_used=\d+',
        'component_cycle_sync_pair_candidates=\d+',
        'component_cycle_sync=([0-9.]+) s',
        'pbctopo_topology_bonds_available=\d+',
        'pbctopo_topology_bond_edges_added=\d+',
        'pbctopo_topology_bond_edges_skipped=\d+',
        'pbctopo_topology_source=(gromacs_topology|gromacs_topology\+geometry_cutoff|geometry_cutoff)',
        'pbctopo_topology_component_source=(topology|topology\+geometry|geometry)'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected pbctopo lazy-graph log pattern: $pattern"
    }
}

$report = Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbctopo_frame_report.csv')
if ($report.status -ne 'certified' -or $report.mode -ne 'component_graph' -or
    $report.component_count -ne '4' -or $report.hard_certified -ne '1') {
    throw 'Tetramer lazy-graph report contract changed'
}

$incrementalMatch = [regex]::Match(
    $logText, 'component_graph_refine_trial_build_incremental=(\d+)')
if (-not $incrementalMatch.Success -or
    [int]$incrementalMatch.Groups[1].Value -le 0) {
    throw 'component_graph_refine_trial_build_incremental did not record any incremental refine trials'
}

$shortlistBuilds = [regex]::Match(
    $logText, 'component_graph_spatial_shortlist_builds=(\d+)')
$shortlistFallbacks = [regex]::Match(
    $logText, 'component_graph_spatial_shortlist_fallbacks=(\d+)')
$shortlistPairs = [regex]::Match(
    $logText, 'component_graph_spatial_shortlist_pairs=(\d+)')
$shortlistAvoided = [regex]::Match(
    $logText, 'component_graph_spatial_shortlist_shift_scores_avoided=(\d+)')
if (-not $shortlistBuilds.Success -or
    [int]$shortlistBuilds.Groups[1].Value -le 0 -or
    -not $shortlistFallbacks.Success -or
    [int]$shortlistFallbacks.Groups[1].Value -ne 0 -or
    -not $shortlistPairs.Success -or
    [int]$shortlistPairs.Groups[1].Value -le 0 -or
    -not $shortlistAvoided.Success -or
    [int]$shortlistAvoided.Groups[1].Value -le 0) {
    throw 'Spatial pair shortlist was not active on the tetramer graph fixture'
}

Write-Host 'PASS pbctopo_component_graph_lazy_smoke'
