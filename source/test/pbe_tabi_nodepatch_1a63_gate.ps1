param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_nodepatch_1a63_gate'),
    [int]$IterationGate = 22,
    [int]$CaseTimeoutSeconds = 600
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$pqr = (Resolve-Path -LiteralPath (Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

function New-Case(
    [string]$Name,
    [string]$Discretization,
    [string]$OperatorAudit = 'on'
) {
    $runDir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $arg = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
        -PqrPath $pqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -Threads 4 `
        -MaxIter 720 -Restart 128 -Tolerance 1.0e-4 -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -ProbeRadius 1.4 -OpeningTheta 0.45 `
        -AuditMode sampled -AuditSampleCount 512 `
        -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -PrimarySurface grid-ses -SesGridSpacing 1.25 -DgSolv $true
    if ($Discretization -eq 'node-patch') {
        Add-Content -LiteralPath $arg -Encoding ASCII -Value @(
            '-pbe-bie-discretization=node-patch',
            "-pbe-bie-node-patch-operator-audit=$OperatorAudit"
        )
    }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0 -or $run.timed_out) {
        throw "$Name execution failed: exit=$($run.exit_code) timeout=$($run.timed_out)"
    }
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[-1]
    return [pscustomobject]@{
        name = $Name
        run = $run
        profile = $profile
        energy = $energy
    }
}

$centroid = New-Case 'centroid_default' 'centroid'
$node = New-Case 'node_patch_opt_in' 'node-patch'
$nodeNoAudit = New-Case 'node_patch_audit_off' 'node-patch' 'off'

if ($centroid.profile.status -ne 'ok' -or $node.profile.status -ne 'ok' -or
    $nodeNoAudit.profile.status -ne 'ok') {
    throw "1A63 status failure: centroid=$($centroid.profile.status) node=$($node.profile.status) no_audit=$($nodeNoAudit.profile.status)"
}
if ($centroid.profile.discretization -ne 'centroid') {
    throw "Default discretization drifted from centroid: $($centroid.profile.discretization)"
}
if ($node.profile.discretization -ne 'node_patch') {
    throw "Node-patch opt-in was not active: $($node.profile.discretization)"
}
if ([int]$node.profile.iterations -gt $IterationGate) {
    throw "Node-patch iteration gate failed: $($node.profile.iterations) > $IterationGate"
}
if ([double]$node.profile.true_residual_rel -gt 1.0e-4) {
    throw "Node-patch true residual gate failed: $($node.profile.true_residual_rel)"
}
if ([int64]$node.profile.matvec_audit_rows -lt 128 -or
    [double]$node.profile.tree_vs_direct_matvec_rel -gt 2.0e-2) {
    throw "Node-patch pre-GMRES operator audit failed: rows=$($node.profile.matvec_audit_rows) rel=$($node.profile.tree_vs_direct_matvec_rel)"
}
if ([int64]$nodeNoAudit.profile.matvec_audit_rows -ne 0) {
    throw "Node-patch audit-off still evaluated direct rows: $($nodeNoAudit.profile.matvec_audit_rows)"
}
if ([int]$nodeNoAudit.profile.iterations -ne [int]$node.profile.iterations -or
    [double]$nodeNoAudit.profile.true_residual_rel -gt 1.0e-4) {
    throw "Node-patch audit-off solve drifted: iterations=$($nodeNoAudit.profile.iterations) residual=$($nodeNoAudit.profile.true_residual_rel)"
}
$auditOffEnergyRel = [math]::Abs(
    ([double]$nodeNoAudit.energy.E_rxn_kJmol - [double]$node.energy.E_rxn_kJmol) /
    [double]$node.energy.E_rxn_kJmol)
if ($auditOffEnergyRel -gt 1.0e-12) {
    throw "Node-patch audit-off changed the converged energy: $auditOffEnergyRel"
}
$auditOffSummary = Get-Content -LiteralPath (
    Join-Path (Join-Path $WorkRoot 'node_patch_audit_off') 'pbe_summary.dat') -Raw
if ($auditOffSummary -notmatch '(?m)^tabi_node_patch_operator_audit off\r?$') {
    throw 'Node-patch audit-off was not recorded in pbe_summary.dat.'
}
if ($node.profile.near_kernel_policy -ne 'packed-avx2-node-patch-self-flat-pv') {
    throw "Node-patch AVX2/self policy drifted: $($node.profile.near_kernel_policy)"
}
if ([double]$node.profile.energy_projection_seconds -le 0.0 -or
    $node.profile.energy_projection_kernel -notmatch '^dynamic-atom-tile4-(avx2-exact-exp|scalar)$') {
    throw "Node-patch energy projection telemetry mismatch: seconds=$($node.profile.energy_projection_seconds) kernel=$($node.profile.energy_projection_kernel)"
}
if ([int64]$node.profile.unknown_count -ge [int64]$centroid.profile.unknown_count) {
    throw 'Node-patch did not reduce the coupled unknown count.'
}
$energyRel = [math]::Abs(
    ([double]$node.energy.E_rxn_kJmol - [double]$centroid.energy.E_rxn_kJmol) /
    [double]$centroid.energy.E_rxn_kJmol)
if ($energyRel -gt 8.0e-2) {
    throw "Node-patch energy drift exceeded 8%: $energyRel"
}
if ([double]$node.profile.total_seconds -ge [double]$centroid.profile.total_seconds) {
    throw "Node-patch was not faster: node=$($node.profile.total_seconds) centroid=$($centroid.profile.total_seconds)"
}

$summary = @(
    [pscustomobject]@{
        mode='centroid'; panels=$centroid.profile.panels
        unknowns=$centroid.profile.unknown_count; iterations=$centroid.profile.iterations
        true_residual=$centroid.profile.true_residual_rel
        operator_audit_rel=$centroid.profile.tree_vs_direct_matvec_rel
        total_seconds=$centroid.profile.total_seconds
        peak_rss_mb=$centroid.run.peak_working_set_mb
        energy_kj_mol=$centroid.energy.E_rxn_kJmol
    },
    [pscustomobject]@{
        mode='node-patch'; panels=$node.profile.panels
        unknowns=$node.profile.unknown_count; iterations=$node.profile.iterations
        true_residual=$node.profile.true_residual_rel
        operator_audit_rel=$node.profile.tree_vs_direct_matvec_rel
        total_seconds=$node.profile.total_seconds
        peak_rss_mb=$node.run.peak_working_set_mb
        energy_kj_mol=$node.energy.E_rxn_kJmol
    },
    [pscustomobject]@{
        mode='node-patch-audit-off'; panels=$nodeNoAudit.profile.panels
        unknowns=$nodeNoAudit.profile.unknown_count; iterations=$nodeNoAudit.profile.iterations
        true_residual=$nodeNoAudit.profile.true_residual_rel
        operator_audit_rel=$nodeNoAudit.profile.tree_vs_direct_matvec_rel
        total_seconds=$nodeNoAudit.profile.total_seconds
        peak_rss_mb=$nodeNoAudit.run.peak_working_set_mb
        energy_kj_mol=$nodeNoAudit.energy.E_rxn_kJmol
    }
)
$summary | Export-Csv -LiteralPath (Join-Path $WorkRoot 'nodepatch_1a63_gate.csv') -NoTypeInformation
Write-Host ("PBE TABI node-patch 1A63 gate PASS: iterations={0} audit={1} rows={2} node_s={3} no_audit_s={4} centroid_s={5}" -f `
    $node.profile.iterations, $node.profile.tree_vs_direct_matvec_rel,
    $node.profile.matvec_audit_rows, $node.profile.total_seconds,
    $nodeNoAudit.profile.total_seconds, $centroid.profile.total_seconds)
