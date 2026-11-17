[CmdletBinding()]
param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_ses_patch_ray_1a63'),
    [ValidateRange(1,256)][int]$Threads = 4,
    [ValidateRange(0.25,4.0)][double]$SesGridSpacing = 0.75,
    [ValidateRange(60,3600)][int]$CaseTimeoutSeconds = 1200,
    [switch]$SkipLegacy,
    [ValidateSet('diagnostic','enforce')]
    [string]$LegacyEnergyComparison = 'diagnostic',
    [double]$ExternalReferenceEnergyKjMol = [double]::NaN,
    [string]$ExternalReferenceLabel = 'user-supplied',
    [ValidateRange(0.0001,1.0)]
    [double]$ExternalReferenceRelativeTolerance = 0.02
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Fixture = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_charmm_ph7.pqr'
$LibRoot = Join-Path $RepoRoot 'example\1\lib_ff'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
if (-not (Test-Path -LiteralPath $Fixture)) { throw "Missing 1A63 fixture: $Fixture" }
if (-not (Test-Path -LiteralPath (Join-Path $LibRoot 'ffAMBER99SBnonbonded.dat'))) {
    throw "Missing force-field tables: $LibRoot"
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$culture = [Globalization.CultureInfo]::InvariantCulture
function F([object]$Value) { return [double]::Parse([string]$Value,$culture) }
function TP([string]$Path) { return ($Path -replace '\\','/') }

function Invoke-Engine([string]$Engine) {
    $runDir = Join-Path $WorkRoot $Engine
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $pdb = Join-Path $runDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdb
    $spacing = $SesGridSpacing.ToString('G15',$culture)
    $lib = (TP ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\'))) + '/'
    $args = @(
        "-f=$(TP $pdb)", "-common_topo=$(TP (Write-PbeDriverTopology $pdb))", '-1pdb=true', '-ff=AMBER99SB', '-g=1&1',
        "-path-a=$lib", "-path-n=$lib", '-pbe=bie',
        "-pbe-test-pqr=$(TP $Fixture)", '-pbe-surface=vdw',
        '-pbe-bie-primary-surface=analytic-ses',
        "-pbe-bie-analytic-ses-engine=$Engine",
        "-pbe-bie-ses-grid-spacing=$spacing",
        '-pbe-bie-discretization=node-patch',
        '-pbe-bie-fast-operator=p1', '-pbe-bie-fmm-engine=cartesian',
        '-pbe-bie-tree-leaf-size=64', '-pbe-bie-opening-theta=0.6',
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off', '-pbe-bie-precond=auto',
        '-pbe-bie-gmres-restart=128', '-pbe-max-iter=720',
        '-pbe-rms=1e-4', '-pbe-dphi=1e-4', '-pbe-int-diel=1',
        '-pbe-ext-diel=80', '-pbe-salt=0.15', '-pbe-temp=298.15',
        '-pbe-probe=1.4', '-pbe-timing=detail',
        '-pbe-bie-true-residual-probe=restart',
        '-pbe-bie-true-residual-probe-scope=sampled',
        '-pbe-bie-true-residual-sample-count=512',
        '-pbe-bie-matvec-audit=off', "-nt=$Threads"
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $args -Encoding ASCII
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out -or $run.exit_code -ne 0) {
        $stderr = Get-Content -Raw -LiteralPath (Join-Path $runDir 'stderr.log') -ErrorAction SilentlyContinue
        throw "1A63 $Engine failed (exit=$($run.exit_code), timeout=$($run.timed_out)): $stderr"
    }
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[-1]
    $quality = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_surface_quality.csv'))[-1]
    $residual = F $profile.true_residual_rel
    if ($energy.status -ne 'ok' -or $residual -gt 1e-4) {
        throw "1A63 $Engine convergence failed: status=$($energy.status) residual=$residual"
    }
    if ($quality.calderon_trace_space_status -ne 'pass' -or
        [int64]$quality.calderon_trace_boundary_edges -ne 0 -or
        [int64]$quality.calderon_trace_nonmanifold_edges -ne 0) {
        throw "1A63 $Engine surface topology failed: status=$($quality.calderon_trace_space_status) boundary=$($quality.calderon_trace_boundary_edges) nonmanifold=$($quality.calderon_trace_nonmanifold_edges)"
    }
    if ($Engine -eq 'patch-ray-experimental' -and
        $quality.analytic_audit_status -ne 'pass') {
        throw "1A63 $Engine analytic surface audit failed: status=$($quality.analytic_audit_status)"
    }
    [pscustomobject]@{
        engine = $Engine
        status = $energy.status
        panels = [int64]$profile.panels
        mesh_triangles = [int64]$profile.mesh_triangles
        iterations = [int64]$profile.iterations
        true_residual = $residual
        energy_kj_mol = F $energy.E_rxn_kJmol
        surface_area_A2 = F $quality.total_area
        enclosed_volume_A3 = F $quality.enclosed_volume
        topology_status = $quality.calderon_trace_space_status
        boundary_edges = [int64]$quality.calderon_trace_boundary_edges
        nonmanifold_edges = [int64]$quality.calderon_trace_nonmanifold_edges
        analytic_audit_status = $quality.analytic_audit_status
        surface_seconds = F $profile.surface_build_seconds
        gmres_seconds = F $profile.gmres_seconds
        wall_seconds = $run.wall_seconds
        peak_rss_mb = $run.peak_working_set_mb
        run_dir = $runDir
    }
}

$rows = [Collections.Generic.List[object]]::new()
if (-not $SkipLegacy) {
    $rows.Add((Invoke-Engine 'current'))
}
$experimental = Invoke-Engine 'patch-ray-experimental'
$rows.Add($experimental)

$legacyAbsoluteKj = [double]::NaN
$legacyRelative = [double]::NaN
$legacyAllowedKj = [double]::NaN
if (-not $SkipLegacy) {
    $legacy = @($rows | Where-Object engine -eq 'current')[0]
    $legacyAbsoluteKj = [Math]::Abs($experimental.energy_kj_mol - $legacy.energy_kj_mol)
    $legacyRelative = $legacyAbsoluteKj / [Math]::Max(1e-12,[Math]::Abs($legacy.energy_kj_mol))
    $legacyAllowedKj = [Math]::Max(0.2 * 4.184, 0.002 * [Math]::Abs($legacy.energy_kj_mol))
    if ($LegacyEnergyComparison -eq 'enforce' -and $legacyAbsoluteKj -gt $legacyAllowedKj) {
        throw "1A63 legacy energy parity failed: abs_kJ=$legacyAbsoluteKj rel=$legacyRelative allowed_kJ=$legacyAllowedKj"
    }
}

$externalAbsoluteKj = [double]::NaN
$externalRelative = [double]::NaN
if (-not [double]::IsNaN($ExternalReferenceEnergyKjMol)) {
    if ([string]::IsNullOrWhiteSpace($ExternalReferenceLabel)) {
        throw 'ExternalReferenceLabel must identify a supplied energy reference.'
    }
    $externalAbsoluteKj = [Math]::Abs($experimental.energy_kj_mol - $ExternalReferenceEnergyKjMol)
    $externalRelative = $externalAbsoluteKj / [Math]::Max(1e-12,[Math]::Abs($ExternalReferenceEnergyKjMol))
    if ($externalRelative -gt $ExternalReferenceRelativeTolerance) {
        throw "1A63 external-reference energy failed: abs_kJ=$externalAbsoluteKj rel=$externalRelative allowed_rel=$ExternalReferenceRelativeTolerance"
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (Join-Path $WorkRoot 'comparison.csv')
$legacyRole = if ($SkipLegacy) { 'not_run' } else { $LegacyEnergyComparison }
[pscustomobject]@{
    geometry_solver_status = 'pass'
    legacy_comparison_role = $legacyRole
    legacy_energy_abs_diff_kj_mol = $legacyAbsoluteKj
    legacy_energy_relative_diff = $legacyRelative
    legacy_allowed_kj_mol = $legacyAllowedKj
    external_reference_label = if ([double]::IsNaN($ExternalReferenceEnergyKjMol)) { 'not_run' } else { $ExternalReferenceLabel }
    external_reference_energy_kj_mol = $ExternalReferenceEnergyKjMol
    external_energy_abs_diff_kj_mol = $externalAbsoluteKj
    external_energy_relative_diff = $externalRelative
    external_allowed_relative_diff = $ExternalReferenceRelativeTolerance
    status = if ([double]::IsNaN($ExternalReferenceEnergyKjMol)) {
        'pass_geometry_solver_energy_reference_not_run'
    } else {
        'pass_geometry_solver_external_energy_reference'
    }
} | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (Join-Path $WorkRoot 'gate.csv')
$rows | Format-Table engine,status,panels,iterations,true_residual,energy_kj_mol,surface_seconds,gmres_seconds,wall_seconds,peak_rss_mb
Write-Host "pbe_tabi_ses_patch_ray_1a63_gate=pass geometry_solver=pass legacy_role=$legacyRole legacy_energy_relative_diff=$legacyRelative external_reference=$ExternalReferenceLabel external_energy_relative_diff=$externalRelative"
