[CmdletBinding()]
param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_ses_patch_ray_6lzg_gate'),
    [string]$ReferenceRoot = '',
    [ValidateRange(1,256)][int]$Threads = 4,
    [ValidateRange(120,3600)][int]$CaseTimeoutSeconds = 1200
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Fixture = Join-Path $RepoRoot 'test\pbe_inputs_public_6lzg\6lzg_complex_charmm_ph7.pqr'
$LibRoot = Join-Path $RepoRoot 'example\1\lib_ff'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
if (-not (Test-Path -LiteralPath $Fixture)) { throw "Missing 6LZG fixture: $Fixture" }
if (-not [string]::IsNullOrWhiteSpace($ReferenceRoot)) {
    $ReferenceRoot = [IO.Path]::GetFullPath($ReferenceRoot)
    foreach ($name in 'pbe_energy.csv','pbe_tabi_profile.csv') {
        if (-not (Test-Path -LiteralPath (Join-Path $ReferenceRoot $name))) {
            throw "Missing 6LZG reference artifact: $(Join-Path $ReferenceRoot $name)"
        }
    }
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$culture = [Globalization.CultureInfo]::InvariantCulture
function F([object]$Value) { return [double]::Parse([string]$Value,$culture) }
function TP([string]$Path) { return ($Path -replace '\\','/') }
$lib = (TP ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\'))) + '/'

function Invoke-SesCase([string]$Name, [string]$Engine) {
    $runDir = Join-Path $WorkRoot $Name
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $pdb = Join-Path $runDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdb
    $argumentLines = @(
        "-f=$(TP $pdb)", "-common_topo=$(TP (Write-PbeDriverTopology $pdb))", '-1pdb=true', '-ff=AMBER99SB', '-g=1&1',
        "-path-a=$lib", "-path-n=$lib", '-pbe=bie',
        "-pbe-test-pqr=$(TP $Fixture)", '-pbe-surface=vdw',
        '-pbe-max-iter=720', '-pbe-rms=1e-4', '-pbe-dphi=1e-4',
        '-pbe-int-diel=1', '-pbe-ext-diel=80', '-pbe-salt=0.15',
        '-pbe-temp=298.15', '-pbe-probe=1.4', '-pbe-dgsolv=false',
        '-pbe-timing=detail', '-pbe-bie-primary-surface=analytic-ses',
        "-pbe-bie-analytic-ses-engine=$Engine",
        '-pbe-bie-ses-grid-spacing=1.25', '-pbe-bie-discretization=node-patch',
        '-pbe-bie-node-patch-operator-audit=off', '-pbe-bie-sum=treecode',
        '-pbe-bie-formulation=jcp2013-treecode',
        '-pbe-bie-fmm-engine=cartesian', '-pbe-bie-tree-leaf-size=64',
        '-pbe-bie-opening-theta=0.45', '-pbe-bie-taylor-order=4',
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off',
        '-pbe-bie-precond=fmm-nearfield-ilu0-coarse',
        '-pbe-bie-precond-block-panel-cap=32',
        '-pbe-bie-precond-overlap-panel-cap=64',
        '-pbe-bie-precond-ilu-neighbor-cap=64',
        '-pbe-bie-precond-ilu-drop-tol=0',
        '-pbe-bie-precond-ilu-regularization=adaptive',
        '-pbe-bie-precond-ilu-regularization-max=0.01',
        '-pbe-bie-precond-memory-budget-mb=1024',
        '-pbe-bie-precond-coarse-aggregate-cap=8',
        '-pbe-bie-precond-coarse-verify=off',
        '-pbe-bie-gmres-restart=128',
        '-pbe-bie-true-residual-probe=restart',
        '-pbe-bie-true-residual-probe-scope=sampled',
        '-pbe-bie-true-residual-sample-count=512',
        '-pbe-bie-matvec-audit=sampled',
        '-pbe-bie-matvec-audit-sample-count=512', "-nt=$Threads"
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argumentLines -Encoding ASCII
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out -or $run.exit_code -ne 0) {
        $stderr = Get-Content -Raw -LiteralPath (Join-Path $runDir 'stderr.log') -ErrorAction SilentlyContinue
        throw "6LZG $Engine run failed (exit=$($run.exit_code), timeout=$($run.timed_out)): $stderr"
    }
    return [pscustomobject]@{ run_dir = $runDir; process = $run }
}

if ([string]::IsNullOrWhiteSpace($ReferenceRoot)) {
    $referenceRun = Invoke-SesCase -Name 'current' -Engine 'current'
    $referenceDir = $referenceRun.run_dir
} else {
    $referenceRun = $null
    $referenceDir = $ReferenceRoot
}
$experimentalRun = Invoke-SesCase -Name 'patch-ray-experimental' -Engine 'patch-ray-experimental'
$experimentalDir = $experimentalRun.run_dir
$referenceEnergy = @(Import-Csv -LiteralPath (Join-Path $referenceDir 'pbe_energy.csv'))[-1]
$referenceProfile = @(Import-Csv -LiteralPath (Join-Path $referenceDir 'pbe_tabi_profile.csv'))[-1]
$energy = @(Import-Csv -LiteralPath (Join-Path $experimentalDir 'pbe_energy.csv'))[-1]
$profile = @(Import-Csv -LiteralPath (Join-Path $experimentalDir 'pbe_tabi_profile.csv'))[-1]
$absoluteKj = [Math]::Abs((F $energy.E_rxn_kJmol) - (F $referenceEnergy.E_rxn_kJmol))
$relative = $absoluteKj / [Math]::Max(1e-12,[Math]::Abs((F $referenceEnergy.E_rxn_kJmol)))
$allowedKj = [Math]::Max(0.2 * 4.184, 0.002 * [Math]::Abs((F $referenceEnergy.E_rxn_kJmol)))
if ($energy.status -ne 'ok' -or (F $profile.true_residual_rel) -gt 1e-4 -or
    $absoluteKj -gt $allowedKj -or
    [int64]$profile.panels -ne [int64]$referenceProfile.panels) {
    throw "6LZG gate failed: status=$($energy.status) residual=$($profile.true_residual_rel) panels=$($profile.panels)/$($referenceProfile.panels) abs_kJ=$absoluteKj allowed_kJ=$allowedKj"
}
[pscustomobject]@{
    status = 'pass'
    reference_energy_kj_mol = F $referenceEnergy.E_rxn_kJmol
    experimental_energy_kj_mol = F $energy.E_rxn_kJmol
    energy_abs_diff_kj_mol = $absoluteKj
    energy_relative_diff = $relative
    panels = [int64]$profile.panels
    reference_iterations = [int64]$referenceProfile.iterations
    experimental_iterations = [int64]$profile.iterations
    true_residual = F $profile.true_residual_rel
    reference_surface_seconds = F $referenceProfile.surface_build_seconds
    experimental_surface_seconds = F $profile.surface_build_seconds
    reference_gmres_seconds = F $referenceProfile.gmres_seconds
    experimental_gmres_seconds = F $profile.gmres_seconds
    reference_wall_seconds = if ($null -eq $referenceRun) { [double]::NaN } else { $referenceRun.process.wall_seconds }
    experimental_wall_seconds = $experimentalRun.process.wall_seconds
    experimental_peak_rss_mb = $experimentalRun.process.peak_working_set_mb
} | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (Join-Path $WorkRoot 'gate.csv')
Import-Csv (Join-Path $WorkRoot 'gate.csv') | Format-List
Write-Host "pbe_tabi_ses_patch_ray_6lzg_gate=pass"
