param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_bie_cartesian_full_fmm_1a63'),
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 1800,
    [ValidateRange(0.25, 4.0)][double]$SesGridSpacing = 1.25,
    [ValidateRange(0.05, 1.0)][double]$FullTheta = 0.65,
    [ValidateRange(4, 2048)][int]$FullLeafSize = 80,
    [switch]$AutomaticFullTree,
    [switch]$OnlyFull,
    [switch]$IncludeOrder1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$Fixture = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_charmm_ph7.pqr'
$ExpectedFixtureSha256 = '5058896fba0c2bfde704afe2ba23d74e2aafcea1f1a95fac67ce643e9264203b'
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing TITAN executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $Fixture)) { throw "Missing 1A63 fixture: $Fixture" }
$fixtureSha256 = (Get-FileHash -LiteralPath $Fixture -Algorithm SHA256).Hash.ToLowerInvariant()
if ($fixtureSha256 -ne $ExpectedFixtureSha256) {
    throw "Unexpected 1A63 CHARMM fixture SHA-256: $fixtureSha256"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'example\1\lib_ff\ffAMBER99SBnonbonded.dat'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
$LibRoot = Join-Path $RuntimeRepoRoot 'example\1\lib_ff'
if (-not (Test-Path -LiteralPath (Join-Path $LibRoot 'ffAMBER99SBnonbonded.dat'))) {
    throw "Missing AMBER99SB force-field tables: $LibRoot"
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}
function Test-Finite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Invoke-1A63Case([string]$Variant, [string]$Engine, [int]$Order,
        [double]$Theta) {
    $runDir = Join-Path $WorkRoot $Variant
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
    $thetaText = $Theta.ToString('G15', $InvariantCulture)
    $spacingText = $SesGridSpacing.ToString('G15', $InvariantCulture)
    $leafSize = if ($Engine -eq 'cartesian-full-experimental') {
        $FullLeafSize
    } else {
        64
    }
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb=true',
        # This force field parses the dummy driver only. The BIE atoms, charges,
        # and radii come entirely from the PDB2PQR/CHARMM fixture below.
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=bie',
        "-pbe-test-pqr=$(Convert-ToTitanPath $Fixture)",
        '-pbe-surface=vdw',
        '-pbe-bie-primary-surface=analytic-ses',
        "-pbe-bie-ses-grid-spacing=$spacingText",
        '-pbe-bie-discretization=node-patch',
        '-pbe-max-iter=720',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        '-pbe-probe=1.4',
        '-pbe-timing=detail',
        '-pbe-bie-fast-operator=p1',
        "-pbe-bie-fmm-engine=$Engine",
        '-pbe-bie-surface-points=8',
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off',
        '-pbe-bie-precond=fmm-nearfield-ilu0-coarse',
        '-pbe-bie-precond-ilu-neighbor-cap=64',
        '-pbe-bie-precond-ilu-drop-tol=0',
        '-pbe-bie-precond-ilu-regularization=adaptive',
        '-pbe-bie-precond-coarse-aggregate-cap=8',
        '-pbe-bie-precond-coarse-base=diagonal',
        '-pbe-bie-precond-coarse-verify=off',
        '-pbe-bie-gmres-restart=128',
        '-pbe-bie-true-residual-probe=restart',
        '-pbe-bie-true-residual-probe-scope=sampled',
        '-pbe-bie-true-residual-sample-count=512',
        '-pbe-bie-matvec-audit=sampled-detail',
        '-pbe-bie-matvec-audit-sample-count=512',
        "-nt=$Threads"
    )
    if ($Engine -eq 'cartesian-full-experimental') {
        $argLines += "-pbe-bie-fmm-cartesian-order=$Order"
    }
    if ($Engine -ne 'cartesian-full-experimental' -or -not $AutomaticFullTree) {
        $argLines += "-pbe-bie-tree-leaf-size=$leafSize"
        $argLines += "-pbe-bie-opening-theta=$thetaText"
    }
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) { throw "1A63 $Variant timed out after $CaseTimeoutSeconds seconds." }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw -LiteralPath (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw -LiteralPath (Join-Path $runDir 'stderr.log'))
        throw "1A63 $Variant failed with exit $($run.exit_code):`n$log"
    }
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    if ($profile.surface_generator -ne 'titan_analytic_patch_ray_ses' -or
        $profile.discretization -ne 'node_patch') {
        throw "1A63 $Variant surface contract failed: generator=$($profile.surface_generator) discretization=$($profile.discretization)"
    }
    $trueResidual = Convert-ToDouble $profile.true_residual_rel
    if ($energy.status -ne 'ok' -or -not (Test-Finite $trueResidual) -or
        $trueResidual -gt 1e-4) {
        throw "1A63 $Variant convergence gate failed: status=$($energy.status) residual=$trueResidual"
    }
    $audit = Convert-ToDouble $profile.tree_vs_direct_matvec_rel
    if ($Engine -eq 'cartesian-full-experimental' -and
        (-not (Test-Finite $audit) -or $audit -gt 1e-4)) {
        throw "1A63 $Variant direct matvec audit failed: $audit"
    }
    return [pscustomobject]@{
        variant = $Variant
        fixture_force_field = 'PDB2PQR-CHARMM'
        fixture_sha256 = $fixtureSha256
        engine = $Engine
        order = $Order
        theta = $Theta
        leaf_size = $leafSize
        surface = $profile.surface_generator
        discretization = $profile.discretization
        ses_grid_spacing_A = $SesGridSpacing
        status = $energy.status
        panels = [int64]$profile.panels
        iterations = [int64]$profile.iterations
        true_residual = $trueResidual
        direct_matvec_relative_l2 = $audit
        energy_kj_mol = Convert-ToDouble $energy.E_rxn_kJmol
        wall_seconds = $run.wall_seconds
        peak_rss_mb = $run.peak_working_set_mb
        prepare_seconds = Convert-ToDouble $profile.fmm_local_prepare_seconds
        gmres_seconds = Convert-ToDouble $profile.gmres_seconds
        surface_build_seconds = Convert-ToDouble $profile.surface_build_seconds
        matvec_seconds_total = Convert-ToDouble $profile.matvec_seconds_total
        matvec_seconds_mean = Convert-ToDouble $profile.matvec_seconds_avg
        nodes = [int64]$profile.bie_cartesian_full_fmm_nodes
        far_pairs = [int64]$profile.bie_cartesian_full_fmm_far_pairs
        direct_pairs = [int64]$profile.bie_cartesian_full_fmm_direct_pairs
        plan_bytes = [int64]$profile.bie_cartesian_full_fmm_plan_bytes
        workspace_bytes = [int64]$profile.bie_cartesian_full_fmm_workspace_bytes
        run_dir = $runDir
    }
}

$rows = New-Object Collections.Generic.List[object]
if (-not $OnlyFull) {
    $rows.Add((Invoke-1A63Case 'production_p1_m2p' 'cartesian' 0 0.6)) | Out-Null
    $rows.Add((Invoke-1A63Case 'production_p1_m2p_theta015' 'cartesian' 0 0.15)) | Out-Null
}
$rows.Add((Invoke-1A63Case 'full_order4' 'cartesian-full-experimental' 4 $FullTheta)) | Out-Null
if ($IncludeOrder1) {
    $rows.Add((Invoke-1A63Case 'full_order1' 'cartesian-full-experimental' 1 0.05)) | Out-Null
}

$matrix = $rows.ToArray()
$full4 = @($matrix | Where-Object variant -eq 'full_order4')[0]
if ($OnlyFull) {
    $baselineReferenceValid = $true
    $energyTolerance = 0.0
    $energyDelta = 0.0
    $speedup = 0.0
    $decision = 'full_only_accuracy_pass'
} else {
    $baseline = @($matrix | Where-Object variant -eq 'production_p1_m2p_theta015')[0]
    $baselineReferenceValid =
        (Test-Finite ([double]$baseline.direct_matvec_relative_l2)) -and
        [double]$baseline.direct_matvec_relative_l2 -le 1e-4
    $energyTolerance = [Math]::Max(0.2 * 4.184,
        0.002 * [Math]::Abs([double]$baseline.energy_kj_mol))
    $energyDelta = [Math]::Abs([double]$full4.energy_kj_mol -
        [double]$baseline.energy_kj_mol)
    $speedup = [double]$baseline.matvec_seconds_mean /
        [Math]::Max(1e-300, [double]$full4.matvec_seconds_mean)
    $decision = if (-not $baselineReferenceValid) {
        'legacy_reference_invalid_keep_experimental'
    } elseif ($energyDelta -gt $energyTolerance) {
        'energy_gate_failed'
    } elseif ($speedup -lt 1.3) {
        'accuracy_pass_keep_experimental_runtime_gap'
    } else {
        'accuracy_and_speed_gate_pass'
    }
}

$matrix | Export-Csv -LiteralPath (Join-Path $WorkRoot 'matrix.csv') `
    -NoTypeInformation -Encoding ASCII
[pscustomobject]@{
    decision = $decision
    energy_delta_kj_mol = $energyDelta
    energy_tolerance_kj_mol = $energyTolerance
    full4_vs_legacy_matvec_speedup = $speedup
    legacy_reference_direct_audit_valid = $baselineReferenceValid
} | Export-Csv -LiteralPath (Join-Path $WorkRoot 'decision.csv') `
    -NoTypeInformation -Encoding ASCII
if ($decision -eq 'energy_gate_failed') {
    throw "1A63 full order 4 energy delta $energyDelta exceeds $energyTolerance kJ/mol."
}
Write-Host "PASS 1A63 BIE Cartesian full FMM decision=$decision speedup=$speedup energy_delta_kJ=$energyDelta"
$matrix | Format-Table -AutoSize
