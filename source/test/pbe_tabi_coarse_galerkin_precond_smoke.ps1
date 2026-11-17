param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_coarse_galerkin_precond_smoke'),
    [switch]$IncludeBudgetFallback
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI coarse Galerkin smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Invoke-CoarseCase(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$SurfacePoints,
    [string]$CoarseBase,
    [int]$MemoryBudgetMb,
    [string]$Preconditioner = 'coarse-galerkin-probe',
    [int]$TreeLeafSize = 8,
    [int]$MaxIter = 30,
    [double]$RmsTolerance = 1.0e-4,
    [ValidateSet('off','full-small')]
    [string]$CoarseVerify = 'full-small'
) {
    New-CleanDir $RunDir
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libRoot = Join-Path $RepoRoot 'lib_ff'
    if (-not (Test-Path -LiteralPath (Join-Path $libRoot 'ffAMBER99SBnonbonded.dat'))) {
        $libRoot = Join-Path $RepoRoot 'example\1\lib_ff'
    }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $libRoot).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        "-common_topo=$(Convert-ToTitanPath (Write-PbeDriverTopology $pdbPath))",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-salt=0.15',
        "-pbe-max-iter=$MaxIter",
        "-pbe-rms=$($RmsTolerance.ToString('G17', [System.Globalization.CultureInfo]::InvariantCulture))",
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-tabi-tree-preset=stage3u-safe',
        '-pbe-tabi-fmm-m2l=single-matvec-analytic',
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-pbe-tabi-tree-leaf-size=$TreeLeafSize",
        "-pbe-tabi-precond=$Preconditioner",
        '-pbe-tabi-precond-coarse-aggregate-cap=8',
        "-pbe-tabi-precond-coarse-base=$CoarseBase",
        "-pbe-tabi-precond-coarse-verify=$CoarseVerify",
        "-pbe-tabi-precond-memory-budget-mb=$MemoryBudgetMb",
        '-pbe-tabi-precond-block-near=buffer-aware',
        '-pbe-tabi-precond-block-panel-cap=32',
        '-pbe-tabi-gmres-restart=8',
        '-pbe-tabi-progress=detail',
        '-pbe-tabi-progress-interval=1',
        '-nt=2'
    ) -Encoding ASCII
    Push-Location $RunDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
    } finally {
        Pop-Location
    }
    return @(Import-Csv -LiteralPath (Join-Path $RunDir 'pbe_timing.csv'))[0]
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot
$spherePqr = Join-Path $PSScriptRoot 'pbe_inputs_local\sphere.pqr'
if (-not (Test-Path -LiteralPath $spherePqr)) {
    $fixtureDir = Join-Path $WorkRoot 'fixtures'
    New-CleanDir $fixtureDir
    $spherePqr = Join-Path $fixtureDir 'sphere.pqr'
    Set-Content -LiteralPath $spherePqr -Value @(
        'ATOM      1 CA   GLY A   1      0.0000  0.0000  0.0000  1.0000  2.0000',
        'END'
    ) -Encoding ASCII
}

$activeDir = Join-Path $WorkRoot 'active'
$timing = Invoke-CoarseCase -RunDir $activeDir -PqrPath $spherePqr -SurfacePoints 48 -CoarseBase 'diagonal' -MemoryBudgetMb 256
foreach ($col in @(
    'tabi_preconditioner_coarse_aggregate_cap',
    'tabi_preconditioner_coarse_aggregates',
    'tabi_preconditioner_coarse_dimension',
    'tabi_preconditioner_coarse_base',
    'tabi_preconditioner_coarse_matrix_bytes',
    'tabi_preconditioner_coarse_pivot_ratio',
    'tabi_preconditioner_coarse_status',
    'tabi_preconditioner_coarse_setup_matvec_calls',
    'tabi_preconditioner_coarse_setup_seconds',
    'tabi_preconditioner_coarse_moments_setup_seconds',
    'tabi_preconditioner_coarse_far_setup_seconds',
    'tabi_preconditioner_coarse_audit_setup_seconds',
    'tabi_preconditioner_coarse_apply_seconds',
    'tabi_preconditioner_coarse_apply_calls',
    'tabi_preconditioner_coarse_rp_identity_error_max',
    'tabi_preconditioner_coarse_verify_rel_error_max'
)) {
    if (-not ($timing.PSObject.Properties.Name -contains $col)) {
        throw "Missing timing column $col"
    }
}
if ([string]$timing.tabi_preconditioner_requested -ne 'coarse-galerkin-probe') {
    throw 'Expected coarse Galerkin preconditioner request.'
}
if ([string]$timing.tabi_preconditioner_effective -ne 'coarse_galerkin_probe_diagonal') {
    throw 'Expected active diagonal-base coarse Galerkin preconditioner.'
}
if ([string]$timing.tabi_preconditioner_coarse_status -ne 'active') {
    throw 'Expected active coarse status.'
}
if ([int]$timing.tabi_preconditioner_coarse_aggregates -le 0 -or
    [int]$timing.tabi_preconditioner_coarse_dimension -le 0) {
    throw 'Expected a non-empty coarse space.'
}
if ([int]$timing.tabi_preconditioner_coarse_dimension -gt 16) {
    throw 'Coarse dimension exceeded the configured coupled aggregate cap.'
}
if ([double]$timing.tabi_preconditioner_coarse_matrix_bytes -le 0 -or
    [double]$timing.tabi_preconditioner_coarse_setup_matvec_calls -lt [int]$timing.tabi_preconditioner_coarse_dimension) {
    throw 'Expected exact Galerkin setup work.'
}
if ([double]$timing.tabi_preconditioner_coarse_rp_identity_error_max -gt 1.0e-12 -or
    [double]$timing.tabi_preconditioner_coarse_verify_rel_error_max -gt 1.0e-12) {
    throw 'Expected exact coarse identity/Galerkin verification.'
}
if ([double]$timing.tabi_preconditioner_coarse_apply_calls -le 0) {
    throw 'Expected coarse preconditioner applications.'
}

$balancedDir = Join-Path $WorkRoot 'balanced_ilu0'
$balanced = Invoke-CoarseCase -RunDir $balancedDir -PqrPath $spherePqr `
    -SurfacePoints 48 -CoarseBase 'diagonal' -MemoryBudgetMb 256 `
    -Preconditioner 'fmm-nearfield-ilu0-coarse' -CoarseVerify 'off'
if ([string]$balanced.tabi_preconditioner_requested -ne 'fmm-nearfield-ilu0-coarse') {
    throw 'Expected opt-in balanced ILU0 plus coarse request.'
}
if ([string]$balanced.tabi_preconditioner_effective -notlike
    'fmm_nearfield_ilu0_coarse_galerkin_balanced_active_balanced_*') {
    throw "Expected active balanced ILU0 plus coarse preconditioner, got $($balanced.tabi_preconditioner_effective)."
}
if ([string]$balanced.tabi_preconditioner_coarse_status -notin @(
        'active_balanced_hybrid_monopole',
        'active_balanced_exact_fallback') -or
    [string]$balanced.tabi_preconditioner_coarse_base -ne 'fmm-nearfield-ilu0') {
    throw 'Expected active balanced coarse correction with near-field ILU0 base.'
}
if ([int]$balanced.tabi_preconditioner_coarse_apply_calls -le 0) {
    throw 'Expected retained hybrid operator columns and balanced applications.'
}
if ([string]$balanced.tabi_preconditioner_coarse_status -eq
        'active_balanced_hybrid_monopole' -and
    ([int]$balanced.tabi_preconditioner_coarse_setup_matvec_calls -lt 1 -or
     [int]$balanced.tabi_preconditioner_coarse_setup_matvec_calls -ge
        [int]$balanced.tabi_preconditioner_coarse_dimension)) {
    throw 'Hybrid coarse setup must audit with fewer full matvecs than its dimension.'
}
if ([double]$balanced.tabi_preconditioner_estimated_bytes -le
    [double]$balanced.tabi_preconditioner_coarse_matrix_bytes) {
    throw 'Expected the memory estimate to include ILU and retained A*P columns.'
}
foreach ($col in @(
    'tabi_preconditioner_graph_setup_seconds',
    'tabi_preconditioner_near_assembly_setup_seconds',
    'tabi_preconditioner_ilu_factor_setup_seconds',
    'tabi_preconditioner_coarse_moments_setup_seconds',
    'tabi_preconditioner_coarse_far_setup_seconds',
    'tabi_preconditioner_coarse_audit_setup_seconds'
)) {
    if (-not ($balanced.PSObject.Properties.Name -contains $col) -or
        [double]$balanced.$col -le 0.0) {
        throw "Expected positive setup telemetry in $col"
    }
}
if ([string]$balanced.tabi_preconditioner_graph_cache_reused -ne 'true') {
    throw 'Expected ILU setup to reuse the normalized graph from budget selection.'
}

$exactFallback = Invoke-CoarseCase `
    -RunDir (Join-Path $WorkRoot 'balanced_exact_fallback') `
    -PqrPath $spherePqr -SurfacePoints 48 -CoarseBase 'diagonal' `
    -MemoryBudgetMb 256 -Preconditioner 'fmm-nearfield-ilu0-coarse' `
    -CoarseVerify 'full-small'
if ([string]$exactFallback.tabi_preconditioner_coarse_status -ne
        'active_balanced_exact_fallback' -or
    [int]$exactFallback.tabi_preconditioner_coarse_setup_matvec_calls -le
        [int]$exactFallback.tabi_preconditioner_coarse_dimension) {
    throw 'Full-small must exercise the audited exact Galerkin fallback.'
}

$autoDir = Join-Path $WorkRoot 'auto_unchanged'
$auto = Invoke-CoarseCase -RunDir $autoDir -PqrPath $spherePqr `
    -SurfacePoints 48 -CoarseBase 'diagonal' -MemoryBudgetMb 256 `
    -Preconditioner 'auto'
if ([string]$auto.tabi_preconditioner_requested -ne 'fmm-nearfield-ilu0' -or
    [string]$auto.tabi_preconditioner_effective -notlike 'fmm_nearfield_ilu0*' -or
    [string]$auto.tabi_preconditioner_effective -like '*coarse*' -or
    [string]$auto.tabi_preconditioner_coarse_status -ne 'off') {
    throw 'Auto must retain the existing single-level near-field ILU0 pipeline.'
}

$progress = @(Import-Csv -LiteralPath (Join-Path $activeDir 'pbe_tabi_longrun_progress.csv'))[-1]
foreach ($col in @(
    'preconditioner_coarse_status',
    'preconditioner_coarse_dimension',
    'preconditioner_coarse_setup_matvec_calls',
    'preconditioner_graph_setup_seconds',
    'preconditioner_graph_cache_reused',
    'preconditioner_near_assembly_setup_seconds',
    'preconditioner_ilu_factor_setup_seconds',
    'preconditioner_coarse_moments_setup_seconds',
    'preconditioner_coarse_far_setup_seconds',
    'preconditioner_coarse_audit_setup_seconds',
    'preconditioner_coarse_apply_calls'
)) {
    if (-not ($progress.PSObject.Properties.Name -contains $col)) {
        throw "Missing progress column $col"
    }
}

if ($IncludeBudgetFallback) {
    $compactPqr = Join-Path $PSScriptRoot 'pbe_inputs_local\compact-mini-protein.pqr'
    if (-not (Test-Path -LiteralPath $compactPqr)) { throw "Missing compact PQR fixture: $compactPqr" }
    $fallback = Invoke-CoarseCase -RunDir (Join-Path $WorkRoot 'budget_fallback') -PqrPath $compactPqr -SurfacePoints 256 -CoarseBase 'block-budgeted' -MemoryBudgetMb 16 -TreeLeafSize 128 -MaxIter 1 -RmsTolerance 2.0
    if ([string]$fallback.tabi_preconditioner_coarse_status -ne 'budget_fallback' -or
        [string]$fallback.tabi_preconditioner_effective -ne 'diagonal_coarse_budget_cap') {
        throw 'Expected coarse plus block base to fallback before allocation when over budget.'
    }
}

Write-Host "PBE TABI coarse Galerkin preconditioner smoke PASS: $WorkRoot"
