[CmdletBinding()]
param(
    [string]$TitanExe = '',
    [string]$WorkRoot = '',
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 900
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$TestRoot = [IO.Path]::GetFullPath($PSScriptRoot)
$AllowedWorkRootBase = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
if ([string]::IsNullOrWhiteSpace($TitanExe)) {
    $TitanExe = Join-Path $RepoRoot 'titan.exe'
}
if ([string]::IsNullOrWhiteSpace($WorkRoot)) {
    $WorkRoot = Join-Path $PSScriptRoot 'tmp\pbe_tabi_cartesian_retry_runtime'
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    $text = ([string]$Value).Trim()
    if ($text -match '^(?i:nan)$') { return [double]::NaN }
    if ($text -match '^\+(?i:inf|infinity)$') {
        return [double]::PositiveInfinity
    }
    if ($text -match '^-(?i:inf|infinity)$') {
        return [double]::NegativeInfinity
    }
    return [double]::Parse($text, $InvariantCulture)
}

function Test-Finite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Assert-RelativeClose([double]$Actual, [double]$Expected,
                              [double]$Tolerance, [string]$Label) {
    $scale = [Math]::Max(1.0, [Math]::Abs($Expected))
    if ([Math]::Abs($Actual - $Expected) / $scale -gt $Tolerance) {
        throw "$Label changed: actual=$Actual expected=$Expected tolerance=$Tolerance"
    }
}

function Get-PathForComparison([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    $trimChars = [char[]]@([IO.Path]::DirectorySeparatorChar,
        [IO.Path]::AltDirectorySeparatorChar)
    $driveRoot = [IO.Path]::GetPathRoot($full)
    if ($full.Equals($driveRoot, [StringComparison]::OrdinalIgnoreCase)) {
        return $driveRoot
    }
    return $full.TrimEnd($trimChars)
}

function Assert-SafeWorkRoot([string]$Path) {
    $full = Get-PathForComparison $Path
    $repo = Get-PathForComparison $RepoRoot
    $test = Get-PathForComparison $TestRoot
    $allowedBase = Get-PathForComparison $AllowedWorkRootBase
    $driveRoot = Get-PathForComparison ([IO.Path]::GetPathRoot($full))
    if ($full.Equals($repo, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot equal to repository root: $full"
    }
    if ($full.Equals($test, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot equal to test root: $full"
    }
    if ($full.Equals($driveRoot, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot drive root: $full"
    }
    if ($full.Equals($allowedBase, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot equal to the protected test/tmp base: $full"
    }
    $allowedPrefix = $allowedBase +
        [IO.Path]::DirectorySeparatorChar
    if (-not $full.StartsWith($allowedPrefix,
                              [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot outside protected test/tmp base: $full"
    }
    return $full
}

$WorkRoot = Assert-SafeWorkRoot $WorkRoot

function Assert-UnderWorkRoot([string]$Path) {
    $full = Get-PathForComparison $Path
    $rootBase = Get-PathForComparison $WorkRoot
    $root = $rootBase +
        [IO.Path]::DirectorySeparatorChar
    if (-not ($full.Equals($rootBase, [StringComparison]::OrdinalIgnoreCase) -or
            $full.StartsWith($root, [StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside resolved Cartesian retry smoke root: $full"
    }
}

function New-CleanDirectory([string]$Path) {
    Assert-UnderWorkRoot $Path
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-LibRoot {
    $runtimeRoot = $RepoRoot
    if (-not (Test-Path -LiteralPath (Join-Path $runtimeRoot 'lib_ff'))) {
        $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
        if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
            $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
        }
        $runtimeRoot = Split-Path -Parent $commonGitDir
    }
    $candidate = Join-Path $runtimeRoot 'lib_ff'
    if (-not (Test-Path -LiteralPath (Join-Path $candidate 'ffAMBER99SBnonbonded.dat'))) {
        $candidate = Join-Path $runtimeRoot 'example\1\lib_ff'
    }
    if (-not (Test-Path -LiteralPath (Join-Path $candidate 'ffAMBER99SBnonbonded.dat'))) {
        throw "Missing AMBER99SB force-field directory below $runtimeRoot"
    }
    return $candidate
}

function Write-DriverPdb([string]$RunDir) {
    @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    ) | Set-Content -LiteralPath (Join-Path $RunDir 'driver.pdb') -Encoding ASCII
}

function Write-OneAtomPqr([string]$RunDir) {
    $path = Join-Path $RunDir 'one_atom.pqr'
    @(
        'ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.0000  1.5500',
        'END'
    ) | Set-Content -LiteralPath $path -Encoding ASCII
    return $path
}

function Get-BaseArgs(
    [string]$RunDir,
    [string]$PqrPath,
    [string]$Mode,
    [double]$Theta = 0.20,
    [int]$MaxIter = 720,
    [string]$PrimarySurface = 'analytic-ses',
    [int]$SurfacePoints = 96
) {
    $driver = Convert-ToTitanPath (Join-Path $RunDir 'driver.pdb')
    $pqr = Convert-ToTitanPath $PqrPath
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $script:LibRoot).Path.TrimEnd('\\') + '\\')
    $thetaText = $Theta.ToString('G17', $InvariantCulture)
    $lines = @(
        "-f=$driver",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$lib",
        "-path-n=$lib",
        '-pbe=bie',
        "-pbe-test-pqr=$pqr",
        '-pbe-surface=vdw',
        "-pbe-bie-primary-surface=$PrimarySurface",
        '-pbe-bie-ses-grid-spacing=1.25',
        '-pbe-bie-geometry=triangle',
        '-pbe-bie-discretization=node-patch',
        '-pbe-bie-node-patch-operator-audit=on',
        '-pbe-bie-sum=treecode',
        # This selects the established JCP-treecode input route.  The actual
        # operator remains the explicit Cartesian full-FMM engine below.
        '-pbe-bie-fast-operator=p1',
        '-pbe-bie-formulation=jcp2013-treecode',
        '-pbe-bie-derivatives=analytic',
        '-pbe-bie-fmm-m2l=single-matvec-analytic',
        '-pbe-bie-fmm-local=off',
        '-pbe-bie-fmm-engine=cartesian-full-experimental',
        '-pbe-bie-fmm-cartesian-order=1',
        "-pbe-bie-fmm-cartesian-error-audit=$Mode",
        '-pbe-bie-tree-leaf-size=64',
        "-pbe-bie-opening-theta=$thetaText",
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off',
        '-pbe-bie-precond=fmm-nearfield-ilu0',
        '-pbe-bie-precond-ilu-neighbor-cap=64',
        '-pbe-bie-precond-ilu-drop-tol=0',
        '-pbe-bie-precond-ilu-regularization=adaptive',
        '-pbe-bie-precond-ilu-regularization-max=0.01',
        '-pbe-bie-precond-memory-budget-mb=256',
        "-pbe-max-iter=$MaxIter",
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        '-pbe-probe=1.4',
        '-pbe-timing=detail',
        '-pbe-bie-gmres-restart=128',
        '-pbe-bie-true-residual-probe=restart',
        '-pbe-bie-true-residual-probe-scope=sampled',
        '-pbe-bie-true-residual-sample-count=512',
        '-pbe-bie-matvec-audit=sampled-detail',
        '-pbe-bie-matvec-audit-sample-count=512',
        '-pbe-bie-surface-output=off',
        '-pbe-bie-progress=detail',
        '-pbe-bie-progress-interval=1',
        "-pbe-bie-surface-points=$SurfacePoints",
        "-nt=$Threads"
    )
    return ,$lines
}

function Replace-Arg([string[]]$ArgumentLines, [string]$Prefix,
                     [string]$Replacement) {
    $found = $false
    $updated = foreach ($arg in $ArgumentLines) {
        if ($arg.StartsWith($Prefix, [StringComparison]::Ordinal)) {
            $found = $true
            $Replacement
        } else {
            $arg
        }
    }
    if (-not $found) { $updated += $Replacement }
    return ,$updated
}

function Invoke-Case([string]$Name, [string[]]$ArgumentLines,
                     [bool]$ExpectSuccess = $true) {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDirectory $runDir
    Write-DriverPdb $runDir
    $argPath = Join-Path $runDir 'arg.in'
    $ArgumentLines | Set-Content -LiteralPath $argPath -Encoding ASCII
    $psi = New-Object Diagnostics.ProcessStartInfo
    $psi.FileName = $TitanExe
    $psi.WorkingDirectory = $runDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = '-a "' + ($argPath -replace '"', '\"') + '"'
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $psi
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
        $process.Kill($true)
        throw "Cartesian retry smoke timed out: $Name"
    }
    $process.WaitForExit()
    $stdout = $stdoutTask.GetAwaiter().GetResult()
    $stderr = $stderrTask.GetAwaiter().GetResult()
    $exitCode = $process.ExitCode
    $process.Dispose()
    Set-Content -LiteralPath (Join-Path $runDir 'stdout.log') -Value $stdout -Encoding UTF8
    Set-Content -LiteralPath (Join-Path $runDir 'stderr.log') -Value $stderr -Encoding UTF8
    if ($ExpectSuccess -and $exitCode -ne 0) {
        throw "Cartesian retry smoke case ${Name} failed with exit ${exitCode}:`n$stdout`n$stderr"
    }
    if (-not $ExpectSuccess -and $exitCode -eq 0) {
        throw "Cartesian retry invalid scope unexpectedly succeeded: $Name"
    }
    return [pscustomobject]@{
        Name = $Name
        RunDir = $runDir
        ExitCode = $exitCode
        Stdout = $stdout
        Stderr = $stderr
    }
}

function Get-RunProfile($Case) {
    $path = Join-Path $Case.RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing TABI profile: $path" }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -ne 1) { throw "Expected exactly one TABI profile row: $path" }
    return $rows[0]
}

function Get-RunEnergy($Case) {
    $path = Join-Path $Case.RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing PBE energy output: $path" }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -ne 1) { throw "Expected exactly one PBE energy row: $path" }
    return $rows[0]
}

function Assert-ExactCoverage($Profile, [string]$Prefix) {
    $expected = [int64]$Profile."${Prefix}direct_pairs_expected"
    $applied = [int64]$Profile."${Prefix}direct_pairs_applied"
    if ($expected -le 0 -or $applied -ne $expected) {
        throw "Retry telemetry direct coverage is not exact for ${Prefix}: expected=$expected applied=$applied"
    }
}

function Assert-RetryProgress($Case, [int]$Attempts) {
    $canonical = Join-Path $Case.RunDir 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $canonical)) {
        throw "Retry canonical long-run progress is missing: $canonical"
    }
    for ($attempt = 1; $attempt -le $Attempts; ++$attempt) {
        $copy = Join-Path $Case.RunDir "pbe_tabi_longrun_progress_retry_attempt${attempt}.csv"
        if (-not (Test-Path -LiteralPath $copy)) {
            throw "Retry attempt progress is missing: $copy"
        }
        $lines = [IO.File]::ReadAllLines($copy)
        if ($lines.Count -lt 2 -or $lines[0] -notmatch '^elapsed_seconds,phase,') {
            throw "Retry attempt progress is not the canonical schema: $copy"
        }
    }
    $finalCopy = Join-Path $Case.RunDir "pbe_tabi_longrun_progress_retry_attempt${Attempts}.csv"
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $canonical).Hash -ne
        (Get-FileHash -Algorithm SHA256 -LiteralPath $finalCopy).Hash) {
        throw 'Retry canonical progress does not end at the final attempt.'
    }
    $trace = Join-Path $Case.RunDir 'pbe_tabi_fmm_cartesian_retry.csv'
    if (-not (Test-Path -LiteralPath $trace)) { throw "Retry trace is missing: $trace" }
    $traceRows = @(Import-Csv -LiteralPath $trace)
    if ($traceRows.Count -lt $Attempts -or
        $traceRows[0].PSObject.Properties.Name -notcontains 'identity_reason' -or
        $traceRows[0].PSObject.Properties.Name -notcontains 'classification_reason') {
        throw 'Retry trace schema or event stream is incomplete.'
    }
    return $traceRows
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing TITAN executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$script:LibRoot = Get-LibRoot
New-CleanDirectory $WorkRoot

# Parser/validator reject every calibrated-mode scope violation before solve.
$invalidRoot = Join-Path $WorkRoot 'invalid_template'
New-CleanDirectory $invalidRoot
Write-DriverPdb $invalidRoot
$invalidPqr = Write-OneAtomPqr $invalidRoot
$validTemplate = Get-BaseArgs $invalidRoot $invalidPqr 'channel-shadow-retry-global-theta015' 0.20 400 'sphere-patch' 8
$invalidCases = @(
    @{ name = 'jcp_treecode_route_direct'; prefix = '-pbe-bie-fast-operator='; value = '-pbe-bie-fast-operator=direct'; diagnostic = 'Cartesian full FMM requires the P1 JCP treecode' },
    @{ name = 'engine_not_full'; prefix = '-pbe-bie-fmm-engine='; value = '-pbe-bie-fmm-engine=cartesian'; diagnostic = 'fmm-cartesian-order is only valid' },
    @{ name = 'order4'; prefix = '-pbe-bie-fmm-cartesian-order='; value = '-pbe-bie-fmm-cartesian-order=4'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'nodepatch_off'; prefix = '-pbe-bie-discretization='; value = '-pbe-bie-discretization=centroid'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'salt_zero'; prefix = '-pbe-salt='; value = '-pbe-salt=0'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'audit_off'; prefix = '-pbe-bie-node-patch-operator-audit='; value = '-pbe-bie-node-patch-operator-audit=off'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'surface_only'; prefix = '-pbe-bie-surface-only='; value = '-pbe-bie-surface-only=true'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'precond_auto'; prefix = '-pbe-bie-precond='; value = '-pbe-bie-precond=auto'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'precond_diagonal'; prefix = '-pbe-bie-precond='; value = '-pbe-bie-precond=diagonal'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'precond_coarse'; prefix = '-pbe-bie-precond='; value = '-pbe-bie-precond=fmm-nearfield-ilu0-coarse'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'precond_distance_knn'; prefix = '-pbe-bie-precond='; value = '-pbe-bie-precond=distance-knn-ilu0'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'leaf80'; prefix = '-pbe-bie-tree-leaf-size='; value = '-pbe-bie-tree-leaf-size=80'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'theta025'; prefix = '-pbe-bie-opening-theta='; value = '-pbe-bie-opening-theta=0.25'; diagnostic = 'channel-shadow-retry-global-theta015' },
    @{ name = 'p4_not_off'; prefix = '-pbe-bie-fmm-local='; value = '-pbe-bie-fmm-local=single-matvec-p4'; diagnostic = 'Cartesian full FMM is a separate engine' }
)
foreach ($invalid in $invalidCases) {
    $caseArgs = Replace-Arg $validTemplate $invalid.prefix $invalid.value
    $case = Invoke-Case (Join-Path 'invalid' $invalid.name) $caseArgs $false
    $combined = $case.Stdout + "`n" + $case.Stderr
    if ($combined -notmatch $invalid.diagnostic) {
        throw "Retry invalid-scope diagnostic was not explicit for $($invalid.name):`n$combined"
    }
}

# A one-atom node-patch surface is intentionally all-direct at leaf=64.  It
# proves the A1-accept path does not perturb the ordinary channel-shadow solve.
$smallTemplateDir = Join-Path $WorkRoot 'small_template'
New-CleanDirectory $smallTemplateDir
Write-DriverPdb $smallTemplateDir
$smallPqr = Write-OneAtomPqr $smallTemplateDir
$smallShadowArgs = Get-BaseArgs $smallTemplateDir $smallPqr 'channel-shadow' 0.20 400 'sphere-patch' 8
$smallRetryArgs = Get-BaseArgs $smallTemplateDir $smallPqr 'channel-shadow-retry-global-theta015' 0.20 400 'sphere-patch' 8
$smallShadow = Invoke-Case 'small_channel_shadow' $smallShadowArgs
$smallRetry = Invoke-Case 'small_retry_a1_accept' $smallRetryArgs
$smallShadowProfile = Get-RunProfile $smallShadow
$smallRetryProfile = Get-RunProfile $smallRetry
$smallShadowEnergy = Get-RunEnergy $smallShadow
$smallRetryEnergy = Get-RunEnergy $smallRetry
if ([int64]$smallRetryProfile.bie_cartesian_full_fmm_far_pairs -ne 0) {
    throw 'Small retry A1-pass case was not all-direct.'
}
if ($smallRetryProfile.fmm_cartesian_retry_mode -ne 'channel-shadow-retry-global-theta015' -or
    $smallRetryProfile.fmm_cartesian_retry_enabled -ne 'true' -or
    [int]$smallRetryProfile.fmm_cartesian_retry_attempts_started -ne 1 -or
    [int]$smallRetryProfile.fmm_cartesian_retry_attempts_completed -ne 1 -or
    [int]$smallRetryProfile.fmm_cartesian_retry_accepted_attempt -ne 1 -or
    $smallRetryProfile.fmm_cartesian_retry_attempt1_classification -ne 'accept' -or
    $smallRetryProfile.fmm_cartesian_retry_attempt2_classification -ne 'off') {
    throw 'Small retry A1 acceptance telemetry is invalid.'
}
foreach ($field in @('iterations', 'matvec_calls', 'fmm_local_gmres_matvec_calls',
        'bie_cartesian_full_fmm_direct_panel_pairs_applied_total')) {
    if ([int64]$smallShadowProfile.$field -ne [int64]$smallRetryProfile.$field) {
        throw "Small A1 acceptance changed counted solve telemetry for $field"
    }
}
Assert-RelativeClose (Convert-ToDouble $smallRetryEnergy.E_rxn_kJmol) (Convert-ToDouble $smallShadowEnergy.E_rxn_kJmol) 1.0e-12 'small retry energy'
Assert-RelativeClose (Convert-ToDouble $smallRetryProfile.fmm_cartesian_solution_audit_rel) (Convert-ToDouble $smallShadowProfile.fmm_cartesian_solution_audit_rel) 1.0e-12 'small retry solution audit'
Assert-ExactCoverage $smallRetryProfile 'fmm_cartesian_retry_attempt1_solution_audit_'
$smallTrace = Assert-RetryProgress $smallRetry 1
if ($smallTrace.event -notcontains 'attempt_accepted') {
    throw 'Small retry A1-pass trace omitted the acceptance event.'
}

# 1SH1 is the accepted natural discriminator: O1 theta=.20 triggers only from
# the post-solve solution audit, then a fresh theta=.15 solve passes.  Compare
# the selected A2 result with an independent ordinary channel-shadow theta=.15 run.
$pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
if (-not (Test-Path -LiteralPath $pqr1sh1)) { throw "Missing 1SH1 fixture: $pqr1sh1" }
$naturalTemplateDir = Join-Path $WorkRoot 'natural_template'
New-CleanDirectory $naturalTemplateDir
Write-DriverPdb $naturalTemplateDir
$retry1sh1Args = Get-BaseArgs $naturalTemplateDir $pqr1sh1 'channel-shadow-retry-global-theta015' 0.20 720 'analytic-ses' 96
$manual15Args = Get-BaseArgs $naturalTemplateDir $pqr1sh1 'channel-shadow' 0.15 720 'analytic-ses' 96
$retry1sh1 = Invoke-Case '1sh1_retry_a1_to_a2' $retry1sh1Args
$manual15 = Invoke-Case '1sh1_manual_theta015_channel_shadow' $manual15Args
$retryProfile = Get-RunProfile $retry1sh1
$manualProfile = Get-RunProfile $manual15
$retryEnergy = Get-RunEnergy $retry1sh1
$manualEnergy = Get-RunEnergy $manual15
if ($retryEnergy.status -ne 'ok' -or $retryProfile.fmm_cartesian_retry_mode -ne 'channel-shadow-retry-global-theta015' -or
    [int]$retryProfile.fmm_cartesian_retry_attempts_started -ne 2 -or
    [int]$retryProfile.fmm_cartesian_retry_attempts_completed -ne 2 -or
    [int]$retryProfile.fmm_cartesian_retry_accepted_attempt -ne 2 -or
    $retryProfile.fmm_cartesian_retry_attempt1_classification -ne 'retry_trigger' -or
    $retryProfile.fmm_cartesian_retry_attempt2_classification -ne 'accept' -or
    $retryProfile.fmm_cartesian_retry_geometry_identity_valid -ne 'true' -or
    $retryProfile.fmm_cartesian_retry_tree_identity_valid -ne 'true' -or
    $retryProfile.fmm_cartesian_retry_policy_transition_valid -ne 'true') {
    throw '1SH1 retry ladder did not reach the required A1->A2 acceptance state.'
}
if ([Math]::Abs((Convert-ToDouble $retryProfile.fmm_cartesian_retry_attempt1_theta) - 0.20) -gt 1e-12 -or
    [Math]::Abs((Convert-ToDouble $retryProfile.fmm_cartesian_retry_attempt2_theta) - 0.15) -gt 1e-12 -or
    (Convert-ToDouble $retryProfile.fmm_cartesian_retry_attempt1_solution_audit_rel) -le 1.0e-4 -or
    (Convert-ToDouble $retryProfile.fmm_cartesian_retry_attempt2_solution_audit_rel) -gt 1.0e-4) {
    throw '1SH1 retry solution-audit gates are inconsistent with the calibrated ladder.'
}
foreach ($prefix in @('fmm_cartesian_retry_attempt1_', 'fmm_cartesian_retry_attempt2_')) {
    Assert-ExactCoverage $retryProfile "${prefix}solution_audit_"
    if ([int64]$retryProfile."${prefix}direct_pairs_expected" -le 0 -or
        [int64]$retryProfile."${prefix}direct_pairs_applied_last" -ne
            [int64]$retryProfile."${prefix}direct_pairs_expected") {
        throw "1SH1 retry full-FMM direct coverage is invalid for $prefix"
    }
    if ((Convert-ToDouble $retryProfile."${prefix}preconditioner_setup_seconds") -le 0.0 -or
        $retryProfile."${prefix}preconditioner_effective" -ne 'fmm_nearfield_ilu0') {
        throw "1SH1 retry did not rebuild a plain ILU0 preconditioner for $prefix"
    }
}
if ([uint64]$retryProfile.fmm_cartesian_retry_attempt1_geometry_hash -eq 0 -or
    [uint64]$retryProfile.fmm_cartesian_retry_attempt1_geometry_hash -ne
        [uint64]$retryProfile.fmm_cartesian_retry_attempt2_geometry_hash -or
    [uint64]$retryProfile.fmm_cartesian_retry_attempt1_tree_hash -ne
        [uint64]$retryProfile.fmm_cartesian_retry_attempt2_tree_hash -or
    [uint64]$retryProfile.fmm_cartesian_retry_attempt1_policy_hash -eq
        [uint64]$retryProfile.fmm_cartesian_retry_attempt2_policy_hash) {
    throw '1SH1 retry identity hashes violate the fixed-geometry/theta-policy contract.'
}
foreach ($field in @('iterations', 'matvec_calls', 'fmm_local_gmres_matvec_calls',
        'bie_cartesian_full_fmm_direct_panel_pairs_applied_total')) {
    if ([int64]$retryProfile.$field -ne [int64]$manualProfile.$field) {
        throw "Accepted retry A2 differs from manual theta=.15 for $field"
    }
}
foreach ($field in @('bie_cartesian_full_fmm_tree_hash',
        'bie_cartesian_full_fmm_worklist_hash',
        'bie_cartesian_full_fmm_policy_hash')) {
    if ([uint64]$retryProfile.$field -ne [uint64]$manualProfile.$field) {
        throw "Accepted retry A2 differs from manual theta=.15 for $field"
    }
}
Assert-RelativeClose (Convert-ToDouble $retryEnergy.E_rxn_kJmol) (Convert-ToDouble $manualEnergy.E_rxn_kJmol) 1.0e-12 'accepted retry/manual theta=.15 energy'
Assert-RelativeClose (Convert-ToDouble $retryProfile.true_residual_rel) (Convert-ToDouble $manualProfile.true_residual_rel) 1.0e-12 'accepted retry/manual true residual'
Assert-RelativeClose (Convert-ToDouble $retryProfile.fmm_cartesian_solution_audit_rel) (Convert-ToDouble $manualProfile.fmm_cartesian_solution_audit_rel) 1.0e-12 'accepted retry/manual solution audit'
$retryTrace = Assert-RetryProgress $retry1sh1 2
foreach ($event in @('retry_triggered', 'cross_attempt_identity', 'attempt_accepted')) {
    if ($retryTrace.event -notcontains $event) { throw "1SH1 retry trace omits $event" }
}

# A valid retry scope with a deliberately insufficient GMRES budget must fail
# closed after A1: no A2, no accepted energy, and retained diagnostics only.
$failureArgs = Get-BaseArgs $naturalTemplateDir $pqr1sh1 'channel-shadow-retry-global-theta015' 0.20 1 'analytic-ses' 96
$failure = Invoke-Case '1sh1_retry_a1_invalid' $failureArgs
$failureProfile = Get-RunProfile $failure
$failureEnergy = Get-RunEnergy $failure
if ($failureEnergy.status -ne 'failed' -or
    $failureEnergy.reaction_energy_source -ne 'not_computed_cartesian_retry_failed' -or
    [int]$failureProfile.fmm_cartesian_retry_attempts_started -ne 1 -or
    [int]$failureProfile.fmm_cartesian_retry_attempts_completed -ne 1 -or
    [int]$failureProfile.fmm_cartesian_retry_accepted_attempt -ne 0 -or
    $failureProfile.fmm_cartesian_retry_attempt1_classification -ne 'fail_closed' -or
    $failureProfile.fmm_cartesian_retry_attempt2_classification -ne 'off') {
    throw 'A1 invalid retry path did not fail closed without launching A2.'
}
if (Test-Finite (Convert-ToDouble $failureEnergy.E_rxn_kJmol)) {
    throw 'A1 invalid retry path leaked an unaccepted energy.'
}
$failureTrace = Assert-RetryProgress $failure 1
if (@($failureTrace | Where-Object { $_.attempt -eq '2' }).Count -ne 0) {
    throw 'A1 invalid retry path launched an unexpected A2.'
}

Write-Host ("PASS Cartesian retry runtime: small A1 accepted; 1SH1 A1=.20 -> A2=.15; invalid scopes={0}" -f $invalidCases.Count)
