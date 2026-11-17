[CmdletBinding()]
param(
    [string]$TitanExe = '',
    [string]$WorkRoot = '',
    [string[]]$FixtureKeys = @('1sh1', 'barnase', '1a63'),
    [ValidateRange(60, 3600)]
    [int]$CaseTimeoutSeconds = 1800
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest

$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$TestRoot = [IO.Path]::GetFullPath($PSScriptRoot)
$AllowedWorkRootBase = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
$ApprovedTitanSha256 = '6A3D6B605E93F27572FFD4831369F028F568ED0D25F41EF740D842034BA3F90C'
[int64]$ApprovedTitanBytes = 22775296
if ([string]::IsNullOrWhiteSpace($TitanExe)) {
    $TitanExe = Join-Path $RepoRoot 'titan.exe'
}
if ([string]::IsNullOrWhiteSpace($WorkRoot)) {
    $WorkRoot = Join-Path $AllowedWorkRootBase 'p2b2d_cartesian_retry_validation_20260802'
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
        throw "Refusing -WorkRoot equal to protected test/tmp base: $full"
    }
    $prefix = $allowedBase + [IO.Path]::DirectorySeparatorChar
    if (-not $full.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing -WorkRoot outside protected test/tmp base: $full"
    }
    return $full
}

function Assert-ApprovedTitanExe([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not (Test-Path -LiteralPath $full)) {
        throw "Missing TITAN executable: $full"
    }
    $item = Get-Item -LiteralPath $full
    $sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $full).Hash.ToUpperInvariant()
    if ($item.Length -ne $ApprovedTitanBytes -or $sha256 -ne $ApprovedTitanSha256) {
        throw "Refusing unapproved TITAN executable: bytes=$($item.Length) sha256=$sha256 expected_bytes=$ApprovedTitanBytes expected_sha256=$ApprovedTitanSha256"
    }
    return [pscustomobject]@{
        path = $full
        bytes = [int64]$item.Length
        last_write_utc = $item.LastWriteTimeUtc
        sha256 = $sha256
    }
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Get-TextValue([object]$Row, [string]$Name, [string]$Default = '') {
    if ($null -eq $Row) { return $Default }
    $property = $Row.PSObject.Properties[$Name]
    if ($null -eq $property -or $null -eq $property.Value) { return $Default }
    return [string]$property.Value
}

function Get-DoubleValue([object]$Row, [string]$Name,
                         [double]$Default = [double]::NaN) {
    $text = (Get-TextValue $Row $Name '').Trim()
    if ($text -match '^(?i:nan)$') { return [double]::NaN }
    if ($text -match '^\+(?i:inf|infinity)$') { return [double]::PositiveInfinity }
    if ($text -match '^-(?i:inf|infinity)$') { return [double]::NegativeInfinity }
    $parsed = 0.0
    if ([double]::TryParse($text, [Globalization.NumberStyles]::Float,
            $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Int64Value([object]$Row, [string]$Name, [int64]$Default = 0) {
    $parsed = 0L
    if ([int64]::TryParse((Get-TextValue $Row $Name ''), [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-BoolValue([object]$Row, [string]$Name, [bool]$Default = $false) {
    $text = (Get-TextValue $Row $Name '').Trim().ToLowerInvariant()
    if ($text -eq 'true') { return $true }
    if ($text -eq 'false') { return $false }
    return $Default
}

function Test-FiniteValue([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Format-Double([double]$Value) {
    if (-not (Test-FiniteValue $Value)) { return 'NA' }
    return $Value.ToString('G17', $InvariantCulture)
}

function Test-StrictClose([double]$Actual, [double]$Expected,
                           [double]$Tolerance = 1.0e-12) {
    if (-not (Test-FiniteValue $Actual) -or -not (Test-FiniteValue $Expected)) {
        return $false
    }
    $scale = [Math]::Max(1.0, [Math]::Abs($Expected))
    return [Math]::Abs($Actual - $Expected) / $scale -le $Tolerance
}

function Test-ProfileSixSignificantMatch([double]$TraceValue, [double]$ProfileValue) {
    if (-not (Test-FiniteValue $TraceValue) -or -not (Test-FiniteValue $ProfileValue)) {
        return $false
    }
    # pbe_tabi_profile.csv currently exports these snapshot scalars at six significant digits.
    # Compare that exact exported representation rather than inventing precision the profile does not retain.
    return $TraceValue.ToString('G6', $InvariantCulture) -eq
        $ProfileValue.ToString('G6', $InvariantCulture)
}

function Get-TextSha256([string]$Text) {
    $sha = [Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [Text.Encoding]::UTF8.GetBytes($Text)
        return ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-', '').ToLowerInvariant()
    } finally {
        $sha.Dispose()
    }
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

function New-RunDirectory([string]$Path) {
    if (Test-Path -LiteralPath $Path) {
        throw "Refusing to overwrite P2B2D run directory: $Path"
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function New-ArgFile([string]$RunDir, [object]$Fixture, [string]$Mode,
                     [double]$Theta) {
    Write-DriverPdb $RunDir
    $driver = Convert-ToTitanPath (Join-Path $RunDir 'driver.pdb')
    $driverTopology = Convert-ToTitanPath (Write-PbeDriverTopology (Join-Path $RunDir 'driver.pdb'))
    $pqr = Convert-ToTitanPath (Join-Path $RepoRoot $Fixture.pqr_relative)
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $script:LibRoot).Path.TrimEnd('\\') + '\\')
    $thetaText = $Theta.ToString('G17', $InvariantCulture)
    $lines = @(
        "-f=$driver",
        "-common_topo=$driverTopology",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$lib",
        "-path-n=$lib",
        '-pbe=bie',
        "-pbe-test-pqr=$pqr",
        '-pbe-surface=vdw',
        '-pbe-bie-primary-surface=analytic-ses',
        '-pbe-bie-ses-grid-spacing=1.25',
        '-pbe-bie-geometry=triangle',
        '-pbe-bie-discretization=node-patch',
        '-pbe-bie-node-patch-operator-audit=on',
        '-pbe-bie-sum=treecode',
        # This is the established JCP-treecode CLI route only; the operator is
        # explicitly pinned below to Cartesian full FMM, not public P1.
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
        '-pbe-max-iter=720',
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
        '-pbe-bie-surface-points=96',
        '-nt=4'
    )
    $text = $lines -join "`n"
    if ($text -match 'gmres-p4|fmm-cartesian-order=4') {
        throw 'P2B2D generated arguments breach O1-only scope.'
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $lines -Encoding ASCII
    return [pscustomobject]@{ path = $argPath; text = $text }
}

function Invoke-SolverCase([string]$Name, [string]$RunDir, [object]$Arg) {
    $stdout = ''
    $stderr = ''
    $errorText = ''
    $timedOut = $false
    $exitCode = -1
    $peak = 0L
    $watch = [Diagnostics.Stopwatch]::StartNew()
    $process = $null
    try {
        $psi = New-Object Diagnostics.ProcessStartInfo
        $psi.FileName = $TitanExe
        $psi.WorkingDirectory = $RunDir
        $psi.UseShellExecute = $false
        $psi.CreateNoWindow = $true
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.Arguments = '-a "' + ($Arg.path -replace '"', '\"') + '"'
        $process = New-Object Diagnostics.Process
        $process.StartInfo = $psi
        [void]$process.Start()
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
            $timedOut = $true
            $process.Kill($true)
        }
        $process.WaitForExit()
        $stdout = $stdoutTask.GetAwaiter().GetResult()
        $stderr = $stderrTask.GetAwaiter().GetResult()
        $exitCode = $process.ExitCode
        $peak = [int64]$process.PeakWorkingSet64
    } catch {
        $errorText = $_.Exception.ToString()
    } finally {
        $watch.Stop()
        if ($null -ne $process) { $process.Dispose() }
    }
    Set-Content -LiteralPath (Join-Path $RunDir 'stdout.log') -Value $stdout -Encoding UTF8
    Set-Content -LiteralPath (Join-Path $RunDir 'stderr.log') -Value $stderr -Encoding UTF8
    if (-not [string]::IsNullOrWhiteSpace($errorText)) {
        Set-Content -LiteralPath (Join-Path $RunDir 'harness_exception.txt') -Value $errorText -Encoding UTF8
    }
    return [pscustomobject]@{
        name = $Name; run_dir = $RunDir; exit_code = $exitCode; timed_out = $timedOut
        wall_seconds = $watch.Elapsed.TotalSeconds; process_peak_rss_bytes_diagnostic = $peak
        error = $errorText
    }
}

function Get-CsvSingleRow([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    try {
        $rows = @(Import-Csv -LiteralPath $Path)
        if ($rows.Count -ne 1) { return $null }
        return $rows[0]
    } catch {
        return $null
    }
}

function Get-CsvLastRow([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    try {
        $rows = @(Import-Csv -LiteralPath $Path)
        if ($rows.Count -eq 0) { return $null }
        return $rows[-1]
    } catch {
        return $null
    }
}

function Read-RequiredSingleCsvArtifact([string]$Path, [string]$Kind) {
    $result = [ordered]@{
        kind = $Kind
        path = $Path
        exists = $false
        row_count = 0
        row = $null
        error = ''
        ok = $false
    }
    if (-not (Test-Path -LiteralPath $Path)) {
        $result.error = 'missing'
        return [pscustomobject]$result
    }
    $result.exists = $true
    try {
        $rows = @(Import-Csv -LiteralPath $Path)
        $result.row_count = $rows.Count
        if ($rows.Count -ne 1) {
            $result.error = "expected_exactly_one_row_got_$($rows.Count)"
            return [pscustomobject]$result
        }
        $result.row = $rows[0]
        $result.ok = $true
    } catch {
        $result.error = 'import_failed'
    }
    return [pscustomobject]$result
}

function Test-RunExecutionArtifacts([object]$Case, [object]$ProfileArtifact,
                                    [object]$EnergyArtifact, [object]$MemoryArtifact) {
    $profile = $ProfileArtifact.row
    $energy = $EnergyArtifact.row
    $memory = $MemoryArtifact.row
    $energyValue = Get-DoubleValue $energy 'E_rxn_kJmol'
    $energySource = Get-TextValue $energy 'reaction_energy_source'
    $memoryPeak = Get-Int64Value $memory 'peak_working_set_bytes'
    $processOk = $Case.exit_code -eq 0 -and -not $Case.timed_out -and
        [string]::IsNullOrWhiteSpace($Case.error)
    $profileOk = $ProfileArtifact.ok -and (Get-TextValue $profile 'status') -eq 'ok'
    $energyOk = $EnergyArtifact.ok -and (Get-TextValue $energy 'status') -eq 'ok' -and
        (Test-FiniteValue $energyValue) -and -not [string]::IsNullOrWhiteSpace($energySource) -and
        $energySource -ne 'not_computed_cartesian_retry_failed'
    $memoryOk = $MemoryArtifact.ok -and (Get-TextValue $memory 'status') -eq 'ok' -and
        $memoryPeak -gt 0
    $ok = $processOk -and $profileOk -and $energyOk -and $memoryOk
    $reason = if ($ok) { 'ok' } else {
        "process=$processOk profile=$profileOk energy=$energyOk memory=$memoryOk; profile_artifact=$($ProfileArtifact.error); energy_artifact=$($EnergyArtifact.error); memory_artifact=$($MemoryArtifact.error)"
    }
    return [pscustomobject]@{
        ok = $ok; reason = $reason
        process_ok = $processOk; profile_ok = $profileOk; energy_ok = $energyOk; memory_ok = $memoryOk
        energy_value = $energyValue; energy_source = $energySource; memory_peak = $memoryPeak
    }
}

function Test-DirectCoverage([object]$Profile, [string]$Prefix = '') {
    if ($Prefix -eq '') {
        $expected = Get-Int64Value $Profile 'bie_cartesian_full_fmm_direct_panel_pairs_expected'
        $appliedLast = Get-Int64Value $Profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_last'
        $appliedTotal = Get-Int64Value $Profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_total'
    } else {
        $expected = Get-Int64Value $Profile "${Prefix}direct_pairs_expected"
        $appliedLast = Get-Int64Value $Profile "${Prefix}direct_pairs_applied_last"
        $appliedTotal = Get-Int64Value $Profile "${Prefix}direct_pairs_applied_total"
    }
    return [pscustomobject]@{
        expected = $expected; applied_last = $appliedLast; applied_total = $appliedTotal
        ok = $expected -gt 0 -and $appliedLast -eq $expected -and
            $appliedTotal -ge $expected -and ($appliedTotal % $expected) -eq 0
    }
}

function Test-SolutionCoverage([object]$Profile, [string]$Prefix = '') {
    $auditPrefix = if ($Prefix -eq '') { 'fmm_cartesian_' } else { $Prefix }
    $expected = Get-Int64Value $Profile "${auditPrefix}solution_audit_direct_pairs_expected"
    $applied = Get-Int64Value $Profile "${auditPrefix}solution_audit_direct_pairs_applied"
    $targets = Get-Int64Value $Profile "${auditPrefix}solution_audit_targets"
    $rows = Get-Int64Value $Profile "${auditPrefix}solution_audit_rows"
    $preHash = Get-TextValue $Profile "${auditPrefix}pre_audit_target_hash"
    $targetHash = Get-TextValue $Profile "${auditPrefix}solution_audit_target_hash"
    $rowsOk = $targets -gt 0 -and $targets -le ([int64]::MaxValue / 2L) -and
        $rows -eq (2L * $targets)
    return [pscustomobject]@{
        expected = $expected; applied = $applied; targets = $targets; rows = $rows
        pre_hash = $preHash; target_hash = $targetHash
        ok = $expected -gt 0 -and $applied -eq $expected -and $rowsOk -and
            $preHash -ne '' -and $preHash -eq $targetHash
    }
}

function Test-PlainIlu0([object]$Profile, [string]$Prefix = '') {
    $panelField = if ($Prefix -eq '') { 'panels' } else { "${Prefix}panel_count" }
    $panelCount = Get-Int64Value $Profile $panelField
    $requested = Get-TextValue $Profile "${Prefix}preconditioner_requested"
    $effective = Get-TextValue $Profile "${Prefix}preconditioner_effective"
    $quality = Get-TextValue $Profile "${Prefix}preconditioner_quality_status"
    $coarse = Get-TextValue $Profile "${Prefix}preconditioner_coarse_status"
    $selected = Get-DoubleValue $Profile "${Prefix}preconditioner_ilu_regularization_selected"
    $shifted = Get-Int64Value $Profile "${Prefix}preconditioner_ilu_regularization_shifted_rows"
    $failureStage = Get-TextValue $Profile "${Prefix}preconditioner_ilu_unshifted_failure_stage"
    $validPanels = Get-Int64Value $Profile "${Prefix}preconditioner_valid_panels"
    $ill = Get-Int64Value $Profile "${Prefix}preconditioner_ill_conditioned_blocks"
    $rejected = Get-Int64Value $Profile "${Prefix}preconditioner_quality_rejected_blocks"
    return [pscustomobject]@{
        requested = $requested; effective = $effective; quality = $quality; coarse = $coarse
        selected = $selected; shifted = $shifted; failure_stage = $failureStage
        valid_panels = $validPanels; ill = $ill; rejected = $rejected
        ok = $requested -eq 'fmm-nearfield-ilu0' -and
            $effective -eq 'fmm_nearfield_ilu0' -and $quality -eq 'pass' -and
            $coarse -eq 'off' -and (Test-FiniteValue $selected) -and $selected -eq 0.0 -and
            $shifted -eq 0 -and $failureStage -eq 'off' -and
            $panelCount -gt 0 -and $validPanels -eq $panelCount -and
            $ill -eq 0 -and $rejected -eq 0
    }
}

function Test-TopTechnical([object]$Profile, [double]$ExpectedTheta,
                            [string]$ExpectedMode) {
    if ($null -eq $Profile) {
        return [pscustomobject]@{ ok = $false; reason = 'profile_missing' }
    }
    $direct = Test-DirectCoverage $Profile
    $solution = Test-SolutionCoverage $Profile
    $ilu = Test-PlainIlu0 $Profile
    $residual = Get-DoubleValue $Profile 'true_residual_rel'
    $audit = Get-DoubleValue $Profile 'fmm_cartesian_solution_audit_rel'
    $routeOk = (Get-TextValue $Profile 'status') -eq 'ok' -and
        (Get-BoolValue $Profile 'bie_cartesian_full_fmm_active') -and
        (Get-BoolValue $Profile 'fmm_m2l_finite_outputs') -and
        (Get-TextValue $Profile 'bie_cartesian_full_fmm_requested') -eq 'cartesian-full-experimental' -and
        (Get-TextValue $Profile 'bie_cartesian_full_fmm_effective') -eq 'cartesian-full-experimental' -and
        (Get-Int64Value $Profile 'bie_cartesian_full_fmm_order') -eq 1 -and
        (Get-Int64Value $Profile 'bie_cartesian_full_fmm_leaf_size') -eq 64 -and
        (Test-StrictClose (Get-DoubleValue $Profile 'bie_cartesian_full_fmm_opening_theta') $ExpectedTheta) -and
        (Get-TextValue $Profile 'fmm_cartesian_pre_audit_mode') -eq 'combined+channel-shadow' -and
        (Get-BoolValue $Profile 'fmm_cartesian_solution_audit_enabled') -and
        (Get-BoolValue $Profile 'fmm_cartesian_solution_audit_finite') -and
        (Get-BoolValue $Profile 'fmm_cartesian_solution_audit_relative_available') -and
        (Get-TextValue $Profile 'fmm_cartesian_solution_audit_status') -eq 'relative' -and
        -not (Get-BoolValue $Profile 'fallback_used')
    $numericOk = (Test-FiniteValue $residual) -and $residual -ge 0.0 -and
        $residual -le 1.0e-4 -and (Test-FiniteValue $audit) -and $audit -ge 0.0
    $modeOk = if ($ExpectedMode -eq 'channel-shadow') {
        (Get-TextValue $Profile 'fmm_cartesian_retry_mode') -eq 'off'
    } else {
        (Get-TextValue $Profile 'fmm_cartesian_retry_mode') -eq $ExpectedMode
    }
    $ok = $routeOk -and $numericOk -and $direct.ok -and $solution.ok -and $ilu.ok -and $modeOk
    $reason = if ($ok) { 'ok' } else {
        "route=$routeOk numeric=$numericOk direct=$($direct.ok) solution=$($solution.ok) ilu=$($ilu.ok) mode=$modeOk"
    }
    return [pscustomobject]@{ ok = $ok; reason = $reason; direct = $direct; solution = $solution; ilu = $ilu }
}

function Test-RetryAttempt([object]$Profile, [int]$Attempt, [double]$ExpectedTheta) {
    if ($null -eq $Profile) {
        return [pscustomobject]@{ ok = $false; reason = 'profile_missing' }
    }
    $prefix = "fmm_cartesian_retry_attempt${Attempt}_"
    $direct = Test-DirectCoverage $Profile $prefix
    $solution = Test-SolutionCoverage $Profile $prefix
    $ilu = Test-PlainIlu0 $Profile $prefix
    $residual = Get-DoubleValue $Profile "${prefix}true_residual_rel"
    $audit = Get-DoubleValue $Profile "${prefix}solution_audit_rel"
    $routeOk = (Get-TextValue $Profile "${prefix}status") -eq 'ok' -and
        (Get-BoolValue $Profile "${prefix}converged") -and
        (Get-BoolValue $Profile "${prefix}finite") -and
        (Get-BoolValue $Profile "${prefix}cartesian_full_active") -and
        (Get-BoolValue $Profile "${prefix}cartesian_full_finite") -and
        (Get-TextValue $Profile "${prefix}cartesian_full_requested") -eq 'cartesian-full-experimental' -and
        (Get-TextValue $Profile "${prefix}cartesian_full_effective") -eq 'cartesian-full-experimental' -and
        (Get-Int64Value $Profile "${prefix}cartesian_full_order") -eq 1 -and
        (Get-Int64Value $Profile "${prefix}cartesian_full_leaf_size") -eq 64 -and
        (Test-StrictClose (Get-DoubleValue $Profile "${prefix}theta") $ExpectedTheta) -and
        (Get-BoolValue $Profile "${prefix}coupled") -and
        (Get-TextValue $Profile "${prefix}discretization") -eq 'node_patch' -and
        (Get-BoolValue $Profile "${prefix}node_patch_operator_audit_enabled") -and
        (Test-FiniteValue (Get-DoubleValue $Profile "${prefix}salt_molar")) -and
        (Get-DoubleValue $Profile "${prefix}salt_molar") -gt 0.0 -and
        -not (Get-BoolValue $Profile "${prefix}fallback_used")
    $numericOk = (Test-FiniteValue $residual) -and $residual -ge 0.0 -and
        $residual -le 1.0e-4 -and (Test-FiniteValue $audit) -and $audit -ge 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $Profile "${prefix}solution_audit_abs_l2")) -and
        (Get-DoubleValue $Profile "${prefix}solution_audit_abs_l2") -ge 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $Profile "${prefix}solution_audit_direct_l2")) -and
        (Get-DoubleValue $Profile "${prefix}solution_audit_direct_l2") -gt 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $Profile "${prefix}solution_audit_seconds")) -and
        (Get-DoubleValue $Profile "${prefix}solution_audit_seconds") -ge 0.0
    $rebuildOk = (Get-DoubleValue $Profile "${prefix}preconditioner_setup_seconds") -gt 0.0 -and
        (Get-DoubleValue $Profile "${prefix}fmm_prepare_seconds") -gt 0.0 -and
        (Get-DoubleValue $Profile "${prefix}gmres_seconds") -gt 0.0 -and
        (Get-DoubleValue $Profile "${prefix}wall_seconds") -gt 0.0
    $accountedMemory = Get-Int64Value $Profile "${prefix}accounted_memory_bytes"
    $rssBefore = Get-Int64Value $Profile "${prefix}rss_before_bytes"
    $rssAfter = Get-Int64Value $Profile "${prefix}rss_after_bytes"
    $memoryOk = $accountedMemory -gt 0 -and $rssBefore -gt 0 -and $rssAfter -gt 0
    $ok = $routeOk -and $numericOk -and $direct.ok -and $solution.ok -and $ilu.ok -and $rebuildOk -and $memoryOk
    $reason = if ($ok) { 'ok' } else {
        "route=$routeOk numeric=$numericOk direct=$($direct.ok) solution=$($solution.ok) ilu=$($ilu.ok) rebuild=$rebuildOk memory=$memoryOk"
    }
    return [pscustomobject]@{
        ok = $ok; reason = $reason; direct = $direct; solution = $solution; ilu = $ilu
        residual = $residual; audit = $audit; rebuild_ok = $rebuildOk
        accounted_memory_bytes = $accountedMemory; rss_before_bytes = $rssBefore; rss_after_bytes = $rssAfter
    }
}

function Test-TraceAttemptSnapshot([object]$TraceRow, [object]$Profile,
                                   [int]$Attempt, [double]$ExpectedTheta,
                                   [string]$ExpectedClassification) {
    $reasons = [Collections.Generic.List[string]]::new()
    if ($null -eq $TraceRow -or $null -eq $Profile) {
        return [pscustomobject]@{ ok = $false; reason = 'trace_or_profile_missing' }
    }
    $prefix = "fmm_cartesian_retry_attempt${Attempt}_"
    if ((Get-TextValue $TraceRow 'status') -ne 'ok') { $reasons.Add('status') }
    if (-not (Get-BoolValue $TraceRow 'converged')) { $reasons.Add('converged') }
    if ((Get-TextValue $TraceRow 'classification') -ne $ExpectedClassification -or
        (Get-TextValue $Profile "${prefix}classification") -ne $ExpectedClassification) {
        $reasons.Add('classification')
    }
    if ((Get-TextValue $TraceRow 'solution_audit_status') -ne
        (Get-TextValue $Profile "${prefix}solution_audit_status")) {
        $reasons.Add('solution_audit_status')
    }
    $traceTheta = Get-DoubleValue $TraceRow 'theta'
    $profileTheta = Get-DoubleValue $Profile "${prefix}theta"
    if (-not (Test-StrictClose $traceTheta $ExpectedTheta 0.0) -or
        -not (Test-ProfileSixSignificantMatch $traceTheta $profileTheta)) {
        $reasons.Add('theta_mismatch')
    }
    foreach ($field in @('solution_audit_rel', 'true_residual_rel')) {
        $traceValue = Get-DoubleValue $TraceRow $field
        $profileValue = Get-DoubleValue $Profile "$prefix$field"
        if (-not (Test-ProfileSixSignificantMatch $traceValue $profileValue)) {
            $reasons.Add("${field}_profile_mismatch")
        }
    }
    if ((Get-Int64Value $TraceRow 'panels') -ne (Get-Int64Value $Profile "${prefix}panel_count")) {
        $reasons.Add('panels_mismatch')
    }
    foreach ($field in @('geometry_hash', 'tree_hash', 'worklist_hash', 'policy_hash')) {
        $actual = Get-TextValue $TraceRow $field
        $expected = Get-TextValue $Profile "$prefix$field"
        if ([string]::IsNullOrWhiteSpace($actual) -or $actual -eq '0' -or $actual -ne $expected) {
            $reasons.Add("${field}_mismatch")
        }
    }
    return [pscustomobject]@{ ok = $reasons.Count -eq 0; reason = ($reasons -join ';') }
}

function Test-RetryArtifacts([string]$RunDir, [object]$Profile) {
    $reasons = [Collections.Generic.List[string]]::new()
    $canonical = Join-Path $RunDir 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $canonical)) { $reasons.Add('canonical_progress_missing') }
    for ($attempt = 1; $attempt -le 2; ++$attempt) {
        $copy = Join-Path $RunDir "pbe_tabi_longrun_progress_retry_attempt${attempt}.csv"
        if (-not (Test-Path -LiteralPath $copy)) {
            $reasons.Add("attempt${attempt}_progress_missing")
        } else {
            $lines = [IO.File]::ReadAllLines($copy)
            if ($lines.Count -lt 2 -or $lines[0] -notmatch '^elapsed_seconds,phase,') {
                $reasons.Add("attempt${attempt}_progress_schema")
            }
        }
    }
    $finalCopy = Join-Path $RunDir 'pbe_tabi_longrun_progress_retry_attempt2.csv'
    if ((Test-Path -LiteralPath $canonical) -and (Test-Path -LiteralPath $finalCopy) -and
        (Get-FileHash -Algorithm SHA256 -LiteralPath $canonical).Hash -ne
        (Get-FileHash -Algorithm SHA256 -LiteralPath $finalCopy).Hash) {
        $reasons.Add('canonical_progress_not_final_attempt')
    }
    $trace = Join-Path $RunDir 'pbe_tabi_fmm_cartesian_retry.csv'
    if (-not (Test-Path -LiteralPath $trace)) {
        $reasons.Add('retry_trace_missing')
    } else {
        try {
            $traceRows = @(Import-Csv -LiteralPath $trace)
        } catch {
            $traceRows = @()
            $reasons.Add('retry_trace_import_failed')
        }
        if ($traceRows.Count -eq 0) {
            $reasons.Add('retry_trace_empty')
        } else {
            foreach ($row in $traceRows) {
                $attempt = Get-Int64Value $row 'attempt' ([int64]-1)
                if ($attempt -ne 1 -and $attempt -ne 2) {
                    $reasons.Add("trace_attempt_out_of_range_$attempt")
                }
            }
            foreach ($attempt in 1..2) {
                $started = @($traceRows | Where-Object {
                    (Get-TextValue $_ 'event') -eq 'attempt_started' -and (Get-Int64Value $_ 'attempt' ([int64]-1)) -eq $attempt
                })
                $completed = @($traceRows | Where-Object {
                    (Get-TextValue $_ 'event') -eq 'attempt_completed' -and (Get-Int64Value $_ 'attempt' ([int64]-1)) -eq $attempt
                })
                if ($started.Count -ne 1) { $reasons.Add("attempt${attempt}_start_count_$($started.Count)") }
                if ($completed.Count -ne 1) {
                    $reasons.Add("attempt${attempt}_completed_count_$($completed.Count)")
                    continue
                }
                $expectedTheta = if ($attempt -eq 1) { 0.20 } else { 0.15 }
                $expectedClass = if ($attempt -eq 1) { 'retry_trigger' } else { 'accept' }
                $snapshot = Test-TraceAttemptSnapshot $completed[0] $Profile $attempt $expectedTheta $expectedClass
                if (-not $snapshot.ok) { $reasons.Add("attempt${attempt}_snapshot_$($snapshot.reason)") }
            }
            $accepted = @($traceRows | Where-Object { (Get-TextValue $_ 'event') -eq 'attempt_accepted' })
            if ($accepted.Count -ne 1) {
                $reasons.Add("accepted_count_$($accepted.Count)")
            } elseif ((Get-Int64Value $accepted[0] 'attempt' ([int64]-1)) -ne 2) {
                $reasons.Add('accepted_attempt_not_2')
            } else {
                $acceptedSnapshot = Test-TraceAttemptSnapshot $accepted[0] $Profile 2 0.15 'accept'
                if (-not $acceptedSnapshot.ok) { $reasons.Add("accepted_snapshot_$($acceptedSnapshot.reason)") }
            }
            foreach ($requiredEvent in @(
                @{ name = 'retry_triggered'; attempt = 1 },
                @{ name = 'cross_attempt_identity'; attempt = 2 })) {
                $count = @($traceRows | Where-Object {
                    (Get-TextValue $_ 'event') -eq $requiredEvent.name -and
                    (Get-Int64Value $_ 'attempt' ([int64]-1)) -eq $requiredEvent.attempt
                }).Count
                if ($count -ne 1) { $reasons.Add("trace_$($requiredEvent.name)_count_$count") }
            }
        }
    }
    return [pscustomobject]@{ ok = $reasons.Count -eq 0; reason = ($reasons -join ';') }
}

function New-RunRecord([object]$Fixture, [string]$Role, [string]$Mode,
                       [double]$RequestedTheta, [object]$Case,
                       [object]$Arg, [object]$Profile, [object]$Energy,
                       [object]$Memory, [object]$ProfileArtifact,
                       [object]$EnergyArtifact, [object]$MemoryArtifact,
                       [object]$Execution) {
    $effectiveTopTheta = if ($Mode -eq 'channel-shadow-retry-global-theta015') {
        0.15
    } else {
        $RequestedTheta
    }
    $top = Test-TopTechnical $Profile $effectiveTopTheta $Mode
    $direct = if ($null -eq $Profile) { $null } else { Test-DirectCoverage $Profile }
    $solution = if ($null -eq $Profile) { $null } else { Test-SolutionCoverage $Profile }
    $argHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $Arg.path).Hash.ToLowerInvariant()
    $technicalOk = $top.ok -and $Execution.ok
    $technicalReason = if ($technicalOk) { 'ok' } else {
        "top=$($top.reason); execution=$($Execution.reason)"
    }
    return [pscustomobject][ordered]@{
        fixture_key = $Fixture.key
        fixture = $Fixture.name
        role = $Role
        audit_mode = $Mode
        requested_theta = Format-Double $RequestedTheta
        run_dir = $Case.run_dir
        arg_sha256 = $argHash
        arg_logical_sha256 = Get-TextSha256 $Arg.text
        arg_bytes = (Get-Item -LiteralPath $Arg.path).Length
        exit_code = $Case.exit_code
        timed_out = $Case.timed_out
        harness_error = $Case.error
        wall_seconds = Format-Double $Case.wall_seconds
        process_peak_rss_bytes_diagnostic = $Case.process_peak_rss_bytes_diagnostic
        execution_artifacts_ok = $Execution.ok
        execution_artifacts_reason = $Execution.reason
        profile_artifact_exists = $ProfileArtifact.exists
        profile_artifact_row_count = $ProfileArtifact.row_count
        profile_artifact_error = $ProfileArtifact.error
        energy_artifact_exists = $EnergyArtifact.exists
        energy_artifact_row_count = $EnergyArtifact.row_count
        energy_artifact_error = $EnergyArtifact.error
        memory_artifact_exists = $MemoryArtifact.exists
        memory_artifact_row_count = $MemoryArtifact.row_count
        memory_artifact_error = $MemoryArtifact.error
        profile_status = Get-TextValue $Profile 'status' 'missing'
        energy_status = Get-TextValue $Energy 'status' 'missing'
        energy_kj_mol = Format-Double (Get-DoubleValue $Energy 'E_rxn_kJmol')
        energy_source = Get-TextValue $Energy 'reaction_energy_source'
        technical_ok = $technicalOk
        technical_reason = $technicalReason
        panels = Get-Int64Value $Profile 'panels'
        iterations = Get-Int64Value $Profile 'iterations'
        true_residual_rel = Format-Double (Get-DoubleValue $Profile 'true_residual_rel')
        gmres_seconds = Format-Double (Get-DoubleValue $Profile 'gmres_seconds')
        matvec_calls = Get-Int64Value $Profile 'matvec_calls'
        matvec_seconds_total = Format-Double (Get-DoubleValue $Profile 'matvec_seconds_total')
        full_fmm_apply_calls = Get-Int64Value $Profile 'fmm_local_gmres_matvec_calls'
        preconditioner_setup_seconds = Format-Double (Get-DoubleValue $Profile 'preconditioner_setup_seconds')
        preconditioner_ilu_factor_setup_seconds = Format-Double (Get-DoubleValue $Profile 'preconditioner_ilu_factor_setup_seconds')
        preconditioner_apply_calls = Get-Int64Value $Profile 'preconditioner_apply_calls'
        preconditioner_apply_seconds = Format-Double (Get-DoubleValue $Profile 'preconditioner_apply_seconds')
        profile_memory_status = Get-TextValue $Memory 'status' 'missing'
        profile_peak_rss_bytes = Get-Int64Value $Memory 'peak_working_set_bytes'
        profile_rss_before_bytes = Get-Int64Value $Memory 'rss_before_bytes'
        profile_rss_after_bytes = Get-Int64Value $Memory 'rss_after_bytes'
        profile_memory_estimated_bytes_total = Get-Int64Value $Memory 'estimated_bytes_total'
        profile_memory_total_seconds = Format-Double (Get-DoubleValue $Memory 'total_seconds')
        profile_memory_solver_seconds = Format-Double (Get-DoubleValue $Memory 'solver_seconds')
        profile_memory_bytes = Get-Int64Value $Profile 'memory_bytes'
        fmm_prepare_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_local_prepare_seconds')
        direct_seconds = Format-Double (Get-DoubleValue $Profile 'bie_cartesian_full_fmm_direct_seconds')
        m2l_seconds = Format-Double (Get-DoubleValue $Profile 'bie_cartesian_full_fmm_m2l_seconds')
        full_active = Get-TextValue $Profile 'bie_cartesian_full_fmm_active'
        full_requested = Get-TextValue $Profile 'bie_cartesian_full_fmm_requested'
        full_effective = Get-TextValue $Profile 'bie_cartesian_full_fmm_effective'
        full_order = Get-Int64Value $Profile 'bie_cartesian_full_fmm_order'
        full_leaf = Get-Int64Value $Profile 'bie_cartesian_full_fmm_leaf_size'
        full_theta = Format-Double (Get-DoubleValue $Profile 'bie_cartesian_full_fmm_opening_theta')
        tree_hash = Get-TextValue $Profile 'bie_cartesian_full_fmm_tree_hash'
        worklist_hash = Get-TextValue $Profile 'bie_cartesian_full_fmm_worklist_hash'
        policy_hash = Get-TextValue $Profile 'bie_cartesian_full_fmm_policy_hash'
        far_leaf_pairs = Get-Int64Value $Profile 'bie_cartesian_full_fmm_far_pairs'
        direct_leaf_pairs = Get-Int64Value $Profile 'bie_cartesian_full_fmm_direct_pairs'
        direct_pairs_expected = if ($null -eq $direct) { 0 } else { $direct.expected }
        direct_pairs_applied_last = if ($null -eq $direct) { 0 } else { $direct.applied_last }
        direct_pairs_applied_total = if ($null -eq $direct) { 0 } else { $direct.applied_total }
        direct_coverage_ok = if ($null -eq $direct) { $false } else { $direct.ok }
        solution_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_solution_audit_rel')
        solution_abs_l2 = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_solution_audit_abs_l2')
        solution_direct_l2 = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_solution_audit_direct_l2')
        solution_targets = if ($null -eq $solution) { 0 } else { $solution.targets }
        solution_rows = if ($null -eq $solution) { 0 } else { $solution.rows }
        solution_target_hash = if ($null -eq $solution) { '' } else { $solution.target_hash }
        solution_pairs_expected = if ($null -eq $solution) { 0 } else { $solution.expected }
        solution_pairs_applied = if ($null -eq $solution) { 0 } else { $solution.applied }
        solution_coverage_ok = if ($null -eq $solution) { $false } else { $solution.ok }
        channel_pairs_expected = Get-Int64Value $Profile 'fmm_cartesian_channel_direct_pairs_expected'
        channel_pairs_applied = Get-Int64Value $Profile 'fmm_cartesian_channel_direct_pairs_applied'
        channel_k1_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_channel_k1_rel')
        channel_k2_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_channel_k2_rel')
        channel_k3_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_channel_k3_rel')
        channel_k4_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_channel_k4_rel')
        retry_mode = Get-TextValue $Profile 'fmm_cartesian_retry_mode'
        retry_attempts_started = Get-Int64Value $Profile 'fmm_cartesian_retry_attempts_started'
        retry_attempts_completed = Get-Int64Value $Profile 'fmm_cartesian_retry_attempts_completed'
        retry_accepted_attempt = Get-Int64Value $Profile 'fmm_cartesian_retry_accepted_attempt'
        retry_trigger_reason = Get-TextValue $Profile 'fmm_cartesian_retry_trigger_reason'
        retry_final_reason = Get-TextValue $Profile 'fmm_cartesian_retry_final_reason'
        retry_geometry_identity_valid = Get-TextValue $Profile 'fmm_cartesian_retry_geometry_identity_valid'
        retry_tree_identity_valid = Get-TextValue $Profile 'fmm_cartesian_retry_tree_identity_valid'
        retry_policy_transition_valid = Get-TextValue $Profile 'fmm_cartesian_retry_policy_transition_valid'
        retry_cumulative_wall_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_cumulative_wall_seconds')
        retry_attempt1_classification = Get-TextValue $Profile 'fmm_cartesian_retry_attempt1_classification'
        retry_attempt1_theta = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_theta')
        retry_attempt1_true_residual_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_true_residual_rel')
        retry_attempt1_solution_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_solution_audit_rel')
        retry_attempt1_geometry_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt1_geometry_hash'
        retry_attempt1_tree_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt1_tree_hash'
        retry_attempt1_worklist_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt1_worklist_hash'
        retry_attempt1_policy_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt1_policy_hash'
        retry_attempt1_direct_expected = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_direct_pairs_expected'
        retry_attempt1_direct_applied_last = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_direct_pairs_applied_last'
        retry_attempt1_solution_expected = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_solution_audit_direct_pairs_expected'
        retry_attempt1_solution_applied = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_solution_audit_direct_pairs_applied'
        retry_attempt1_preconditioner_setup_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_preconditioner_setup_seconds')
        retry_attempt1_gmres_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_gmres_seconds')
        retry_attempt1_wall_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt1_wall_seconds')
        retry_attempt1_accounted_memory_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_accounted_memory_bytes'
        retry_attempt1_rss_before_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_rss_before_bytes'
        retry_attempt1_rss_after_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt1_rss_after_bytes'
        retry_attempt2_classification = Get-TextValue $Profile 'fmm_cartesian_retry_attempt2_classification'
        retry_attempt2_theta = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_theta')
        retry_attempt2_true_residual_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_true_residual_rel')
        retry_attempt2_solution_rel = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_solution_audit_rel')
        retry_attempt2_geometry_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt2_geometry_hash'
        retry_attempt2_tree_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt2_tree_hash'
        retry_attempt2_worklist_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt2_worklist_hash'
        retry_attempt2_policy_hash = Get-TextValue $Profile 'fmm_cartesian_retry_attempt2_policy_hash'
        retry_attempt2_direct_expected = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_direct_pairs_expected'
        retry_attempt2_direct_applied_last = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_direct_pairs_applied_last'
        retry_attempt2_solution_expected = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_solution_audit_direct_pairs_expected'
        retry_attempt2_solution_applied = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_solution_audit_direct_pairs_applied'
        retry_attempt2_preconditioner_setup_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_preconditioner_setup_seconds')
        retry_attempt2_gmres_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_gmres_seconds')
        retry_attempt2_wall_seconds = Format-Double (Get-DoubleValue $Profile 'fmm_cartesian_retry_attempt2_wall_seconds')
        retry_attempt2_accounted_memory_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_accounted_memory_bytes'
        retry_attempt2_rss_before_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_rss_before_bytes'
        retry_attempt2_rss_after_bytes = Get-Int64Value $Profile 'fmm_cartesian_retry_attempt2_rss_after_bytes'
    }
}

function Add-Gate([string]$Fixture, [string]$Gate, [bool]$Passed,
                  [string]$Detail) {
    $script:GateRows.Add([pscustomobject][ordered]@{
        fixture_key = $Fixture; gate = $Gate; passed = $Passed; detail = $Detail
    }) | Out-Null
}

function Add-Comparison([string]$Fixture, [string]$Comparison,
                        [string]$Field, [string]$Kind, [bool]$Passed,
                        [string]$Left, [string]$Right) {
    $script:ComparisonRows.Add([pscustomobject][ordered]@{
        fixture_key = $Fixture; comparison = $Comparison; field = $Field
        kind = $Kind; passed = $Passed; left = $Left; right = $Right
    }) | Out-Null
    return $Passed
}

function Compare-TextField([string]$Fixture, [string]$Comparison,
                           [string]$Field, [string]$Left, [string]$Right) {
    return Add-Comparison $Fixture $Comparison $Field 'exact' ($Left -eq $Right) $Left $Right
}

function Compare-DoubleField([string]$Fixture, [string]$Comparison,
                             [string]$Field, [double]$Left, [double]$Right) {
    return Add-Comparison $Fixture $Comparison $Field 'relative_1e-12' `
        (Test-StrictClose $Left $Right) (Format-Double $Left) (Format-Double $Right)
}

function Compare-A1Snapshot([object]$Fixture, [object]$A, [object]$B) {
    $fixtureKey = $Fixture.key
    $profileA = $A.profile
    $profileB = $B.profile
    $prefix = 'fmm_cartesian_retry_attempt1_'
    $ok = $true
    $textMap = @(
        @{ top = 'status'; attempt = 'status' },
        @{ top = 'bie_cartesian_full_fmm_active'; attempt = 'cartesian_full_active' },
        @{ top = 'bie_cartesian_full_fmm_requested'; attempt = 'cartesian_full_requested' },
        @{ top = 'bie_cartesian_full_fmm_effective'; attempt = 'cartesian_full_effective' },
        @{ top = 'bie_cartesian_full_fmm_tree_hash'; attempt = 'tree_hash' },
        @{ top = 'bie_cartesian_full_fmm_worklist_hash'; attempt = 'worklist_hash' },
        @{ top = 'bie_cartesian_full_fmm_policy_hash'; attempt = 'policy_hash' },
        @{ top = 'fmm_cartesian_pre_audit_target_hash'; attempt = 'pre_audit_target_hash' },
        @{ top = 'fmm_cartesian_solution_audit_enabled'; attempt = 'solution_audit_enabled' },
        @{ top = 'fmm_cartesian_solution_audit_finite'; attempt = 'solution_audit_finite' },
        @{ top = 'fmm_cartesian_solution_audit_relative_available'; attempt = 'solution_audit_relative_available' },
        @{ top = 'fmm_cartesian_solution_audit_status'; attempt = 'solution_audit_status' },
        @{ top = 'fmm_cartesian_solution_audit_target_hash'; attempt = 'solution_audit_target_hash' },
        @{ top = 'fallback_used'; attempt = 'fallback_used' },
        @{ top = 'preconditioner_requested'; attempt = 'preconditioner_requested' },
        @{ top = 'preconditioner_effective'; attempt = 'preconditioner_effective' },
        @{ top = 'preconditioner_coarse_status'; attempt = 'preconditioner_coarse_status' },
        @{ top = 'preconditioner_ilu_unshifted_failure_stage'; attempt = 'preconditioner_ilu_unshifted_failure_stage' },
        @{ top = 'preconditioner_quality_status'; attempt = 'preconditioner_quality_status' }
    )
    foreach ($item in $textMap) {
        if (-not (Compare-TextField $fixtureKey 'A_vs_B_attempt1' $item.top `
                (Get-TextValue $profileA $item.top) (Get-TextValue $profileB "$prefix$($item.attempt)"))) {
            $ok = $false
        }
    }
    $numberMap = @(
        @{ top = 'true_residual_rel'; attempt = 'true_residual_rel' },
        @{ top = 'bie_cartesian_full_fmm_order'; attempt = 'cartesian_full_order' },
        @{ top = 'bie_cartesian_full_fmm_leaf_size'; attempt = 'cartesian_full_leaf_size' },
        @{ top = 'bie_cartesian_full_fmm_opening_theta'; attempt = 'theta' },
        @{ top = 'panels'; attempt = 'panel_count' },
        @{ top = 'bie_cartesian_full_fmm_direct_panel_pairs_expected'; attempt = 'direct_pairs_expected' },
        @{ top = 'bie_cartesian_full_fmm_direct_panel_pairs_applied_last'; attempt = 'direct_pairs_applied_last' },
        @{ top = 'bie_cartesian_full_fmm_direct_panel_pairs_applied_total'; attempt = 'direct_pairs_applied_total' },
        @{ top = 'fmm_cartesian_solution_audit_rel'; attempt = 'solution_audit_rel' },
        @{ top = 'fmm_cartesian_solution_audit_abs_l2'; attempt = 'solution_audit_abs_l2' },
        @{ top = 'fmm_cartesian_solution_audit_direct_l2'; attempt = 'solution_audit_direct_l2' },
        @{ top = 'fmm_cartesian_solution_audit_targets'; attempt = 'solution_audit_targets' },
        @{ top = 'fmm_cartesian_solution_audit_rows'; attempt = 'solution_audit_rows' },
        @{ top = 'fmm_cartesian_solution_audit_direct_pairs_expected'; attempt = 'solution_audit_direct_pairs_expected' },
        @{ top = 'fmm_cartesian_solution_audit_direct_pairs_applied'; attempt = 'solution_audit_direct_pairs_applied' },
        @{ top = 'preconditioner_ilu_regularization_selected'; attempt = 'preconditioner_ilu_regularization_selected' },
        @{ top = 'preconditioner_ilu_regularization_shifted_rows'; attempt = 'preconditioner_ilu_regularization_shifted_rows' },
        @{ top = 'preconditioner_valid_panels'; attempt = 'preconditioner_valid_panels' },
        @{ top = 'preconditioner_ill_conditioned_blocks'; attempt = 'preconditioner_ill_conditioned_blocks' },
        @{ top = 'preconditioner_quality_rejected_blocks'; attempt = 'preconditioner_quality_rejected_blocks' }
    )
    foreach ($item in $numberMap) {
        if (-not (Compare-DoubleField $fixtureKey 'A_vs_B_attempt1' $item.top `
                (Get-DoubleValue $profileA $item.top) (Get-DoubleValue $profileB "$prefix$($item.attempt)"))) {
            $ok = $false
        }
    }
    Add-Comparison $fixtureKey 'A_vs_B_attempt1' 'attempt1_energy_kj_mol' 'not_exposed' $true 'ordinary A available' 'retry attempt snapshot does not expose energy' | Out-Null
    Add-Comparison $fixtureKey 'A_vs_B_attempt1' 'attempt1_geometry_hash' 'not_exposed' $true 'ordinary A has no top-level geometry hash' (Get-TextValue $profileB "${prefix}geometry_hash") | Out-Null
    return $ok
}

function Compare-FinalToIndependent([object]$Fixture, [object]$B, [object]$C,
                                    [string]$Comparison = 'B_final_vs_C') {
    $fixtureKey = $Fixture.key
    $left = $B.profile
    $right = $C.profile
    $ok = $true
    foreach ($field in @(
        'status', 'bie_cartesian_full_fmm_tree_hash', 'bie_cartesian_full_fmm_worklist_hash',
        'bie_cartesian_full_fmm_policy_hash', 'fmm_cartesian_pre_audit_target_hash',
        'fmm_cartesian_solution_audit_target_hash', 'fmm_cartesian_solution_audit_status',
        'preconditioner_effective', 'preconditioner_quality_status')) {
        if (-not (Compare-TextField $fixtureKey $Comparison $field `
                (Get-TextValue $left $field) (Get-TextValue $right $field))) { $ok = $false }
    }
    foreach ($field in @(
        'iterations', 'matvec_calls', 'fmm_local_gmres_matvec_calls',
        'preconditioner_apply_calls', 'bie_cartesian_full_fmm_direct_panel_pairs_expected',
        'bie_cartesian_full_fmm_direct_panel_pairs_applied_last',
        'bie_cartesian_full_fmm_direct_panel_pairs_applied_total',
        'fmm_cartesian_solution_audit_targets', 'fmm_cartesian_solution_audit_rows',
        'fmm_cartesian_solution_audit_direct_pairs_expected',
        'fmm_cartesian_solution_audit_direct_pairs_applied',
        'fmm_cartesian_channel_direct_pairs_expected', 'fmm_cartesian_channel_direct_pairs_applied')) {
        if (-not (Compare-DoubleField $fixtureKey $Comparison $field `
                (Get-DoubleValue $left $field) (Get-DoubleValue $right $field))) { $ok = $false }
    }
    foreach ($field in @(
        'true_residual_rel', 'fmm_cartesian_solution_audit_rel',
        'fmm_cartesian_solution_audit_abs_l2', 'fmm_cartesian_solution_audit_direct_l2',
        'fmm_cartesian_channel_k1_rel', 'fmm_cartesian_channel_k2_rel',
        'fmm_cartesian_channel_k3_rel', 'fmm_cartesian_channel_k4_rel')) {
        if (-not (Compare-DoubleField $fixtureKey $Comparison $field `
                (Get-DoubleValue $left $field) (Get-DoubleValue $right $field))) { $ok = $false }
    }
    if (-not (Compare-DoubleField $fixtureKey $Comparison 'energy_kj_mol' `
            (Get-DoubleValue $B.energy 'E_rxn_kJmol') (Get-DoubleValue $C.energy 'E_rxn_kJmol'))) { $ok = $false }
    return $ok
}

function Get-O4ControlRows {
    $controlRoot = Join-Path $RepoRoot 'test\tmp\p2a_full_fmm_error_control_20260801\runs_explicit'
    foreach ($fixture in $script:SelectedFixtures) {
        $runDir = Join-Path $controlRoot "$($fixture.key)_o4_l64_t0p20"
        $argPath = Join-Path $runDir 'arg.in'
        $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
        $argText = if (Test-Path -LiteralPath $argPath) { Get-Content -LiteralPath $argPath -Raw } else { '' }
        $pqr = Convert-ToTitanPath (Join-Path $RepoRoot $fixture.pqr_relative)
        $protocolMatch = $argText -match [Regex]::Escape("-pbe-test-pqr=$pqr") -and
            $argText -match '(?m)^-pbe-bie-primary-surface=analytic-ses$' -and
            $argText -match '(?m)^-pbe-bie-ses-grid-spacing=1.25$' -and
            $argText -match '(?m)^-pbe-bie-discretization=node-patch$' -and
            $argText -match '(?m)^-pbe-bie-fmm-cartesian-order=4$' -and
            $argText -match '(?m)^-pbe-bie-tree-leaf-size=64$' -and
            $argText -match '(?m)^-pbe-bie-opening-theta=0.2$' -and
            $argText -match '(?m)^-pbe-salt=0.15$' -and $argText -match '(?m)^-nt=4$'
        $profile = Get-CsvSingleRow $profilePath
        $binaryProvenance = $false
        $hashTelemetry = (Get-TextValue $profile 'bie_cartesian_full_fmm_tree_hash') -ne ''
        [pscustomobject][ordered]@{
            fixture_key = $fixture.key; fixture = $fixture.name
            historical_o4_run_dir = $runDir
            arg_exists = Test-Path -LiteralPath $argPath
            profile_exists = Test-Path -LiteralPath $profilePath
            protocol_fields_match = $protocolMatch
            binary_provenance_available = $binaryProvenance
            current_hash_telemetry_available = $hashTelemetry
            assessment = if ($protocolMatch -and $binaryProvenance -and $hashTelemetry) {
                'matched_external_control_available_not_used_for_retry'
            } else {
                'control_gap_not_used_for_retry'
            }
            reason = 'historical O4 raw args exist, but no pinned binary provenance and no current identity telemetry; excluded from P2B2D retry gates'
        }
    }
}

$TitanExe = [IO.Path]::GetFullPath($TitanExe)
$script:ApprovedTitan = Assert-ApprovedTitanExe $TitanExe
$WorkRoot = Assert-SafeWorkRoot $WorkRoot
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite existing P2B2D output root: $WorkRoot"
}

$script:Fixtures = @(
    [pscustomobject]@{
        key = '1sh1'; name = '1SH1'
        pqr_relative = 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
        pqr_sha256 = 'c661ed6012607d92bbab2d43b970ae1ec954a16d0694ce60f7037a5cd5c9bd6e'
        source_relative = 'test\pbe_inputs_public\1sh1_source.pdb'
        source_sha256 = '2f6b7bc602aca77f8c3d918250bbdec4e321498dd5fd4e02123e286869915530'
    },
    [pscustomobject]@{
        key = 'barnase'; name = 'barnase'
        pqr_relative = 'test\pbe_inputs_public\barnase_full_amber_ph7.pqr'
        pqr_sha256 = '5a8136a5cc07ab472fcfa0e817d9ecb16d2d705b2bb88ac026ea7fbe16885002'
        source_relative = 'test\pbe_inputs_public\barnase_source.pdb'
        source_sha256 = '5c059cc01374299f45898fe23e2a8be51bf48e35ef13c682a850c425b6d95cd9'
    },
    [pscustomobject]@{
        key = '1a63'; name = '1A63'
        pqr_relative = 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
        pqr_sha256 = 'a26d0c3e1fa32fde448a3571a3ab720e4fc8a2bb40d5989f77731ee678681c7a'
        source_relative = 'test\pbe_inputs_public\1a63_source.pdb'
        source_sha256 = '3040aaff38e1ce4dcd3d1123b861bec0548a7b42065c1a0fcaa878b091be257e'
    }
)

$validFixtureKeys = @('1sh1', 'barnase', '1a63')
$normalizedFixtureKeys = [Collections.Generic.List[string]]::new()
foreach ($fixtureKey in @($FixtureKeys)) {
    $normalized = ([string]$fixtureKey).Trim().ToLowerInvariant()
    if ([string]::IsNullOrWhiteSpace($normalized)) {
        throw 'P2B2D fixture selection contains an empty key.'
    }
    $normalizedFixtureKeys.Add($normalized) | Out-Null
}
if ($normalizedFixtureKeys.Count -eq 0) {
    throw 'P2B2D fixture selection must not be empty.'
}
if (@($normalizedFixtureKeys | Sort-Object -Unique).Count -ne $normalizedFixtureKeys.Count) {
    throw 'P2B2D fixture selection must not contain duplicate keys.'
}
$FixtureKeys = @($normalizedFixtureKeys)
foreach ($fixtureKey in $FixtureKeys) {
    if ($fixtureKey -notin $validFixtureKeys) {
        throw "Unsupported P2B2D fixture key: $fixtureKey"
    }
}
$script:SelectedFixtures = @($script:Fixtures | Where-Object {
    $FixtureKeys -contains $_.key
})
if ($script:SelectedFixtures.Count -ne $FixtureKeys.Count) {
    throw 'P2B2D fixture selection did not resolve every requested fixture.'
}

$fixtureRows = [Collections.Generic.List[object]]::new()
foreach ($fixture in $script:SelectedFixtures) {
    foreach ($kind in @('pqr', 'source')) {
        $relative = $fixture."${kind}_relative"
        $expected = $fixture."${kind}_sha256"
        $path = Join-Path $RepoRoot $relative
        if (-not (Test-Path -LiteralPath $path)) { throw "Missing $kind fixture: $path" }
        $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant()
        if ($actual -ne $expected) { throw "Fixture hash mismatch for $($fixture.name) ${kind}: $actual" }
        $fixtureRows.Add([pscustomobject][ordered]@{
            fixture_key = $fixture.key; fixture = $fixture.name; kind = $kind
            relative_path = $relative; sha256 = $actual; bytes = (Get-Item -LiteralPath $path).Length
        }) | Out-Null
    }
}
$script:LibRoot = Get-LibRoot

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$RunsPath = Join-Path $WorkRoot 'p2b2d_runs.csv'
$GatesPath = Join-Path $WorkRoot 'p2b2d_gate_results.csv'
$ComparisonPath = Join-Path $WorkRoot 'p2b2d_comparisons.csv'
$PerformancePath = Join-Path $WorkRoot 'p2b2d_performance.csv'
$SummaryPath = Join-Path $WorkRoot 'p2b2d_summary.csv'
$MetadataPath = Join-Path $WorkRoot 'p2b2d_metadata.csv'
$FixturePath = Join-Path $WorkRoot 'p2b2d_fixture_hashes.csv'
$O4Path = Join-Path $WorkRoot 'p2b2d_o4_control_assessment.csv'

$fixtureRows | Export-Csv -LiteralPath $FixturePath -NoTypeInformation -Encoding ASCII
@([pscustomobject][ordered]@{
    scope = 'P2B2D three-fixture validation only; O1 retry safeguard; no P3; no O4 runtime fallback'
    output_root = $WorkRoot
    titan_exe = $TitanExe
    titan_exe_bytes = $script:ApprovedTitan.bytes
    titan_exe_last_write_utc = $script:ApprovedTitan.last_write_utc.ToString('o')
    titan_exe_sha256 = $script:ApprovedTitan.sha256
    git_head = (& git -C $RepoRoot rev-parse HEAD).Trim()
    engine = 'cartesian-full-experimental'
    order = 1
    leaf = 64
    theta_ladder = '0.20 -> 0.15 exactly once'
    salt_molar = '0.15'
    nt = 4
    preconditioner = 'explicit fmm-nearfield-ilu0; unshifted selected state required'
    retry_mode = 'channel-shadow-retry-global-theta015'
    source_provenance = 'P2B1/P2A AMBER analytic-SES node-patch protocol; P2B2C required JCP route flag and retry progress telemetry'
}) | Export-Csv -LiteralPath $MetadataPath -NoTypeInformation -Encoding ASCII

$script:RunRows = [Collections.Generic.List[object]]::new()
$script:GateRows = [Collections.Generic.List[object]]::new()
$script:ComparisonRows = [Collections.Generic.List[object]]::new()

function Save-PartialResults {
    @($script:RunRows) | Export-Csv -LiteralPath $RunsPath -NoTypeInformation -Encoding ASCII
}

function Run-Scenario([object]$Fixture, [string]$FixtureRoot, [string]$Role,
                       [string]$Mode, [double]$Theta) {
    $runDir = Join-Path $FixtureRoot $Role
    New-RunDirectory $runDir
    $arg = New-ArgFile $runDir $Fixture $Mode $Theta
    Write-Host "run: $($Fixture.name) $Role"
    $case = Invoke-SolverCase "$($Fixture.key)_$Role" $runDir $arg
    $profileArtifact = Read-RequiredSingleCsvArtifact (Join-Path $runDir 'pbe_tabi_profile.csv') 'profile'
    $energyArtifact = Read-RequiredSingleCsvArtifact (Join-Path $runDir 'pbe_energy.csv') 'energy'
    $memoryArtifact = Read-RequiredSingleCsvArtifact (Join-Path $runDir 'pbe_memory_profile.csv') 'memory'
    $profile = $profileArtifact.row
    $energy = $energyArtifact.row
    $memory = $memoryArtifact.row
    $execution = Test-RunExecutionArtifacts $case $profileArtifact $energyArtifact $memoryArtifact
    $record = New-RunRecord $Fixture $Role $Mode $Theta $case $arg $profile $energy $memory `
        $profileArtifact $energyArtifact $memoryArtifact $execution
    $script:RunRows.Add($record) | Out-Null
    Save-PartialResults
    return [pscustomobject]@{
        case = $case; arg = $arg; profile = $profile; energy = $energy; memory = $memory
        profile_artifact = $profileArtifact; energy_artifact = $energyArtifact; memory_artifact = $memoryArtifact
        execution = $execution; record = $record
    }
}

function Evaluate-Fixture([object]$Fixture, [object]$A, [object]$B, [object]$C) {
    $key = $Fixture.key
    $aTechnical = Test-TopTechnical $A.profile 0.20 'channel-shadow'
    $cTechnical = Test-TopTechnical $C.profile 0.15 'channel-shadow'
    Add-Gate $key 'ordinary_A_execution_artifacts_valid' $A.execution.ok $A.execution.reason
    Add-Gate $key 'retry_execution_artifacts_valid' $B.execution.ok $B.execution.reason
    Add-Gate $key 'ordinary_C_execution_artifacts_valid' $C.execution.ok $C.execution.reason
    Add-Gate $key 'ordinary_A_technical_valid' ($A.execution.ok -and $aTechnical.ok) $aTechnical.reason
    Add-Gate $key 'ordinary_C_technical_valid' ($C.execution.ok -and $cTechnical.ok) $cTechnical.reason

    $bProfile = $B.profile
    $attempt1 = Test-RetryAttempt $bProfile 1 0.20
    $attempt2 = Test-RetryAttempt $bProfile 2 0.15
    $retryTopOk = $B.execution.ok -and $null -ne $bProfile -and (Get-TextValue $bProfile 'fmm_cartesian_retry_mode') -eq 'channel-shadow-retry-global-theta015' -and
        (Get-Int64Value $bProfile 'fmm_cartesian_retry_attempts_started') -eq 2 -and
        (Get-Int64Value $bProfile 'fmm_cartesian_retry_attempts_completed') -eq 2 -and
        (Get-Int64Value $bProfile 'fmm_cartesian_retry_accepted_attempt') -eq 2 -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_classification') -eq 'retry_trigger' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_classification') -eq 'accept'
    Add-Gate $key 'retry_state_machine_exactly_two_attempts' $retryTopOk 'mode/attempt count/accepted attempt/classifications'
    $a1Natural = $B.execution.ok -and $attempt1.ok -and $attempt1.audit -gt 1.0e-4
    $a2Accepted = $B.execution.ok -and $attempt2.ok -and $attempt2.audit -le 1.0e-4
    Add-Gate $key 'retry_attempt1_natural_trigger' $a1Natural "audit=$(Format-Double $attempt1.audit); $($attempt1.reason)"
    Add-Gate $key 'retry_attempt2_accept' $a2Accepted "audit=$(Format-Double $attempt2.audit); $($attempt2.reason)"

    $identity = $B.execution.ok -and $null -ne $bProfile -and
        (Get-BoolValue $bProfile 'fmm_cartesian_retry_geometry_identity_valid') -and
        (Get-BoolValue $bProfile 'fmm_cartesian_retry_tree_identity_valid') -and
        (Get-BoolValue $bProfile 'fmm_cartesian_retry_policy_transition_valid') -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_geometry_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_geometry_hash') -eq (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_geometry_hash') -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_tree_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_tree_hash') -eq (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_tree_hash') -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_worklist_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_worklist_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_worklist_hash') -ne (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_worklist_hash') -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_policy_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_policy_hash') -ne '0' -and
        (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt1_policy_hash') -ne (Get-TextValue $bProfile 'fmm_cartesian_retry_attempt2_policy_hash')
    Add-Gate $key 'retry_identity_and_theta_policy_transition' $identity 'geometry/tree equal; worklist/policy nonzero and changed with theta'
    $artifacts = Test-RetryArtifacts $B.case.run_dir $bProfile
    Add-Gate $key 'retry_trace_and_progress_complete' $artifacts.ok $artifacts.reason

    $a1Parity = Compare-A1Snapshot $Fixture $A $B
    Add-Gate $key 'ordinary_A_matches_retry_attempt1_common_telemetry' $a1Parity 'strict common-field parity; energy/geometry absent from ordinary top-level telemetry are recorded separately'
    $finalParity = Compare-FinalToIndependent $Fixture $B $C
    Add-Gate $key 'accepted_retry_A2_matches_independent_C' $finalParity 'energy/audit/residual/hashes/deterministic counters at relative 1e-12 or exact'
}

$results = @{}
foreach ($fixture in $script:SelectedFixtures) {
    $fixtureRoot = Join-Path $WorkRoot $fixture.key
    New-RunDirectory $fixtureRoot
    $a = Run-Scenario $fixture $fixtureRoot 'A_channel_shadow_theta020' 'channel-shadow' 0.20
    $b = Run-Scenario $fixture $fixtureRoot 'B_retry_theta020_to015' 'channel-shadow-retry-global-theta015' 0.20
    $c = Run-Scenario $fixture $fixtureRoot 'C_channel_shadow_theta015' 'channel-shadow' 0.15
    $results[$fixture.key] = [pscustomobject]@{ A = $a; B = $b; C = $c }
    Evaluate-Fixture $fixture $a $b $c
}

if ($FixtureKeys -contains '1sh1') {
    $repeatFixture = @($script:SelectedFixtures | Where-Object { $_.key -eq '1sh1' })[0]
    $repeatRoot = Join-Path $WorkRoot '1sh1_retry_repeat'
    New-RunDirectory $repeatRoot
    $repeat = Run-Scenario $repeatFixture $repeatRoot 'B_retry_theta020_to015_repeat' 'channel-shadow-retry-global-theta015' 0.20
    $primary = $results['1sh1'].B
    $repeatProfile = $repeat.profile
    Add-Gate '1sh1' 'retry_repeat_execution_artifacts_valid' $repeat.execution.ok $repeat.execution.reason
    $repeatAttempt1 = Test-RetryAttempt $repeatProfile 1 0.20
    $repeatAttempt2 = Test-RetryAttempt $repeatProfile 2 0.15
    $repeatBasic = $repeat.execution.ok -and $repeatAttempt1.ok -and $repeatAttempt2.ok -and
        $repeatAttempt1.audit -gt 1.0e-4 -and $repeatAttempt2.audit -le 1.0e-4 -and
        $null -ne $repeatProfile -and
        (Get-TextValue $repeatProfile 'fmm_cartesian_retry_attempt1_classification') -eq 'retry_trigger' -and
        (Get-TextValue $repeatProfile 'fmm_cartesian_retry_attempt2_classification') -eq 'accept' -and
        (Get-Int64Value $repeatProfile 'fmm_cartesian_retry_accepted_attempt') -eq 2
    Add-Gate '1sh1' 'retry_repeat_state_machine' $repeatBasic 'independent root has the same bounded A1->A2 state machine'
    $repeatArtifacts = Test-RetryArtifacts $repeat.case.run_dir $repeatProfile
    Add-Gate '1sh1' 'retry_repeat_trace_and_progress_complete' $repeatArtifacts.ok $repeatArtifacts.reason
    $repeatParity = Compare-FinalToIndependent $repeatFixture $primary $repeat '1SH1_retry_repeat_final'
    foreach ($suffix in @(
        'retry_attempt1_classification', 'retry_attempt1_solution_audit_rel',
        'retry_attempt1_geometry_hash', 'retry_attempt1_tree_hash',
        'retry_attempt1_worklist_hash', 'retry_attempt1_policy_hash',
        'retry_attempt2_classification', 'retry_attempt2_solution_audit_rel',
        'retry_attempt2_geometry_hash', 'retry_attempt2_tree_hash',
        'retry_attempt2_worklist_hash', 'retry_attempt2_policy_hash')) {
        $left = Get-TextValue $primary.profile "fmm_cartesian_$suffix"
        $right = Get-TextValue $repeat.profile "fmm_cartesian_$suffix"
        if (-not (Compare-TextField '1sh1' '1SH1_retry_repeat_attempts' $suffix $left $right)) { $repeatParity = $false }
    }
    Add-Gate '1sh1' 'retry_repeat_scientific_hash_counter_parity' $repeatParity 'same final and per-attempt scientific values/hashes/counters; timing/RSS excluded'
}

$performanceRows = [Collections.Generic.List[object]]::new()
foreach ($fixture in $script:SelectedFixtures) {
    $item = $results[$fixture.key]
    $aWall = $item.A.case.wall_seconds
    $bWall = $item.B.case.wall_seconds
    $cWall = $item.C.case.wall_seconds
    $performanceRows.Add([pscustomobject][ordered]@{
        fixture_key = $fixture.key; fixture = $fixture.name
        A_wall_seconds = Format-Double $aWall
        B_retry_total_wall_seconds = Format-Double $bWall
        B_attempt1_wall_seconds = Format-Double (Get-DoubleValue $item.B.profile 'fmm_cartesian_retry_attempt1_wall_seconds')
        B_attempt2_wall_seconds = Format-Double (Get-DoubleValue $item.B.profile 'fmm_cartesian_retry_attempt2_wall_seconds')
        C_wall_seconds = Format-Double $cWall
        retry_overhead_vs_A_fraction = Format-Double (($bWall - $aWall) / [Math]::Max($aWall, 1.0e-30))
        retry_overhead_vs_C_fraction = Format-Double (($bWall - $cWall) / [Math]::Max($cWall, 1.0e-30))
        A_profile_peak_working_set_bytes = $item.A.record.profile_peak_rss_bytes
        B_profile_peak_working_set_bytes = $item.B.record.profile_peak_rss_bytes
        C_profile_peak_working_set_bytes = $item.C.record.profile_peak_rss_bytes
        A_profile_memory_bytes = $item.A.record.profile_memory_bytes
        B_profile_memory_bytes = $item.B.record.profile_memory_bytes
        C_profile_memory_bytes = $item.C.record.profile_memory_bytes
        A_profile_estimated_bytes_total = $item.A.record.profile_memory_estimated_bytes_total
        B_profile_estimated_bytes_total = $item.B.record.profile_memory_estimated_bytes_total
        C_profile_estimated_bytes_total = $item.C.record.profile_memory_estimated_bytes_total
        A_profile_rss_before_bytes = $item.A.record.profile_rss_before_bytes
        A_profile_rss_after_bytes = $item.A.record.profile_rss_after_bytes
        B_profile_rss_before_bytes = $item.B.record.profile_rss_before_bytes
        B_profile_rss_after_bytes = $item.B.record.profile_rss_after_bytes
        C_profile_rss_before_bytes = $item.C.record.profile_rss_before_bytes
        C_profile_rss_after_bytes = $item.C.record.profile_rss_after_bytes
        A_process_peak_rss_bytes_diagnostic = $item.A.record.process_peak_rss_bytes_diagnostic
        B_process_peak_rss_bytes_diagnostic = $item.B.record.process_peak_rss_bytes_diagnostic
        C_process_peak_rss_bytes_diagnostic = $item.C.record.process_peak_rss_bytes_diagnostic
        B_attempt1_accounted_memory_bytes = $item.B.record.retry_attempt1_accounted_memory_bytes
        B_attempt1_rss_before_bytes = $item.B.record.retry_attempt1_rss_before_bytes
        B_attempt1_rss_after_bytes = $item.B.record.retry_attempt1_rss_after_bytes
        B_attempt2_accounted_memory_bytes = $item.B.record.retry_attempt2_accounted_memory_bytes
        B_attempt2_rss_before_bytes = $item.B.record.retry_attempt2_rss_before_bytes
        B_attempt2_rss_after_bytes = $item.B.record.retry_attempt2_rss_after_bytes
        A_iterations = $item.A.record.iterations
        B_final_iterations = $item.B.record.iterations
        C_iterations = $item.C.record.iterations
        A_matvec_calls = $item.A.record.matvec_calls
        B_final_matvec_calls = $item.B.record.matvec_calls
        C_matvec_calls = $item.C.record.matvec_calls
        A_direct_pairs_total = $item.A.record.direct_pairs_applied_total
        B_final_direct_pairs_total = $item.B.record.direct_pairs_applied_total
        C_direct_pairs_total = $item.C.record.direct_pairs_applied_total
        note = 'single sequential observations only; no speedup claim'
    }) | Out-Null
}

$script:RunRows | Export-Csv -LiteralPath $RunsPath -NoTypeInformation -Encoding ASCII
$script:GateRows | Export-Csv -LiteralPath $GatesPath -NoTypeInformation -Encoding ASCII
$script:ComparisonRows | Export-Csv -LiteralPath $ComparisonPath -NoTypeInformation -Encoding ASCII
$performanceRows | Export-Csv -LiteralPath $PerformancePath -NoTypeInformation -Encoding ASCII
@(Get-O4ControlRows) | Export-Csv -LiteralPath $O4Path -NoTypeInformation -Encoding ASCII

$failed = @($script:GateRows | Where-Object { -not $_.passed })
@([pscustomobject][ordered]@{
    all_gates_passed = $failed.Count -eq 0
    gate_count = $script:GateRows.Count
    failed_gate_count = $failed.Count
    run_count = $script:RunRows.Count
    sequential = $true
    o4_runtime_runs = 0
    o4_role = 'historical external control only; not used in retry acceptance'
    a1_energy_snapshot = 'not exposed by production retry attempt telemetry; no telemetry extension made in validation-only slice'
    a1_geometry_ordinary_parity = 'not exposed by ordinary top-level profile; retry internal identity gate remains required'
}) | Export-Csv -LiteralPath $SummaryPath -NoTypeInformation -Encoding ASCII

if ($failed.Count -ne 0) {
    $labels = ($failed | ForEach-Object { "$($_.fixture_key):$($_.gate)" }) -join ', '
    throw "P2B2D validation failed closed: $labels. See $GatesPath"
}

Write-Host "PASS P2B2D three-fixture retry validation: runs=$($script:RunRows.Count) gates=$($script:GateRows.Count) root=$WorkRoot"
