[CmdletBinding()]
param(
    [string]$BaselineExe = '',
    [string]$InspectorExe = '',
    [string]$WorkRoot = '',
    [ValidateRange(1, 50)]
    [int]$MeasuredPairs = 7,
    [ValidateRange(60, 3600)]
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$ValidateOnly
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$AllowedWorkRootBase = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
$EvidenceRoot = Join-Path $AllowedWorkRootBase 'p2b2d_cartesian_retry_three_fixture_final_20260802'
$DefaultWorkRoot = Join-Path $AllowedWorkRootBase 'p3d0_direct_phase_attribution_20260802'
$ApprovedBaselineSha256 = '6A3D6B605E93F27572FFD4831369F028F568ED0D25F41EF740D842034BA3F90C'
[int64]$ApprovedBaselineBytes = 22775296
$ProfileIdentityFields = @(
    'status', 'panels', 'iterations', 'matvec_calls',
    'fmm_local_gmres_matvec_calls', 'jcp_direct_kernel',
    'bie_cartesian_full_fmm_active', 'bie_cartesian_full_fmm_requested',
    'bie_cartesian_full_fmm_effective', 'bie_cartesian_full_fmm_order',
    'bie_cartesian_full_fmm_leaf_size', 'bie_cartesian_full_fmm_opening_theta',
    'bie_cartesian_full_fmm_tree_hash', 'bie_cartesian_full_fmm_worklist_hash',
    'bie_cartesian_full_fmm_policy_hash', 'bie_cartesian_full_fmm_direct_pairs',
    'bie_cartesian_full_fmm_direct_panel_pairs_expected',
    'bie_cartesian_full_fmm_direct_panel_pairs_applied_last',
    'bie_cartesian_full_fmm_direct_panel_pairs_applied_total',
    'preconditioner_effective', 'preconditioner_ilu_regularization_selected',
    'preconditioner_ilu_regularization_shifted_rows',
    'preconditioner_ilu_unshifted_failure_stage', 'fallback_used',
    'fmm_m2l_finite_outputs', 'fmm_local_operator_finite',
    'fmm_cartesian_pre_audit_mode', 'fmm_cartesian_retry_mode',
    'true_residual_rel', 'fmm_cartesian_solution_audit_rel',
    'fmm_cartesian_solution_audit_status',
    'fmm_cartesian_solution_audit_target_hash',
    'fmm_cartesian_solution_audit_direct_pairs_expected',
    'fmm_cartesian_solution_audit_direct_pairs_applied'
)

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
    $base = Get-PathForComparison $AllowedWorkRootBase
    $repo = Get-PathForComparison $RepoRoot
    $testRoot = Get-PathForComparison $PSScriptRoot
    $drive = Get-PathForComparison ([IO.Path]::GetPathRoot($full))
    if ($full.Equals($base, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($repo, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($testRoot, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($drive, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing unsafe P3D0 WorkRoot: $full"
    }
    $prefix = $base + [IO.Path]::DirectorySeparatorChar
    if (-not $full.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing P3D0 WorkRoot outside strict test/tmp descendant: $full"
    }
    return $full
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

function Test-FiniteValue([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Format-Double([double]$Value) {
  if (-not (Test-FiniteValue $Value)) { return 'NA' }
  return $Value.ToString('G17', $InvariantCulture)
}

function Get-DoubleText([string]$Text,
                        [double]$Default = [double]::NaN) {
    $parsed = 0.0
    if ([double]::TryParse($Text, [Globalization.NumberStyles]::Float,
            $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-TextSha256([string]$Text) {
    $sha = [Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [Text.Encoding]::UTF8.GetBytes($Text)
        return ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-', '').ToUpperInvariant()
    } finally {
        $sha.Dispose()
    }
}

function Get-FileProvenance([string]$Path, [string]$Label) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not (Test-Path -LiteralPath $full -PathType Leaf)) {
        throw "Missing ${Label}: $full"
    }
    $item = Get-Item -LiteralPath $full
    return [pscustomobject]@{
        path = $full; bytes = [int64]$item.Length
        sha256 = (Get-FileHash -LiteralPath $full -Algorithm SHA256).Hash.ToUpperInvariant()
        timestamp_utc = $item.LastWriteTimeUtc.ToString('o')
    }
}

function Assert-ApprovedBaseline([string]$Path) {
    $provenance = Get-FileProvenance $Path 'approved baseline titan.exe'
    if ($provenance.bytes -ne $ApprovedBaselineBytes -or
        $provenance.sha256 -ne $ApprovedBaselineSha256) {
        throw "P3D0 baseline mismatch: bytes=$($provenance.bytes) sha256=$($provenance.sha256)"
    }
    return $provenance
}

function Get-ArgValue([string[]]$ArgumentLines, [string]$Key) {
    $prefix = "-$Key="
    $found = @($ArgumentLines | Where-Object {
        $_.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)
    })
    if ($found.Count -ne 1) {
        throw "Expected exactly one $prefix line in P3D0 seed arg.in; got $($found.Count)"
    }
    return $found[0].Substring($prefix.Length)
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Read-RequiredSingleCsv([string]$Path, [string]$Kind) {
    $artifact = [ordered]@{ kind=$Kind; path=$Path; ok=$false; row=$null; error='' }
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        $artifact.error = 'missing'
        return [pscustomobject]$artifact
    }
    try {
        $rows = @(Import-Csv -LiteralPath $Path)
        if ($rows.Count -ne 1) {
            $artifact.error = "expected_one_row_got_$($rows.Count)"
            return [pscustomobject]$artifact
        }
        $artifact.row = $rows[0]
        $artifact.ok = $true
    } catch {
        $artifact.error = 'import_failed'
    }
    return [pscustomobject]$artifact
}

function Assert-ApprovedSeed {
    $seedDir = Join-Path (Join-Path $EvidenceRoot '1sh1') 'A_channel_shadow_theta020'
    $argPath = Join-Path $seedDir 'arg.in'
    $driverPath = Join-Path $seedDir 'driver.pdb'
    $profilePath = Join-Path $seedDir 'pbe_tabi_profile.csv'
    foreach ($path in @($argPath, $driverPath, $profilePath)) {
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
            throw "Missing approved P2B2D 1SH1/.20 seed artifact: $path"
        }
    }
    $argumentLines = @([IO.File]::ReadAllLines($argPath))
    $required = @{
        'pbe'='bie'; 'pbe-bie-discretization'='node-patch';
        'pbe-bie-node-patch-operator-audit'='on'; 'pbe-bie-sum'='treecode';
        'pbe-bie-fast-operator'='p1';
        'pbe-bie-fmm-engine'='cartesian-full-experimental';
        'pbe-bie-fmm-cartesian-order'='1';
        'pbe-bie-fmm-cartesian-error-audit'='channel-shadow';
        'pbe-bie-tree-leaf-size'='64';
        'pbe-bie-precond'='fmm-nearfield-ilu0';
        'pbe-salt'='0.15'; 'nt'='4'
    }
    foreach ($key in $required.Keys) {
        if ((Get-ArgValue $argumentLines $key) -ne $required[$key]) {
            throw "P3D0 seed route mismatch: $key"
        }
    }
    $theta = 0.0
    if (-not [double]::TryParse((Get-ArgValue $argumentLines 'pbe-bie-opening-theta'),
            [Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$theta) -or
        [Math]::Abs($theta - 0.20) -gt 1.0e-15) {
        throw 'P3D0 seed theta must be exactly .20'
    }
    $pqrPath = [IO.Path]::GetFullPath((Get-ArgValue $argumentLines 'pbe-test-pqr'))
    if (-not (Test-Path -LiteralPath $pqrPath -PathType Leaf)) {
        throw "Missing seed PQR: $pqrPath"
    }
    return [pscustomobject]@{
        seed_dir=$seedDir; arg_path=$argPath; driver_path=$driverPath; pqr_path=$pqrPath
        argument_lines=$argumentLines
        arg_sha256=(Get-FileHash -LiteralPath $argPath -Algorithm SHA256).Hash.ToUpperInvariant()
        driver_sha256=(Get-FileHash -LiteralPath $driverPath -Algorithm SHA256).Hash.ToUpperInvariant()
        pqr_sha256=(Get-FileHash -LiteralPath $pqrPath -Algorithm SHA256).Hash.ToUpperInvariant()
    }
}

function New-RunDirectory([string]$Path) {
    if (Test-Path -LiteralPath $Path) {
        throw "Refusing to overwrite P3D0 run directory: $Path"
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function New-RunArg([object]$Seed, [string]$RunDir) {
    Copy-Item -LiteralPath $Seed.driver_path -Destination (Join-Path $RunDir 'driver.pdb') -ErrorAction Stop
    $driver = Convert-ToTitanPath (Join-Path $RunDir 'driver.pdb')
    $lines = [Collections.Generic.List[string]]::new()
    $semantic = [Collections.Generic.List[string]]::new()
    foreach ($line in $Seed.argument_lines) {
        if ($line.StartsWith('-f=', [StringComparison]::OrdinalIgnoreCase)) {
            $lines.Add("-f=$driver")
            $semantic.Add('-f=<run-driver.pdb>')
        } else {
            $lines.Add($line)
            $semantic.Add($line)
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    [IO.File]::WriteAllLines($argPath, $lines, [Text.Encoding]::ASCII)
    return [pscustomobject]@{
        path=$argPath
        sha256=(Get-FileHash -LiteralPath $argPath -Algorithm SHA256).Hash.ToUpperInvariant()
        semantic_sha256=(Get-TextSha256 ($semantic -join "`n"))
    }
}

function Invoke-Case([string]$Exe, [string]$RunDir, [string]$Role,
                     [int]$PairIndex, [string]$Phase) {
    $stdout = ''; $stderr = ''; $errorText = ''; $timedOut = $false
    $exitCode = -1; [int64]$processPeak = 0; $process = $null
    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    $sidecarPath = if ($Role -eq 'inspector') {
        Join-Path $RunDir 'p3d0_direct_phase.csv'
    } else { '' }
    try {
        $psi = [Diagnostics.ProcessStartInfo]::new()
        $psi.FileName = $Exe
        $psi.WorkingDirectory = $RunDir
        $psi.UseShellExecute = $false
        $psi.CreateNoWindow = $true
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.Arguments = '-a "arg.in"'
        # Windows PowerShell 5.1 initially enumerates this StringDictionary to
        # $null. Materialize the property without caching it, then manipulate
        # the ProcessStartInfo property directly for each child process.
        [void]$psi.EnvironmentVariables
        $psi.EnvironmentVariables.Remove('TITAN_P3C0_PACKED_PLAN_INSPECT_PATH')
        $psi.EnvironmentVariables.Remove('TITAN_P3D0_DIRECT_PHASE_PATH')
        if ($Role -eq 'inspector') {
            $psi.EnvironmentVariables['TITAN_P3D0_DIRECT_PHASE_PATH'] = $sidecarPath
        }
        $process = [Diagnostics.Process]::new()
        $process.StartInfo = $psi
        [void]$process.Start()
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
            $timedOut = $true
            $process.Kill()
        }
        $process.WaitForExit()
        $stdout = $stdoutTask.GetAwaiter().GetResult()
        $stderr = $stderrTask.GetAwaiter().GetResult()
        $exitCode = $process.ExitCode
        $processPeak = [int64]$process.PeakWorkingSet64
    } catch {
        $errorText = $_.Exception.ToString()
    } finally {
        $stopwatch.Stop()
        if ($null -ne $process) { $process.Dispose() }
    }
    [IO.File]::WriteAllText((Join-Path $RunDir 'stdout.log'), $stdout, [Text.Encoding]::UTF8)
    [IO.File]::WriteAllText((Join-Path $RunDir 'stderr.log'), $stderr, [Text.Encoding]::UTF8)
    if (-not [string]::IsNullOrWhiteSpace($errorText)) {
        [IO.File]::WriteAllText((Join-Path $RunDir 'harness_exception.txt'), $errorText, [Text.Encoding]::UTF8)
    }
    return [pscustomobject]@{
        role=$Role; phase=$Phase; pair_index=$PairIndex; run_dir=$RunDir
        exit_code=$exitCode; timed_out=$timedOut; harness_error=$errorText
        wall_seconds=$stopwatch.Elapsed.TotalSeconds
        process_peak_rss_bytes_diagnostic=$processPeak
        sidecar_path=$sidecarPath
    }
}

function Get-ProfileFingerprint([object]$Profile) {
    $missing = [Collections.Generic.List[string]]::new()
    $lines = [Collections.Generic.List[string]]::new()
    foreach ($field in $ProfileIdentityFields) {
        if ($null -eq $Profile.PSObject.Properties[$field]) { $missing.Add($field) }
        $lines.Add("$field=" + (Get-TextValue $Profile $field '<missing>'))
    }
    return [pscustomobject]@{
        ok=($missing.Count -eq 0); fingerprint=(Get-TextSha256 ($lines -join "`n"))
        missing=($missing -join ';')
    }
}

function Test-Run([object]$Case, [object]$ProfileArtifact,
                  [object]$EnergyArtifact, [object]$MemoryArtifact) {
    $profile = $ProfileArtifact.row
    $energy = $EnergyArtifact.row
    $memory = $MemoryArtifact.row
    $expected = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_expected'
    $appliedLast = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_last'
    $appliedTotal = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_total'
    $matvec = Get-Int64Value $profile 'fmm_local_gmres_matvec_calls'
    $directSeconds = Get-DoubleValue $profile 'bie_cartesian_full_fmm_direct_seconds'
    $energyValue = Get-DoubleValue $energy 'E_rxn_kJmol'
    $peak = Get-Int64Value $memory 'peak_working_set_bytes'
    $execution = $Case.exit_code -eq 0 -and -not $Case.timed_out -and
        [string]::IsNullOrWhiteSpace($Case.harness_error) -and
        $ProfileArtifact.ok -and (Get-TextValue $profile 'status') -eq 'ok' -and
        $EnergyArtifact.ok -and (Get-TextValue $energy 'status') -eq 'ok' -and
        -not [string]::IsNullOrWhiteSpace((Get-TextValue $energy 'reaction_energy_source')) -and
        (Test-FiniteValue $energyValue) -and $MemoryArtifact.ok -and
        (Get-TextValue $memory 'status') -eq 'ok' -and $peak -gt 0
    $route = (Get-TextValue $profile 'jcp_direct_kernel') -eq 'packed-avx2-vector-exp-node-patch-corrected' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_active') -eq 'true' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_requested') -eq 'cartesian-full-experimental' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_effective') -eq 'cartesian-full-experimental' -and
        (Get-Int64Value $profile 'bie_cartesian_full_fmm_order') -eq 1 -and
        (Get-Int64Value $profile 'bie_cartesian_full_fmm_leaf_size') -eq 64 -and
        [Math]::Abs((Get-DoubleValue $profile 'bie_cartesian_full_fmm_opening_theta') - 0.20) -le 1.0e-15 -and
        (Get-TextValue $profile 'fmm_cartesian_pre_audit_mode') -eq 'combined+channel-shadow' -and
        (Get-TextValue $profile 'fmm_cartesian_retry_mode') -eq 'off' -and
        (Get-TextValue $profile 'fallback_used') -eq 'false' -and
        (Get-TextValue $profile 'fmm_m2l_finite_outputs') -eq 'true'
    $coverage = $expected -gt 0 -and $appliedLast -eq $expected -and
        $matvec -gt 0 -and $appliedTotal -eq ($expected * $matvec)
    $ilu = (Get-TextValue $profile 'preconditioner_effective') -eq 'fmm_nearfield_ilu0' -and
        (Get-DoubleValue $profile 'preconditioner_ilu_regularization_selected') -eq 0.0 -and
        (Get-Int64Value $profile 'preconditioner_ilu_regularization_shifted_rows') -eq 0 -and
        (Get-TextValue $profile 'preconditioner_ilu_unshifted_failure_stage') -eq 'off'
    $numeric = (Test-FiniteValue $directSeconds) -and $directSeconds -gt 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $profile 'true_residual_rel')) -and
        (Get-DoubleValue $profile 'true_residual_rel') -ge 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $profile 'fmm_cartesian_solution_audit_rel')) -and
        (Get-DoubleValue $profile 'fmm_cartesian_solution_audit_rel') -ge 0.0
    return [pscustomobject]@{
        ok=($execution -and $route -and $coverage -and $ilu -and $numeric)
        execution_ok=$execution; route_ok=$route; coverage_ok=$coverage; ilu_ok=$ilu; numeric_ok=$numeric
        reason=("execution={0};route={1};coverage={2};ilu={3};numeric={4}" -f $execution,$route,$coverage,$ilu,$numeric)
        expected=$expected; applied_last=$appliedLast; applied_total=$appliedTotal; matvec=$matvec
        direct_seconds=$directSeconds; energy=$energyValue; peak_rss=$peak
    }
}

function Test-InspectorSidecar([object]$Artifact, [object]$Profile,
                               [object]$Correctness) {
    if (-not $Artifact.ok) {
        return [pscustomobject]@{ ok=$false; reason="sidecar_$($Artifact.error)"; row=$null }
    }
    $row = $Artifact.row
    $requiredText = @{
        'status'='ok'; 'direct_kernel'='packed-avx2-vector-exp-node-patch-corrected'
    }
    $route = $true
    foreach ($key in $requiredText.Keys) {
        $route = $route -and ((Get-TextValue $row $key) -eq $requiredText[$key])
    }
    foreach ($pair in @(
            @('tree_hash','bie_cartesian_full_fmm_tree_hash'),
            @('worklist_hash','bie_cartesian_full_fmm_worklist_hash'),
            @('policy_hash','bie_cartesian_full_fmm_policy_hash'),
            @('plan_direct_leaf_pairs','bie_cartesian_full_fmm_direct_pairs'),
            @('plan_direct_panel_pairs','bie_cartesian_full_fmm_direct_panel_pairs_expected'))) {
        $route = $route -and ((Get-TextValue $row $pair[0]) -eq (Get-TextValue $Profile $pair[1]))
    }
    $applyCount = Get-Int64Value $row 'apply_count'
    $expected = Get-Int64Value $row 'direct_pairs_expected_total'
    $applied = Get-Int64Value $row 'direct_pairs_applied_total'
    $perApply = Get-Int64Value $row 'plan_direct_panel_pairs'
    $countOk = $applyCount -gt 0 -and $applyCount -eq $Correctness.matvec -and
        $expected -eq $applied -and $expected -eq ($perApply * $applyCount) -and
        $expected -eq $Correctness.applied_total
    $nonnegative = @(
        'packed_unknown_cache_reorder_seconds', 'packed_worker_scratch_prep_seconds',
        'packed_parallel_region_seconds', 'packed_worker_metric_reduction_seconds',
        'packed_finite_scan_seconds', 'packed_wrapper_residual_seconds',
        'packed_phase_sum_seconds', 'packed_total_seconds',
        'node_patch_correction_seconds', 'adapter_direct_total_seconds',
        'adapter_direct_residual_seconds', 'adapter_coverage_bookkeeping_seconds',
        'adapter_residual_coverage_seconds', 'adapter_direct_phase_sum_seconds'
    )
    $timingOk = $true
    foreach ($name in $nonnegative) {
        $value = Get-DoubleValue $row $name
        $timingOk = $timingOk -and (Test-FiniteValue $value) -and $value -ge 0.0
    }
    $packedTotal = Get-DoubleValue $row 'packed_total_seconds'
    $packedSum = Get-DoubleValue $row 'packed_phase_sum_seconds'
    $packedClosure = Get-DoubleValue $row 'packed_closure_seconds'
    $directTotal = Get-DoubleValue $row 'adapter_direct_total_seconds'
    $directSum = Get-DoubleValue $row 'adapter_direct_phase_sum_seconds'
    $directClosure = Get-DoubleValue $row 'adapter_direct_closure_seconds'
    $packedTolerance = [Math]::Max(1.0e-6 * $applyCount, 1.0e-5 * $packedTotal)
    $directTolerance = [Math]::Max(1.0e-6 * $applyCount, 1.0e-5 * $directTotal)
    # This compares the high-precision test-only sidecar against the
    # production-profile CSV, whose textual direct_seconds value is rounded.
    # It does not relax either internal sidecar closure invariant.
    $profileSerializationTolerance = [Math]::Max(5.0e-6,
        5.0e-6 * [Math]::Abs($Correctness.direct_seconds))
    $closureOk = (Test-FiniteValue $packedClosure) -and (Test-FiniteValue $directClosure) -and
        $packedSum -le ($packedTotal + $packedTolerance) -and
        [Math]::Abs($packedClosure) -le $packedTolerance -and
        $directSum -le ($directTotal + $directTolerance) -and
        [Math]::Abs($directClosure) -le $directTolerance -and
        [Math]::Abs($directTotal - $Correctness.direct_seconds) -le $profileSerializationTolerance
    return [pscustomobject]@{
        ok=($route -and $countOk -and $timingOk -and $closureOk)
        reason=("route={0};count={1};timing={2};closure={3}" -f $route,$countOk,$timingOk,$closureOk)
        row=$row
    }
}

function Get-Median([double[]]$Values) {
    [double[]]$sorted = @($Values | Sort-Object)
    [int]$count = $sorted.Length
    if ($count -eq 0) { return [double]::NaN }
    [int]$middle = [Math]::Floor($count / 2.0)
    if (($count % 2) -eq 1) { return $sorted[$middle] }
    return ($sorted[$middle - 1] + $sorted[$middle]) / 2.0
}

function Get-Mad([double[]]$Values) {
    $median = Get-Median $Values
    $deviations = [Collections.Generic.List[double]]::new()
    foreach ($value in $Values) {
        $deviations.Add([Math]::Abs($value - $median))
    }
    return Get-Median $deviations.ToArray()
}

function Get-Minimum([double[]]$Values) {
    if ($Values.Count -eq 0) { return [double]::NaN }
    return [double](($Values | Measure-Object -Minimum).Minimum)
}

function Get-Maximum([double[]]$Values) {
    if ($Values.Count -eq 0) { return [double]::NaN }
    return [double](($Values | Measure-Object -Maximum).Maximum)
}

function Write-Csv([object[]]$Rows, [string]$Path) {
    if ($Rows.Count -eq 0) {
        [IO.File]::WriteAllText($Path, '', [Text.Encoding]::UTF8)
        return
    }
    $Rows | Export-Csv -LiteralPath $Path -NoTypeInformation -Encoding UTF8
}

if ([string]::IsNullOrWhiteSpace($WorkRoot)) { $WorkRoot = $DefaultWorkRoot }
$WorkRoot = Assert-SafeWorkRoot $WorkRoot
if (-not (Test-Path -LiteralPath $WorkRoot -PathType Container)) {
    throw "P3D0 WorkRoot must be pre-created by the isolated inspector build: $WorkRoot"
}
if ([string]::IsNullOrWhiteSpace($BaselineExe)) { $BaselineExe = Join-Path $RepoRoot 'titan.exe' }
if ([string]::IsNullOrWhiteSpace($InspectorExe)) { $InspectorExe = Join-Path $WorkRoot 'titan_p3d0_timing_inspector.exe' }
$baseline = Assert-ApprovedBaseline $BaselineExe
$inspector = Get-FileProvenance $InspectorExe 'P3D0 test-only inspector binary'
if ($inspector.sha256 -eq $baseline.sha256) {
    throw 'P3D0 inspector binary is unexpectedly byte-identical to the baseline'
}
if ($MeasuredPairs -lt 7) { throw 'P3D0 pilot requires at least seven measured pairs' }
if (-not (Test-Path -LiteralPath $EvidenceRoot -PathType Container)) {
    throw "Missing approved P2B2D evidence root: $EvidenceRoot"
}
$seed = Assert-ApprovedSeed
if ($ValidateOnly) {
    Write-Output "P3D0_VALIDATE_ONLY_PASS baseline=$($baseline.sha256) inspector=$($inspector.sha256) seed=$($seed.arg_sha256)"
    exit 0
}

$artifactNames = @(
    'p3d0_pilot_provenance.csv', 'p3d0_paired_raw.csv', 'p3d0_pilot_gates.csv',
    'p3d0_phase_summary.csv', 'p3d0_preregistration_next_step.csv'
)
foreach ($name in $artifactNames) {
    if (Test-Path -LiteralPath (Join-Path $WorkRoot $name)) {
        throw "Refusing to overwrite P3D0 artifact: $name"
    }
}
$runsRoot = Join-Path $WorkRoot 'runs'
if (Test-Path -LiteralPath $runsRoot) { throw "Refusing to overwrite P3D0 runs: $runsRoot" }
New-Item -ItemType Directory -Path $runsRoot -Force | Out-Null

@([pscustomobject][ordered]@{
    schema='p3d0-direct-subphase-pilot-v1'; created_utc=[DateTime]::UtcNow.ToString('o')
    baseline_path=$baseline.path; baseline_bytes=$baseline.bytes; baseline_sha256=$baseline.sha256; baseline_timestamp_utc=$baseline.timestamp_utc
    inspector_path=$inspector.path; inspector_bytes=$inspector.bytes; inspector_sha256=$inspector.sha256; inspector_timestamp_utc=$inspector.timestamp_utc
    evidence_root=$EvidenceRoot; seed_dir=$seed.seed_dir; seed_arg_sha256=$seed.arg_sha256; seed_driver_sha256=$seed.driver_sha256; pqr_path=$seed.pqr_path; pqr_sha256=$seed.pqr_sha256
    fixture='1sh1'; theta='0.20'; order=1; leaf_size=64; nt=4; error_audit='channel-shadow'; retry_mode='off'; warmup_pairs=1; measured_pairs=$MeasuredPairs
    timing_contract='test-only broad timers; no clock or atomic in task/source/lane loops; sidecar after solve'
}) | Export-Csv -LiteralPath (Join-Path $WorkRoot 'p3d0_pilot_provenance.csv') -NoTypeInformation -Encoding UTF8

$records = [Collections.Generic.List[object]]::new()
$gates = [Collections.Generic.List[object]]::new()
$pairRows = [Collections.Generic.List[object]]::new()
for ($pair = 0; $pair -le $MeasuredPairs; ++$pair) {
    $phase = if ($pair -eq 0) { 'warmup' } else { 'measured' }
    $roles = if ($pair -eq 0 -or ($pair % 2) -eq 1) { @('baseline','inspector') } else { @('inspector','baseline') }
    $byRole = @{}
    foreach ($role in $roles) {
        $runDir = Join-Path $runsRoot ("{0}_pair_{1:D2}_{2}" -f $phase, $pair, $role)
        New-RunDirectory $runDir
        $argument = New-RunArg $seed $runDir
        $binary = if ($role -eq 'baseline') { $baseline } else { $inspector }
        $case = Invoke-Case $binary.path $runDir $role $pair $phase
        $profileArtifact = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_tabi_profile.csv') 'profile'
        $energyArtifact = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_energy.csv') 'energy'
        $memoryArtifact = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_memory_profile.csv') 'memory'
        $correctness = Test-Run $case $profileArtifact $energyArtifact $memoryArtifact
        $profileFingerprint = if ($profileArtifact.ok) { Get-ProfileFingerprint $profileArtifact.row } else {
            [pscustomobject]@{ ok=$false; fingerprint=''; missing='profile_missing' }
        }
        $energyPath = Join-Path $runDir 'pbe_energy.csv'
        $energySha = if (Test-Path -LiteralPath $energyPath -PathType Leaf) {
            (Get-FileHash -LiteralPath $energyPath -Algorithm SHA256).Hash.ToUpperInvariant()
        } else { '' }
        $sidecarArtifact = if ($role -eq 'inspector') {
            Read-RequiredSingleCsv $case.sidecar_path 'P3D0 sidecar'
        } else { [pscustomobject]@{ ok=$true; row=$null; error='' } }
        $sidecar = if ($role -eq 'inspector') {
            Test-InspectorSidecar $sidecarArtifact $profileArtifact.row $correctness
        } else { [pscustomobject]@{ ok=$true; row=$null; reason='not_applicable' } }
        $profile = $profileArtifact.row
        $sidecarRow = $sidecar.row
        $record = [pscustomobject][ordered]@{
            fixture='1sh1'; theta='0.20'; phase=$phase; pair_index=$pair; role=$role; schedule_order=($roles -join '_then_')
            binary_path=$binary.path; binary_sha256=$binary.sha256; binary_bytes=$binary.bytes
            run_dir=$runDir; arg_path=$argument.path; arg_sha256=$argument.sha256; arg_semantic_sha256=$argument.semantic_sha256
            exit_code=$case.exit_code; timed_out=$case.timed_out; harness_error=$case.harness_error; wall_seconds=Format-Double $case.wall_seconds; process_peak_rss_bytes_diagnostic=$case.process_peak_rss_bytes_diagnostic
            profile_path=$profileArtifact.path; profile_artifact_ok=$profileArtifact.ok; energy_path=$energyArtifact.path; energy_artifact_ok=$energyArtifact.ok; memory_path=$memoryArtifact.path; memory_artifact_ok=$memoryArtifact.ok
            correctness_ok=$correctness.ok; correctness_reason=$correctness.reason; profile_identity_ok=$profileFingerprint.ok; profile_identity_fingerprint=$profileFingerprint.fingerprint; profile_identity_missing=$profileFingerprint.missing
            energy_sha256=$energySha; energy_kj_mol=Format-Double $correctness.energy; direct_seconds=Format-Double $correctness.direct_seconds; direct_panel_pairs_expected=$correctness.expected; direct_panel_pairs_applied_total=$correctness.applied_total; direct_panel_pairs_per_second=Format-Double $(if ($correctness.direct_seconds -gt 0.0) { $correctness.applied_total / $correctness.direct_seconds } else { [double]::NaN })
            panels=Get-TextValue $profile 'panels'; iterations=Get-TextValue $profile 'iterations'; matvec_calls=Get-TextValue $profile 'matvec_calls'; fmm_gmres_matvec_calls=Get-TextValue $profile 'fmm_local_gmres_matvec_calls'; true_residual_rel=Get-TextValue $profile 'true_residual_rel'; solution_audit_rel=Get-TextValue $profile 'fmm_cartesian_solution_audit_rel'
            tree_hash=Get-TextValue $profile 'bie_cartesian_full_fmm_tree_hash'; worklist_hash=Get-TextValue $profile 'bie_cartesian_full_fmm_worklist_hash'; policy_hash=Get-TextValue $profile 'bie_cartesian_full_fmm_policy_hash'; direct_leaf_pairs=Get-TextValue $profile 'bie_cartesian_full_fmm_direct_pairs'; peak_rss_bytes=$correctness.peak_rss
            sidecar_path=$case.sidecar_path; sidecar_ok=$sidecar.ok; sidecar_reason=$sidecar.reason
            sidecar_apply_count=Get-TextValue $sidecarRow 'apply_count'; sidecar_packed_total_seconds=Get-TextValue $sidecarRow 'packed_total_seconds'; sidecar_parallel_seconds=Get-TextValue $sidecarRow 'packed_parallel_region_seconds'; sidecar_node_patch_seconds=Get-TextValue $sidecarRow 'node_patch_correction_seconds'; sidecar_adapter_direct_seconds=Get-TextValue $sidecarRow 'adapter_direct_total_seconds'; sidecar_packed_closure_seconds=Get-TextValue $sidecarRow 'packed_closure_seconds'; sidecar_adapter_closure_seconds=Get-TextValue $sidecarRow 'adapter_direct_closure_seconds'
        }
        $records.Add($record); $byRole[$role] = $record
        foreach ($gate in @(
                @('execution_and_artifacts',$correctness.execution_ok,$correctness.reason),
                @('route',$correctness.route_ok,$correctness.reason),
                @('coverage',$correctness.coverage_ok,$correctness.reason),
                @('plain_unshifted_ilu0',$correctness.ilu_ok,$correctness.reason),
                @('finite_numeric',$correctness.numeric_ok,$correctness.reason),
                @('sidecar',$sidecar.ok,$sidecar.reason))) {
            $gates.Add([pscustomobject][ordered]@{ fixture='1sh1'; theta='0.20'; phase=$phase; pair_index=$pair; role=$role; gate=$gate[0]; pass=[bool]$gate[1]; detail=$gate[2] })
        }
        Write-Output ("P3D0_RUN phase={0} pair={1} role={2} correct={3} sidecar={4} direct_s={5}" -f $phase,$pair,$role,$correctness.ok,$sidecar.ok,$record.direct_seconds)
    }
    $baselineRecord = $byRole['baseline']
    $inspectorRecord = $byRole['inspector']
    $energyParity = $baselineRecord.energy_sha256 -ne '' -and $baselineRecord.energy_sha256 -eq $inspectorRecord.energy_sha256
    $profileParity = $baselineRecord.profile_identity_ok -and $inspectorRecord.profile_identity_ok -and $baselineRecord.profile_identity_fingerprint -eq $inspectorRecord.profile_identity_fingerprint
    $semanticParity = $baselineRecord.arg_semantic_sha256 -eq $inspectorRecord.arg_semantic_sha256
    $pairExact = $baselineRecord.correctness_ok -and $inspectorRecord.correctness_ok -and $inspectorRecord.sidecar_ok -and $energyParity -and $profileParity -and $semanticParity
    $baselineDirect = Get-DoubleText $baselineRecord.direct_seconds
    $inspectorDirect = Get-DoubleText $inspectorRecord.direct_seconds
    $pairDeltaPercent = if ((Test-FiniteValue $baselineDirect) -and
        (Test-FiniteValue $inspectorDirect) -and $baselineDirect -gt 0.0) {
        100.0 * (($inspectorDirect - $baselineDirect) / $baselineDirect)
    } else { [double]::NaN }
    $pairRows.Add([pscustomobject][ordered]@{
        fixture='1sh1'; theta='0.20'; phase=$phase; pair_index=$pair; schedule_order=($roles -join '_then_'); exact_parity_ok=$pairExact
        baseline_direct_seconds=$baselineRecord.direct_seconds; inspector_direct_seconds=$inspectorRecord.direct_seconds
        direct_delta_percent=Format-Double $pairDeltaPercent
        energy_byte_identical=$energyParity; profile_identity_identical=$profileParity; semantic_args_identical=$semanticParity
        baseline_energy_sha256=$baselineRecord.energy_sha256; inspector_energy_sha256=$inspectorRecord.energy_sha256
        baseline_profile_identity=$baselineRecord.profile_identity_fingerprint; inspector_profile_identity=$inspectorRecord.profile_identity_fingerprint
    })
    foreach ($gate in @(
            @('pbe_energy_byte_identical',$energyParity,"baseline=$($baselineRecord.energy_sha256);inspector=$($inspectorRecord.energy_sha256)"),
            @('profile_identity_exact',$profileParity,"baseline=$($baselineRecord.profile_identity_fingerprint);inspector=$($inspectorRecord.profile_identity_fingerprint)"),
            @('pair_exact',$pairExact,"energy=$energyParity;profile=$profileParity;semantic=$semanticParity"))) {
        $gates.Add([pscustomobject][ordered]@{ fixture='1sh1'; theta='0.20'; phase=$phase; pair_index=$pair; role='pair'; gate=$gate[0]; pass=[bool]$gate[1]; detail=$gate[2] })
    }
}

$measuredPairRecords = [Collections.Generic.List[object]]::new()
$deltaValues = [Collections.Generic.List[double]]::new()
[int]$notExactPairCount = 0
foreach ($pairRecord in $pairRows) {
    if (-not $pairRecord.exact_parity_ok) { ++$notExactPairCount }
    if ($pairRecord.phase -eq 'measured') {
        $measuredPairRecords.Add($pairRecord)
        $delta = Get-DoubleText $pairRecord.direct_delta_percent
        if (Test-FiniteValue $delta) { $deltaValues.Add($delta) }
    }
}
[double[]]$deltas = $deltaValues.ToArray()
$deltaMedian = Get-Median $deltas
$deltaMad = Get-Mad $deltas
$deltaMin = Get-Minimum $deltas
$deltaMax = Get-Maximum $deltas
$allExact = $measuredPairRecords.Count -eq $MeasuredPairs -and $deltas.Count -eq $MeasuredPairs -and $notExactPairCount -eq 0
$overheadMedianOk = (Test-FiniteValue $deltaMedian) -and [Math]::Abs($deltaMedian) -le 3.0
$regressionOk = (Test-FiniteValue $deltaMax) -and $deltaMax -le 5.0
$gates.Add([pscustomobject][ordered]@{ fixture='1sh1'; theta='0.20'; phase='measured'; pair_index=-1; role='aggregate'; gate='instrumentation_median_abs_delta_le_3pct'; pass=$overheadMedianOk; detail=(Format-Double $deltaMedian) })
$gates.Add([pscustomobject][ordered]@{ fixture='1sh1'; theta='0.20'; phase='measured'; pair_index=-1; role='aggregate'; gate='no_unexplained_direct_regression_gt_5pct'; pass=$regressionOk; detail=(Format-Double $deltaMax) })

$phaseRows = [Collections.Generic.List[object]]::new()
$inspectorMeasured = [Collections.Generic.List[object]]::new()
[int]$inspectorSidecarFailureCount = 0
foreach ($record in $records) {
    if ($record.phase -eq 'measured' -and $record.role -eq 'inspector') {
        $inspectorMeasured.Add($record)
        if (-not $record.sidecar_ok) { ++$inspectorSidecarFailureCount }
    }
}
$phaseMap = @(
    @('packed_unknown_cache_reorder','packed_unknown_cache_reorder_seconds','packed_unknown_cache_reorder_fraction'),
    @('packed_worker_scratch_prep','packed_worker_scratch_prep_seconds','packed_worker_scratch_prep_fraction'),
    @('packed_parallel_region','packed_parallel_region_seconds','packed_parallel_region_fraction'),
    @('packed_worker_metric_reduction','packed_worker_metric_reduction_seconds','packed_worker_metric_reduction_fraction'),
    @('packed_finite_scan','packed_finite_scan_seconds','packed_finite_scan_fraction'),
    @('packed_wrapper_residual','packed_wrapper_residual_seconds','packed_wrapper_residual_fraction'),
    @('node_patch_correction','node_patch_correction_seconds','node_patch_correction_fraction'),
    @('adapter_direct_residual','adapter_direct_residual_seconds','adapter_direct_residual_fraction'),
    @('adapter_coverage_bookkeeping','adapter_coverage_bookkeeping_seconds','adapter_coverage_bookkeeping_fraction')
)
if ($inspectorMeasured.Count -eq $MeasuredPairs -and $inspectorSidecarFailureCount -eq 0) {
    foreach ($phase in $phaseMap) {
        $secondsValues = [Collections.Generic.List[double]]::new()
        $fractionValues = [Collections.Generic.List[double]]::new()
        foreach ($inspectorRecord in $inspectorMeasured) {
            $sidecarRow = (Read-RequiredSingleCsv $inspectorRecord.sidecar_path 'P3D0 sidecar summary').row
            $secondsValues.Add((Get-DoubleValue $sidecarRow $phase[1]))
            $fractionValues.Add((Get-DoubleValue $sidecarRow $phase[2]))
        }
        [double[]]$seconds = $secondsValues.ToArray()
        [double[]]$fractions = $fractionValues.ToArray()
        $phaseRows.Add([pscustomobject][ordered]@{
            fixture='1sh1'; theta='0.20'; phase=$phase[0]; measured_repeats=$seconds.Count
            seconds_median=Format-Double (Get-Median $seconds); seconds_mad=Format-Double (Get-Mad $seconds); seconds_min=Format-Double (Get-Minimum $seconds); seconds_max=Format-Double (Get-Maximum $seconds)
            fraction_of_adapter_direct_median=Format-Double (Get-Median $fractions); fraction_mad=Format-Double (Get-Mad $fractions); fraction_min=Format-Double (Get-Minimum $fractions); fraction_max=Format-Double (Get-Maximum $fractions)
        })
    }
} else {
    $phaseRows.Add([pscustomobject][ordered]@{
        fixture='1sh1'; theta='0.20'; phase='sidecar_unavailable'; measured_repeats=$inspectorMeasured.Count
        seconds_median='NA'; seconds_mad='NA'; seconds_min='NA'; seconds_max='NA'
        fraction_of_adapter_direct_median='NA'; fraction_mad='NA'; fraction_min='NA'; fraction_max='NA'
    })
}
$phaseRows.Add([pscustomobject][ordered]@{
    fixture='1sh1'; theta='0.20'; phase='instrumentation_overhead_paired_direct'; measured_repeats=$deltas.Count
    seconds_median='NA'; seconds_mad='NA'; seconds_min='NA'; seconds_max='NA'
    fraction_of_adapter_direct_median=Format-Double ($deltaMedian / 100.0); fraction_mad=Format-Double ($deltaMad / 100.0); fraction_min=Format-Double ($deltaMin / 100.0); fraction_max=Format-Double ($deltaMax / 100.0)
})

Write-Csv @($records) (Join-Path $WorkRoot 'p3d0_paired_raw.csv')
Write-Csv @($gates) (Join-Path $WorkRoot 'p3d0_pilot_gates.csv')
Write-Csv @($phaseRows) (Join-Path $WorkRoot 'p3d0_phase_summary.csv')
$decision = if ($allExact -and $overheadMedianOk -and $regressionOk) { 'P3D0A_PILOT_PASS' } else { 'P3D0A_PILOT_FAIL' }
$next = if ($decision -eq 'P3D0A_PILOT_PASS') {
    'Review before any P3D0B matrix; no optimization is authorized by this pilot.'
} else {
    'Fail closed: no full matrix or kernel candidate; preserve raw artifacts and review timing/parity failure.'
}
Write-Csv @([pscustomobject][ordered]@{
    decision=$decision; next_step=$next; exact_pair_parity=$allExact; median_delta_percent=Format-Double $deltaMedian; mad_delta_percent=Format-Double $deltaMad; min_delta_percent=Format-Double $deltaMin; max_delta_percent=Format-Double $deltaMax; median_overhead_gate=$overheadMedianOk; max_regression_gate=$regressionOk; completed_utc=[DateTime]::UtcNow.ToString('o')
}) (Join-Path $WorkRoot 'p3d0_preregistration_next_step.csv')
if ($decision -ne 'P3D0A_PILOT_PASS') {
    throw "P3D0 pilot failed; see $(Join-Path $WorkRoot 'p3d0_pilot_gates.csv')"
}
Write-Output "P3D0A_PILOT_PASS root=$WorkRoot median_direct_delta_pct=$(Format-Double $deltaMedian)"
