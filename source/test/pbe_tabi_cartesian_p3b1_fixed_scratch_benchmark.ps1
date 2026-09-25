[CmdletBinding()]
param(
    [string]$BaselineExe = '',
    [string]$CandidateExe = '',
    [string]$WorkRoot = '',
    [string[]]$FixtureKeys = @('1sh1', 'barnase', '1a63'),
    [double[]]$Thetas = @(0.20, 0.15),
    [ValidateRange(1, 20)]
    [int]$WarmupRepeats = 1,
    [ValidateRange(1, 50)]
    [int]$MeasuredRepeats = 7,
    [ValidateRange(60, 3600)]
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$ValidateOnly,
    [switch]$Pilot
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$AllowedWorkRootBase = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
$EvidenceRoot = Join-Path $AllowedWorkRootBase 'p2b2d_cartesian_retry_three_fixture_final_20260802'
$DefaultWorkRoot = Join-Path $AllowedWorkRootBase 'p3b1_fixed_four_lane_scratch_20260802'
$ApprovedBaselineSha256 = '6A3D6B605E93F27572FFD4831369F028F568ED0D25F41EF740D842034BA3F90C'
[int64]$ApprovedBaselineBytes = 22775296
$ScienceOutputNames = @(
    'cal.txt', 'cal.xvg', 'electrostatic.dat', 'total_energy.dat',
    'vanderwaals.dat', 'pbe_energy.csv', 'pbe_summary.dat'
)
$InvariantProfileFields = @(
    'status', 'panels', 'iterations', 'true_residual_rel',
    'jcp_direct_kernel', 'bie_cartesian_full_fmm_active',
    'bie_cartesian_full_fmm_requested', 'bie_cartesian_full_fmm_effective',
    'bie_cartesian_full_fmm_order', 'bie_cartesian_full_fmm_leaf_size',
    'bie_cartesian_full_fmm_opening_theta',
    'bie_cartesian_full_fmm_tree_hash',
    'bie_cartesian_full_fmm_worklist_hash',
    'bie_cartesian_full_fmm_policy_hash',
    'bie_cartesian_full_fmm_direct_pairs',
    'bie_cartesian_full_fmm_direct_panel_pairs_expected',
    'bie_cartesian_full_fmm_direct_panel_pairs_applied_last',
    'bie_cartesian_full_fmm_direct_panel_pairs_applied_total',
    'bie_cartesian_full_fmm_plan_bytes',
    'bie_cartesian_full_fmm_workspace_bytes',
    'bie_cartesian_full_fmm_cache_bytes', 'near_self_pairs',
    'leaf_pair_batch_kernel', 'leaf_pair_batch_vectorized_sources',
    'leaf_pair_batch_vector_lanes', 'matvec_calls',
    'fmm_local_gmres_matvec_calls', 'preconditioner_apply_calls',
    'preconditioner_effective', 'preconditioner_ilu_regularization_selected',
    'preconditioner_ilu_regularization_shifted_rows',
    'preconditioner_ilu_unshifted_failure_stage', 'fallback_used',
    'fmm_m2l_finite_outputs', 'fmm_cartesian_pre_audit_mode',
    'fmm_cartesian_pre_audit_target_hash',
    'fmm_cartesian_solution_audit_enabled',
    'fmm_cartesian_solution_audit_finite',
    'fmm_cartesian_solution_audit_relative_available',
    'fmm_cartesian_solution_audit_status',
    'fmm_cartesian_solution_audit_rel',
    'fmm_cartesian_solution_audit_target_hash',
    'fmm_cartesian_solution_audit_targets',
    'fmm_cartesian_solution_audit_rows',
    'fmm_cartesian_solution_audit_direct_pairs_expected',
    'fmm_cartesian_solution_audit_direct_pairs_applied',
    'fmm_cartesian_retry_mode'
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
    $drive = Get-PathForComparison ([IO.Path]::GetPathRoot($full))
    if ($full.Equals($base, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($repo, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($drive, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing unsafe P3B1 WorkRoot: $full"
    }
    $prefix = $base + [IO.Path]::DirectorySeparatorChar
    if (-not $full.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing P3B1 WorkRoot outside test/tmp: $full"
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

function Get-TextSha256([string]$Text) {
    $sha = [Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [Text.Encoding]::UTF8.GetBytes($Text)
        return ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-', '').ToLowerInvariant()
    } finally {
        $sha.Dispose()
    }
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Read-RequiredSingleCsv([string]$Path, [string]$Kind) {
    $result = [ordered]@{ kind = $Kind; path = $Path; ok = $false; row = $null; error = '' }
    if (-not (Test-Path -LiteralPath $Path)) {
        $result.error = 'missing'
        return [pscustomobject]$result
    }
    try {
        $rows = @(Import-Csv -LiteralPath $Path)
        if ($rows.Count -ne 1) {
            $result.error = "expected_one_row_got_$($rows.Count)"
            return [pscustomobject]$result
        }
        $result.row = $rows[0]
        $result.ok = $true
    } catch {
        $result.error = 'import_failed'
    }
    return [pscustomobject]$result
}

function Assert-ExactBaseline([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not (Test-Path -LiteralPath $full)) { throw "Missing baseline binary: $full" }
    $item = Get-Item -LiteralPath $full
    $hash = (Get-FileHash -LiteralPath $full -Algorithm SHA256).Hash.ToUpperInvariant()
    if ($item.Length -ne $ApprovedBaselineBytes -or $hash -ne $ApprovedBaselineSha256) {
        throw "Baseline binary mismatch: bytes=$($item.Length) sha256=$hash"
    }
    return [pscustomobject]@{ path = $full; bytes = [int64]$item.Length; sha256 = $hash; timestamp_utc = $item.LastWriteTimeUtc.ToString('o') }
}

function Get-BinaryProvenance([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not (Test-Path -LiteralPath $full)) { throw "Missing candidate binary: $full" }
    $item = Get-Item -LiteralPath $full
    return [pscustomobject]@{ path = $full; bytes = [int64]$item.Length; sha256 = (Get-FileHash -LiteralPath $full -Algorithm SHA256).Hash.ToUpperInvariant(); timestamp_utc = $item.LastWriteTimeUtc.ToString('o') }
}

function Get-ArgValue([string[]]$Lines, [string]$Key) {
    $prefix = "-$Key="
    $found = @($Lines | Where-Object { $_.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase) })
    if ($found.Count -ne 1) { throw "Expected exactly one $prefix line in seed arg.in; got $($found.Count)" }
    return $found[0].Substring($prefix.Length)
}

function Assert-Seed([string]$FixtureKey, [double]$Theta) {
    $thetaTag = if ($Theta -eq 0.20) { 'A_channel_shadow_theta020' } else { 'C_channel_shadow_theta015' }
    $seedDir = Join-Path (Join-Path $EvidenceRoot $FixtureKey) $thetaTag
    $seedArg = Join-Path $seedDir 'arg.in'
    $seedDriver = Join-Path $seedDir 'driver.pdb'
    if (-not (Test-Path -LiteralPath $seedArg) -or -not (Test-Path -LiteralPath $seedDriver)) {
        throw "Missing approved P2B2D seed for $FixtureKey theta=$($Theta): $seedDir"
    }
    $lines = @([IO.File]::ReadAllLines($seedArg))
    $checks = @{
        'pbe' = 'bie'; 'pbe-bie-discretization' = 'node-patch';
        'pbe-bie-node-patch-operator-audit' = 'on'; 'pbe-bie-sum' = 'treecode';
        'pbe-bie-fast-operator' = 'p1'; 'pbe-bie-fmm-engine' = 'cartesian-full-experimental';
        'pbe-bie-fmm-cartesian-order' = '1';
        'pbe-bie-fmm-cartesian-error-audit' = 'channel-shadow';
        'pbe-bie-tree-leaf-size' = '64'; 'pbe-bie-precond' = 'fmm-nearfield-ilu0';
        'pbe-salt' = '0.15'; 'nt' = '4'
    }
    foreach ($key in $checks.Keys) {
        if ((Get-ArgValue $lines $key) -ne $checks[$key]) {
            throw "P3B1 seed route mismatch for $($FixtureKey): $key"
        }
    }
    $seedTheta = 0.0
    if (-not [double]::TryParse((Get-ArgValue $lines 'pbe-bie-opening-theta'),
            [Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$seedTheta) -or
        [Math]::Abs($seedTheta - $Theta) -gt 1.0e-15) {
        throw "P3B1 seed theta mismatch for $FixtureKey"
    }
    $pqr = [IO.Path]::GetFullPath((Get-ArgValue $lines 'pbe-test-pqr'))
    if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing seed PQR: $pqr" }
    return [pscustomobject]@{
        fixture_key = $FixtureKey; theta = $Theta; seed_dir = $seedDir; seed_arg = $seedArg
        seed_arg_sha256 = (Get-FileHash -LiteralPath $seedArg -Algorithm SHA256).Hash.ToUpperInvariant()
        seed_driver = $seedDriver; seed_driver_sha256 = (Get-FileHash -LiteralPath $seedDriver -Algorithm SHA256).Hash.ToUpperInvariant()
        pqr = $pqr; pqr_sha256 = (Get-FileHash -LiteralPath $pqr -Algorithm SHA256).Hash.ToUpperInvariant()
        lines = $lines
    }
}

function New-RunDirectory([string]$Path) {
    if (Test-Path -LiteralPath $Path) { throw "Refusing to overwrite P3B1 run: $Path" }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function New-RunArg([object]$Seed, [string]$RunDir) {
    Copy-Item -LiteralPath $Seed.seed_driver -Destination (Join-Path $RunDir 'driver.pdb') -ErrorAction Stop
    $driver = Convert-ToTitanPath (Join-Path $RunDir 'driver.pdb')
    $lines = [Collections.Generic.List[string]]::new()
    $semantic = [Collections.Generic.List[string]]::new()
    foreach ($line in $Seed.lines) {
        if ($line.StartsWith('-f=', [StringComparison]::OrdinalIgnoreCase)) {
            $lines.Add("-f=$driver")
            $semantic.Add('-f=<run-driver.pdb>')
        } else {
            $lines.Add($line)
            $semantic.Add($line)
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $lines -Encoding ASCII
    return [pscustomobject]@{
        path = $argPath; sha256 = (Get-FileHash -LiteralPath $argPath -Algorithm SHA256).Hash.ToUpperInvariant()
        semantic_sha256 = Get-TextSha256 ($semantic -join "`n")
    }
}

function Invoke-Run([string]$Exe, [string]$RunDir, [string]$Name) {
    $stdout = ''; $stderr = ''; $errorText = ''; $timedOut = $false; $exitCode = -1; $processPeak = 0L
    $watch = [Diagnostics.Stopwatch]::StartNew(); $process = $null
    try {
        $psi = New-Object Diagnostics.ProcessStartInfo
        $psi.FileName = $Exe; $psi.WorkingDirectory = $RunDir
        $psi.UseShellExecute = $false; $psi.CreateNoWindow = $true
        $psi.RedirectStandardOutput = $true; $psi.RedirectStandardError = $true
        $psi.Arguments = '-a "arg.in"'
        $process = New-Object Diagnostics.Process; $process.StartInfo = $psi
        [void]$process.Start()
        $stdoutTask = $process.StandardOutput.ReadToEndAsync(); $stderrTask = $process.StandardError.ReadToEndAsync()
        if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
            $timedOut = $true; $process.Kill($true)
        }
        $process.WaitForExit(); $stdout = $stdoutTask.GetAwaiter().GetResult(); $stderr = $stderrTask.GetAwaiter().GetResult()
        $exitCode = $process.ExitCode; $processPeak = [int64]$process.PeakWorkingSet64
    } catch { $errorText = $_.Exception.ToString() }
    finally { $watch.Stop(); if ($null -ne $process) { $process.Dispose() } }
    Set-Content -LiteralPath (Join-Path $RunDir 'stdout.log') -Value $stdout -Encoding UTF8
    Set-Content -LiteralPath (Join-Path $RunDir 'stderr.log') -Value $stderr -Encoding UTF8
    if (-not [string]::IsNullOrWhiteSpace($errorText)) {
        Set-Content -LiteralPath (Join-Path $RunDir 'harness_exception.txt') -Value $errorText -Encoding UTF8
    }
    return [pscustomobject]@{ name = $Name; exit_code = $exitCode; timed_out = $timedOut; harness_error = $errorText; wall_seconds = $watch.Elapsed.TotalSeconds; process_peak_rss_bytes_diagnostic = $processPeak }
}

function Get-OutputFingerprint([string]$RunDir) {
    $lines = [Collections.Generic.List[string]]::new(); $missing = [Collections.Generic.List[string]]::new()
    foreach ($name in $ScienceOutputNames) {
        $path = Join-Path $RunDir $name
        if (-not (Test-Path -LiteralPath $path)) { $missing.Add($name); continue }
        $lines.Add("$name=" + (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToUpperInvariant())
    }
    return [pscustomobject]@{ ok = $missing.Count -eq 0; fingerprint = Get-TextSha256 ($lines -join "`n"); missing = ($missing -join ';') }
}

function Get-ProfileFingerprint([object]$Profile) {
    $missing = [Collections.Generic.List[string]]::new(); $lines = [Collections.Generic.List[string]]::new()
    foreach ($field in $InvariantProfileFields) {
        if ($null -eq $Profile.PSObject.Properties[$field]) { $missing.Add($field) }
        $lines.Add("$field=" + (Get-TextValue $Profile $field '<missing>'))
    }
    return [pscustomobject]@{ ok = $missing.Count -eq 0; fingerprint = Get-TextSha256 ($lines -join "`n"); missing = ($missing -join ';') }
}

function Test-RunCorrectness([object]$Case, [object]$ProfileArtifact, [object]$EnergyArtifact, [object]$MemoryArtifact, [object]$Outputs, [double]$Theta) {
    $profile = $ProfileArtifact.row; $energy = $EnergyArtifact.row; $memory = $MemoryArtifact.row
    $expected = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_expected'
    $appliedLast = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_last'
    $appliedTotal = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_panel_pairs_applied_total'
    $directSeconds = Get-DoubleValue $profile 'bie_cartesian_full_fmm_direct_seconds'
    $energyValue = Get-DoubleValue $energy 'E_rxn_kJmol'
    $memoryPeak = Get-Int64Value $memory 'peak_working_set_bytes'
    $execution = $Case.exit_code -eq 0 -and -not $Case.timed_out -and [string]::IsNullOrWhiteSpace($Case.harness_error) -and
        $ProfileArtifact.ok -and (Get-TextValue $profile 'status') -eq 'ok' -and
        $EnergyArtifact.ok -and (Get-TextValue $energy 'status') -eq 'ok' -and (Test-FiniteValue $energyValue) -and
        $MemoryArtifact.ok -and (Get-TextValue $memory 'status') -eq 'ok' -and $memoryPeak -gt 0
    $route = (Get-TextValue $profile 'jcp_direct_kernel') -eq 'packed-avx2-vector-exp-node-patch-corrected' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_active') -eq 'true' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_requested') -eq 'cartesian-full-experimental' -and
        (Get-TextValue $profile 'bie_cartesian_full_fmm_effective') -eq 'cartesian-full-experimental' -and
        (Get-Int64Value $profile 'bie_cartesian_full_fmm_order') -eq 1 -and
        (Get-Int64Value $profile 'bie_cartesian_full_fmm_leaf_size') -eq 64 -and
        [Math]::Abs((Get-DoubleValue $profile 'bie_cartesian_full_fmm_opening_theta') - $Theta) -le 1.0e-15 -and
        (Get-TextValue $profile 'fmm_cartesian_pre_audit_mode') -eq 'combined+channel-shadow' -and
        (Get-TextValue $profile 'fmm_cartesian_retry_mode') -eq 'off' -and
        (Get-TextValue $profile 'fallback_used') -eq 'false' -and
        (Get-TextValue $profile 'fmm_m2l_finite_outputs') -eq 'true'
    $coverage = $expected -gt 0 -and $appliedLast -eq $expected -and $appliedTotal -ge $expected -and ($appliedTotal % $expected) -eq 0
    $ilu = (Get-TextValue $profile 'preconditioner_effective') -eq 'fmm_nearfield_ilu0' -and
        (Get-DoubleValue $profile 'preconditioner_ilu_regularization_selected') -eq 0.0 -and
        (Get-Int64Value $profile 'preconditioner_ilu_regularization_shifted_rows') -eq 0 -and
        (Get-TextValue $profile 'preconditioner_ilu_unshifted_failure_stage') -eq 'off'
    $numeric = (Test-FiniteValue $directSeconds) -and $directSeconds -gt 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $profile 'true_residual_rel')) -and
        (Get-DoubleValue $profile 'true_residual_rel') -ge 0.0 -and
        (Test-FiniteValue (Get-DoubleValue $profile 'fmm_cartesian_solution_audit_rel')) -and
        (Get-DoubleValue $profile 'fmm_cartesian_solution_audit_rel') -ge 0.0
    $ok = $execution -and $route -and $coverage -and $ilu -and $numeric -and $Outputs.ok
    return [pscustomobject]@{
        ok = $ok; reason = "execution=$execution route=$route coverage=$coverage ilu=$ilu numeric=$numeric outputs=$($Outputs.ok)"
        expected = $expected; applied_last = $appliedLast; applied_total = $appliedTotal
        direct_seconds = $directSeconds; energy = $energyValue; peak_rss = $memoryPeak
    }
}

function New-Record([object]$Seed, [string]$Role, [string]$Phase, [int]$PairIndex,
                    [object]$Binary, [object]$Case, [object]$Arg,
                    [object]$ProfileArtifact, [object]$EnergyArtifact,
                    [object]$MemoryArtifact, [object]$Output, [object]$ProfileFingerprint,
                    [object]$Correctness) {
    $profile = $ProfileArtifact.row
    $throughput = if ($Correctness.direct_seconds -gt 0.0) { [double]$Correctness.applied_total / $Correctness.direct_seconds } else { [double]::NaN }
    return [pscustomobject][ordered]@{
        fixture_key = $Seed.fixture_key; theta = Format-Double $Seed.theta; phase = $Phase; pair_index = $PairIndex; role = $Role
        binary_path = $Binary.path; binary_sha256 = $Binary.sha256; binary_bytes = $Binary.bytes
        run_dir = $Case.run_dir; arg_sha256 = $Arg.sha256; arg_semantic_sha256 = $Arg.semantic_sha256
        seed_arg_sha256 = $Seed.seed_arg_sha256; pqr_path = $Seed.pqr; pqr_sha256 = $Seed.pqr_sha256; driver_sha256 = $Seed.seed_driver_sha256
        exit_code = $Case.exit_code; timed_out = $Case.timed_out; harness_error = $Case.harness_error
        wall_seconds = Format-Double $Case.wall_seconds; process_peak_rss_bytes_diagnostic = $Case.process_peak_rss_bytes_diagnostic
        correctness_ok = $Correctness.ok; correctness_reason = $Correctness.reason
        output_fingerprint = $Output.fingerprint; output_files_ok = $Output.ok; output_missing = $Output.missing
        profile_invariant_fingerprint = $ProfileFingerprint.fingerprint; profile_fields_ok = $ProfileFingerprint.ok; profile_fields_missing = $ProfileFingerprint.missing
        energy_kj_mol = Format-Double $Correctness.energy; true_residual_rel = Format-Double (Get-DoubleValue $profile 'true_residual_rel')
        direct_seconds = Format-Double $Correctness.direct_seconds; direct_pairs_expected = $Correctness.expected; direct_pairs_applied_last = $Correctness.applied_last; direct_pairs_applied_total = $Correctness.applied_total; direct_panel_pairs_per_second = Format-Double $throughput
        direct_leaf_pairs = Get-Int64Value $profile 'bie_cartesian_full_fmm_direct_pairs'; far_leaf_pairs = Get-Int64Value $profile 'bie_cartesian_full_fmm_far_pairs'
        tree_hash = Get-TextValue $profile 'bie_cartesian_full_fmm_tree_hash'; worklist_hash = Get-TextValue $profile 'bie_cartesian_full_fmm_worklist_hash'; policy_hash = Get-TextValue $profile 'bie_cartesian_full_fmm_policy_hash'
        plan_bytes = Get-Int64Value $profile 'bie_cartesian_full_fmm_plan_bytes'; workspace_bytes = Get-Int64Value $profile 'bie_cartesian_full_fmm_workspace_bytes'; cache_bytes = Get-Int64Value $profile 'bie_cartesian_full_fmm_cache_bytes'
        jcp_direct_kernel = Get-TextValue $profile 'jcp_direct_kernel'; near_self_pairs = Get-Int64Value $profile 'near_self_pairs'; leaf_pair_batch_kernel = Get-TextValue $profile 'leaf_pair_batch_kernel'; leaf_pair_batch_vectorized_sources = Get-Int64Value $profile 'leaf_pair_batch_vectorized_sources'; leaf_pair_batch_vector_lanes = Get-Int64Value $profile 'leaf_pair_batch_vector_lanes'
        node_patch_audit_mode = Get-TextValue $profile 'fmm_cartesian_pre_audit_mode'; channel_shadow_seconds = Format-Double (Get-DoubleValue $profile 'fmm_cartesian_channel_shadow_seconds'); solution_audit_seconds = Format-Double (Get-DoubleValue $profile 'fmm_cartesian_solution_audit_seconds')
        matvec_calls = Get-Int64Value $profile 'matvec_calls'; matvec_seconds_total = Format-Double (Get-DoubleValue $profile 'matvec_seconds_total'); gmres_seconds = Format-Double (Get-DoubleValue $profile 'gmres_seconds'); iterations = Get-Int64Value $profile 'iterations'
        profile_peak_rss_bytes = $Correctness.peak_rss; profile_estimated_bytes_total = Get-Int64Value $MemoryArtifact.row 'estimated_bytes_total'
    }
}

function Get-Median([double[]]$Values) {
    $sorted = @($Values | Sort-Object); $count = $sorted.Count
    if ($count -eq 0) { return [double]::NaN }
    if (($count % 2) -eq 1) { return [double]$sorted[[int]($count / 2)] }
    return (([double]$sorted[$count / 2 - 1]) + ([double]$sorted[$count / 2])) / 2.0
}

function Get-Mad([double[]]$Values) {
    $median = Get-Median $Values
    $deviations = [double[]]@($Values | ForEach-Object { [Math]::Abs($_ - $median) })
    return Get-Median $deviations
}

function Write-Csv([object[]]$Rows, [string]$Path) {
    if ($Rows.Count -eq 0) { Set-Content -LiteralPath $Path -Value '' -Encoding UTF8; return }
    $Rows | Export-Csv -LiteralPath $Path -NoTypeInformation -Encoding UTF8
}

if ([string]::IsNullOrWhiteSpace($WorkRoot)) { $WorkRoot = $DefaultWorkRoot }
$WorkRoot = Assert-SafeWorkRoot $WorkRoot
if (-not (Test-Path -LiteralPath $WorkRoot)) { throw "P3B1 WorkRoot must be created before benchmark: $WorkRoot" }
if ([string]::IsNullOrWhiteSpace($BaselineExe)) { $BaselineExe = Join-Path $WorkRoot 'titan_baseline.exe' }
if ([string]::IsNullOrWhiteSpace($CandidateExe)) { $CandidateExe = Join-Path $WorkRoot 'titan_candidate.exe' }

$normalizedKeys = [Collections.Generic.List[string]]::new(); $seenKeys = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
foreach ($rawKey in $FixtureKeys) {
    if ([string]::IsNullOrWhiteSpace($rawKey)) { throw 'P3B1 FixtureKeys may not contain an empty entry' }
    $key = $rawKey.Trim().ToLowerInvariant()
    if ($key -notin @('1sh1', 'barnase', '1a63')) { throw "Unsupported P3B1 fixture key: $key" }
    if (-not $seenKeys.Add($key)) { throw "Duplicate P3B1 fixture key: $key" }
    $normalizedKeys.Add($key)
}
if ($normalizedKeys.Count -eq 0) { throw 'P3B1 requires at least one fixture key' }
$normalizedThetas = [Collections.Generic.List[double]]::new(); $seenThetas = [Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
foreach ($theta in $Thetas) {
    if (-not (Test-FiniteValue $theta) -or
        ([Math]::Abs($theta - 0.20) -gt 1.0e-15 -and [Math]::Abs($theta - 0.15) -gt 1.0e-15)) {
        throw "P3B1 theta must be exactly .20 or .15: $theta"
    }
    $canonicalTheta = if ([Math]::Abs($theta - 0.20) -le 1.0e-15) { 0.20 } else { 0.15 }
    $thetaKey = $canonicalTheta.ToString('G17', $InvariantCulture)
    if (-not $seenThetas.Add($thetaKey)) { throw "Duplicate P3B1 theta: $thetaKey" }
    $normalizedThetas.Add($canonicalTheta)
}
if ($normalizedThetas.Count -eq 0) { throw 'P3B1 requires at least one theta' }
if ($Pilot) {
    if ($normalizedKeys.Count -ne 1 -or $normalizedKeys[0] -ne '1sh1' -or
        $normalizedThetas.Count -ne 1 -or [Math]::Abs($normalizedThetas[0] - 0.20) -gt 1.0e-15 -or
        $WarmupRepeats -ne 1 -or $MeasuredRepeats -ne 1) {
        throw 'P3B1 -Pilot is restricted to 1SH1 theta=.20 with one warmup and one measured pair'
    }
} elseif ($MeasuredRepeats -lt 7) {
    throw 'P3B1 full benchmark requires at least seven measured paired repeats'
}

$baseline = Assert-ExactBaseline $BaselineExe
$candidate = Get-BinaryProvenance $CandidateExe
if ($candidate.sha256 -eq $baseline.sha256) { throw 'P3B1 candidate binary is byte-identical to baseline' }
if (-not (Test-Path -LiteralPath $EvidenceRoot)) { throw "Missing P2B2D evidence root: $EvidenceRoot" }

$seeds = [Collections.Generic.List[object]]::new()
foreach ($key in $normalizedKeys) { foreach ($theta in $normalizedThetas) { $seeds.Add((Assert-Seed $key $theta)) } }
if ($ValidateOnly) {
    Write-Output "P3B1_VALIDATE_ONLY_PASS baseline=$($baseline.sha256) candidate=$($candidate.sha256) configs=$($seeds.Count)"
    exit 0
}

$artifactNames = @('p3b1_metadata.json', 'p3b1_raw_runs.csv', 'p3b1_correctness.csv', 'p3b1_summary.csv', 'p3b1_decision.csv')
foreach ($name in $artifactNames) {
    if (Test-Path -LiteralPath (Join-Path $WorkRoot $name)) { throw "Refusing to overwrite P3B1 artifact: $name" }
}
$metadata = [ordered]@{
    schema = 'p3b1-fixed-four-lane-packed-scratch-v1'; created_utc = [DateTime]::UtcNow.ToString('o')
    baseline = $baseline; candidate = $candidate; evidence_root = $EvidenceRoot
    warmup_repeats = $WarmupRepeats; measured_repeats = $MeasuredRepeats; timeout_seconds = $CaseTimeoutSeconds; pilot = [bool]$Pilot
    preregistration = [ordered]@{
        direct_speedup_threshold = 0.10; required_fixture_count = 2; max_regression = 0.03; max_rss_increase = 0.05
        fixture_claim_rule = 'both theta020 and theta015 must reach the direct-speed threshold and exceed their MAD noise floor'
        correctness_rule = 'all runs must have exact selected scientific-output and invariant-profile fingerprints within every baseline/candidate pair'
    }
    seeds = @($seeds | ForEach-Object { [ordered]@{ fixture_key=$_.fixture_key; theta=$_.theta; seed_arg_sha256=$_.seed_arg_sha256; pqr=$_.pqr; pqr_sha256=$_.pqr_sha256; driver_sha256=$_.seed_driver_sha256 } })
}
$metadata | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $WorkRoot 'p3b1_metadata.json') -Encoding UTF8

$records = [Collections.Generic.List[object]]::new(); $comparisons = [Collections.Generic.List[object]]::new()
foreach ($seed in $seeds) {
    $thetaTag = if ($seed.theta -eq 0.20) { 'theta020' } else { 'theta015' }
    $configRoot = Join-Path (Join-Path $WorkRoot 'runs') (Join-Path $seed.fixture_key $thetaTag)
    for ($pair = 1; $pair -le ($WarmupRepeats + $MeasuredRepeats); ++$pair) {
        $phase = if ($pair -le $WarmupRepeats) { 'warmup' } else { 'measured' }
        $roles = if (($pair % 2) -eq 1) { @('baseline', 'candidate') } else { @('candidate', 'baseline') }
        $pairRecords = @{}
        foreach ($role in $roles) {
            $binary = if ($role -eq 'baseline') { $baseline } else { $candidate }
            $runDir = Join-Path $configRoot ("{0}_{1:D2}_{2}" -f $phase, $pair, $role)
            New-RunDirectory $runDir
            $arg = New-RunArg $seed $runDir
            $case = Invoke-Run $binary.path $runDir ("$($seed.fixture_key)-$thetaTag-$phase-$pair-$role")
            $profile = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_tabi_profile.csv') 'profile'
            $energy = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_energy.csv') 'energy'
            $memory = Read-RequiredSingleCsv (Join-Path $runDir 'pbe_memory_profile.csv') 'memory'
            $output = Get-OutputFingerprint $runDir
            $profileFingerprint = if ($profile.ok) { Get-ProfileFingerprint $profile.row } else { [pscustomobject]@{ ok=$false; fingerprint=''; missing='profile_missing' } }
            $correctness = Test-RunCorrectness $case $profile $energy $memory $output $seed.theta
            $case | Add-Member -NotePropertyName run_dir -NotePropertyValue $runDir
            $record = New-Record $seed $role $phase $pair $binary $case $arg $profile $energy $memory $output $profileFingerprint $correctness
            $records.Add($record); $pairRecords[$role] = $record
            Write-Output ("P3B1_RUN fixture={0} theta={1} phase={2} pair={3} role={4} correct={5} direct_s={6}" -f $seed.fixture_key, $thetaTag, $phase, $pair, $role, $record.correctness_ok, $record.direct_seconds)
        }
        $baseRecord = $pairRecords['baseline']; $candidateRecord = $pairRecords['candidate']
        $exact = $baseRecord.correctness_ok -and $candidateRecord.correctness_ok -and
            $baseRecord.output_fingerprint -eq $candidateRecord.output_fingerprint -and
            $baseRecord.profile_invariant_fingerprint -eq $candidateRecord.profile_invariant_fingerprint -and
            $baseRecord.arg_semantic_sha256 -eq $candidateRecord.arg_semantic_sha256
        $comparisons.Add([pscustomobject][ordered]@{
            fixture_key=$seed.fixture_key; theta=Format-Double $seed.theta; phase=$phase; pair_index=$pair; exact_parity_ok=$exact
            parity_reason=("baseline_correct={0};candidate_correct={1};output_equal={2};profile_equal={3};semantic_args_equal={4}" -f $baseRecord.correctness_ok,$candidateRecord.correctness_ok,($baseRecord.output_fingerprint -eq $candidateRecord.output_fingerprint),($baseRecord.profile_invariant_fingerprint -eq $candidateRecord.profile_invariant_fingerprint),($baseRecord.arg_semantic_sha256 -eq $candidateRecord.arg_semantic_sha256))
            baseline_direct_seconds=$baseRecord.direct_seconds; candidate_direct_seconds=$candidateRecord.direct_seconds
            baseline_pairs_per_second=$baseRecord.direct_panel_pairs_per_second; candidate_pairs_per_second=$candidateRecord.direct_panel_pairs_per_second
            baseline_peak_rss_bytes=$baseRecord.profile_peak_rss_bytes; candidate_peak_rss_bytes=$candidateRecord.profile_peak_rss_bytes
            baseline_output_fingerprint=$baseRecord.output_fingerprint; candidate_output_fingerprint=$candidateRecord.output_fingerprint
            baseline_profile_fingerprint=$baseRecord.profile_invariant_fingerprint; candidate_profile_fingerprint=$candidateRecord.profile_invariant_fingerprint
        })
    }
}

$summaries = [Collections.Generic.List[object]]::new(); $fixtureClaims = [Collections.Generic.List[object]]::new()
foreach ($key in $normalizedKeys) {
    $thetaClaims = [Collections.Generic.List[object]]::new()
    foreach ($theta in $normalizedThetas) {
        $thetaText = Format-Double $theta
        $baseRows = @($records | Where-Object { $_.fixture_key -eq $key -and $_.theta -eq $thetaText -and $_.phase -eq 'measured' -and $_.role -eq 'baseline' })
        $candidateRows = @($records | Where-Object { $_.fixture_key -eq $key -and $_.theta -eq $thetaText -and $_.phase -eq 'measured' -and $_.role -eq 'candidate' })
        $pairRows = @($comparisons | Where-Object { $_.fixture_key -eq $key -and $_.theta -eq $thetaText -and $_.phase -eq 'measured' })
        [double[]]$baseTimes = @($baseRows | ForEach-Object { [double]$_.direct_seconds })
        [double[]]$candidateTimes = @($candidateRows | ForEach-Object { [double]$_.direct_seconds })
        [double[]]$baseRss = @($baseRows | ForEach-Object { [double]$_.profile_peak_rss_bytes })
        [double[]]$candidateRss = @($candidateRows | ForEach-Object { [double]$_.profile_peak_rss_bytes })
        $baseMedian = Get-Median $baseTimes; $candidateMedian = Get-Median $candidateTimes
        $baseMad = Get-Mad $baseTimes; $candidateMad = Get-Mad $candidateTimes
        $speedup = if ($candidateMedian -gt 0.0) { $baseMedian / $candidateMedian - 1.0 } else { [double]::NaN }
        $signal = [Math]::Abs($baseMedian - $candidateMedian); $noise = [Math]::Max($baseMad, $candidateMad)
        $pairDeltas = [double[]]@($pairRows | ForEach-Object { [double]$_.candidate_direct_seconds - [double]$_.baseline_direct_seconds })
        $parityOk = $pairRows.Count -eq $MeasuredRepeats -and @($pairRows | Where-Object { -not $_.exact_parity_ok }).Count -eq 0
        $rssOk = (Get-Median $candidateRss) -le 1.05 * (Get-Median $baseRss)
        $speedOk = (Test-FiniteValue $speedup) -and $speedup -ge 0.10 -and $signal -gt $noise
        $regressionOk = (Test-FiniteValue $speedup) -and $speedup -ge -0.03
        $summary = [pscustomobject][ordered]@{
            fixture_key=$key; theta=$thetaText; measured_pairs=$pairRows.Count; parity_ok=$parityOk
            baseline_direct_seconds_median=Format-Double $baseMedian; baseline_direct_seconds_mad=Format-Double $baseMad; baseline_direct_seconds_min=Format-Double ($baseTimes | Measure-Object -Minimum | Select-Object -ExpandProperty Minimum); baseline_direct_seconds_max=Format-Double ($baseTimes | Measure-Object -Maximum | Select-Object -ExpandProperty Maximum)
            candidate_direct_seconds_median=Format-Double $candidateMedian; candidate_direct_seconds_mad=Format-Double $candidateMad; candidate_direct_seconds_min=Format-Double ($candidateTimes | Measure-Object -Minimum | Select-Object -ExpandProperty Minimum); candidate_direct_seconds_max=Format-Double ($candidateTimes | Measure-Object -Maximum | Select-Object -ExpandProperty Maximum)
            paired_delta_seconds_median=Format-Double (Get-Median $pairDeltas); direct_speedup_fraction=Format-Double $speedup; signal_seconds=Format-Double $signal; noise_mad_seconds=Format-Double $noise; speed_gate_ok=$speedOk; regression_gate_ok=$regressionOk
            baseline_pairs_per_second_median=Format-Double (Get-Median ([double[]]@($baseRows | ForEach-Object { [double]$_.direct_panel_pairs_per_second }))); candidate_pairs_per_second_median=Format-Double (Get-Median ([double[]]@($candidateRows | ForEach-Object { [double]$_.direct_panel_pairs_per_second })))
            baseline_peak_rss_bytes_median=[int64](Get-Median $baseRss); candidate_peak_rss_bytes_median=[int64](Get-Median $candidateRss); rss_gate_ok=$rssOk
        }
        $summaries.Add($summary); $thetaClaims.Add($summary)
    }
    $fixtureClaims.Add([pscustomobject]@{ fixture_key=$key; claim_ok=(@($thetaClaims | Where-Object { $_.speed_gate_ok }).Count -eq $normalizedThetas.Count); all_regression_ok=(@($thetaClaims | Where-Object { -not $_.regression_gate_ok }).Count -eq 0); all_rss_ok=(@($thetaClaims | Where-Object { -not $_.rss_gate_ok }).Count -eq 0); all_parity_ok=(@($thetaClaims | Where-Object { -not $_.parity_ok }).Count -eq 0) })
}
$allPairsExact = @($comparisons | Where-Object { -not $_.exact_parity_ok }).Count -eq 0
$fixturesSpeed = @($fixtureClaims | Where-Object { $_.claim_ok }).Count
$noRegression = @($summaries | Where-Object { -not $_.regression_gate_ok }).Count -eq 0
$rssOk = @($summaries | Where-Object { -not $_.rss_gate_ok }).Count -eq 0
$fixtureParity = @($fixtureClaims | Where-Object { -not $_.all_parity_ok }).Count -eq 0
$go = $allPairsExact -and $fixtureParity -and $fixturesSpeed -ge 2 -and $noRegression -and $rssOk
$decision = if ($Pilot) { if ($allPairsExact -and $fixtureParity) { 'PILOT_PASS' } else { 'PILOT_FAIL' } } elseif ($go) { 'GO_P3B1' } else { 'NO_GO_P3B1' }
$reason = if ($Pilot) { "pilot exact=$allPairsExact fixture_parity=$fixtureParity" } elseif ($go) { 'go_preregistered_direct_gate' } else { "no_go exact=$allPairsExact fixture_parity=$fixtureParity speed_fixtures=$fixturesSpeed no_regression=$noRegression rss=$rssOk" }

Write-Csv @($records) (Join-Path $WorkRoot 'p3b1_raw_runs.csv')
Write-Csv @($comparisons) (Join-Path $WorkRoot 'p3b1_correctness.csv')
Write-Csv @($summaries) (Join-Path $WorkRoot 'p3b1_summary.csv')
Write-Csv @([pscustomobject][ordered]@{ decision = $decision; reason = $reason; direct_speed_fixture_count = $fixturesSpeed; all_exact_parity = $allPairsExact; no_regression = $noRegression; rss_ok = $rssOk; baseline_sha256 = $baseline.sha256; candidate_sha256 = $candidate.sha256; completed_utc = [DateTime]::UtcNow.ToString('o') }) (Join-Path $WorkRoot 'p3b1_decision.csv')
Write-Output "P3B1_DECISION=$decision reason=$reason"
