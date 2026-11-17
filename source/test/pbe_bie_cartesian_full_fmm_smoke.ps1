param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_bie_cartesian_full_fmm'),
    [int]$Threads = 1,
    [int]$CaseTimeoutSeconds = 900,
    [ValidateSet(1, 4)][int[]]$Orders = @(1, 4),
    [double[]]$SaltMolar = @(0.0, 0.15),
    [ValidateSet('off', 'channel-shadow')]
    [string]$ErrorAudit = 'off',
    # node-patch builds curved dual cells, which main solves only with the
    # classic treecode (curved_node_patch_requires_direct_or_classic_treecode);
    # curved cells in the Cartesian FMM belong to the bie-ses route.
    [ValidateSet('centroid', 'node-patch')]
    [string[]]$Discretizations = @('centroid')
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_driver_topology.ps1')
Set-StrictMode -Version Latest
$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Test-Finite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Get-CsvFieldCount([string]$Line) {
    if ($Line.Length -eq 0) { return 0 }
    $count = 1
    $inQuotes = $false
    for ($index = 0; $index -lt $Line.Length; ++$index) {
        if ($Line[$index] -eq '"') {
            if ($inQuotes -and $index + 1 -lt $Line.Length -and
                $Line[$index + 1] -eq '"') {
                ++$index
            } else {
                $inQuotes = -not $inQuotes
            }
        } elseif ($Line[$index] -eq ',' -and -not $inQuotes) {
            ++$count
        }
    }
    if ($inQuotes) { throw 'Unterminated quoted CSV field in TABI profile.' }
    return $count
}

function Assert-UnderWorkRoot([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    $root = [IO.Path]::GetFullPath($WorkRoot).TrimEnd('\') + '\'
    if (-not ($full + '\').StartsWith($root,
            [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside Cartesian full FMM smoke root: $full"
    }
}

function New-CleanDirectory([string]$Path) {
    Assert-UnderWorkRoot $Path
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing TITAN executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
$libRoot = Join-Path $RuntimeRepoRoot 'lib_ff'
if (-not (Test-Path -LiteralPath (Join-Path $libRoot 'ffAMBER99SBnonbonded.dat'))) {
    $libRoot = Join-Path $RuntimeRepoRoot 'example\1\lib_ff'
}
if (-not (Test-Path -LiteralPath (Join-Path $libRoot 'ffAMBER99SBnonbonded.dat'))) {
    throw "Missing force-field directory: $libRoot"
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$rows = New-Object Collections.Generic.List[object]
foreach ($order in $Orders) {
foreach ($salt in $SaltMolar) {
foreach ($discretization in $Discretizations) {
    if (-not (Test-Finite $salt) -or $salt -lt 0.0) {
        throw "Invalid salt concentration: $salt"
    }
    if ($ErrorAudit -eq 'channel-shadow' -and
        ($order -ne 1 -or $salt -le 0.0 -or $discretization -ne 'node-patch')) {
        throw 'Cartesian channel-shadow smoke requires order 1, coupled salt, and node-patch.'
    }
    $openingTheta = if ($order -eq 1) { 0.05 } else { 0.55 }
    $saltTag = $salt.ToString('G15', $InvariantCulture).Replace('.', 'p')
    $discretizationTag = $discretization.Replace('-', '')
    $runDir = Join-Path $WorkRoot "order${order}_salt${saltTag}_${discretizationTag}"
    New-CleanDirectory $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    ) | Set-Content -LiteralPath $pdbPath -Encoding ASCII
    $pqrPath = Join-Path $runDir 'four_atom.pqr'
    @(
        'ATOM      1  C   GLY A   1      -6.000   0.000   0.000  1.0000  1.7000',
        'ATOM      2  O   GLY A   1      -2.000   0.000   0.000 -1.0000  1.5000',
        'ATOM      3  N   GLY A   2       2.000   0.000   0.000  1.0000  1.5500',
        'ATOM      4  O   GLY A   2       6.000   0.000   0.000 -1.0000  1.5000',
        'END'
    ) | Set-Content -LiteralPath $pqrPath -Encoding ASCII

    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $libRoot).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
    @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        "-common_topo=$(Convert-ToTitanPath (Write-PbeDriverTopology $pdbPath))",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=bie',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-max-iter=400',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        "-pbe-salt=$($salt.ToString('G15', $InvariantCulture))",
        '-pbe-temp=298.15',
        '-pbe-probe=1.4',
        '-pbe-timing=detail',
        '-pbe-bie-fast-operator=p1',
        "-pbe-bie-discretization=$discretization",
        '-pbe-bie-fmm-engine=cartesian-full-experimental',
        "-pbe-bie-fmm-cartesian-order=$order",
        "-pbe-bie-fmm-cartesian-error-audit=$ErrorAudit",
        '-pbe-bie-surface-points=8',
        '-pbe-bie-tree-leaf-size=8',
        "-pbe-bie-opening-theta=$($openingTheta.ToString('G15', $InvariantCulture))",
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off',
        '-pbe-bie-precond=diagonal',
        '-pbe-bie-gmres-restart=64',
        '-pbe-bie-matvec-audit=full-small',
        '-pbe-bie-matvec-audit-sample-count=512',
        "-nt=$Threads"
    ) | Set-Content -LiteralPath $argPath -Encoding ASCII

    $stdoutPath = Join-Path $runDir 'stdout.log'
    $stderrPath = Join-Path $runDir 'stderr.log'
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
        throw "Cartesian full FMM order $order timed out."
    }
    $process.WaitForExit()
    $stdout = $stdoutTask.GetAwaiter().GetResult()
    $stderr = $stderrTask.GetAwaiter().GetResult()
    Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding UTF8
    Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding UTF8
    $exitCode = $process.ExitCode
    $process.Dispose()
    if ($exitCode -ne 0) {
        throw "Cartesian full FMM order $order failed with exit ${exitCode}:`n$stdout`n$stderr"
    }

    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    if ($energy.status -ne 'ok') {
        throw "Cartesian full FMM order $order status=$($energy.status)"
    }
    if ($profile.fmm_local_gmres_active -ne 'true') {
        throw "Cartesian full FMM order $order did not use the prepared GMRES operator."
    }
    if ($profile.fmm_m2l_mode -ne 'cartesian-full-experimental') {
        throw "Cartesian full FMM order $order effective mode=$($profile.fmm_m2l_mode)"
    }
    if ([int]$profile.fmm_m2l_expansion_order -ne $order) {
        throw "Cartesian full FMM requested order $order but reported $($profile.fmm_m2l_expansion_order)."
    }
    if ($profile.bie_cartesian_full_fmm_active -ne 'true' -or
        $profile.bie_cartesian_full_fmm_requested -ne 'cartesian-full-experimental' -or
        $profile.bie_cartesian_full_fmm_effective -ne 'cartesian-full-experimental') {
        throw "Cartesian full FMM order $order requested/effective telemetry is inconsistent."
    }
    foreach ($hashColumn in @('bie_cartesian_full_fmm_tree_hash',
            'bie_cartesian_full_fmm_worklist_hash',
            'bie_cartesian_full_fmm_policy_hash')) {
        if ([uint64]$profile.$hashColumn -eq 0) {
            throw "Cartesian full FMM order $order missing $hashColumn."
        }
    }
    $effectiveLeaf = [int]$profile.bie_cartesian_full_fmm_leaf_size
    $effectiveTheta = Convert-ToDouble $profile.bie_cartesian_full_fmm_opening_theta
    if ($effectiveLeaf -ne 8 -or [Math]::Abs($effectiveTheta - $openingTheta) -gt 1e-12) {
        throw "Cartesian full FMM order $order effective tree telemetry drifted: leaf=$effectiveLeaf theta=$effectiveTheta"
    }
    $expectedDirectPanels = [int64]$profile.bie_cartesian_full_fmm_direct_panel_pairs_expected
    $appliedDirectPanels = [int64]$profile.bie_cartesian_full_fmm_direct_panel_pairs_applied_last
    $totalAppliedDirectPanels = [int64]$profile.bie_cartesian_full_fmm_direct_panel_pairs_applied_total
    if ($expectedDirectPanels -le 0 -or $appliedDirectPanels -ne $expectedDirectPanels -or
        $totalAppliedDirectPanels -lt $expectedDirectPanels -or
        ($totalAppliedDirectPanels % $expectedDirectPanels) -ne 0) {
        throw "Cartesian full FMM order $order direct coverage telemetry is inconsistent: expected=$expectedDirectPanels last=$appliedDirectPanels total=$totalAppliedDirectPanels"
    }
    $cacheEntries = [int64]$profile.bie_cartesian_full_fmm_cache_entries
    $cacheBytes = [int64]$profile.bie_cartesian_full_fmm_cache_bytes
    $planBytes = [int64]$profile.bie_cartesian_full_fmm_plan_bytes
    if ($cacheEntries -gt 0 -and ($cacheBytes -le 0 -or
            $cacheBytes -gt $planBytes)) {
        throw "Cartesian full FMM order $order cache byte telemetry is inconsistent."
    }
    foreach ($parityColumn in @('bie_cartesian_full_fmm_m2m_parity_error',
            'bie_cartesian_full_fmm_l2l_parity_error')) {
        $parity = Convert-ToDouble $profile.$parityColumn
        if (-not (Test-Finite $parity) -or $parity -gt 1e-11) {
            throw "Cartesian full FMM order $order parity gate failed for ${parityColumn}: $parity"
        }
    }
    $audit = if ($discretization -eq 'node-patch') {
        if ($profile.fmm_cartesian_pre_audit_enabled -ne 'true' -or
            $profile.fmm_cartesian_pre_audit_finite -ne 'true' -or
            [int64]$profile.fmm_cartesian_pre_audit_rows -le 0 -or
            [uint64]$profile.fmm_cartesian_pre_audit_probe_hash -eq 0 -or
            [uint64]$profile.fmm_cartesian_pre_audit_target_hash -eq 0) {
            throw "Cartesian full FMM node-patch pre-audit telemetry is incomplete."
        }
        Convert-ToDouble $profile.fmm_cartesian_pre_audit_aggregate_rel
    } else {
        Convert-ToDouble $profile.tree_vs_direct_matvec_rel
    }
    if (-not (Test-Finite $audit) -or $audit -gt 1e-4) {
        throw "Cartesian full FMM order $order direct matvec audit failed: $audit"
    }
    if ($discretization -eq 'node-patch') {
        # Compatibility contract: node-patch's legacy audit remains the
        # fixed pre-GMRES probe, even when channel-shadow adds a distinct
        # post-GMRES solution audit.
        $legacyAudit = Convert-ToDouble $profile.tree_vs_direct_matvec_rel
        if ([Math]::Abs($legacyAudit - $audit) /
                [Math]::Max(1.0, [Math]::Abs($audit)) -gt 1.0e-12 -or
            [int64]$profile.matvec_audit_rows -ne
                [int64]$profile.fmm_cartesian_pre_audit_rows -or
            [Math]::Abs((Convert-ToDouble $profile.matvec_audit_seconds) -
                (Convert-ToDouble $profile.fmm_cartesian_pre_audit_seconds)) -gt 1.0e-12) {
            throw "Cartesian node-patch legacy pre-audit compatibility drifted."
        }
    }
    if ($ErrorAudit -eq 'off') {
        $profileLines = [IO.File]::ReadAllLines((Join-Path $runDir 'pbe_tabi_profile.csv'))
        if ($profileLines.Count -ne 2 -or
            (Get-CsvFieldCount $profileLines[0]) -ne
                (Get-CsvFieldCount $profileLines[1])) {
            throw 'Cartesian retry foundation profile header/value column count mismatch.'
        }
        foreach ($retryColumn in @(
                'fmm_cartesian_retry_mode',
                'fmm_cartesian_retry_attempts_started',
                'fmm_cartesian_retry_attempts_completed',
                'fmm_cartesian_retry_attempt1_classification',
                'fmm_cartesian_retry_attempt1_cartesian_full_requested',
                'fmm_cartesian_retry_attempt1_solution_audit_direct_l2',
                'fmm_cartesian_retry_attempt2_classification')) {
            if ($profile.PSObject.Properties.Name -notcontains $retryColumn) {
                throw "Cartesian retry foundation profile column is missing: $retryColumn"
            }
        }
        if ($profile.fmm_cartesian_retry_mode -ne 'off' -or
            $profile.fmm_cartesian_retry_enabled -ne 'false' -or
            [int64]$profile.fmm_cartesian_retry_attempts_started -ne 0 -or
            [int64]$profile.fmm_cartesian_retry_attempts_completed -ne 0 -or
            [int]$profile.fmm_cartesian_retry_accepted_attempt -ne 0 -or
            $profile.fmm_cartesian_retry_trigger_reason -ne 'off' -or
            $profile.fmm_cartesian_retry_final_reason -ne 'off' -or
            $profile.fmm_cartesian_retry_attempt1_classification -ne 'off' -or
            $profile.fmm_cartesian_retry_attempt2_classification -ne 'off') {
            throw 'Cartesian retry foundation default-off profile state drifted.'
        }
        foreach ($attemptIndex in 1..2) {
            $prefix = "fmm_cartesian_retry_attempt${attemptIndex}_"
            if ($profile."${prefix}cartesian_full_active" -ne 'false' -or
                $profile."${prefix}cartesian_full_finite" -ne 'false' -or
                $profile."${prefix}cartesian_full_requested" -ne 'off' -or
                $profile."${prefix}cartesian_full_effective" -ne 'off' -or
                [int]$profile."${prefix}cartesian_full_order" -ne 0 -or
                [int]$profile."${prefix}cartesian_full_leaf_size" -ne 0 -or
                $profile."${prefix}solution_audit_status" -ne 'off' -or
                [double]$profile."${prefix}solution_audit_abs_l2" -ne 0.0 -or
                [double]$profile."${prefix}solution_audit_direct_l2" -ne 0.0 -or
                [double]$profile."${prefix}solution_audit_seconds" -ne 0.0) {
                throw "Cartesian retry foundation attempt $attemptIndex default profile state drifted."
            }
        }
        if ($discretization -eq 'node-patch' -and
            $profile.fmm_cartesian_pre_audit_mode -ne 'combined') {
            throw "Cartesian default pre-audit mode drifted: $($profile.fmm_cartesian_pre_audit_mode)"
        }
        if ($profile.fmm_cartesian_solution_audit_enabled -ne 'false' -or
            $profile.fmm_cartesian_solution_audit_status -ne 'off' -or
            [int64]$profile.fmm_cartesian_solution_audit_targets -ne 0 -or
            [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_expected -ne 0 -or
            [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_applied -ne 0) {
            throw 'Cartesian default-off unexpectedly enabled solution audit telemetry.'
        }
        foreach ($channel in 1..4) {
            if ($profile."fmm_cartesian_channel_k${channel}_status" -ne 'off') {
                throw "Cartesian default channel $channel was unexpectedly enabled."
            }
        }
    } else {
        if ($profile.fmm_cartesian_pre_audit_mode -ne 'combined+channel-shadow' -or
            [uint64]$profile.fmm_cartesian_channel_probe_hash -eq 0 -or
            [int64]$profile.fmm_cartesian_channel_direct_pairs_expected -le 0 -or
            [int64]$profile.fmm_cartesian_channel_direct_pairs_applied -ne
                [int64]$profile.fmm_cartesian_channel_direct_pairs_expected -or
            (Convert-ToDouble $profile.fmm_cartesian_channel_shadow_seconds) -le 0.0) {
            throw 'Cartesian channel-shadow telemetry is incomplete.'
        }
        if ($profile.fmm_cartesian_solution_audit_enabled -ne 'true' -or
            $profile.fmm_cartesian_solution_audit_finite -ne 'true' -or
            $profile.fmm_cartesian_solution_audit_status -notin @('relative', 'absolute_only') -or
            [int64]$profile.fmm_cartesian_solution_audit_targets -ne
                [int64]$profile.fmm_cartesian_pre_audit_targets -or
            [int64]$profile.fmm_cartesian_solution_audit_rows -ne
                2 * [int64]$profile.fmm_cartesian_solution_audit_targets -or
            [uint64]$profile.fmm_cartesian_solution_audit_target_hash -eq 0 -or
            [uint64]$profile.fmm_cartesian_solution_audit_target_hash -ne
                [uint64]$profile.fmm_cartesian_pre_audit_target_hash -or
            [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_expected -le 0 -or
            [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_applied -ne
                [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_expected -or
            (Convert-ToDouble $profile.fmm_cartesian_solution_audit_seconds) -le 0.0 -or
            (-not (Test-Finite (Convert-ToDouble $profile.fmm_cartesian_solution_audit_abs_l2))) -or
            (-not (Test-Finite (Convert-ToDouble $profile.fmm_cartesian_solution_audit_direct_l2)))) {
            throw 'Cartesian sampled solution-audit telemetry is incomplete.'
        }
        if ($profile.fmm_cartesian_solution_audit_status -eq 'relative' -and
            ((-not (Test-Finite (Convert-ToDouble $profile.fmm_cartesian_solution_audit_rel))) -or
             $profile.fmm_cartesian_solution_audit_relative_available -ne 'true')) {
            throw 'Cartesian sampled solution-audit relative metric is invalid.'
        }
        if ($profile.fmm_cartesian_solution_audit_status -eq 'absolute_only' -and
            $profile.fmm_cartesian_solution_audit_relative_available -ne 'false') {
            throw 'Cartesian sampled solution-audit denominator state is inconsistent.'
        }
        foreach ($channel in 1..4) {
            $prefix = "fmm_cartesian_channel_k${channel}_"
            if ($profile."${prefix}finite" -ne 'true' -or
                [int64]$profile."${prefix}samples" -ne
                    [int64]$profile.fmm_cartesian_pre_audit_targets -or
                $profile."${prefix}status" -notin @('relative', 'absolute_only')) {
                throw "Cartesian channel-shadow K$channel telemetry is invalid."
            }
            foreach ($suffix in @('abs_l2', 'rms', 'reference_l2')) {
                if (-not (Test-Finite (Convert-ToDouble $profile."${prefix}${suffix}"))) {
                    throw "Cartesian channel-shadow K$channel $suffix is non-finite."
                }
            }
            $relativeError = $profile."${prefix}rel"
            if ($profile."${prefix}status" -eq 'relative' -and
                ((-not (Test-Finite (Convert-ToDouble $relativeError))) -or
                 (Convert-ToDouble $relativeError) -gt 1.0e-4)) {
                throw "Cartesian channel-shadow K$channel mapping/diagonal audit failed: rel=$relativeError"
            }
        }
    }
    $trueResidual = Convert-ToDouble $profile.true_residual_rel
    if (-not (Test-Finite $trueResidual) -or $trueResidual -gt 1e-4) {
        throw "Cartesian full FMM order $order true residual failed: $trueResidual"
    }
    $rows.Add([pscustomobject]@{
        order = $order
        error_audit = $ErrorAudit
        salt_molar = $salt
        discretization = $discretization
        status = $energy.status
        panels = [int64]$profile.panels
        iterations = [int64]$profile.iterations
        true_residual = $trueResidual
        direct_matvec_relative_l2 = $audit
        legacy_pre_audit_rel = Convert-ToDouble $profile.tree_vs_direct_matvec_rel
        legacy_pre_audit_rows = [int64]$profile.matvec_audit_rows
        pre_audit_rel = Convert-ToDouble $profile.fmm_cartesian_pre_audit_aggregate_rel
        pre_audit_probe_hash = [uint64]$profile.fmm_cartesian_pre_audit_probe_hash
        pre_audit_target_hash = [uint64]$profile.fmm_cartesian_pre_audit_target_hash
        solution_audit_enabled = $profile.fmm_cartesian_solution_audit_enabled
        solution_audit_status = $profile.fmm_cartesian_solution_audit_status
        solution_audit_rel = Convert-ToDouble $profile.fmm_cartesian_solution_audit_rel
        solution_audit_seconds = Convert-ToDouble $profile.fmm_cartesian_solution_audit_seconds
        solution_audit_target_hash = [uint64]$profile.fmm_cartesian_solution_audit_target_hash
        solution_audit_rows = [int64]$profile.fmm_cartesian_solution_audit_rows
        solution_audit_pairs_expected = [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_expected
        solution_audit_pairs_applied = [int64]$profile.fmm_cartesian_solution_audit_direct_pairs_applied
        channel_k1_rel = Convert-ToDouble $profile.fmm_cartesian_channel_k1_rel
        channel_k2_rel = Convert-ToDouble $profile.fmm_cartesian_channel_k2_rel
        channel_k3_rel = Convert-ToDouble $profile.fmm_cartesian_channel_k3_rel
        channel_k4_rel = Convert-ToDouble $profile.fmm_cartesian_channel_k4_rel
        tree_hash = [uint64]$profile.bie_cartesian_full_fmm_tree_hash
        worklist_hash = [uint64]$profile.bie_cartesian_full_fmm_worklist_hash
        policy_hash = [uint64]$profile.bie_cartesian_full_fmm_policy_hash
        matvec_calls = [int64]$profile.matvec_calls
        fmm_gmres_apply_calls = [int64]$profile.fmm_local_gmres_matvec_calls
        direct_panel_pairs_applied_total = [int64]$profile.bie_cartesian_full_fmm_direct_panel_pairs_applied_total
        energy_kj_mol = Convert-ToDouble $energy.E_rxn_kJmol
        plan_seconds = Convert-ToDouble $profile.fmm_local_prepare_seconds
        matvec_seconds = Convert-ToDouble $profile.fmm_local_gmres_apply_seconds_total
        run_dir = $runDir
    }) | Out-Null
}
}
}

$rows.ToArray() | Export-Csv -LiteralPath (Join-Path $WorkRoot 'matrix.csv') `
    -NoTypeInformation -Encoding ASCII
Write-Host "PASS BIE Cartesian full FMM orders $($Orders -join ', ') for salt $($SaltMolar -join ', ') and discretizations $($Discretizations -join ', '): $WorkRoot"
$rows.ToArray() | Format-Table -AutoSize
