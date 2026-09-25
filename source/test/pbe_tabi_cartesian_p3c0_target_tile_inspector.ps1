[CmdletBinding()]
param(
    [string]$InspectorExe = '',
    [string]$WorkRoot = '',
    [string[]]$FixtureKeys = @('1sh1', 'barnase', '1a63'),
    [double[]]$Thetas = @(0.20, 0.15)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$AllowedWorkRootBase = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
$EvidenceRoot = Join-Path $AllowedWorkRootBase 'p2b2d_cartesian_retry_three_fixture_final_20260802'
if ([string]::IsNullOrWhiteSpace($WorkRoot)) {
    $WorkRoot = Join-Path $AllowedWorkRootBase 'p3c0_target_tile_hoist_audit_20260802'
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
    $base = Get-PathForComparison $AllowedWorkRootBase
    $repo = Get-PathForComparison $RepoRoot
    $drive = Get-PathForComparison ([IO.Path]::GetPathRoot($full))
    if ($full.Equals($base, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($repo, [StringComparison]::OrdinalIgnoreCase) -or
        $full.Equals($drive, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing unsafe P3C0 WorkRoot: $full"
    }
    $prefix = $base + [IO.Path]::DirectorySeparatorChar
    if (-not $full.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing P3C0 WorkRoot outside test/tmp: $full"
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
    $parsed = 0.0
    if ([double]::TryParse((Get-TextValue $Row $Name ''),
            [Globalization.NumberStyles]::Float, $InvariantCulture,
            [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-ArgValue([string[]]$ArgumentLines, [string]$Key) {
    $prefix = "-$Key="
    $found = @($ArgumentLines | Where-Object {
        $_.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)
    })
    if ($found.Count -ne 1) {
        throw "Expected exactly one $prefix line in P3C0 seed arg.in; got $($found.Count)"
    }
    return $found[0].Substring($prefix.Length)
}

function Assert-FixtureInput([string]$FixtureKey, [double]$Theta) {
    $thetaTag = if ([Math]::Abs($Theta - 0.20) -le 1.0e-15) {
        'A_channel_shadow_theta020'
    } elseif ([Math]::Abs($Theta - 0.15) -le 1.0e-15) {
        'C_channel_shadow_theta015'
    } else {
        throw "P3C0 supports only theta .20 or .15; got $Theta"
    }
    $seedDir = Join-Path (Join-Path $EvidenceRoot $FixtureKey) $thetaTag
    $seedArg = Join-Path $seedDir 'arg.in'
    $seedProfile = Join-Path $seedDir 'pbe_tabi_profile.csv'
    $seedDriver = Join-Path $seedDir 'driver.pdb'
    foreach ($required in @($seedArg, $seedProfile, $seedDriver)) {
        if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
            throw "Missing approved P2B2D seed artifact: $required"
        }
    }
    $argumentLines = @([IO.File]::ReadAllLines($seedArg))
    $checks = @{
        'pbe' = 'bie'; 'pbe-bie-discretization' = 'node-patch';
        'pbe-bie-node-patch-operator-audit' = 'on'; 'pbe-bie-sum' = 'treecode';
        'pbe-bie-fast-operator' = 'p1';
        'pbe-bie-fmm-engine' = 'cartesian-full-experimental';
        'pbe-bie-fmm-cartesian-order' = '1';
        'pbe-bie-fmm-cartesian-error-audit' = 'channel-shadow';
        'pbe-bie-tree-leaf-size' = '64';
        'pbe-bie-precond' = 'fmm-nearfield-ilu0';
        'pbe-salt' = '0.15'; 'nt' = '4'
    }
    foreach ($key in $checks.Keys) {
        if ((Get-ArgValue $argumentLines $key) -ne $checks[$key]) {
            throw "P3C0 seed route mismatch for $FixtureKey theta=${Theta}: $key"
        }
    }
    $seedTheta = 0.0
    if (-not [double]::TryParse((Get-ArgValue $argumentLines 'pbe-bie-opening-theta'),
            [Globalization.NumberStyles]::Float, $InvariantCulture,
            [ref]$seedTheta) -or [Math]::Abs($seedTheta - $Theta) -gt 1.0e-15) {
        throw "P3C0 seed theta mismatch for $FixtureKey"
    }
    return [pscustomobject]@{
        fixture_key = $FixtureKey; theta = $Theta; seed_dir = $seedDir
        seed_arg = $seedArg; seed_profile = $seedProfile; seed_driver = $seedDriver
        argument_lines = $argumentLines
    }
}

function New-RunArg([object]$Seed, [string]$RunDir) {
    Copy-Item -LiteralPath $Seed.seed_driver -Destination (Join-Path $RunDir 'driver.pdb') -ErrorAction Stop
    $driver = ((Join-Path $RunDir 'driver.pdb') -replace '\\', '/')
    $lines = [Collections.Generic.List[string]]::new()
    foreach ($line in $Seed.argument_lines) {
        if ($line.StartsWith('-f=', [StringComparison]::OrdinalIgnoreCase)) {
            $lines.Add("-f=$driver")
        } else {
            $lines.Add($line)
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    [IO.File]::WriteAllLines($argPath, $lines, [Text.Encoding]::ASCII)
    return $argPath
}

function Read-OneRow([string]$Path, [string]$Label) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Missing ${Label}: $Path"
    }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -ne 1) {
        throw "Expected one $Label row; got $($rows.Count): $Path"
    }
    return $rows[0]
}

function Assert-ProfileIdentity([object]$Reference, [object]$Profile,
                                [object]$Inspection, [double]$Theta) {
    $required = @{
        'status' = 'ok';
        'jcp_direct_kernel' = 'packed-avx2-vector-exp-node-patch-corrected';
        'bie_cartesian_full_fmm_active' = 'true';
        'bie_cartesian_full_fmm_requested' = 'cartesian-full-experimental';
        'bie_cartesian_full_fmm_effective' = 'cartesian-full-experimental';
        'bie_cartesian_full_fmm_order' = '1';
        'bie_cartesian_full_fmm_leaf_size' = '64'
    }
    foreach ($name in $required.Keys) {
        if ((Get-TextValue $Profile $name) -ne $required[$name]) {
            throw "P3C0 profile route mismatch: $name"
        }
    }
    if ([Math]::Abs((Get-DoubleValue $Profile 'bie_cartesian_full_fmm_opening_theta') - $Theta) -gt 1.0e-15) {
        throw 'P3C0 profile theta mismatch'
    }
    foreach ($name in @('panels', 'bie_cartesian_full_fmm_tree_hash',
            'bie_cartesian_full_fmm_worklist_hash',
            'bie_cartesian_full_fmm_policy_hash',
            'bie_cartesian_full_fmm_direct_pairs',
            'bie_cartesian_full_fmm_direct_panel_pairs_expected')) {
        if ((Get-TextValue $Profile $name) -ne (Get-TextValue $Reference $name)) {
            throw "P3C0 baseline identity mismatch: $name"
        }
    }
    if ((Get-TextValue $Inspection 'available') -ne 'true') {
        throw 'P3C0 test-only inspector reports unavailable packed plan'
    }
    foreach ($pair in @(
            @('total_pairs', 'bie_cartesian_full_fmm_direct_panel_pairs_expected'),
            @('plan_direct_panel_pairs', 'bie_cartesian_full_fmm_direct_panel_pairs_expected'),
            @('plan_direct_pairs', 'bie_cartesian_full_fmm_direct_pairs'),
            @('plan_tree_hash', 'bie_cartesian_full_fmm_tree_hash'),
            @('plan_worklist_hash', 'bie_cartesian_full_fmm_worklist_hash'))) {
        if ((Get-TextValue $Inspection $pair[0]) -ne (Get-TextValue $Profile $pair[1])) {
            throw "P3C0 inspector/profile mismatch: $($pair[0])"
        }
    }
    if ([Math]::Abs((Get-DoubleValue $Inspection 'theta') - $Theta) -gt 1.0e-15 -or
        (Get-TextValue $Inspection 'leaf_size') -ne '64') {
        throw 'P3C0 inspector configuration mismatch'
    }
}

$WorkRoot = Assert-SafeWorkRoot $WorkRoot
if ([string]::IsNullOrWhiteSpace($InspectorExe)) {
    $InspectorExe = Join-Path $WorkRoot 'titan_p3c0_inspector.exe'
}
$InspectorExe = [IO.Path]::GetFullPath($InspectorExe)
if (-not (Test-Path -LiteralPath $InspectorExe -PathType Leaf)) {
    throw "Missing P3C0 test-only inspector binary: $InspectorExe"
}
if ($FixtureKeys.Count -eq 0 -or @($FixtureKeys | Where-Object {
        [string]::IsNullOrWhiteSpace($_) }).Count -ne 0) {
    throw 'P3C0 FixtureKeys must be non-empty'
}
if (@($FixtureKeys | Sort-Object -Unique).Count -ne $FixtureKeys.Count) {
    throw 'P3C0 FixtureKeys contains duplicates'
}
foreach ($fixtureKey in $FixtureKeys) {
    if ($fixtureKey -notin @('1sh1', 'barnase', '1a63')) {
        throw "Unsupported P3C0 fixture key: $fixtureKey"
    }
}
if ($Thetas.Count -ne 2 -or @($Thetas | Sort-Object) -join ',' -ne '0.15,0.2') {
    throw 'P3C0 Thetas must be exactly 0.20 and 0.15'
}

$runsRoot = Join-Path $WorkRoot 'runs'
$rawPath = Join-Path $WorkRoot 'p3c0_plan_inspection_raw.csv'
$provenancePath = Join-Path $WorkRoot 'p3c0_inspector_provenance.csv'
if ((Test-Path -LiteralPath $runsRoot) -or (Test-Path -LiteralPath $rawPath) -or
    (Test-Path -LiteralPath $provenancePath)) {
    throw "Refusing to overwrite existing P3C0 inspection outputs under $WorkRoot"
}
New-Item -ItemType Directory -Path $runsRoot -ErrorAction Stop | Out-Null
$binary = Get-Item -LiteralPath $InspectorExe
[pscustomobject]@{
    inspector_path = $binary.FullName; inspector_bytes = [int64]$binary.Length
    inspector_sha256 = (Get-FileHash -LiteralPath $InspectorExe -Algorithm SHA256).Hash.ToUpperInvariant()
    inspector_timestamp_utc = $binary.LastWriteTimeUtc.ToString('o')
    evidence_root = $EvidenceRoot
} | Export-Csv -LiteralPath $provenancePath -NoTypeInformation

$rows = [Collections.Generic.List[object]]::new()
foreach ($fixtureKey in $FixtureKeys) {
    foreach ($theta in $Thetas) {
        $tag = if ([Math]::Abs($theta - 0.20) -le 1.0e-15) { 'theta020' } else { 'theta015' }
        $runDir = Join-Path (Join-Path $runsRoot $fixtureKey) $tag
        $row = [ordered]@{
            fixture_key = $fixtureKey; theta = $theta; run_dir = $runDir
            exit_code = -1; status = 'failed'; error = ''
            seed_profile = ''; profile_path = ''; inspection_path = ''
            panels = ''; tree_hash = ''; worklist_hash = ''; policy_hash = ''
            direct_pairs = ''; direct_panel_pairs = ''; direct_seconds = ''
            target_batches = ''; full_four_lane_batches = ''; tail_batches = ''
            full_four_lane_task_count = ''; tail_task_count = ''
            target_lanes = ''; task_count = ''; source_span_count = ''
            source_spans_per_task = ''; tasks_per_batch_mean = ''
            tasks_per_batch_min = ''; tasks_per_batch_p50 = ''
            tasks_per_batch_p95 = ''; tasks_per_batch_max = ''
            source_panels_per_task_mean = ''; source_panels_per_task_min = ''
            source_panels_per_task_p50 = ''; source_panels_per_task_p95 = ''
            source_panels_per_task_max = ''; metadata_bytes = ''; total_pairs = ''
        }
        try {
            if (Test-Path -LiteralPath $runDir) { throw "Refusing to reuse run directory: $runDir" }
            New-Item -ItemType Directory -Path $runDir -ErrorAction Stop | Out-Null
            $seed = Assert-FixtureInput $fixtureKey $theta
            $reference = Read-OneRow $seed.seed_profile 'approved reference profile'
            $argPath = New-RunArg $seed $runDir
            $inspectionPath = Join-Path $runDir 'p3c0_packed_plan_inspection.csv'
            $previousInspectorPath = $env:TITAN_P3C0_PACKED_PLAN_INSPECT_PATH
            try {
                $env:TITAN_P3C0_PACKED_PLAN_INSPECT_PATH = $inspectionPath
                Push-Location $runDir
                try {
                    & $InspectorExe -a $argPath 1> stdout.log 2> stderr.log
                    $row.exit_code = $LASTEXITCODE
                } finally {
                    Pop-Location
                }
            } finally {
                if ($null -eq $previousInspectorPath) {
                    Remove-Item Env:TITAN_P3C0_PACKED_PLAN_INSPECT_PATH -ErrorAction SilentlyContinue
                } else {
                    $env:TITAN_P3C0_PACKED_PLAN_INSPECT_PATH = $previousInspectorPath
                }
            }
            if ($row.exit_code -ne 0) { throw "inspector process exit=$($row.exit_code)" }
            $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
            $profile = Read-OneRow $profilePath 'P3C0 profile'
            $inspection = Read-OneRow $inspectionPath 'P3C0 inspector output'
            Assert-ProfileIdentity $reference $profile $inspection $theta
            $row.status = 'ok'
            $row.seed_profile = $seed.seed_profile; $row.profile_path = $profilePath
            $row.inspection_path = $inspectionPath
            foreach ($name in @('panels','bie_cartesian_full_fmm_tree_hash',
                    'bie_cartesian_full_fmm_worklist_hash',
                    'bie_cartesian_full_fmm_policy_hash',
                    'bie_cartesian_full_fmm_direct_pairs',
                    'bie_cartesian_full_fmm_direct_panel_pairs_expected',
                    'bie_cartesian_full_fmm_direct_seconds')) {
                $key = switch ($name) {
                    'bie_cartesian_full_fmm_tree_hash' { 'tree_hash'; break }
                    'bie_cartesian_full_fmm_worklist_hash' { 'worklist_hash'; break }
                    'bie_cartesian_full_fmm_policy_hash' { 'policy_hash'; break }
                    'bie_cartesian_full_fmm_direct_pairs' { 'direct_pairs'; break }
                    'bie_cartesian_full_fmm_direct_panel_pairs_expected' { 'direct_panel_pairs'; break }
                    'bie_cartesian_full_fmm_direct_seconds' { 'direct_seconds'; break }
                    default { 'panels' }
                }
                $row[$key] = Get-TextValue $profile $name
            }
            foreach ($name in @('target_batches','full_four_lane_batches','tail_batches',
                    'full_four_lane_task_count','tail_task_count',
                    'target_lanes','task_count','source_span_count',
                    'source_spans_per_task','tasks_per_batch_mean',
                    'tasks_per_batch_min','tasks_per_batch_p50','tasks_per_batch_p95',
                    'tasks_per_batch_max','source_panels_per_task_mean',
                    'source_panels_per_task_min','source_panels_per_task_p50',
                    'source_panels_per_task_p95','source_panels_per_task_max',
                    'metadata_bytes','total_pairs')) {
                $row[$name] = Get-TextValue $inspection $name
            }
        } catch {
            $row.error = $_.Exception.Message
        }
        $rows.Add([pscustomobject]$row)
    }
}
$rows | Export-Csv -LiteralPath $rawPath -NoTypeInformation
if (@($rows | Where-Object { $_.status -ne 'ok' }).Count -ne 0) {
    throw "P3C0 inspector matrix has failed rows; see $rawPath"
}
Write-Host "P3C0 target-tile inspector PASS: $rawPath"
