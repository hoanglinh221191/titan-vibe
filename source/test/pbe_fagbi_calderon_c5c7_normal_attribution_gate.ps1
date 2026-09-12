[CmdletBinding()]
param(
    [string[]]$Fixtures = @('1sh1', '1a63', 'barnase'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c7_normal_attribution',
    [string]$Compiler = 'clang++',
    [double]$OwnerDistanceSlack = 0.25,
    [int]$CaseTimeoutSeconds = 900,
    [switch]$Quick,
    [switch]$RerunFailed,
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$Fixtures = @($Fixtures | ForEach-Object {
    @([string]$_ -split ',') | ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})
if ($Quick) { $Fixtures = @('1sh1') }
if ($OwnerDistanceSlack -lt 0.0 -or $OwnerDistanceSlack -gt 1.0) {
    throw 'OwnerDistanceSlack must be in 0..1 Angstrom.'
}

$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c7_normal_attribution_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c7_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c7_normal_attribution_matrix.csv'
$MismatchPath = Join-Path $WorkRoot 'calderon_c5c7_normal_mismatches.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c7_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c7_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_CALDERON_C5C7_NORMAL_ATTRIBUTION_REPORT.txt'

function Get-TokenValue([string]$Text, [string]$Name) {
    $pattern = '(?:^|\s)' + [regex]::Escape($Name) + '=([^\s]+)'
    $match = [regex]::Match($Text, $pattern)
    if ($match.Success) { return $match.Groups[1].Value }
    return ''
}

function Resolve-C5C7Fixture([string]$Name, [string]$RunDir) {
    switch ($Name) {
        '1sh1' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr' }
        '1a63' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr' }
        'barnase' {
            $public = Join-Path $RepoRoot 'test\pbe_inputs_public\barnase_full_amber_ph7.pqr'
            if (Test-Path -LiteralPath $public) { return $public }
            return Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $RunDir
        }
        default { throw "Unsupported C5C.7 fixture '$Name'. Use 1sh1, 1a63, or barnase." }
    }
}

function Invoke-C5C7Case(
    [string]$Executable,
    [string]$PqrPath,
    [string]$Policy,
    [double]$Slack,
    [string]$CaseMismatchPath,
    [int]$TimeoutSeconds
) {
    $runDir = Split-Path -Parent $CaseMismatchPath
    $stdoutPath = Join-Path $runDir 'stdout.txt'
    $stderrPath = Join-Path $runDir 'stderr.txt'
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $Executable
    $startInfo.Arguments = '"' + $PqrPath + '" ' + $Policy + ' ' +
        $Slack.ToString('R', [System.Globalization.CultureInfo]::InvariantCulture) +
        ' "' + $CaseMismatchPath + '"'
    $startInfo.WorkingDirectory = $RepoRoot
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $startInfo
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $null = $process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $timedOut = -not $process.WaitForExit($TimeoutSeconds * 1000)
    if ($timedOut) {
        try { $process.Kill() } catch {}
        $process.WaitForExit()
    }
    $timer.Stop()
    $stdout = $stdoutTask.Result
    $stderr = $stderrTask.Result
    Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding ASCII
    Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding ASCII
    return [pscustomobject]@{
        timed_out = $timedOut
        exit_code = if ($timedOut) { -1 } else { $process.ExitCode }
        wall_seconds = $timer.Elapsed.TotalSeconds
        diagnostic = (($stdout + ' ' + $stderr) -replace '[\r\n]+', ' ').Trim()
    }
}

$sources = @(
    (Join-Path $RepoRoot 'test\titan_calderon_c5c7_normal_attribution_smoke.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_bilevel_grid.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_transfer.cpp'),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp')
)
$needsBuild = $Force -or -not (Test-Path -LiteralPath $ExecutablePath)
if (-not $needsBuild) {
    $exeTime = (Get-Item -LiteralPath $ExecutablePath).LastWriteTimeUtc
    $needsBuild = @($sources | Where-Object {
        (Get-Item -LiteralPath $_).LastWriteTimeUtc -gt $exeTime
    }).Count -gt 0
}
if ($needsBuild) {
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5C.7 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources +
        @('-o', $ExecutablePath)
    Write-Host "build: C5C.7 normal attribution with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($LASTEXITCODE -ne 0) { throw "C5C.7 build failed; see $BuildLogPath" }
}

Write-Host 'run: C5C.7 synthetic deterministic comparison smoke'
$smokeOutput = & $ExecutablePath 2>&1
if ($LASTEXITCODE -ne 0 -or (@($smokeOutput) -join ' ') -notmatch
    'calderon_c5c7_normal_attribution_smoke=pass') {
    throw "C5C.7 synthetic smoke failed: $(@($smokeOutput) -join ' ')"
}

$existing = @{}
if ((Test-Path -LiteralPath $MatrixPath) -and -not $Force) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $existing["$($row.fixture)|$($row.policy)"] = $row
    }
}
$rows = New-Object System.Collections.Generic.List[object]
$policies = @('facet-nearest', 'owner-sdf-aware')
foreach ($fixture in $Fixtures) {
    foreach ($policy in $policies) {
        $key = "$fixture|$policy"
        $prior = $existing[$key]
        if ($null -ne $prior -and -not $Force -and
            (-not $RerunFailed -or [string]$prior.run_status -eq 'ok')) {
            Write-Host "resume: $key"
            $rows.Add($prior)
            continue
        }
        $runDir = Join-Path $WorkRoot ($fixture + '_' + ($policy -replace '-', '_'))
        New-Item -ItemType Directory -Force -Path $runDir | Out-Null
        $pqr = Resolve-C5C7Fixture -Name $fixture -RunDir $runDir
        if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing C5C.7 fixture: $pqr" }
        $caseMismatch = Join-Path $runDir 'normal_mismatches.csv'
        Write-Host "run: $key (timeout ${CaseTimeoutSeconds}s)"
        $invocation = Invoke-C5C7Case -Executable $ExecutablePath `
            -PqrPath ([System.IO.Path]::GetFullPath($pqr)) -Policy $policy `
            -Slack $OwnerDistanceSlack -CaseMismatchPath $caseMismatch `
            -TimeoutSeconds $CaseTimeoutSeconds
        $diagnostic = $invocation.diagnostic
        $runStatus = if ($invocation.timed_out) { 'timeout' } `
            elseif ($invocation.exit_code -ne 0) { 'failed' } else { 'ok' }
        $rows.Add([pscustomobject]@{
            fixture = $fixture
            policy = $policy
            slack = $OwnerDistanceSlack
            run_status = $runStatus
            exit_code = $invocation.exit_code
            wall_seconds = $invocation.wall_seconds
            transfer_status = Get-TokenValue -Text $diagnostic -Name 'status'
            input_atoms = Get-TokenValue -Text $diagnostic -Name 'input_atoms'
            primary_panels = Get-TokenValue -Text $diagnostic -Name 'primary_panels'
            auxiliary_vertices = Get-TokenValue -Text $diagnostic -Name 'auxiliary_vertices'
            relative_area_mismatch = Get-TokenValue -Text $diagnostic -Name 'area_rel'
            primary_to_auxiliary_rms_distance = Get-TokenValue -Text $diagnostic -Name 'p2a_rms'
            primary_to_auxiliary_max_distance = Get-TokenValue -Text $diagnostic -Name 'p2a_max'
            auxiliary_to_primary_rms_distance = Get-TokenValue -Text $diagnostic -Name 'a2p_rms'
            auxiliary_to_primary_max_distance = Get-TokenValue -Text $diagnostic -Name 'a2p_max'
            primary_coverage_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_coverage'
            auxiliary_coverage_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_coverage'
            primary_normal_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_normal_fraction'
            auxiliary_normal_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_normal_fraction'
            primary_weighted_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_weighted_fraction'
            auxiliary_weighted_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_weighted_fraction'
            primary_same_owner_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_same_owner'
            auxiliary_same_owner_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_same_owner'
            primary_owner_ambiguity_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_ambiguity'
            auxiliary_owner_ambiguity_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_ambiguity'
            primary_negative_samples = Get-TokenValue -Text $diagnostic -Name 'p2a_negative'
            auxiliary_negative_samples = Get-TokenValue -Text $diagnostic -Name 'a2p_negative'
            primary_positive_alternatives = Get-TokenValue -Text $diagnostic -Name 'p2a_positive_alternative'
            auxiliary_positive_alternatives = Get-TokenValue -Text $diagnostic -Name 'a2p_positive_alternative'
            forward_constant_relative = Get-TokenValue -Text $diagnostic -Name 'forward_constant'
            adjoint_relative = Get-TokenValue -Text $diagnostic -Name 'adjoint'
            transfer_build_ms = Get-TokenValue -Text $diagnostic -Name 'build_ms'
            mismatch_csv = $caseMismatch
            diagnostic = $diagnostic
        })
    }
}

$merged = @{}
foreach ($row in $existing.Values) { $merged["$($row.fixture)|$($row.policy)"] = $row }
foreach ($row in $rows) { $merged["$($row.fixture)|$($row.policy)"] = $row }
$finalRows = @($merged.Values | Sort-Object fixture, policy)
$finalRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$mismatches = New-Object System.Collections.Generic.List[object]
foreach ($row in $finalRows) {
    if (-not [string]::IsNullOrWhiteSpace([string]$row.mismatch_csv) -and
        (Test-Path -LiteralPath $row.mismatch_csv)) {
        foreach ($sample in @(Import-Csv -LiteralPath $row.mismatch_csv)) {
            $mismatches.Add([pscustomobject]@{
                fixture = $row.fixture
                policy = $row.policy
                direction = $sample.direction
                source_index = $sample.source_index
                source_atom_index = $sample.source_atom_index
                matched_index = $sample.matched_index
                matched_atom_index = $sample.matched_atom_index
                distance = $sample.distance
                normal_dot = $sample.normal_dot
                same_owner = $sample.same_owner
                owner_ambiguous = $sample.owner_ambiguous
                positive_alternative = $sample.positive_alternative
            })
        }
    }
}
$mismatches | Export-Csv -LiteralPath $MismatchPath -NoTypeInformation

$selected = @($finalRows | Where-Object { $Fixtures -contains [string]$_.fixture })
$corrected = @($selected | Where-Object { $_.policy -eq 'owner-sdf-aware' })
$runtimeGap = @($selected | Where-Object { $_.run_status -ne 'ok' }).Count -gt 0
$attributionGap = $corrected.Count -ne $Fixtures.Count -or @($corrected | Where-Object {
    [string]::IsNullOrWhiteSpace([string]$_.primary_weighted_alignment_fraction) -or
    [string]::IsNullOrWhiteSpace([string]$_.auxiliary_weighted_alignment_fraction)
}).Count -gt 0
$geometryGap = @($corrected | Where-Object {
    $_.transfer_status -ne 'pass' -or
    [double]$_.primary_coverage_fraction -lt 1.0 -or
    [double]$_.auxiliary_coverage_fraction -lt 1.0 -or
    [double]$_.primary_normal_alignment_fraction -lt 0.90 -or
    [double]$_.auxiliary_normal_alignment_fraction -lt 0.90 -or
    [double]$_.primary_weighted_alignment_fraction -lt 0.90 -or
    [double]$_.auxiliary_weighted_alignment_fraction -lt 0.90
}).Count -gt 0
$distanceGap = $false
foreach ($row in $corrected) {
    $legacy = @($selected | Where-Object {
        $_.fixture -eq $row.fixture -and $_.policy -eq 'facet-nearest'
    } | Select-Object -First 1)
    if ($legacy.Count -eq 1 -and
        ([double]$row.primary_to_auxiliary_max_distance -gt
             [double]$legacy[0].primary_to_auxiliary_max_distance + $OwnerDistanceSlack + 1e-8 -or
         [double]$row.auxiliary_to_primary_max_distance -gt
             [double]$legacy[0].auxiliary_to_primary_max_distance + $OwnerDistanceSlack + 1e-8)) {
        $distanceGap = $true
    }
}
$decision = if ($runtimeGap) {
    'calderon_c5c7_runtime_limited'
} elseif ($attributionGap) {
    'calderon_c5c7_attribution_gap'
} elseif ($geometryGap) {
    'calderon_c5c7_representation_gap'
} elseif ($distanceGap) {
    'calderon_c5c7_distance_tradeoff_gap'
} else {
    'calderon_c5c7_owner_sdf_geometry_pass'
}

[pscustomobject]@{
    decision = $decision
    rows = $selected.Count
    fixtures = $Fixtures.Count
    owner_distance_slack = $OwnerDistanceSlack
    quick = [bool]$Quick
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    "rows=$($selected.Count)",
    "fixtures=$($Fixtures.Count)",
    "owner_distance_slack=$OwnerDistanceSlack",
    "quick=$([bool]$Quick)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.7 Barnase Seam And Normal Attribution Report',
    "Decision: $decision",
    'Scope: diagnostic primary/auxiliary geometry matching only. Sparse transfer stencils, Sinkhorn balancing, solver, GMRES, FMM, and production defaults are unchanged.',
    "Policy: owner/SDF-aware matching considers only candidates within $OwnerDistanceSlack Angstrom of the Euclidean nearest distance, uses interpolated union-SDF vertex normals, and deterministically prefers aligned same-owner correspondences.",
    'Gate: both sample-count and area/P1-mass-weighted normal alignment must remain at least 0.90; distance coverage must remain complete.',
    'Rows:',
    (($selected | ForEach-Object {
        "  $($_.fixture) $($_.policy): run=$($_.run_status) status=$($_.transfer_status) p2a=$($_.primary_normal_alignment_fraction) a2p=$($_.auxiliary_normal_alignment_fraction) p2a_weighted=$($_.primary_weighted_alignment_fraction) a2p_weighted=$($_.auxiliary_weighted_alignment_fraction) p2a_owner=$($_.primary_same_owner_fraction) a2p_owner=$($_.auxiliary_same_owner_fraction) p2a_ambiguity=$($_.primary_owner_ambiguity_fraction) a2p_ambiguity=$($_.auxiliary_owner_ambiguity_fraction) p2a_negative=$($_.primary_negative_samples) a2p_negative=$($_.auxiliary_negative_samples) p2a_max=$($_.primary_to_auxiliary_max_distance) a2p_max=$($_.auxiliary_to_primary_max_distance)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.7 normal-attribution decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Mismatches: $MismatchPath"
if ($decision -ne 'calderon_c5c7_owner_sdf_geometry_pass' -and
    -not $AllowFailedDecision) {
    throw "C5C.7 normal attribution gate did not pass: $decision"
}
