[CmdletBinding()]
param(
    [string[]]$Fixtures = @('1sh1', '1a63', 'barnase'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c6_real_primary_transfer',
    [string]$Compiler = 'clang++',
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
if ($Quick) {
    $Fixtures = @('1sh1')
}
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c6_real_primary_transfer_smoke.exe'
$MassSmokePath = Join-Path $WorkRoot 'titan_calderon_c5c6_mass_compatible_transfer_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c6_real_primary_transfer_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c6_real_primary_transfer_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c6_real_primary_transfer_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c6_real_primary_transfer_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C6_REAL_PRIMARY_TRANSFER_REPORT.txt'

function Get-TokenValue([string]$Text, [string]$Name) {
    $pattern = '(?:^|\s)' + [regex]::Escape($Name) + '=([^\s]+)'
    $match = [regex]::Match($Text, $pattern)
    if ($match.Success) { return $match.Groups[1].Value }
    return ''
}

function Resolve-C5C6Fixture([string]$Name, [string]$RunDir) {
    switch ($Name) {
        '1sh1' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr' }
        '1a63' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr' }
        'barnase' { return Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $RunDir }
        default { throw "Unsupported C5C.6 fixture '$Name'. Use 1sh1, 1a63, or barnase." }
    }
}

function Invoke-C5C6Case([string]$Executable, [string]$PqrPath, [int]$TimeoutSeconds) {
    $stdoutPath = Join-Path (Split-Path -Parent $PqrPath) 'c5c6_stdout.tmp'
    $stderrPath = Join-Path (Split-Path -Parent $PqrPath) 'c5c6_stderr.tmp'
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $Executable
    $startInfo.Arguments = '"' + $PqrPath + '"'
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
    $exitCode = if ($timedOut) { -1 } else { $process.ExitCode }
    return [pscustomobject]@{
        timed_out = $timedOut
        exit_code = $exitCode
        wall_seconds = $timer.Elapsed.TotalSeconds
        diagnostic = (($stdout + ' ' + $stderr) -replace '[\r\n]+', ' ').Trim()
    }
}

if ($Force -or -not (Test-Path -LiteralPath $ExecutablePath) -or
    -not (Test-Path -LiteralPath $MassSmokePath)) {
    $sources = @(
        (Join-Path $RepoRoot 'test\titan_calderon_c5c6_real_primary_transfer_smoke.cpp'),
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
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5C.6 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.6 real-primary transfer diagnostic with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.6 build failed; see $BuildLogPath"
    }
    $massSources = @(
        (Join-Path $RepoRoot 'test\titan_calderon_c5c6_mass_compatible_transfer_smoke.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_bilevel_grid.cpp'),
        (Join-Path $RepoRoot 'src\titan_executor.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_transfer.cpp')
    )
    $massArguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $massSources + @('-o', $MassSmokePath)
    $massBuildOutput = & $Compiler @massArguments 2>&1
    $massBuildExitCode = $LASTEXITCODE
    @($massBuildOutput | ForEach-Object { [string]$_ }) |
        Add-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($massBuildExitCode -ne 0) {
        throw "C5C.6 mass-compatible smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.6 matched/shifted mass-compatible transfer smoke'
$massSmokeOutput = & $MassSmokePath 2>&1
if ($LASTEXITCODE -ne 0 -or (@($massSmokeOutput) -join ' ') -notmatch
    'calderon_c5c6_mass_compatible_transfer_smoke=pass') {
    throw "C5C.6 matched/shifted mass-compatible transfer smoke failed: $(@($massSmokeOutput) -join ' ')"
}
Write-Host 'run: C5C.6 production-primary synthetic smoke'
$primarySmokeOutput = & $ExecutablePath 2>&1
if ($LASTEXITCODE -ne 0 -or (@($primarySmokeOutput) -join ' ') -notmatch
    'calderon_c5c6_real_primary_transfer_smoke=pass') {
    throw "C5C.6 production-primary synthetic smoke failed: $(@($primarySmokeOutput) -join ' ')"
}

$existing = @{}
if ((Test-Path -LiteralPath $MatrixPath) -and -not $Force) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $existing[[string]$row.fixture] = $row
    }
}
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in $Fixtures) {
    $prior = $existing[$fixture]
    if ($null -ne $prior -and -not $Force -and
        (-not $RerunFailed -or [string]$prior.run_status -eq 'ok')) {
        Write-Host "resume: $fixture"
        $rows.Add($prior)
        continue
    }
    $runDir = Join-Path $WorkRoot $fixture
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $pqr = Resolve-C5C6Fixture -Name $fixture -RunDir $runDir
    if (-not (Test-Path -LiteralPath $pqr)) {
        throw "Missing C5C.6 fixture: $pqr"
    }
    Write-Host "run: $fixture C5C.6 transfer (timeout ${CaseTimeoutSeconds}s)"
    $invocation = Invoke-C5C6Case -Executable $ExecutablePath `
        -PqrPath ([System.IO.Path]::GetFullPath($pqr)) -TimeoutSeconds $CaseTimeoutSeconds
    $diagnostic = $invocation.diagnostic
    $runStatus = if ($invocation.timed_out) { 'timeout' } elseif ($invocation.exit_code -ne 0) { 'failed' } else { 'ok' }
    $rows.Add([pscustomobject]@{
        fixture = $fixture
        run_status = $runStatus
        exit_code = $invocation.exit_code
        timed_out = $invocation.timed_out
        wall_seconds = $invocation.wall_seconds
        transfer_status = Get-TokenValue -Text $diagnostic -Name 'status'
        geometry_status = Get-TokenValue -Text $diagnostic -Name 'geometry_status'
        input_atoms = Get-TokenValue -Text $diagnostic -Name 'input_atoms'
        primary_panels = Get-TokenValue -Text $diagnostic -Name 'primary_panels'
        auxiliary_vertices = Get-TokenValue -Text $diagnostic -Name 'auxiliary_vertices'
        support_edges = Get-TokenValue -Text $diagnostic -Name 'support_edges'
        support_k = Get-TokenValue -Text $diagnostic -Name 'support_k'
        support_components = Get-TokenValue -Text $diagnostic -Name 'support_components'
        repair_edges = Get-TokenValue -Text $diagnostic -Name 'repair_edges'
        balance_iterations = Get-TokenValue -Text $diagnostic -Name 'balance_iterations'
        balance_relative = Get-TokenValue -Text $diagnostic -Name 'balance_relative'
        primary_mass_scale = Get-TokenValue -Text $diagnostic -Name 'mass_scale'
        component_scale_min = Get-TokenValue -Text $diagnostic -Name 'component_scale_min'
        component_scale_max = Get-TokenValue -Text $diagnostic -Name 'component_scale_max'
        forward_constant_relative = Get-TokenValue -Text $diagnostic -Name 'forward_constant'
        reverse_constant_relative = Get-TokenValue -Text $diagnostic -Name 'reverse_constant'
        adjoint_relative = Get-TokenValue -Text $diagnostic -Name 'adjoint'
        relative_area_mismatch = Get-TokenValue -Text $diagnostic -Name 'area_rel'
        primary_to_auxiliary_rms_distance = Get-TokenValue -Text $diagnostic -Name 'p2a_rms'
        primary_to_auxiliary_max_distance = Get-TokenValue -Text $diagnostic -Name 'p2a_max'
        auxiliary_to_primary_rms_distance = Get-TokenValue -Text $diagnostic -Name 'a2p_rms'
        auxiliary_to_primary_max_distance = Get-TokenValue -Text $diagnostic -Name 'a2p_max'
        primary_coverage_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_coverage'
        auxiliary_coverage_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_coverage'
        primary_normal_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'p2a_normal_fraction'
        auxiliary_normal_alignment_fraction = Get-TokenValue -Text $diagnostic -Name 'a2p_normal_fraction'
        transfer_build_ms = Get-TokenValue -Text $diagnostic -Name 'build_ms'
        diagnostic = $diagnostic
    })
}

$merged = @{}
foreach ($row in $existing.Values) { $merged[[string]$row.fixture] = $row }
foreach ($row in $rows) { $merged[[string]$row.fixture] = $row }
$finalRows = @($merged.Values | Sort-Object fixture)
foreach ($row in $finalRows) {
    if ($row.PSObject.Properties.Name -notcontains 'support_k') {
        $inferredSupportK = if ([double]$row.auxiliary_vertices -gt 0) {
            [int][math]::Round([double]$row.support_edges / [double]$row.auxiliary_vertices)
        } else { 0 }
        $row | Add-Member -NotePropertyName support_k -NotePropertyValue $inferredSupportK
    }
}
$finalRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$decisionRows = @($finalRows | Where-Object { $Fixtures -contains [string]$_.fixture })

$runtimeGap = @($decisionRows | Where-Object { $_.run_status -ne 'ok' }).Count -gt 0
$massGap = @($decisionRows | Where-Object {
    [double]$_.balance_relative -gt 1e-9 -or
    [double]$_.forward_constant_relative -gt 1e-8 -or
    [double]$_.reverse_constant_relative -gt 1e-8 -or
    [double]$_.adjoint_relative -gt 1e-8
}).Count -gt 0
$geometryGap = @($decisionRows | Where-Object {
    $_.geometry_status -ne 'pass' -or [double]$_.primary_coverage_fraction -lt 1.0 -or
    [double]$_.auxiliary_coverage_fraction -lt 1.0 -or
    [double]$_.primary_normal_alignment_fraction -lt 0.90 -or
    [double]$_.auxiliary_normal_alignment_fraction -lt 0.90
}).Count -gt 0
$decision = if ($runtimeGap) {
    'calderon_c5c6_runtime_limited'
} elseif ($massGap) {
    'calderon_c5c6_mass_transfer_gap'
} elseif ($geometryGap) {
    'calderon_c5c6_real_primary_geometry_gap'
} else {
    'calderon_c5c6_mass_compatible_real_primary_pass'
}
[pscustomobject]@{
    decision = $decision
    rows = $decisionRows.Count
    quick = [bool]$Quick
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision", "rows=$($decisionRows.Count)", "quick=$([bool]$Quick)") |
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.6 Mass-Compatible Real-Primary Transfer Report',
    "Decision: $decision",
    'Scope: diagnostic transfer between production TABI panel surfaces and watertight implicit-VDW auxiliary P1 spaces. No solver, preconditioner, GMRES, or fallback policy is changed.',
    'Method: sparse local Sinkhorn transport with auxiliary P1 lumped-mass row marginals and component-scaled primary panel-area column marginals. The paired maps preserve constants in both directions and are exact Riesz adjoints.',
    'Geometry gate: distances, area mismatch, full sampled coverage, and the fraction of nonnegative normal alignments are reported independently of mass balancing.',
    'Rows:',
    (($finalRows | ForEach-Object {
        "  $($_.fixture): run=$($_.run_status) transfer=$($_.transfer_status) geometry=$($_.geometry_status) atoms=$($_.input_atoms) primary=$($_.primary_panels) auxiliary=$($_.auxiliary_vertices) support_k=$($_.support_k) balance=$($_.balance_relative) forward_const=$($_.forward_constant_relative) reverse_const=$($_.reverse_constant_relative) adjoint=$($_.adjoint_relative) area_rel=$($_.relative_area_mismatch) p2a_max=$($_.primary_to_auxiliary_max_distance) a2p_max=$($_.auxiliary_to_primary_max_distance) p2a_normal_fraction=$($_.primary_normal_alignment_fraction) a2p_normal_fraction=$($_.auxiliary_normal_alignment_fraction) ms=$($_.transfer_build_ms)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.6 real-primary-transfer decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c6_mass_compatible_real_primary_pass' -and -not $AllowFailedDecision) {
    throw "C5C.6 real-primary transfer gate did not pass: $decision"
}
