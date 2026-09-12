[CmdletBinding()]
param(
    [string[]]$Fixtures = @('1sh1', '1a63', 'barnase'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5d1_prepared_apply',
    [string]$Compiler = 'clang++',
    [int]$CaseTimeoutSeconds = 1800,
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
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5d1_prepared_apply_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5d1_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5d1_prepared_apply_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5d1_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5d1_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_CALDERON_C5D1_PREPARED_APPLY_REPORT.txt'

function Get-TokenValue([string]$Text, [string]$Name) {
    $match = [regex]::Match($Text, '(?:^|\s)' + [regex]::Escape($Name) + '=([^\s]+)')
    if ($match.Success) { return $match.Groups[1].Value }
    return ''
}

function Resolve-C5D1Fixture([string]$Name, [string]$RunDir) {
    switch ($Name) {
        '1sh1' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr' }
        '1a63' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr' }
        'barnase' {
            $public = Join-Path $RepoRoot 'test\pbe_inputs_public\barnase_full_amber_ph7.pqr'
            if (Test-Path -LiteralPath $public) { return $public }
            return Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $RunDir
        }
        default { throw "Unsupported C5D.1 fixture '$Name'." }
    }
}

function Invoke-C5D1Case([string]$Executable, [string]$PqrPath,
                          [string]$RunDir, [int]$TimeoutSeconds) {
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $Executable
    $startInfo.Arguments = if ([string]::IsNullOrWhiteSpace($PqrPath)) { '' } else { '"' + $PqrPath + '"' }
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
    Set-Content -LiteralPath (Join-Path $RunDir 'stdout.txt') -Value $stdout -Encoding ASCII
    Set-Content -LiteralPath (Join-Path $RunDir 'stderr.txt') -Value $stderr -Encoding ASCII
    return [pscustomobject]@{
        timed_out = $timedOut
        exit_code = if ($timedOut) { -1 } else { $process.ExitCode }
        wall_seconds = $timer.Elapsed.TotalSeconds
        diagnostic = (($stdout + ' ' + $stderr) -replace '[\r\n]+', ' ').Trim()
    }
}

$sources = @(
    'test\titan_calderon_c5d1_prepared_apply_smoke.cpp',
    'src\tabi\tabi_calderon_trace_space.cpp',
    'src\tabi\tabi_calderon_aux_surface.cpp',
    'src\tabi\tabi_analytic_ses.cpp',
    'src\tabi\tabi_ses_quartic.cpp',
    'src\tabi\tabi_ses_patch_intersections.cpp',
    'src\tabi\tabi_ses_patch_provenance.cpp',
    'src\tabi\tabi_ses_patch_ray.cpp',
    'src\tabi\tabi_ses_bilevel_grid.cpp',
    'src\tabi\tabi_calderon_transfer.cpp',
    'src\tabi\tabi_calderon_weak_form.cpp',
    'src\tabi\tabi_calderon_prepared.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\titan_surface.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\titan_executor.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp'
)
$sourcePaths = @($sources | ForEach-Object { Join-Path $RepoRoot $_ })
$needsBuild = $Force -or -not (Test-Path -LiteralPath $ExecutablePath)
if (-not $needsBuild) {
    $exeTime = (Get-Item -LiteralPath $ExecutablePath).LastWriteTimeUtc
    $needsBuild = @($sourcePaths | Where-Object {
        (Get-Item -LiteralPath $_).LastWriteTimeUtc -gt $exeTime
    }).Count -gt 0
}
if ($needsBuild) {
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5D.1 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++23','-O2','-mavx2','-mfma',
        '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',"-I$RepoRoot\src") +
        $sourcePaths + @('-o',$ExecutablePath)
    Write-Host "build: C5D.1 prepared Calderon apply with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($LASTEXITCODE -ne 0) { throw "C5D.1 build failed; see $BuildLogPath" }
}

$caseNames = @('synthetic') + $Fixtures
$existing = @{}
if ((Test-Path -LiteralPath $MatrixPath) -and -not $Force) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $existing[[string]$row.fixture] = $row
    }
}
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in $caseNames) {
    $prior = $existing[$fixture]
    if ($null -ne $prior -and -not $Force -and
        (-not $RerunFailed -or [string]$prior.run_status -eq 'ok')) {
        Write-Host "resume: $fixture"
        $rows.Add($prior)
        continue
    }
    $runDir = Join-Path $WorkRoot $fixture
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $pqr = if ($fixture -eq 'synthetic') { '' } else {
        [System.IO.Path]::GetFullPath((Resolve-C5D1Fixture -Name $fixture -RunDir $runDir))
    }
    Write-Host "run: $fixture C5D.1 prepared apply (timeout ${CaseTimeoutSeconds}s)"
    $invocation = Invoke-C5D1Case -Executable $ExecutablePath -PqrPath $pqr `
        -RunDir $runDir -TimeoutSeconds $CaseTimeoutSeconds
    $diagnostic = $invocation.diagnostic
    $runStatus = if ($invocation.timed_out) { 'timeout' } `
        elseif ($invocation.exit_code -ne 0) { 'failed' } else { 'ok' }
    $rows.Add([pscustomobject]@{
        fixture = $fixture
        run_status = $runStatus
        exit_code = $invocation.exit_code
        wall_seconds = $invocation.wall_seconds
        status = Get-TokenValue $diagnostic 'status'
        geometry_status = Get-TokenValue $diagnostic 'geometry_status'
        input_atoms = Get-TokenValue $diagnostic 'input_atoms'
        primary_panels = Get-TokenValue $diagnostic 'primary_panels'
        auxiliary_vertices = Get-TokenValue $diagnostic 'auxiliary_vertices'
        fmm_quadrature_panels = Get-TokenValue $diagnostic 'fmm_quadrature_panels'
        prepare_build_count = Get-TokenValue $diagnostic 'prepare_build_count'
        context_build_count = Get-TokenValue $diagnostic 'context_build_count'
        apply_count = Get-TokenValue $diagnostic 'apply_count'
        prepare_seconds = Get-TokenValue $diagnostic 'prepare_seconds'
        transfer_seconds = Get-TokenValue $diagnostic 'transfer_seconds'
        fmm_setup_seconds = Get-TokenValue $diagnostic 'fmm_setup_seconds'
        first_apply_seconds = Get-TokenValue $diagnostic 'first_apply_seconds'
        second_apply_seconds = Get-TokenValue $diagnostic 'second_apply_seconds'
        pairing_seconds = Get-TokenValue $diagnostic 'pairing_seconds'
        pairing_iterations = Get-TokenValue $diagnostic 'pairing_iterations'
        near_pairs = Get-TokenValue $diagnostic 'near_pairs'
        far_pairs = Get-TokenValue $diagnostic 'far_pairs'
        estimated_bytes = Get-TokenValue $diagnostic 'estimated_bytes'
        dense_status = Get-TokenValue $diagnostic 'dense_status'
        dense_relative = Get-TokenValue $diagnostic 'dense_relative'
        dense_phi_relative = Get-TokenValue $diagnostic 'dense_phi_relative'
        dense_phin_relative = Get-TokenValue $diagnostic 'dense_phin_relative'
        diagnostic = $diagnostic
    })
}

$merged = @{}
foreach ($row in $existing.Values) { $merged[[string]$row.fixture] = $row }
foreach ($row in $rows) { $merged[[string]$row.fixture] = $row }
$finalRows = @($merged.Values | Sort-Object fixture)
$finalRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$selected = @($finalRows | Where-Object { $caseNames -contains [string]$_.fixture })
$synthetic = @($selected | Where-Object { $_.fixture -eq 'synthetic' } | Select-Object -First 1)
$runtimeGap = @($selected | Where-Object { $_.run_status -ne 'ok' }).Count -gt 0
$lifecycleGap = @($selected | Where-Object {
    $_.status -ne 'prepared' -or $_.geometry_status -ne 'pass' -or
    [int]$_.prepare_build_count -ne 1 -or [int]$_.context_build_count -ne 1 -or
    [int]$_.apply_count -lt 2 -or [int64]$_.far_pairs -le 0
}).Count -gt 0
$parityGap = $synthetic.Count -ne 1 -or $synthetic[0].dense_status -ne 'pass' -or
    [double]$synthetic[0].dense_relative -gt 5e-2 -or
    [double]$synthetic[0].dense_phi_relative -gt 1e-1 -or
    [double]$synthetic[0].dense_phin_relative -gt 1e-1
$decision = if ($runtimeGap) {
    'calderon_c5d1_runtime_limited'
} elseif ($lifecycleGap) {
    'calderon_c5d1_lifecycle_gap'
} elseif ($parityGap) {
    'calderon_c5d1_dense_parity_gap'
} else {
    'calderon_c5d1_prepared_apply_pass'
}
[pscustomobject]@{
    decision = $decision
    rows = $selected.Count
    quick = [bool]$Quick
    dense_relative_gate = 0.05
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision","rows=$($selected.Count)","quick=$([bool]$Quick)") |
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5D.1 FMM-Native Prepared Apply Report',
    "Decision: $decision",
    'Scope: diagnostic prepared auxiliary action only. GMRES, requested/effective preconditioner selection, fallback policy, and production defaults are unchanged.',
    'Method: six barycentric dual patches per auxiliary triangle provide mixed P1/dual-P0 source and target quadrature. A frozen analytical FMM context evaluates the JCP action; a prepared sparse Riesz pairing maps the weak action back to strong coordinates.',
    'Dense gate: the synthetic sphere candidate must agree with the exact C5C.4 dense mixed action within 5% aggregate relative L2 and 10% per coupled row.',
    'Rows:',
    (($selected | ForEach-Object {
        "  $($_.fixture): run=$($_.run_status) status=$($_.status) geometry=$($_.geometry_status) primary=$($_.primary_panels) auxiliary=$($_.auxiliary_vertices) quadrature=$($_.fmm_quadrature_panels) builds=$($_.prepare_build_count)/$($_.context_build_count) applies=$($_.apply_count) near=$($_.near_pairs) far=$($_.far_pairs) prepare_s=$($_.prepare_seconds) apply_s=$($_.second_apply_seconds) pairing_iter=$($_.pairing_iterations) dense=$($_.dense_relative) phi=$($_.dense_phi_relative) phin=$($_.dense_phin_relative)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5D.1 prepared-apply decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5d1_prepared_apply_pass' -and
    -not $AllowFailedDecision) {
    throw "C5D.1 prepared apply gate did not pass: $decision"
}
