[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string[]]$Fixtures = @('1sh1', '1a63', 'barnase'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5_trace',
    [int]$CaseTimeoutSeconds = 900,
    [ValidateSet('diagonal','fmm-nearfield-ilu0')]
    [string]$Preconditioner = 'diagonal',
    [ValidateRange(1,64)]
    [int]$MaxIter = 8,
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
    @([string]$_ -split ',') | ForEach-Object { $_.Trim() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})

function Resolve-C5Fixture([string]$Name, [string]$RunDir) {
    switch ($Name.ToLowerInvariant()) {
        '1sh1' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr' }
        '1a63' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr' }
        'barnase' { return Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $RunDir }
        default { throw "Unsupported C5 fixture '$Name'. Use 1sh1, 1a63, or barnase." }
    }
}

if ($Quick) {
    $Fixtures = @('1sh1')
}
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Titan executable not found: $TitanExe"
}

$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$MatrixPath = Join-Path $WorkRoot 'calderon_c5_trace_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5_trace_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5_trace_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5_TRACE_SPACE_REPORT.txt'

$existing = @{}
if ((Test-Path -LiteralPath $MatrixPath) -and -not $Force) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $existing[[string]$row.fixture] = $row
    }
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in $Fixtures) {
    $runDir = Join-Path $WorkRoot $fixture.ToLowerInvariant()
    $prior = $existing[$fixture]
    if ($null -ne $prior -and -not $Force -and
        (-not $RerunFailed -or [string]$prior.run_status -eq 'ok')) {
        Write-Host "resume: $fixture"
        $rows.Add($prior)
        continue
    }

    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $pqr = Resolve-C5Fixture -Name $fixture -RunDir $runDir
    if ([string]::IsNullOrWhiteSpace($pqr) -or -not (Test-Path -LiteralPath $pqr)) {
        throw "Missing PQR fixture for ${fixture}: $pqr"
    }
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter $MaxIter -Restart 32 -Tolerance 1e-6 -Threads 4 `
        -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditMode off `
        -Preconditioner $Preconditioner -PreconditionerIluNeighborCap 64 `
        -NearIntegration analytic_singular -TaylorOrder 4 -CalderonTraceAudit $true
    $badArgs = @((Get-Content -LiteralPath $argPath) | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badArgs.Count -ne 0) {
        throw "C5 arg.in contains non equals-syntax options: $($badArgs -join '; ')"
    }

    Write-Host "run: $fixture (timeout ${CaseTimeoutSeconds}s)"
    $invocation = Invoke-FmmCliTitanCase -TitanExe ([System.IO.Path]::GetFullPath($TitanExe)) `
        -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $qualityPath = Join-Path $runDir 'pbe_tabi_surface_quality.csv'
    $quality = $null
    if (Test-Path -LiteralPath $qualityPath) {
        $quality = @(Import-Csv -LiteralPath $qualityPath) | Select-Object -Last 1
    }

    $runStatus = if ($invocation.timed_out) { 'timeout' } elseif ($invocation.exit_code -ne 0) { 'failed' } else { 'ok' }
    $traceStatus = if ($null -eq $quality) { 'missing_surface_quality' } else { [string]$quality.calderon_trace_space_status }
    $rows.Add([pscustomobject]@{
        fixture = $fixture
        run_status = $runStatus
        exit_code = $invocation.exit_code
        timed_out = $invocation.timed_out
        wall_seconds = $invocation.wall_seconds
        peak_rss_mb = $invocation.peak_working_set_mb
        trace_status = $traceStatus
        trace_vertices = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_vertices }
        trace_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_edges }
        trace_boundary_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_boundary_edges }
        trace_orientation_mismatches = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_orientation_mismatches }
        trace_nonmanifold_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_nonmanifold_edges }
        trace_unreferenced_vertices = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_unreferenced_vertices }
        panels = if ($null -eq $quality) { 0 } else { [int64]$quality.panel_count }
    })
}

$merged = @{}
foreach ($row in $existing.Values) {
    $merged[[string]$row.fixture] = $row
}
foreach ($row in $rows) {
    $merged[[string]$row.fixture] = $row
}
$finalRows = @($merged.Values | Sort-Object fixture)
$finalRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$allTracePass = @($finalRows | Where-Object { $_.trace_status -eq 'pass' }).Count -eq $finalRows.Count
$anyRuntimeFailure = @($finalRows | Where-Object { $_.trace_status -eq 'missing_surface_quality' }).Count -gt 0
$decision = if ($allTracePass) {
    'calderon_c5_trace_space_ready'
} elseif ($anyRuntimeFailure) {
    'calderon_c5_trace_runtime_limited'
} else {
    'calderon_c5_trace_mesh_gap'
}
[pscustomobject]@{ decision = $decision; rows = $finalRows.Count; quick = [bool]$Quick } |
    Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision", "rows=$($finalRows.Count)", "quick=$([bool]$Quick)") |
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5 Trace-Space Report',
    "Decision: $decision",
    'Scope: opt-in closed-surface P1/P0 trace-space readiness only; no Calderon operator or fallback is active.',
    "Protocol: production TABI triangle mesh with FMM/M2L and preconditioner=$Preconditioner; only mesh telemetry gates this diagnostic.",
    'Rows:',
    (($finalRows | ForEach-Object {
        "  $($_.fixture): run=$($_.run_status) trace=$($_.trace_status) vertices=$($_.trace_vertices) edges=$($_.trace_edges) boundary=$($_.trace_boundary_edges) orientation=$($_.trace_orientation_mismatches) nonmanifold=$($_.trace_nonmanifold_edges)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5 trace-space decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5_trace_space_ready' -and -not $AllowFailedDecision) {
    throw "C5 trace-space gate did not pass: $decision"
}
