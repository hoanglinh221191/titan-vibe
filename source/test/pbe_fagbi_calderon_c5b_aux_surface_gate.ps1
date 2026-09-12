[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string[]]$Fixtures = @('1sh1', '1a63', 'barnase'),
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5b_aux',
    [ValidateRange(0.5, 2.5)]
    [double]$GridSpacing = 1.25,
    [int]$CaseTimeoutSeconds = 900,
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

function Resolve-C5BFixture([string]$Name, [string]$RunDir) {
    switch ($Name.ToLowerInvariant()) {
        '1sh1' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr' }
        '1a63' { return Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr' }
        'barnase' { return Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $RunDir }
        default { throw "Unsupported C5B fixture '$Name'. Use 1sh1, 1a63, or barnase." }
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
$MatrixPath = Join-Path $WorkRoot 'calderon_c5b_aux_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5b_aux_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5b_aux_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5B_AUX_SURFACE_REPORT.txt'

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
    $pqr = Resolve-C5BFixture -Name $fixture -RunDir $runDir
    if ([string]::IsNullOrWhiteSpace($pqr) -or -not (Test-Path -LiteralPath $pqr)) {
        throw "Missing PQR fixture for ${fixture}: $pqr"
    }
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter $MaxIter -Restart 32 -Tolerance 1e-6 -Threads 4 `
        -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditMode off `
        -Preconditioner diagonal -PreconditionerIluNeighborCap 64 `
        -NearIntegration analytic_singular -TaylorOrder 4 -CalderonTraceAudit $true `
        -CalderonAuxSurface implicit-vdw -CalderonAuxGridSpacing $GridSpacing
    $badArgs = @((Get-Content -LiteralPath $argPath) | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badArgs.Count -ne 0) {
        throw "C5B arg.in contains non equals-syntax options: $($badArgs -join '; ')"
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
    $solverStatus = if ($null -eq $quality) { 'missing_surface_quality' } else { [string]$quality.status }
    $primaryStatus = if ($null -eq $quality) { 'missing_surface_quality' } else { [string]$quality.calderon_trace_space_status }
    $auxiliaryStatus = if ($null -eq $quality) { 'missing_surface_quality' } else { [string]$quality.calderon_aux_trace_status }
    $rows.Add([pscustomobject]@{
        fixture = $fixture
        run_status = $runStatus
        solver_status = $solverStatus
        exit_code = $invocation.exit_code
        timed_out = $invocation.timed_out
        wall_seconds = $invocation.wall_seconds
        peak_rss_mb = $invocation.peak_working_set_mb
        grid_spacing_a = $GridSpacing
        primary_trace_status = $primaryStatus
        primary_panels = if ($null -eq $quality) { 0 } else { [int64]$quality.panel_count }
        primary_boundary_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_trace_boundary_edges }
        auxiliary_mode = if ($null -eq $quality) { 'missing' } else { [string]$quality.calderon_aux_surface_mode }
        auxiliary_trace_status = $auxiliaryStatus
        auxiliary_panels = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_panels }
        auxiliary_vertices = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_vertices }
        auxiliary_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_edges }
        auxiliary_boundary_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_boundary_edges }
        auxiliary_orientation_mismatches = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_orientation_mismatches }
        auxiliary_nonmanifold_edges = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_nonmanifold_edges }
        auxiliary_unreferenced_vertices = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_unreferenced_vertices }
        auxiliary_grid_vertices = if ($null -eq $quality) { 0 } else { [int64]$quality.calderon_aux_grid_vertices }
        auxiliary_area_a2 = if ($null -eq $quality) { 0.0 } else { [double]$quality.calderon_aux_area_a2 }
        auxiliary_build_ms = if ($null -eq $quality) { 0.0 } else { [double]$quality.calderon_aux_build_ms }
    })
}

$merged = @{}
foreach ($row in $existing.Values) { $merged[[string]$row.fixture] = $row }
foreach ($row in $rows) { $merged[[string]$row.fixture] = $row }
$finalRows = @($merged.Values | Sort-Object fixture)
$finalRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$hasRuntimeIssue = @($finalRows | Where-Object {
    $_.run_status -ne 'ok' -or $_.auxiliary_trace_status -eq 'missing_surface_quality'
}).Count -gt 0
$allPrimaryOpen = @($finalRows | Where-Object {
    $_.primary_trace_status -eq 'open_surface' -and [int64]$_.primary_boundary_edges -gt 0
}).Count -eq $finalRows.Count
$allAuxPass = @($finalRows | Where-Object {
    $_.auxiliary_mode -eq 'implicit-vdw' -and $_.auxiliary_trace_status -eq 'pass' -and
    [int64]$_.auxiliary_boundary_edges -eq 0 -and
    [int64]$_.auxiliary_orientation_mismatches -eq 0 -and
    [int64]$_.auxiliary_nonmanifold_edges -eq 0 -and
    [int64]$_.auxiliary_unreferenced_vertices -eq 0 -and
    [int64]$_.auxiliary_panels -gt 0
}).Count -eq $finalRows.Count
$decision = if ($hasRuntimeIssue) {
    'calderon_c5b_aux_runtime_limited'
} elseif ($allAuxPass) {
    'calderon_c5b_watertight_aux_surface_ready'
} elseif ($allPrimaryOpen) {
    'calderon_c5b_aux_surface_quality_gap'
} else {
    'calderon_c5b_primary_trace_unexpected'
}

[pscustomobject]@{
    decision = $decision
    rows = $finalRows.Count
    quick = [bool]$Quick
    auxiliary_grid_spacing_a = $GridSpacing
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision", "rows=$($finalRows.Count)", "quick=$([bool]$Quick)",
  "auxiliary_grid_spacing_a=$GridSpacing") |
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5B Auxiliary Surface Report',
    "Decision: $decision",
    'Scope: opt-in auxiliary watertight trace-space diagnostic only; the primary TABI panel surface and solver operator are unchanged.',
    "Protocol: implicit VDW union-of-spheres surface, marching tetrahedra, vertex welding, grid spacing=$GridSpacing Angstrom.",
    'Rows:',
    (($finalRows | ForEach-Object {
        "  $($_.fixture): run=$($_.run_status) solver=$($_.solver_status) primary=$($_.primary_trace_status) primary_boundary=$($_.primary_boundary_edges) aux=$($_.auxiliary_trace_status) aux_panels=$($_.auxiliary_panels) aux_boundary=$($_.auxiliary_boundary_edges) aux_orientation=$($_.auxiliary_orientation_mismatches) aux_nonmanifold=$($_.auxiliary_nonmanifold_edges) aux_ms=$($_.auxiliary_build_ms)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5B auxiliary-surface decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5b_watertight_aux_surface_ready' -and -not $AllowFailedDecision) {
    throw "C5B auxiliary-surface gate did not pass: $decision"
}
