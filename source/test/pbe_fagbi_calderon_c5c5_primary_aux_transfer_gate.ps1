[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c5_primary_aux_transfer',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c5_primary_aux_transfer_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c5_primary_aux_transfer_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c5_primary_aux_transfer_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c5_primary_aux_transfer_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c5_primary_aux_transfer_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C5_PRIMARY_AUX_TRANSFER_REPORT.txt'

function Get-TokenValue([string]$Text, [string]$Name) {
    $pattern = '(?:^|\s)' + [regex]::Escape($Name) + '=([^\s]+)'
    $match = [regex]::Match($Text, $pattern)
    if ($match.Success) {
        return $match.Groups[1].Value
    }
    return ''
}

if ($Force -or -not (Test-Path -LiteralPath $ExecutablePath)) {
    $sources = @(
        (Join-Path $RepoRoot 'test\titan_calderon_c5c5_primary_aux_transfer_smoke.cpp'),
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
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5C.5 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.5 primary/auxiliary transfer smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.5 smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.5 primary/auxiliary trace-transfer calibration'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c5_primary_aux_transfer_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'matched_and_shifted_synthetic_primary_auxiliary_transfer'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    matched_status = Get-TokenValue -Text $diagnostic -Name 'matched_status'
    primary_panels = Get-TokenValue -Text $diagnostic -Name 'primary_panels'
    auxiliary_vertices = Get-TokenValue -Text $diagnostic -Name 'auxiliary_vertices'
    primary_to_auxiliary_max_distance = Get-TokenValue -Text $diagnostic -Name 'p2a_max'
    auxiliary_to_primary_max_distance = Get-TokenValue -Text $diagnostic -Name 'a2p_max'
    relative_area_mismatch = Get-TokenValue -Text $diagnostic -Name 'area_rel'
    forward_constant_relative = Get-TokenValue -Text $diagnostic -Name 'forward_constant'
    adjoint_relative = Get-TokenValue -Text $diagnostic -Name 'adjoint'
    adjoint_reverse_constant_relative = Get-TokenValue -Text $diagnostic -Name 'reverse_constant'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c5_primary_aux_transfer_pass'
} elseif ($diagnostic -match 'geometry coverage') {
    'calderon_c5c5_geometry_transfer_gap'
} else {
    'calderon_c5c5_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'synthetic_primary_to_auxiliary_trace_transfer_calibration'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=synthetic_primary_to_auxiliary_trace_transfer_calibration',
    "smoke_exit_code=$runExitCode",
    "forward_constant_relative=$($row.forward_constant_relative)",
    "adjoint_relative=$($row.adjoint_relative)",
    "primary_to_auxiliary_max_distance=$($row.primary_to_auxiliary_max_distance)",
    "auxiliary_to_primary_max_distance=$($row.auxiliary_to_primary_max_distance)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.5 Primary/Auxiliary Transfer Report',
    "Decision: $decision",
    'Scope: sparse panel-to-P1 transfer calibration. The reverse transfer is the explicit area-mass adjoint of the forward interpolation. This does not alter TABI, FMM/M2L, ILU0, GMRES, or fallback policy.',
    'Protocol: a matched synthetic closed surface must preserve forward constants and Riesz adjointness; a deliberately displaced primary panel surface must be rejected by the geometry-coverage gate.',
    'Interpretation: this validates transfer algebra and fail-closed geometry detection. It does not validate compatibility of the open production TABI surface with the implicit VDW auxiliary surface on proteins.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.5 primary/auxiliary-transfer decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c5_primary_aux_transfer_pass' -and -not $AllowFailedDecision) {
    throw "C5C.5 primary/auxiliary-transfer gate did not pass: $decision"
}
