[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c1_mixed_trace',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c1_mixed_trace_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c1_mixed_trace_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c1_mixed_trace_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c1_mixed_trace_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c1_mixed_trace_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C1_MIXED_TRACE_REPORT.txt'

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
        (Join-Path $RepoRoot 'test\titan_calderon_c5c1_mixed_trace_smoke.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_bilevel_grid.cpp'),
        (Join-Path $RepoRoot 'src\titan_executor.cpp'),
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_weak_form.cpp')
    )
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5C.1 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.1 mixed-trace smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.1 smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.1 dense P1/dual-P0 mixed-trace calibration'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c1_mixed_trace_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'synthetic_closed_auxiliary_sphere'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    status = Get-TokenValue -Text $diagnostic -Name 'status'
    p1_dofs = Get-TokenValue -Text $diagnostic -Name 'p1'
    dual_p0_dofs = Get-TokenValue -Text $diagnostic -Name 'dual_p0'
    dual_p0_patches = Get-TokenValue -Text $diagnostic -Name 'patches'
    dual_p0_partition_relative = Get-TokenValue -Text $diagnostic -Name 'partition'
    pairing_pivot_ratio = Get-TokenValue -Text $diagnostic -Name 'pivot'
    pairing_inverse_relative = Get-TokenValue -Text $diagnostic -Name 'inverse'
    k1_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'k1_sym'
    k2_k3_negative_epsilon_closure_relative = Get-TokenValue -Text $diagnostic -Name 'k23_negative_eps'
    k2_k3_positive_epsilon_closure_relative = Get-TokenValue -Text $diagnostic -Name 'k23_positive_eps'
    k2_k3_separated_nonclosure_relative = Get-TokenValue -Text $diagnostic -Name 'k23_nonclosure'
    k4_separated_pair_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'k4_separated_sym'
    zero_screening_kernel_absolute = Get-TokenValue -Text $diagnostic -Name 'zero'
    derivative_block_relation = Get-TokenValue -Text $diagnostic -Name 'k23_relation'
    screened_hypersingular_status = Get-TokenValue -Text $diagnostic -Name 'k4_status'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c1_mixed_trace_pairing_pass'
} elseif ($diagnostic -match 'mixed_riesz_pairing_gap') {
    'calderon_c5c1_mixed_riesz_pairing_gap'
} elseif ($diagnostic -match 'jcp_kernel_parity_gap') {
    'calderon_c5c1_jcp_kernel_parity_gap'
} else {
    'calderon_c5c1_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'p1_primal_dual_p0_riesz_and_jcp_kernel_calibration_only'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=p1_primal_dual_p0_riesz_and_jcp_kernel_calibration_only',
    "smoke_exit_code=$runExitCode",
    "pairing_pivot_ratio=$($row.pairing_pivot_ratio)",
    "k1_symmetry_relative=$($row.k1_symmetry_relative)",
    "k4_separated_pair_symmetry_relative=$($row.k4_separated_pair_symmetry_relative)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.1 Mixed Trace Report',
    "Decision: $decision",
    'Scope: dense P1/primal plus discontinuous P0/barycentric-dual trace-pairing calibration on a synthetic closed auxiliary sphere. It does not alter TABI, FMM/M2L, ILU0, GMRES, or fallback policy.',
    'Protocol: vertex-dual P0 patches form a square P1/dual-P0 Riesz pairing. The gate checks surface partition, pairing inversion, K1 reciprocity, separated-pair K4 reciprocity, and the epsilon=1/kappa=0 zero-screening limit.',
    'JCP K2 and K3 are reported as separate dielectric-weighted blocks. Their distinct Laplace/Yukawa coefficients do not admit a scalar-epsilon adjoint identity; the reported nonclosure is diagnostic, not a failed symmetry claim.',
    'Screened hypersingular K4 weak regularization remains a separate C5C.2 requirement before any fallback lifecycle work.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.1 mixed-trace decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c1_mixed_trace_pairing_pass' -and -not $AllowFailedDecision) {
    throw "C5C.1 mixed-trace gate did not pass: $decision"
}
