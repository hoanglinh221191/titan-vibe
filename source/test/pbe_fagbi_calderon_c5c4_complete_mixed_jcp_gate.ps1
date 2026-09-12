[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c4_complete_mixed_jcp',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c4_complete_mixed_jcp_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c4_complete_mixed_jcp_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c4_complete_mixed_jcp_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c4_complete_mixed_jcp_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c4_complete_mixed_jcp_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C4_COMPLETE_MIXED_JCP_REPORT.txt'

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
        (Join-Path $RepoRoot 'test\titan_calderon_c5c4_complete_mixed_jcp_smoke.cpp'),
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
        throw "C5C.4 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.4 complete mixed-JCP smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.4 smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.4 complete mixed JCP global direct-action calibration'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c4_complete_mixed_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'closed_auxiliary_sphere_complete_mixed_jcp_global_direct_action'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    coarse_status = Get-TokenValue -Text $diagnostic -Name 'coarse_status'
    coarse_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'coarse_p1'
    coarse_strong_system_relative = Get-TokenValue -Text $diagnostic -Name 'coarse_system'
    fine_status = Get-TokenValue -Text $diagnostic -Name 'fine_status'
    fine_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'fine_p1'
    fine_dual_p0_dofs = Get-TokenValue -Text $diagnostic -Name 'fine_dual'
    fine_dual_p0_subtriangles = Get-TokenValue -Text $diagnostic -Name 'fine_patches'
    fine_pairing_pivot_ratio = Get-TokenValue -Text $diagnostic -Name 'fine_pairing_pivot'
    fine_pairing_inverse_relative = Get-TokenValue -Text $diagnostic -Name 'fine_pairing_inverse'
    fine_direct_k4_weak_relative = Get-TokenValue -Text $diagnostic -Name 'fine_k4'
    fine_strong_interaction_relative = Get-TokenValue -Text $diagnostic -Name 'fine_interaction'
    fine_strong_system_relative = Get-TokenValue -Text $diagnostic -Name 'fine_system'
    fine_deterministic_probe_relative = Get-TokenValue -Text $diagnostic -Name 'fine_probe'
    lhs_phi = Get-TokenValue -Text $diagnostic -Name 'lhs_phi'
    lhs_phin = Get-TokenValue -Text $diagnostic -Name 'lhs_phin'
    weak_form_sign = Get-TokenValue -Text $diagnostic -Name 'weak_sign'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c4_complete_mixed_jcp_pass'
} elseif ($diagnostic -match 'global direct oracle|strong action') {
    'calderon_c5c4_complete_mixed_direct_action_gap'
} else {
    'calderon_c5c4_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'complete_mixed_jcp_closed_surface_global_direct_calibration'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=complete_mixed_jcp_closed_surface_global_direct_calibration',
    "smoke_exit_code=$runExitCode",
    "fine_direct_k4_weak_relative=$($row.fine_direct_k4_weak_relative)",
    "fine_strong_interaction_relative=$($row.fine_strong_interaction_relative)",
    "fine_strong_system_relative=$($row.fine_strong_system_relative)",
    "fine_deterministic_probe_relative=$($row.fine_deterministic_probe_relative)",
    "weak_form_sign=$($row.weak_form_sign)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.4 Complete Mixed-JCP Report',
    "Decision: $decision",
    'Scope: dense closed-surface calibration of the complete mixed P1/dual-P0 JCP interaction block. It is diagnostic only and does not alter TABI, FMM/M2L, ILU0, GMRES, fallback policy, or defaults.',
    'Riesz map: diag(B, B^T), where B is the explicit dual-P0/P1 pairing. The complete operator follows the current solver convention A = jump diagonal - K.',
    'Protocol: K1/K2/K3 use the production JCP kernel coefficients. K4 compares the positive screened Maue weak form with a global direct action using Duffy self integration and subdivided adjacent-pair quadrature.',
    'Gate: Riesz pairing must be invertible, the K4 and full strong interaction/system actions must agree with the dense direct reference within 0.025, and the h-refined sphere may not worsen by more than 5 percent.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.4 complete mixed-JCP decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c4_complete_mixed_jcp_pass' -and -not $AllowFailedDecision) {
    throw "C5C.4 complete mixed-JCP gate did not pass: $decision"
}
