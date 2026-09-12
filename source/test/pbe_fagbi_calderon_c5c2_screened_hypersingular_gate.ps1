[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c2_screened_hypersingular',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c2_screened_hypersingular_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c2_screened_hypersingular_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c2_screened_hypersingular_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c2_screened_hypersingular_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c2_screened_hypersingular_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C2_SCREENED_HYPERSINGULAR_REPORT.txt'

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
        (Join-Path $RepoRoot 'test\titan_calderon_c5c2_screened_hypersingular_smoke.cpp'),
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
        throw "C5C.2 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.2 screened-hypersingular smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.2 smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.2 screened K4 weak-form calibration'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c2_screened_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'synthetic_closed_auxiliary_sphere'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    status = Get-TokenValue -Text $diagnostic -Name 'status'
    p1_dofs = Get-TokenValue -Text $diagnostic -Name 'p1'
    triangles = Get-TokenValue -Text $diagnostic -Name 'triangles'
    separated_source_evaluations = Get-TokenValue -Text $diagnostic -Name 'separated_evaluations'
    weak_pre_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'weak_pre_sym'
    weak_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'weak_sym'
    negative_sign_pairwise_diagnostic = Get-TokenValue -Text $diagnostic -Name 'negative_rel'
    positive_sign_pairwise_diagnostic = Get-TokenValue -Text $diagnostic -Name 'positive_rel'
    weak_form_sign = Get-TokenValue -Text $diagnostic -Name 'sign'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c2_screened_weak_form_structural_pass'
} elseif ($diagnostic -match 'weak_symmetry_gap') {
    'calderon_c5c2_weak_symmetry_gap'
} elseif ($diagnostic -match 'screened_weak_parity_gap') {
    'calderon_c5c2_pairwise_nonoracle_detected'
} else {
    'calderon_c5c2_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'screened_k4_maue_weak_form_structural_calibration_only'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=screened_k4_maue_weak_form_structural_calibration_only',
    "smoke_exit_code=$runExitCode",
    "weak_symmetry_relative=$($row.weak_symmetry_relative)",
    "weak_form_sign=$($row.weak_form_sign)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.2 Screened Hypersingular Report',
    "Decision: $decision",
    'Scope: dense P1 screened K4 weak-form structural calibration on a synthetic closed auxiliary sphere. It does not alter TABI, FMM/M2L, ILU0, GMRES, or fallback policy.',
    'Protocol: the modified-Helmholtz Maue bilinear form uses the G_kappa-G_0 surface-curl term plus the kappa-squared G_kappa normal term. The positive sign matches TITAN JCP K4 because the production mixed-normal helper already carries the hypersingular sign convention.',
    'The direct separated-triangle comparison is retained only as a non-oracle diagnostic: Maue integration by parts is a closed-surface identity and cannot be validated pairwise after element-boundary terms are discarded.',
    'A global direct-action or analytic-sphere oracle remains required in C5C.3 before a complete mixed JCP operator or fallback lifecycle can be claimed.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.2 screened-hypersingular decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c2_screened_weak_form_structural_pass' -and -not $AllowFailedDecision) {
    throw "C5C.2 screened-hypersingular gate did not pass: $decision"
}
