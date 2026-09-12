[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c3_global_oracle',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c3_global_oracle_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c3_global_oracle_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c3_global_oracle_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c3_global_oracle_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c3_global_oracle_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C3_GLOBAL_ORACLE_REPORT.txt'

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
        (Join-Path $RepoRoot 'test\titan_calderon_c5c3_global_oracle_smoke.cpp'),
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
        throw "C5C.3 compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C.3 global-oracle smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C.3 smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C.3 screened K4 global direct-action oracle'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c3_global_oracle_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'closed_auxiliary_sphere_global_direct_action'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    coarse_status = Get-TokenValue -Text $diagnostic -Name 'coarse_status'
    coarse_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'coarse_p1'
    coarse_global_direct_weak_relative = Get-TokenValue -Text $diagnostic -Name 'coarse_global'
    fine_status = Get-TokenValue -Text $diagnostic -Name 'fine_status'
    fine_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'fine_p1'
    fine_direct_pre_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'fine_direct_pre_sym'
    fine_direct_symmetry_relative = Get-TokenValue -Text $diagnostic -Name 'fine_direct_sym'
    fine_global_negative_sign_relative = Get-TokenValue -Text $diagnostic -Name 'fine_global_negative'
    fine_global_positive_sign_relative = Get-TokenValue -Text $diagnostic -Name 'fine_global_positive'
    fine_global_direct_weak_relative = Get-TokenValue -Text $diagnostic -Name 'fine_global'
    weak_form_sign = Get-TokenValue -Text $diagnostic -Name 'weak_sign'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c3_global_direct_oracle_pass'
} elseif ($diagnostic -match 'global direct action gate') {
    'calderon_c5c3_global_direct_action_gap'
} else {
    'calderon_c5c3_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'closed_surface_screened_k4_global_direct_action_oracle'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=closed_surface_screened_k4_global_direct_action_oracle',
    "smoke_exit_code=$runExitCode",
    "coarse_global_direct_weak_relative=$($row.coarse_global_direct_weak_relative)",
    "fine_global_direct_weak_relative=$($row.fine_global_direct_weak_relative)",
    "weak_form_sign=$($row.weak_form_sign)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C.3 Global Oracle Report',
    "Decision: $decision",
    'Scope: global dense direct-action calibration for screened JCP K4 on synthetic closed auxiliary spheres. It does not alter TABI, FMM/M2L, ILU0, GMRES, or fallback policy.',
    'Protocol: direct K4 uses Duffy integration for self panels and recursively subdivided quadrature for adjacent panels. Both dense direct and Maue weak matrices are symmetrized only after recording their source-only adaptive assembly asymmetry.',
    'Gate: the positive Maue sign must outperform the negative sign, global direct/weak relative discrepancy must be at most 0.25, and it must not worsen by more than 5 percent under h-refinement.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C.3 global-oracle decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c3_global_direct_oracle_pass' -and -not $AllowFailedDecision) {
    throw "C5C.3 global-oracle gate did not pass: $decision"
}
