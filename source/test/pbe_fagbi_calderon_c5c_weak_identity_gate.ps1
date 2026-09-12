[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5c_weak_identity',
    [string]$Compiler = 'clang++',
    [switch]$Force,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$ExecutablePath = Join-Path $WorkRoot 'titan_calderon_c5c_weak_identity_smoke.exe'
$BuildLogPath = Join-Path $WorkRoot 'calderon_c5c_weak_identity_build.log'
$MatrixPath = Join-Path $WorkRoot 'calderon_c5c_weak_identity_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'calderon_c5c_weak_identity_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'calderon_c5c_weak_identity_summary.dat'
$ReportPath = Join-Path $WorkRoot 'TITAN_CALDERON_C5C_WEAK_IDENTITY_REPORT.txt'

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
        (Join-Path $RepoRoot 'test\titan_calderon_c5c_weak_identity_smoke.cpp'),
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
        throw "C5C compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++20', '-O2', "-I$RepoRoot\src") + $sources + @('-o', $ExecutablePath)
    Write-Host "build: C5C weak identity smoke with $Compiler"
    $buildOutput = & $Compiler @arguments 2>&1
    $buildExitCode = $LASTEXITCODE
    @($buildOutput | ForEach-Object { [string]$_ }) |
        Set-Content -LiteralPath $BuildLogPath -Encoding ASCII
    if ($buildExitCode -ne 0) {
        throw "C5C smoke build failed; see $BuildLogPath"
    }
}

Write-Host 'run: C5C dense weak-form Laplace identity smoke'
$runOutput = & $ExecutablePath 2>&1
$runExitCode = $LASTEXITCODE
$diagnostic = (@($runOutput | ForEach-Object { [string]$_ }) -join ' ').Trim()
$smokePassed = $runExitCode -eq 0 -and $diagnostic -match 'calderon_c5c_weak_identity_smoke=pass'

$row = [pscustomobject]@{
    diagnostic_kind = 'synthetic_closed_auxiliary_sphere'
    run_status = if ($smokePassed) { 'ok' } else { 'failed' }
    exit_code = $runExitCode
    coarse_status = Get-TokenValue -Text $diagnostic -Name 'coarse_status'
    fine_status = Get-TokenValue -Text $diagnostic -Name 'fine_status'
    coarse_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'coarse_p1'
    fine_p1_dofs = Get-TokenValue -Text $diagnostic -Name 'fine_p1'
    coarse_strong_identity_relative = Get-TokenValue -Text $diagnostic -Name 'coarse_strong'
    fine_strong_identity_relative = Get-TokenValue -Text $diagnostic -Name 'fine_strong_identity'
    fine_weak_identity_mass_scaled_relative = Get-TokenValue -Text $diagnostic -Name 'fine_weak_identity_mass_scaled'
    fine_weak_identity_raw_frobenius_relative = Get-TokenValue -Text $diagnostic -Name 'fine_weak_identity_raw'
    diagnostic = $diagnostic
}
$row | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation

$decision = if ($smokePassed) {
    'calderon_c5c_laplace_identity_pass'
} elseif ($diagnostic -match 'identity_gap') {
    'calderon_c5c_identity_gap'
} else {
    'calderon_c5c_harness_gap'
}
[pscustomobject]@{
    decision = $decision
    scope = 'dense_laplace_p1_galerkin_calibration_only'
    production_solver_changed = $false
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    'scope=dense_laplace_p1_galerkin_calibration_only',
    "smoke_exit_code=$runExitCode",
    "fine_strong_identity_relative=$($row.fine_strong_identity_relative)",
    "fine_weak_identity_raw_frobenius_relative=$($row.fine_weak_identity_raw_frobenius_relative)"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN Calderon C5C Weak-Form Identity Report',
    "Decision: $decision",
    'Scope: dense Laplace P1-Galerkin calibration on a synthetic closed auxiliary sphere. It does not alter TABI, FMM/M2L, ILU0, GMRES, or fallback policy.',
    'Protocol: consistent P1 mass, Duffy self integration, recursively subdivided adjacent source pairs, weak Maue hypersingular form, and h-refinement from 74 to 164 P1 degrees of freedom.',
    'Gate: the mass-scaled strong-form Calderon identity must pass; the raw weak coefficient Frobenius value is retained only as a basis-dependent assembly diagnostic.',
    "Result: $diagnostic"
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "C5C weak-identity decision: $decision"
Write-Host "Matrix: $MatrixPath"
if ($decision -ne 'calderon_c5c_laplace_identity_pass' -and -not $AllowFailedDecision) {
    throw "C5C weak-identity gate did not pass: $decision"
}
