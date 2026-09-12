param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$solver = Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp'
$harness = Join-Path $RepoRoot 'test\pbe_fagbi_fmm_e21_kifmm_gmres_smoke.ps1'
$report = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E21_KIFMM_GMRES_REPORT.txt'

foreach ($path in @($solver, $harness, $report)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing E21.1 contract input: $path"
    }
}

$solverText = Get-Content -Raw -LiteralPath $solver
$harnessText = Get-Content -Raw -LiteralPath $harness
$reportText = Get-Content -Raw -LiteralPath $report

foreach ($needle in @(
    'kifmm_iteration_invariants.flush()',
    'pbe_tabi_kifmm_gmres_progress.csv',
    'prepare_start',
    'prepare_done',
    'preconditioner_start',
    'preconditioner_done',
    'gmres_start',
    'matvec_',
    'post_audit_start',
    'post_audit_done')) {
    if (-not $solverText.Contains($needle)) {
        throw "E21.1 solver contract missing '$needle'"
    }
}

foreach ($needle in @(
    '[int]$MaxIter',
    '[int]$Restart',
    'timeout_phase',
    'preconditioner_graph_source',
    'run_completed',
    'iterations_observed',
    'not_observed_due_timeout',
    'pbe_tabi_kifmm_gmres_progress.csv')) {
    if (-not $harnessText.Contains($needle)) {
        throw "E21.1 harness contract missing '$needle'"
    }
}

if (-not $reportText.Contains('runtime-limited before a complete convergence/audit record was produced')) {
    throw 'E21.1 report wording must avoid convergence claims after timeout.'
}

Write-Host 'PASS pbe_fagbi_fmm_e211_gmres_smoke_contract'
