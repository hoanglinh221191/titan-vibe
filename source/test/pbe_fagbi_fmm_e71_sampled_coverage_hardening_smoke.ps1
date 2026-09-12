param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$Source = Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp'
$Report = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E71_SAMPLED_COVERAGE_HARDENING_REPORT.txt'

if (-not (Test-Path -LiteralPath $Source)) {
    throw "Sampled panel matvec source not found: $Source"
}

$SourceText = Get-Content -Raw -LiteralPath $Source
$RequiredStatus = @(
    'direct_reference_coverage_gap',
    'worklist_coverage_gap',
    'fallback_coverage_gap',
    'nonproduction_fallback_used'
)

foreach ($status in $RequiredStatus) {
    if ($SourceText -notmatch [regex]::Escape($status)) {
        throw "Missing fail-closed sampled coverage status: $status"
    }
}

$RequiredExpressions = @(
    'production_direct_reference_valid_pairs != expected_pairs',
    'accepted_far_panel_pairs +',
    'direct_fallback_panel_pairs != expected_pairs',
    'production_fallback_valid_pairs !=',
    'nonproduction_fallback_pairs != 0'
)

foreach ($expr in $RequiredExpressions) {
    if ($SourceText -notmatch [regex]::Escape($expr)) {
        throw "Missing sampled coverage invariant expression: $expr"
    }
}

if (-not (Test-Path -LiteralPath $Report)) {
    throw "E7.1 coverage hardening report not found: $Report"
}

$ReportText = Get-Content -Raw -LiteralPath $Report
foreach ($status in $RequiredStatus) {
    if ($ReportText -notmatch [regex]::Escape($status)) {
        throw "E7.1 report does not document status: $status"
    }
}

Write-Host 'PASS pbe_fagbi_fmm_e71_sampled_coverage_hardening_smoke'
