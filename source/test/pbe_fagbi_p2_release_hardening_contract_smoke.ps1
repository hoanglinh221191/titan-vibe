Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSCommandPath)

function Read-RepoFile {
    param([Parameter(Mandatory = $true)][string]$RelativePath)
    $path = Join-Path $RepoRoot $RelativePath
    if (-not (Test-Path -LiteralPath $path)) {
        throw "missing required file: $RelativePath"
    }
    return Get-Content -Raw -LiteralPath $path
}

function Assert-Contains {
    param(
        [Parameter(Mandatory = $true)][string]$Haystack,
        [Parameter(Mandatory = $true)][string]$Needle,
        [Parameter(Mandatory = $true)][string]$Label
    )
    if ($Haystack -notlike "*$Needle*") {
        throw "expected $Label to contain: $Needle"
    }
}

$header = Read-RepoFile "src/tabi/tabi_fmm_m2l_prototype.h"
Assert-Contains $header "using FmmM2lConfig = PrototypeConfig;" "FMM/M2L header"
Assert-Contains $header "using SurfacePanelView = PrototypePanel;" "FMM/M2L header"
Assert-Contains $header "using JcpPanel = PrototypePanel;" "FMM/M2L header"
Assert-Contains $header "using FmmM2lContext = FmmM2lAnalyticContext;" "FMM/M2L header"

$parser = Read-RepoFile "src/cli/parse_arguments.cpp"
Assert-Contains $parser '"rescan"' "argument parser"
Assert-Contains $parser '"cached-metrics"' "argument parser"
Assert-Contains $parser '"calibrated-tail-v1"' "argument parser"
Assert-Contains $parser '"calibrated-estimator-v1"' "argument parser"

$arguments = Read-RepoFile "src/arguments.cpp"
Assert-Contains $arguments 'pbe.tabi_fmm_local_plan_builder = "legacy-rescan";' "argument normalization"
Assert-Contains $arguments 'pbe.tabi_fmm_local_plan_builder = "metric-cache";' "argument normalization"
Assert-Contains $arguments 'pbe.tabi_convergence_calibration = "stage3j";' "argument normalization"
Assert-Contains $arguments 'pbe.tabi_convergence_estimator_model = "stage3l";' "argument normalization"

$help = Read-RepoFile "src/help.cpp"
Assert-Contains $help "-pbe-bie-fmm-local-plan-builder=rescan|cached-metrics" "help text"
Assert-Contains $help "-pbe-bie-convergence-calibration=calibrated-tail-v1" "help text"
Assert-Contains $help "-pbe-bie-convergence-estimator-model=calibrated-estimator-v1" "help text"

$diagnostics = Read-RepoFile "docs/developer_diagnostics.md"
Assert-Contains $diagnostics "Release-facing names" "developer diagnostics"
Assert-Contains $diagnostics "historical" "developer diagnostics"
Assert-Contains $diagnostics "cached-metrics" "developer diagnostics"

$auditMap = Read-RepoFile "docs/fagbi_source_audit_map.md"
Assert-Contains $auditMap "tabi_solver.cpp" "source audit map"
Assert-Contains $auditMap "tabi_fmm_m2l_preconditioner" "source audit map"
Assert-Contains $auditMap "tabi_fmm_m2l_prototype" "source audit map"
Assert-Contains $auditMap "pbe_tabi_audit_profile.csv" "source audit map"
Assert-Contains $auditMap "internal/future" "source audit map"

"p2 release hardening contract ok"
