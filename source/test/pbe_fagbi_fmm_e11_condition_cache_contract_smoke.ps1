Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
$solveHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_solve.h"
$liteHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_lite.h"
$sampledHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.h"
$smoke = Join-Path $repo "test\titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp"
$harness = Join-Path $repo "test\pbe_fagbi_fmm_e11_condition_cache_audit_gate.ps1"
$report = Join-Path $repo "TITAN_FMM_ENGINE_E11_CONDITION_CACHE_AUDIT_REPORT.txt"

foreach ($path in @($solveHeader, $liteHeader, $sampledHeader, $smoke, $harness, $report)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing required E11 file: $path"
  }
}

$headers = (Get-Content -Raw -LiteralPath $solveHeader) + "`n" +
           (Get-Content -Raw -LiteralPath $liteHeader) + "`n" +
           (Get-Content -Raw -LiteralPath $sampledHeader)
foreach ($token in @(
  "KifmmConditionPolicy",
  "condition_policy",
  "condition_gate_enabled",
  "unique_equivalent_half_widths",
  "unique_equivalent_kappa_h_values",
  "equivalent_cache_hit_rate"
)) {
  if ($headers -notmatch [regex]::Escape($token)) {
    throw "missing E11 condition/cache contract token: $token"
  }
}

$smokeText = Get-Content -Raw -LiteralPath $smoke
foreach ($token in @(
  "--condition-policy=",
  "condition_policy=",
  "unique_half_widths=",
  "unique_kappa_h_values=",
  "cache_hit_rate="
)) {
  if ($smokeText -notmatch [regex]::Escape($token)) {
    throw "missing E11 smoke output token: $token"
  }
}

$harnessText = Get-Content -Raw -LiteralPath $harness
foreach ($token in @(
  "svd-full-debug",
  "svd-values",
  "qr-rank-only",
  "fmm_e11_condition_split_pass",
  "fmm_e11_cache_geometry_hostile"
)) {
  if ($harnessText -notmatch [regex]::Escape($token)) {
    throw "missing E11 harness token: $token"
  }
}

Write-Host "FMM E11 condition/cache contract smoke passed"
