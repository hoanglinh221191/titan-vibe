Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
$header = Join-Path $repo "src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.h"
$fieldsHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_fields.h"
$solveHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_solve.h"
$smoke = Join-Path $repo "test\titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp"
$harness = Join-Path $repo "test\pbe_fagbi_fmm_e10_factor_reuse_gate.ps1"
$report = Join-Path $repo "TITAN_FMM_ENGINE_E10_FACTOR_REUSE_REPORT.txt"

foreach ($path in @($header, $fieldsHeader, $solveHeader, $smoke, $harness, $report)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing required E10 file: $path"
  }
}

$sampledText = Get-Content -Raw -LiteralPath $header
foreach ($token in @(
  "equivalent_matrix_build_seconds",
  "equivalent_condition_seconds",
  "equivalent_factor_seconds",
  "equivalent_rhs_build_seconds",
  "equivalent_solve_seconds",
  "equivalent_factor_cache_hits",
  "equivalent_factor_cache_misses"
)) {
  if ($sampledText -notmatch [regex]::Escape($token)) {
    throw "missing sampled-result telemetry field: $token"
  }
}

$fieldsText = (Get-Content -Raw -LiteralPath $fieldsHeader) + "`n" + (Get-Content -Raw -LiteralPath $solveHeader)
foreach ($token in @(
  "KifmmEquivalentSystemCache",
  "KifmmEquivalentSystemCacheStats",
  "equivalent_factor_cache_hits",
  "equivalent_rhs_build_seconds"
)) {
  if ($fieldsText -notmatch [regex]::Escape($token)) {
    throw "missing field/cache contract token: $token"
  }
}

$smokeText = Get-Content -Raw -LiteralPath $smoke
foreach ($token in @(
  "--e10-factor-profile",
  "equivalent_factor_cache_hits=",
  "equivalent_factor_cache_misses=",
  "equivalent_condition_seconds=",
  "equivalent_rhs_build_seconds=",
  "equivalent_solve_seconds="
)) {
  if ($smokeText -notmatch [regex]::Escape($token)) {
    throw "missing E10 smoke output token: $token"
  }
}

$harnessText = Get-Content -Raw -LiteralPath $harness
foreach ($token in @(
  "1SH1_uniform_128_q6_cached",
  "1SH1_uniform_512_q6_cached",
  "fmm_e10_factor_reuse_pass",
  "fmm_e10_speed_gap",
  "equivalent_factor_cache_hits"
)) {
  if ($harnessText -notmatch [regex]::Escape($token)) {
    throw "missing E10 harness token: $token"
  }
}

Write-Host "FMM E10 factor-reuse contract smoke passed"
