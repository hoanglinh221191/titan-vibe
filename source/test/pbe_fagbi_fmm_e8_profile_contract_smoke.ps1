param()

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = Split-Path -Parent $PSScriptRoot
$header = Join-Path $repoRoot "src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.h"
$cpp = Join-Path $repoRoot "src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.cpp"
$harness = Join-Path $repoRoot "test/pbe_fagbi_fmm_e8_sampled_profile_gate.ps1"
$report = Join-Path $repoRoot "TITAN_FMM_ENGINE_E8_SAMPLED_PROFILE_REPORT.txt"

foreach ($path in @($header, $cpp)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing required source file: $path"
  }
}

$headerText = Get-Content -Raw -LiteralPath $header
$cppText = Get-Content -Raw -LiteralPath $cpp

$requiredHeaderFields = @(
  "direct_reference_seconds",
  "tree_build_seconds",
  "worklist_build_seconds",
  "equivalent_build_solve_seconds",
  "far_eval_seconds",
  "fallback_eval_seconds",
  "total_seconds",
  "unique_far_source_nodes",
  "cached_source_nodes",
  "source_tree_bytes",
  "target_tree_bytes",
  "worklist_bytes",
  "equivalent_charge_storage_bytes",
  "cached_channel_object_estimated_bytes",
  "aggregate_abs_rms",
  "k1_abs_rms",
  "k2_abs_rms",
  "k3_abs_rms",
  "k4_abs_rms",
  "k1_reference_rms",
  "k2_reference_rms",
  "k3_reference_rms",
  "k4_reference_rms",
  "max_equivalent_rank",
  "min_equivalent_sigma_min",
  "max_equivalent_sigma_max",
  "max_equivalent_rhs_norm",
  "max_equivalent_charge_norm",
  "max_equivalent_charge_to_rhs_norm"
)

foreach ($field in $requiredHeaderFields) {
  if ($headerText -notmatch "\b$([regex]::Escape($field))\b") {
    throw "missing E8 sampled profile field in header: $field"
  }
}

$requiredCppTokens = @(
  "std::chrono::steady_clock",
  "tree_memory_bytes",
  "channel_cache_memory_bytes",
  "abs_rms_rows",
  "reference_rms_channel"
)
foreach ($token in $requiredCppTokens) {
  if ($cppText -notmatch [regex]::Escape($token)) {
    throw "missing E8 sampled profile implementation token: $token"
  }
}

foreach ($path in @($harness, $report)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing E8 artifact contract file: $path"
  }
}

$tracked = @(git -C $repoRoot ls-files |
  Where-Object { $_ -match 'src/t/?abi/tabi_fmm_kifmm_(sampled_panel_matvec|real_panel_cloud|lite)\.(h|cpp)$' })
$expectedTracked = @(
  'src/tabi/tabi_fmm_kifmm_lite.cpp',
  'src/tabi/tabi_fmm_kifmm_lite.h',
  'src/tabi/tabi_fmm_kifmm_real_panel_cloud.cpp',
  'src/tabi/tabi_fmm_kifmm_real_panel_cloud.h',
  'src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.cpp',
  'src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.h'
)
$unexpected = @($tracked | Where-Object { $expectedTracked -notcontains $_ })
$missing = @($expectedTracked | Where-Object { $tracked -notcontains $_ })
if ($unexpected.Count -ne 0 -or $missing.Count -ne 0) {
  throw "canonical KIFMM-lite source/header set mismatch; unexpected=[$($unexpected -join ',')]; missing=[$($missing -join ',')]"
}

$harnessText = Get-Content -Raw -LiteralPath $harness
foreach ($token in @(
    "fmm_e8_uniform_quick_profile_pass_full_adversarial_pending",
    "quick_uniform")) {
  if ($harnessText -notmatch [regex]::Escape($token)) {
    throw "missing E8 decision-scope token in harness: $token"
  }
}

$smokePath = Join-Path $repoRoot "test/titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp"
$smokeText = Get-Content -Raw -LiteralPath $smokePath
if ($smokeText -notmatch "uniform_hash_q4_diagnostic") {
  throw "missing q4 diagnostic row in sampled panel smoke"
}

Write-Host "FMM E8 profile contract smoke passed"
