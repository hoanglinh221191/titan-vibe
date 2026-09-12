param()

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = Split-Path -Parent $PSScriptRoot
$header = Join-Path $repoRoot "src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.h"
$cpp = Join-Path $repoRoot "src/tabi/tabi_fmm_kifmm_sampled_panel_matvec.cpp"
$smoke = Join-Path $repoRoot "test/titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp"
$harness = Join-Path $repoRoot "test/pbe_fagbi_fmm_e9_sampled_scaling_profile_gate.ps1"
$report = Join-Path $repoRoot "TITAN_FMM_ENGINE_E9_SAMPLED_SCALING_PROFILE_REPORT.txt"

foreach ($path in @($header, $cpp, $smoke)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing required file: $path"
  }
}

$headerText = Get-Content -Raw -LiteralPath $header
$cppText = Get-Content -Raw -LiteralPath $cpp
$smokeText = Get-Content -Raw -LiteralPath $smoke

$requiredResultFields = @(
  "panel_geometry_hash",
  "unknown_vector_hash",
  "sampled_target_hash",
  "direct_reference_seconds",
  "equivalent_build_solve_seconds",
  "unique_far_source_nodes",
  "max_equivalent_rank",
  "min_equivalent_sigma_min",
  "max_equivalent_sigma_max",
  "max_equivalent_charge_to_rhs_norm"
)
foreach ($field in $requiredResultFields) {
  if ($headerText -notmatch "\b$([regex]::Escape($field))\b") {
    throw "missing E9 result field in header: $field"
  }
}

foreach ($token in @(
    "panel_geometry_hash(",
    "unknown_vector_hash(",
    "sampled_target_hash(")) {
  if ($cppText -notmatch [regex]::Escape($token)) {
    throw "missing E9 hash implementation token: $token"
  }
}

foreach ($token in @(
    "--e9-profile",
    "--theta=",
    "--order=",
    "--sample=",
    "panel_geometry_hash=",
    "unknown_vector_hash=",
    "sampled_target_hash=")) {
  if ($smokeText -notmatch [regex]::Escape($token)) {
    throw "missing E9 sampled smoke token: $token"
  }
}

$tracked = @(git -C $repoRoot ls-files |
  Where-Object { $_ -match 'src/t?abi/tabi_fmm_kifmm_(sampled_panel_matvec|real_panel_cloud|lite)\.(h|cpp)$' })
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

foreach ($path in @($harness, $report)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing E9 artifact contract file: $path"
  }
}

$harnessText = Get-Content -Raw -LiteralPath $harness
$requiredHarnessTokens = @(
  "1SH1_uniform_512_q6",
  "1SH1_uniform_512_q8",
  "1SH1_inward_normal_groove_q6",
  "1SH1_closest_opposing_surface_q6",
  "1SH1_high_aspect_q6",
  "1A63_uniform_512_q6",
  "1A63_high_aspect_q6",
  "q4_q6_q8_trend",
  "theta0p4",
  "theta0p6",
  "theta0p7",
  "fmm_e9_adversarial_scaling_profile_pass"
)
foreach ($token in $requiredHarnessTokens) {
  if ($harnessText -notmatch [regex]::Escape($token)) {
    throw "missing E9 harness token: $token"
  }
}

Write-Host "FMM E9 profile contract smoke passed"
