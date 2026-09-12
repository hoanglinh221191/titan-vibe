param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e11_condition_cache')
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent $PSScriptRoot
$liteHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_lite.h"
$sampledHeader = Join-Path $repo "src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.h"
$harness = Join-Path $repo "test\pbe_fagbi_fmm_e11_condition_cache_audit_gate.ps1"
$report = Join-Path $repo "TITAN_FMM_ENGINE_E11_CONDITION_CACHE_AUDIT_REPORT.txt"

foreach ($path in @($liteHeader, $sampledHeader, $harness, $report)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "missing required E11.1 file: $path"
  }
}

$liteText = Get-Content -Raw -LiteralPath $liteHeader
if ($liteText -notmatch 'condition_policy\{"qr-rank-only"\}') {
  throw 'KifmmLiteConfig default condition_policy must be qr-rank-only'
}

$sampledText = Get-Content -Raw -LiteralPath $sampledHeader
foreach ($token in @("condition_computed", "condition_gate_enabled")) {
  if ($sampledText -notmatch [regex]::Escape($token)) {
    throw "missing sampled condition reporting token: $token"
  }
}

$harnessText = Get-Content -Raw -LiteralPath $harness
foreach ($token in @(
  "condition_computed",
  "fmm_e11_condition_split_pass_cache_geometry_pending",
  "fmm_e11_condition_split_pass_factor_cache_poor"
)) {
  if ($harnessText -notmatch [regex]::Escape($token)) {
    throw "missing harness/reporting token: $token"
  }
}

$decision = Join-Path $WorkRoot 'fmm_e11_decision.csv'
$matrix = Join-Path $WorkRoot 'fmm_e11_matrix.csv'
if ((Test-Path -LiteralPath $decision) -and (Test-Path -LiteralPath $matrix)) {
  $decisionRow = Import-Csv -LiteralPath $decision | Select-Object -First 1
  if ($decisionRow.decision -eq 'fmm_e11_condition_split_pass') {
    throw 'E11.1 decision wording must include cache-geometry pending/factor-cache-poor scope'
  }
  $qr = Import-Csv -LiteralPath $matrix | Where-Object { $_.condition_policy -eq 'qr-rank-only' } | Select-Object -First 1
  if ($null -ne $qr) {
    if ($qr.condition_computed -ne 'false') {
      throw 'qr-rank-only row must report condition_computed=false'
    }
    if ([double]::Parse($qr.max_condition, [Globalization.CultureInfo]::InvariantCulture) -ne -1.0) {
      throw 'qr-rank-only row must report max_condition=-1 sentinel'
    }
  }
}

Write-Host "FMM E11.1 condition reporting contract smoke passed"
