param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e2541_closeout'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$E254WorkRoot = Join-Path $WorkRoot 'e254'

$kifmmSourceFiles = @(Get-ChildItem -LiteralPath (Join-Path $RepoRoot 'src\tabi') -File -Filter 'tabi_fmm_kifmm*' | ForEach-Object { $_.Name })
$stalePatterns = @('*copy*','*stage*','*old*','*bak*','*tmp*')
foreach ($pattern in $stalePatterns) {
  $matches = @($kifmmSourceFiles | Where-Object { $_ -like $pattern })
  if ($matches.Count -gt 0) {
    throw "stale KIFMM source/header variants found: $($matches -join ', ')"
  }
}

& (Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e254_jcp_accuracy_repair_gate.ps1') -WorkRoot $E254WorkRoot -Cxx $Cxx -AllowFailedDecision

$decisionCsv = Join-Path $E254WorkRoot 'fmm_e254_decision.csv'
$matrixCsv = Join-Path $E254WorkRoot 'fmm_e254_matrix.csv'
$decision = @(Import-Csv -LiteralPath $decisionCsv)
$rows = @(Import-Csv -LiteralPath $matrixCsv)

if ($decision.Count -ne 1) {
  throw "expected one decision row, found $($decision.Count)"
}
if ($decision[0].decision -ne 'fmm_e254_synthetic_k1k4_q8_strict_accuracy_pass_optional_high_order_gap') {
  throw "unexpected E25.4.1 decision label: $($decision[0].decision)"
}

$requiredColumns = @(
  'flat_execution_status',
  'flat_accuracy_status',
  'flat_accuracy_failure',
  'flat_status_legacy'
)
foreach ($column in $requiredColumns) {
  if ($null -eq $rows[0].PSObject.Properties[$column]) {
    throw "missing closeout matrix column: $column"
  }
}

$optionalRows = @($rows | Where-Object { $_.optional_order_sweep -eq 'True' -or $_.optional_order_sweep -eq 'true' })
if ($optionalRows.Count -lt 2) {
  throw "expected optional high-order diagnostic rows, found $($optionalRows.Count)"
}
foreach ($row in $optionalRows) {
  if ($row.status -eq 'process_failed' -or $row.process_exit_code -ne '0') {
    throw "optional high-order row did not fail closed cleanly: $($row.row_id), status=$($row.status), exit=$($row.process_exit_code)"
  }
  if ($row.execution_status -ne 'fail' -or $row.accuracy_failure -notlike 'unsupported*') {
    throw "optional high-order row lacks explicit unsupported status: $($row.row_id), failure=$($row.accuracy_failure)"
  }
}

Write-Host 'FMM-E25.4.1 closeout smoke passed'
