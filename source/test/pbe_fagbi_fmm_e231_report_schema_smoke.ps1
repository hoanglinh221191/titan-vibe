param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e23_staged_gmres'),
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$matrixCsv = Join-Path $WorkRoot 'fmm_e23_matrix.csv'
$decisionCsv = Join-Path $WorkRoot 'fmm_e23_decision.csv'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E23_STAGED_GMRES_REPORT.txt'

foreach ($path in @($matrixCsv, $decisionCsv, $reportPath)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "Missing E23 artifact required for schema smoke: $path"
  }
}

$decision = @(Import-Csv -LiteralPath $decisionCsv)[0]
if ($decision.decision -ne 'fmm_e23_1sh1_restart32_gmres_pass_1a63_probe_pending') {
  throw "Unexpected E23 closeout decision: $($decision.decision)"
}

$rows = @(Import-Csv -LiteralPath $matrixCsv)
if ($rows.Count -eq 0) { throw 'E23 matrix is empty.' }

$requiredColumns = @(
  'execution_status',
  'gmres_status',
  'operator_status',
  'postprocess_status',
  'audit_status',
  'audit_interpreted'
)
$columns = @($rows[0].PSObject.Properties.Name)
foreach ($column in $requiredColumns) {
  if ($columns -notcontains $column) {
    throw "E23 matrix missing closeout schema column: $column"
  }
}

$restart = @($rows | Where-Object { $_.stage -eq 'restart32' -and $_.fixture -eq '1sh1' })
if ($restart.Count -ne 1) { throw 'Missing 1SH1 restart32 row.' }
if ($restart[0].gmres_status -ne 'converged') { throw "1SH1 restart32 gmres_status is $($restart[0].gmres_status)" }
if ($restart[0].audit_status -ne 'pass') { throw "1SH1 restart32 audit_status is $($restart[0].audit_status)" }

$oneA63 = @($rows | Where-Object { $_.fixture -eq '1a63' })
if ($oneA63.Count -gt 0) {
  foreach ($row in $oneA63) {
    if ($row.gmres_status -ne 'max_iter_nonconverged') {
      throw "1A63 row $($row.stage) unexpectedly reports gmres_status=$($row.gmres_status)"
    }
    if ($row.audit_interpreted -ne 'false') {
      throw "1A63 nonconverged row $($row.stage) should not interpret post-solve audit."
    }
  }
}

$report = Get-Content -Raw -LiteralPath $reportPath
foreach ($needle in @(
    'Decision: fmm_e23_1sh1_restart32_gmres_pass_1a63_probe_pending',
    '1A63 probe and short rows are not interpreted as convergence or post-solve audit evidence',
    'not a general KIFMM-GMRES backend pass')) {
  if (-not $report.Contains($needle)) {
    throw "E23 report missing closeout wording: $needle"
  }
}

Write-Host 'FMM-E23.1 report/schema smoke: pass'
