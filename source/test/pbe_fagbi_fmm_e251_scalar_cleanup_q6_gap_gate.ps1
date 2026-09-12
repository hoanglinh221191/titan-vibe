param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e251_scalar_q6_gap'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$SmokeScript = Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e250_hierarchical_kifmm_smoke.ps1'
$SmokeWorkRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $SmokeWorkRoot 'titan_fmm_kifmm_hierarchy_smoke.exe'

if (Test-Path -LiteralPath $WorkRoot) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

& $SmokeScript -WorkRoot $SmokeWorkRoot -Cxx $Cxx
if ($LASTEXITCODE -ne 0) { throw 'E25.1 smoke build failed' }

function Parse-E251Line([string]$Line) {
  if (-not $Line.StartsWith('E250_RESULT ')) {
    throw "Unexpected E25.1 line: $Line"
  }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring('E250_RESULT '.Length).Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function Run-Case([string]$CaseId, [int]$Order, [double]$Theta, [double]$CheckScale, [double]$EquivScale, [string]$Kernel, [double]$Kappa) {
  $args = @(
    "--case=$CaseId",
    "--order=$Order",
    "--theta=$Theta",
    "--upward-check-scale=$CheckScale",
    "--equiv-scale=$EquivScale",
    "--kernel=$Kernel",
    "--kappa=$Kappa"
  )
  $output = & $ExePath @args 2>&1
  $line = @($output | Where-Object { $_ -like 'E250_RESULT *' } | Select-Object -Last 1)
  if ($line.Count -ne 1) {
    throw "Case $CaseId order $Order theta $Theta check $CheckScale emitted no E250_RESULT. Output: $($output -join [Environment]::NewLine)"
  }
  $row = Parse-E251Line ([string]$line[0])
  $kappaLabel = $Kappa.ToString('0.#####').Replace('.','p')
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue ("{0}_q{1}_theta{2}_check{3}_{4}_kappa{5}" -f $CaseId,$Order,($Theta.ToString('0.##').Replace('.','p')),($CheckScale.ToString('0.##').Replace('.','p')),$Kernel,$kappaLabel)
  $row | Add-Member -NotePropertyName requested_order -NotePropertyValue $Order
  $row | Add-Member -NotePropertyName requested_theta -NotePropertyValue $Theta
  $row | Add-Member -NotePropertyName requested_check_scale -NotePropertyValue $CheckScale
  $row | Add-Member -NotePropertyName requested_equiv_scale -NotePropertyValue $EquivScale
  $timingStatus = if ([double]$row.p2m_seconds -le [double]$row.total_seconds -and [double]$row.m2m_seconds -eq 0.0) { 'pass' } else { 'fail' }
  $widthStatus = if ([int]$row.unique_half_widths -eq [int]$row.observed_level_count) { 'pass' } else { 'fail' }
  $coverageStatus = if ([uint64]$row.coverage -eq ([uint64]$row.source_count * [uint64]$row.target_count)) { 'pass' } else { 'fail' }
  $accuracyStatus = if ([double]$row.aggregate_rms -le 1.0e-6) { 'pass' } else { 'fail' }
  $farStatus = if ([double]$row.far_fraction -gt 0.0) { 'pass' } else { 'fail' }
  $row | Add-Member -NotePropertyName timing_status -NotePropertyValue $timingStatus
  $row | Add-Member -NotePropertyName width_status -NotePropertyValue $widthStatus
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $row | Add-Member -NotePropertyName accuracy_status -NotePropertyValue $accuracyStatus
  $row | Add-Member -NotePropertyName far_status -NotePropertyValue $farStatus
  $row
}

$rows = New-Object System.Collections.Generic.List[object]
$thetaSet = if ($Quick) { @(0.6, 0.4) } else { @(0.3, 0.4, 0.5, 0.6) }
$checkSet = if ($Quick) { @(2.4, 3.0) } else { @(2.4, 3.0, 3.5) }

foreach ($theta in $thetaSet) {
  foreach ($check in $checkSet) {
    $rows.Add((Run-Case 'scalar_laplace_separated' 6 $theta $check 1.05 'laplace' 0.0))
  }
}
$rows.Add((Run-Case 'scalar_laplace_separated' 8 0.6 2.4 1.05 'laplace' 0.0))
foreach ($kappa in @(0.05, 0.12746, 0.25)) {
  $rows.Add((Run-Case 'scalar_yukawa_mixed' 8 0.6 2.4 1.05 'yukawa' $kappa))
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e251_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e251_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e251_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E251_SCALAR_Q6_GAP_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$badTiming = @($rows | Where-Object { $_.timing_status -ne 'pass' })
$badCoverage = @($rows | Where-Object { $_.coverage_status -ne 'pass' })
$q6Rows = @($rows | Where-Object { [int]$_.requested_order -eq 6 })
$q8Rows = @($rows | Where-Object { [int]$_.requested_order -eq 8 })
$q6Pass = @($rows | Where-Object { [int]$_.requested_order -eq 6 -and $_.accuracy_status -eq 'pass' -and $_.far_status -eq 'pass' })
$q8Pass = @($rows | Where-Object { [int]$_.requested_order -eq 8 -and $_.accuracy_status -eq 'pass' -and $_.far_status -eq 'pass' })

$q6Rms = @($q6Rows | ForEach-Object { [double]$_.aggregate_rms })
$q8Rms = @($q8Rows | ForEach-Object { [double]$_.aggregate_rms })
$q6MinRms = if ($q6Rms.Count -gt 0) { ($q6Rms | Measure-Object -Minimum).Minimum } else { [double]::NaN }
$q6MaxRms = if ($q6Rms.Count -gt 0) { ($q6Rms | Measure-Object -Maximum).Maximum } else { [double]::NaN }
$q8MinRms = if ($q8Rms.Count -gt 0) { ($q8Rms | Measure-Object -Minimum).Minimum } else { [double]::NaN }
$q8MaxRms = if ($q8Rms.Count -gt 0) { ($q8Rms | Measure-Object -Maximum).Maximum } else { [double]::NaN }

$decision = if ($badTiming.Count -gt 0) {
  'fmm_e251_timing_metadata_gap'
} elseif ($badCoverage.Count -gt 0) {
  'fmm_e251_coverage_gap'
} elseif ($q6Pass.Count -gt 0) {
  'fmm_e251_scalar_q6_pass'
} elseif ($q8Pass.Count -gt 0) {
  'fmm_e251_scalar_q6_gap_q8_selected'
} else {
  'fmm_e251_scalar_acceptance_gap'
}

[pscustomobject]@{
  decision = $decision
  rows = $rows.Count
  q6_pass_rows = $q6Pass.Count
  q6_rows = $q6Rows.Count
  q6_min_rms = $q6MinRms
  q6_max_rms = $q6MaxRms
  q8_pass_rows = $q8Pass.Count
  q8_rows = $q8Rows.Count
  q8_min_rms = $q8MinRms
  q8_max_rms = $q8MaxRms
  timing_gap_rows = $badTiming.Count
  coverage_gap_rows = $badCoverage.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.1 Scalar Cleanup And q6 Gap Report',
  '',
  "Decision: $decision",
  '',
  'E25.1 fixes E25.0 scalar hierarchy metadata, adds source-check-surface acceptance accounting, and diagnoses whether q6 remains the scalar bottleneck after the corrected guard.',
  'It remains scalar-only and diagnostic-only: no JCP real panels, GMRES, preconditioner, manuscript, or default path is changed.',
  '',
  ('q6 RMS range: {0} to {1}; q8 RMS range: {2} to {3}.' -f $q6MinRms,$q6MaxRms,$q8MinRms,$q8MaxRms),
  'The corrected separated fixture keeps nonzero far work while ensuring accepted pairs pass the source check-surface guard.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}: order={1}, theta={2}, check={3}, kernel={4}, rms={5}, far_fraction={6}, timing={7}, coverage={8}, mac={9}, check_reject={10}' -f
    $row.row_id,$row.order,$row.theta,$row.requested_check_scale,$row.kernel,$row.aggregate_rms,$row.far_fraction,$row.timing_status,$row.coverage_status,$row.mac_accepted_pairs,$row.check_guard_rejected_pairs)
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.1 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e251_scalar_q6_pass','fmm_e251_scalar_q6_gap_q8_selected') -and -not $AllowFailedDecision) {
  throw "FMM-E25.1 gate did not pass: $decision"
}
