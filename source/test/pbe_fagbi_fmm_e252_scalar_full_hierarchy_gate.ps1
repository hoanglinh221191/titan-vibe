param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e252_scalar_full_hierarchy'),
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
if ($LASTEXITCODE -ne 0) { throw 'E25.2 smoke build failed' }

function Parse-E252Line([string]$Line) {
  if (-not $Line.StartsWith('E250_RESULT ')) {
    throw "Unexpected E25.2 line: $Line"
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

function Ensure-Field([object]$Row, [string]$Name) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue ''
  }
}

function Get-Field([object]$Row, [string]$Name) {
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop) { return $null }
  [string]$prop.Value
}

function Test-FiniteNonnegativeField([object]$Row, [string]$Name) {
  $value = Get-Field $Row $Name
  if ([string]::IsNullOrWhiteSpace($value)) { return $false }
  try {
    $parsed = [double]$value
    return (-not [double]::IsNaN($parsed) -and -not [double]::IsInfinity($parsed) -and $parsed -ge 0.0)
  } catch {
    return $false
  }
}

function Test-NearlyEqualField([object]$Row, [string]$Name, [double]$Expected) {
  $value = Get-Field $Row $Name
  if ([string]::IsNullOrWhiteSpace($value)) { return $false }
  try {
    $parsed = [double]$value
    return ([math]::Abs($parsed - $Expected) -le (1.0e-12 * [math]::Max(1.0, [math]::Abs($Expected))))
  } catch {
    return $false
  }
}

function Run-FullCase(
  [string]$CaseId,
  [int]$Order,
  [double]$Theta,
  [double]$UpwardCheckScale,
  [double]$EquivScale,
  [double]$LocalEquivScale,
  [double]$LocalCheckScale,
  [string]$Kernel,
  [double]$Kappa
) {
  $args = @(
    '--full-hierarchy',
    "--case=$CaseId",
    "--order=$Order",
    "--theta=$Theta",
    "--upward-check-scale=$UpwardCheckScale",
    "--equiv-scale=$EquivScale",
    "--local-equiv-scale=$LocalEquivScale",
    "--local-check-scale=$LocalCheckScale",
    "--kernel=$Kernel",
    "--kappa=$Kappa"
  )
  $output = & $ExePath @args 2>&1
  $line = @($output | Where-Object { $_ -like 'E250_RESULT *' } | Select-Object -Last 1)
  if ($line.Count -ne 1) {
    throw "Case $CaseId order $Order emitted no E250_RESULT. Output: $($output -join [Environment]::NewLine)"
  }
  $row = Parse-E252Line ([string]$line[0])
  foreach ($requiredField in @(
    'upward_check_scale',
    'local_equiv_scale',
    'local_check_scale',
    'condition_policy',
    'source_box_policy',
    'l2l_seconds',
    'l2p_seconds'
  )) {
    Ensure-Field $row $requiredField
  }
  $rowId = ('{0}_q{1}_theta{2}_upcheck{3}_localequiv{4}_{5}_kappa{6}' -f
    $CaseId,
    $Order,
    ($Theta.ToString('0.##').Replace('.','p')),
    ($UpwardCheckScale.ToString('0.##').Replace('.','p')),
    ($LocalEquivScale.ToString('0.##').Replace('.','p')),
    $Kernel,
    ($Kappa.ToString('0.#####').Replace('.','p')))
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
  $row | Add-Member -NotePropertyName requested_order -NotePropertyValue $Order
  $row | Add-Member -NotePropertyName requested_upward_check_scale -NotePropertyValue $UpwardCheckScale
  $row | Add-Member -NotePropertyName requested_local_equiv_scale -NotePropertyValue $LocalEquivScale
  $row | Add-Member -NotePropertyName requested_local_check_scale -NotePropertyValue $LocalCheckScale
  $expectedCoverage = [uint64]$row.source_count * [uint64]$row.target_count
  $coverageStatus = if ([uint64]$row.coverage -eq $expectedCoverage) { 'pass' } else { 'fail' }
  $accuracyStatus = if ([double]$row.aggregate_rms -le 1.0e-6) { 'pass' } else { 'fail' }
  $rmsValue = [double]$row.aggregate_rms
  $finiteStatus = if ($row.status -eq 'ok' -and -not [double]::IsNaN($rmsValue) -and -not [double]::IsInfinity($rmsValue)) { 'pass' } else { 'fail' }
  $farStatus = if ($CaseId -like '*all_near*') {
    if ([uint64]$row.far_pairs -eq 0) { 'pass' } else { 'fail' }
  } else {
    if ([uint64]$row.far_pairs -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  }
  $stressStatus = if ($CaseId -like '*stress*') {
    if ([int]$row.m2m_translations -gt 0 -and [int]$row.l2l_translations -gt 0 -and [int]$row.l2p_target_count -eq [int]$row.target_count) { 'pass' } else { 'fail' }
  } else { 'not_applicable' }
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $row | Add-Member -NotePropertyName accuracy_status -NotePropertyValue $accuracyStatus
  $row | Add-Member -NotePropertyName finite_status -NotePropertyValue $finiteStatus
  $row | Add-Member -NotePropertyName far_status -NotePropertyValue $farStatus
  $row | Add-Member -NotePropertyName stress_status -NotePropertyValue $stressStatus
  $timingStatus = 'pass'
  foreach ($timingField in @(
    'p2m_seconds',
    'm2m_seconds',
    'm2l_seconds',
    'local_solve_seconds',
    'l2l_seconds',
    'l2p_seconds',
    'direct_seconds',
    'oracle_seconds',
    'total_seconds'
  )) {
    if (-not (Test-FiniteNonnegativeField $row $timingField)) {
      $timingStatus = 'fail'
      break
    }
  }
  $policyStatus = if (
    (Get-Field $row 'condition_policy') -eq 'qr-rank-only' -and
    (Get-Field $row 'source_box_policy') -eq 'strict-octree' -and
    (Test-NearlyEqualField $row 'upward_check_scale' $UpwardCheckScale) -and
    (Test-NearlyEqualField $row 'local_equiv_scale' $LocalEquivScale) -and
    (Test-NearlyEqualField $row 'local_check_scale' $LocalCheckScale)
  ) { 'pass' } else { 'fail' }
  $row | Add-Member -NotePropertyName timing_status -NotePropertyValue $timingStatus
  $row | Add-Member -NotePropertyName policy_status -NotePropertyValue $policyStatus
  $row
}

$rows = New-Object System.Collections.Generic.List[object]
$rows.Add((Run-FullCase 'scalar_laplace_all_near_q6' 6 1e-12 3.0 1.05 3.0 1.05 'laplace' 0.0))
$rows.Add((Run-FullCase 'scalar_laplace_separated_q6_check3_local3' 6 0.6 3.0 1.05 3.0 1.05 'laplace' 0.0))
$rows.Add((Run-FullCase 'scalar_laplace_separated_q8_check2p4_local2p4' 8 0.6 2.4 1.05 2.4 1.05 'laplace' 0.0))
$rows.Add((Run-FullCase 'scalar_yukawa_mixed_q6_kappa0p12746_check3' 6 0.6 3.0 1.05 3.0 1.05 'yukawa' 0.12746))
$rows.Add((Run-FullCase 'scalar_yukawa_mixed_q8_kappa0p12746_check2p4' 8 0.6 2.4 1.05 2.4 1.05 'yukawa' 0.12746))
if (-not $Quick) {
  $rows.Add((Run-FullCase 'scalar_laplace_stress_q6_check3' 6 0.6 3.0 1.05 3.0 1.05 'laplace' 0.0))
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e252_scalar_matrix.csv'
$TreeCsv = Join-Path $WorkRoot 'fmm_e252_tree.csv'
$TranslationsCsv = Join-Path $WorkRoot 'fmm_e252_translations.csv'
$TimingCsv = Join-Path $WorkRoot 'fmm_e252_timing.csv'
$PolicyCsv = Join-Path $WorkRoot 'fmm_e252_policy.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e252_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e252_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E252_SCALAR_FULL_HIERARCHY_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows | Select-Object row_id,node_count,leaf_count,observed_level_count,unique_half_widths,max_depth,source_count,target_count,coverage,coverage_status | Export-Csv -LiteralPath $TreeCsv -NoTypeInformation
$rows | Select-Object row_id,p2m_leaf_solves,m2m_translations,m2m_parent_solves,m2l_translations,local_node_solves,l2l_translations,l2p_target_count,local_rhs_adds,upward_rhs_adds,far_pairs,direct_pairs,far_fraction,stress_status | Export-Csv -LiteralPath $TranslationsCsv -NoTypeInformation
$rows | Select-Object row_id,p2m_seconds,m2m_seconds,m2l_seconds,local_solve_seconds,l2l_seconds,l2p_seconds,direct_seconds,oracle_seconds,total_seconds | Export-Csv -LiteralPath $TimingCsv -NoTypeInformation
$rows | Select-Object row_id,order,requested_order,theta,upward_check_scale,requested_upward_check_scale,local_equiv_scale,requested_local_equiv_scale,local_check_scale,requested_local_check_scale,condition_policy,source_box_policy,policy_status,timing_status | Export-Csv -LiteralPath $PolicyCsv -NoTypeInformation

$badFinite = @($rows | Where-Object { $_.finite_status -ne 'pass' })
$badCoverage = @($rows | Where-Object { $_.coverage_status -ne 'pass' })
$badAccuracy = @($rows | Where-Object { $_.accuracy_status -ne 'pass' })
$badFar = @($rows | Where-Object { $_.far_status -ne 'pass' })
$badStress = @($rows | Where-Object { $_.stress_status -eq 'fail' })
$badTiming = @($rows | Where-Object { $_.timing_status -ne 'pass' })
$badPolicy = @($rows | Where-Object { $_.policy_status -ne 'pass' })

$decision = if ($badTiming.Count -gt 0) {
  'fmm_e252_timing_export_gap'
} elseif ($badPolicy.Count -gt 0) {
  'fmm_e252_policy_metadata_gap'
} elseif ($badFinite.Count -gt 0 -or $badCoverage.Count -gt 0) {
  'fmm_e252_hierarchy_structural_gap'
} elseif ($badAccuracy.Count -gt 0) {
  'fmm_e252_scalar_accuracy_gap'
} elseif ($badFar.Count -gt 0 -or $badStress.Count -gt 0) {
  if ($Quick) { 'fmm_e252_scalar_leaf_hierarchy_pass_m2m_l2l_pending' } else { 'fmm_e252_hierarchy_structural_gap' }
} elseif ($Quick) {
  'fmm_e252_scalar_leaf_hierarchy_pass_m2m_l2l_pending'
} else {
  'fmm_e252_scalar_full_hierarchy_synthetic_pass'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  finite_gap_rows = $badFinite.Count
  coverage_gap_rows = $badCoverage.Count
  accuracy_gap_rows = $badAccuracy.Count
  far_gap_rows = $badFar.Count
  stress_gap_rows = $badStress.Count
  timing_gap_rows = $badTiming.Count
  policy_gap_rows = $badPolicy.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "policy=$PolicyCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.2 Scalar Full Hierarchy Report',
  '',
  "Decision: $decision",
  '',
  'E25.2 tests a synthetic scalar P2M-M2M-M2L-local-solve-L2L-L2P hierarchy on Laplace/Yukawa cases only.',
  'It remains diagnostic-only: no JCP channels, real panels, proteins, GMRES, preconditioner, manuscript, or default path is changed.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}: status={1}, order={2}, kernel={3}, rms={4}, far_fraction={5}, p2m={6}, m2m={7}, m2l={8}, local_solves={9}, l2l={10}, l2p_targets={11}, l2l_seconds={12}, l2p_seconds={13}, condition_policy={14}, source_box_policy={15}, accuracy={16}, timing={17}, policy={18}, structural={19}/{20}' -f
    $row.row_id,$row.status,$row.order,$row.kernel,$row.aggregate_rms,$row.far_fraction,$row.p2m_leaf_solves,$row.m2m_translations,$row.m2l_translations,$row.local_node_solves,$row.l2l_translations,$row.l2p_target_count,$row.l2l_seconds,$row.l2p_seconds,$row.condition_policy,$row.source_box_policy,$row.accuracy_status,$row.timing_status,$row.policy_status,$row.far_status,$row.stress_status)
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.2 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e252_scalar_full_hierarchy_synthetic_pass','fmm_e252_scalar_leaf_hierarchy_pass_m2m_l2l_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E25.2 gate did not pass: $decision"
}
