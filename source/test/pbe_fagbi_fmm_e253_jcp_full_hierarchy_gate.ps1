param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e253_jcp_full_hierarchy'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_jcp_smoke.exe'

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) {
      throw "Requested C++ compiler not found: $Requested"
    }
    return $Requested
  }
  foreach ($candidate in @(
      'C:\msys64\ucrt64\bin\clang++.exe',
      'C:\msys64\mingw64\bin\clang++.exe',
      'C:\msys64\ucrt64\bin\g++.exe',
      'C:\msys64\mingw64\bin\g++.exe')) {
    if (Test-Path -LiteralPath $candidate) { return $candidate }
  }
  $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  $cmd = Get-Command g++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  throw 'No C++ compiler found for FMM-E25.3 hierarchy JCP smoke.'
}

function Parse-E253Line([string]$Line) {
  if (-not $Line.StartsWith('E253_RESULT ')) {
    throw "Unexpected E25.3 line: $Line"
  }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring('E253_RESULT '.Length).Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function Test-FiniteNonnegative([object]$Row, [string]$Name) {
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) {
    return $false
  }
  try {
    $value = [double]$prop.Value
    return (-not [double]::IsNaN($value) -and -not [double]::IsInfinity($value) -and $value -ge 0.0)
  } catch {
    return $false
  }
}

function Run-Case(
  [string]$CaseId,
  [int]$Order,
  [double]$Theta,
  [double]$UpwardCheckScale,
  [double]$EquivScale,
  [double]$LocalEquivScale,
  [double]$LocalCheckScale
) {
  $args = @(
    "--case=$CaseId",
    "--order=$Order",
    "--theta=$Theta",
    "--upward-check-scale=$UpwardCheckScale",
    "--equiv-scale=$EquivScale",
    "--local-equiv-scale=$LocalEquivScale",
    "--local-check-scale=$LocalCheckScale"
  )
  $output = & $ExePath @args 2>&1
  $line = @($output | Where-Object { $_ -like 'E253_RESULT *' } | Select-Object -Last 1)
  if ($line.Count -ne 1) {
    throw "Case $CaseId order $Order emitted no E253_RESULT. Output: $($output -join [Environment]::NewLine)"
  }
  $row = Parse-E253Line ([string]$line[0])
  $rowId = ('{0}_q{1}_theta{2}_upcheck{3}_localequiv{4}' -f
    $CaseId,
    $Order,
    ($Theta.ToString('0.##').Replace('.','p')),
    ($UpwardCheckScale.ToString('0.##').Replace('.','p')),
    ($LocalEquivScale.ToString('0.##').Replace('.','p')))
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId

  $expectedCoverage = [uint64]$row.source_count * [uint64]$row.target_count
  $coverageStatus = if ([uint64]$row.coverage -eq $expectedCoverage) { 'pass' } else { 'fail' }
  $aggregate = [double]$row.aggregate_rms
  $maxChannel = [double]$row.max_channel_rms
  $accuracyStatus = if ($CaseId -like '*all_near*') {
    if ($aggregate -le 1.0e-12 -and $maxChannel -le 1.0e-12) { 'pass' } else { 'fail' }
  } else {
    if ($aggregate -le 1.0e-4 -and $maxChannel -le 5.0e-4) { 'pass' } else { 'fail' }
  }
  $finiteStatus = if ($row.status -eq 'ok' -and -not [double]::IsNaN($aggregate) -and -not [double]::IsInfinity($aggregate)) { 'pass' } else { 'fail' }
  $farStatus = if ($CaseId -like '*all_near*') {
    if ([uint64]$row.far_pairs -eq 0 -and [int]$row.m2l_translations -eq 0) { 'pass' } else { 'fail' }
  } elseif ($CaseId -like '*mixed*') {
    if ([uint64]$row.far_pairs -gt 0 -and [uint64]$row.direct_pairs -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  } else {
    if ([uint64]$row.far_pairs -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  }
  $stressStatus = if ($CaseId -like '*stress*') {
    if ([int]$row.m2m_translations -gt 0 -and [int]$row.l2l_translations -gt 0 -and [int]$row.l2p_target_count -eq [int]$row.target_count) { 'pass' } else { 'fail' }
  } else { 'not_applicable' }
  $timingStatus = 'pass'
  foreach ($field in @('p2m_seconds','m2m_seconds','m2l_seconds','local_solve_seconds','l2l_seconds','l2p_seconds','direct_seconds','oracle_seconds','total_seconds')) {
    if (-not (Test-FiniteNonnegative $row $field)) {
      $timingStatus = 'fail'
      break
    }
  }
  $policyStatus = if ($row.condition_policy -eq 'qr-rank-only' -and $row.source_box_policy -eq 'strict-octree') { 'pass' } else { 'fail' }
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $row | Add-Member -NotePropertyName accuracy_status -NotePropertyValue $accuracyStatus
  $row | Add-Member -NotePropertyName finite_status -NotePropertyValue $finiteStatus
  $row | Add-Member -NotePropertyName far_status -NotePropertyValue $farStatus
  $row | Add-Member -NotePropertyName stress_status -NotePropertyValue $stressStatus
  $row | Add-Member -NotePropertyName timing_status -NotePropertyValue $timingStatus
  $row | Add-Member -NotePropertyName policy_status -NotePropertyValue $policyStatus
  $row
}

if (Test-Path -LiteralPath $WorkRoot) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null

$Compiler = Resolve-Cxx $Cxx
$CompileArgs = @(
  '-std=c++23',
  '-O2',
  '-Isrc',
  '-Isrc\eigen',
  'test\titan_fmm_kifmm_hierarchy_jcp_smoke.cpp',
  'src\tabi\tabi_fmm_kifmm_surface.cpp',
  'src\tabi\tabi_fmm_kifmm_kernels.cpp',
  'src\tabi\tabi_fmm_kifmm_solve.cpp',
  'src\tabi\tabi_fmm_kifmm_fields.cpp',
  'src\tabi\tabi_fmm_kifmm_hierarchy.cpp',
  'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp',
    'src\tabi\tabi_fmm_kifmm_cartesian_bridge.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_treecode_series.cpp',
  '-o',
  $ExePath
)

Push-Location $RepoRoot
try {
  & $Compiler @CompileArgs
  if ($LASTEXITCODE -ne 0) {
    throw 'failed to build titan_fmm_kifmm_hierarchy_jcp_smoke'
  }
} finally {
  Pop-Location
}

$rows = New-Object System.Collections.Generic.List[object]
$rows.Add((Run-Case 'jcp_all_near_q6' 6 1e-12 3.0 1.05 3.0 1.05))
$rows.Add((Run-Case 'jcp_laplace_yukawa_separated_q6_check3_local3' 6 0.6 3.0 1.05 3.0 1.05))
$rows.Add((Run-Case 'jcp_laplace_yukawa_separated_q8_check2p4_local2p4' 8 0.6 2.4 1.05 2.4 1.05))
$rows.Add((Run-Case 'jcp_mixed_q6_check3' 6 0.6 3.0 1.05 3.0 1.05))
if (-not $Quick) {
  $rows.Add((Run-Case 'jcp_stress_q6_check3' 6 0.6 3.0 1.05 3.0 1.05))
}

$q6 = @($rows | Where-Object { $_.case -eq 'jcp_laplace_yukawa_separated_q6_check3_local3' } | Select-Object -First 1)
$q8 = @($rows | Where-Object { $_.case -eq 'jcp_laplace_yukawa_separated_q8_check2p4_local2p4' } | Select-Object -First 1)
$q8Status = 'not_applicable'
if ($q6.Count -eq 1 -and $q8.Count -eq 1) {
  if ([double]$q8[0].aggregate_rms -le (1.25 * [double]$q6[0].aggregate_rms) -and [double]$q8[0].max_channel_rms -le (1.25 * [double]$q6[0].max_channel_rms)) {
    $q8Status = 'pass'
  } else {
    $q8Status = 'fail'
  }
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e253_matrix.csv'
$TimingCsv = Join-Path $WorkRoot 'fmm_e253_timing.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e253_channel.csv'
$TreeCsv = Join-Path $WorkRoot 'fmm_e253_tree.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e253_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e253_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E253_JCP_FULL_HIERARCHY_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows | Select-Object row_id,k1_rms,k2_rms,k3_rms,k4_rms,max_channel_rms,max_channel_name,aggregate_rms,accuracy_status | Export-Csv -LiteralPath $ChannelCsv -NoTypeInformation
$rows | Select-Object row_id,p2m_seconds,m2m_seconds,m2l_seconds,local_solve_seconds,l2l_seconds,l2p_seconds,direct_seconds,oracle_seconds,total_seconds,timing_status | Export-Csv -LiteralPath $TimingCsv -NoTypeInformation
$rows | Select-Object row_id,node_count,leaf_count,observed_level_count,unique_half_widths,max_depth,source_count,target_count,coverage,coverage_status,far_pairs,direct_pairs,far_fraction,far_status,stress_status | Export-Csv -LiteralPath $TreeCsv -NoTypeInformation

$badFinite = @($rows | Where-Object { $_.finite_status -ne 'pass' })
$badCoverage = @($rows | Where-Object { $_.coverage_status -ne 'pass' })
$badAccuracy = @($rows | Where-Object { $_.accuracy_status -ne 'pass' })
$badFar = @($rows | Where-Object { $_.far_status -ne 'pass' })
$badStress = @($rows | Where-Object { $_.stress_status -eq 'fail' })
$badTiming = @($rows | Where-Object { $_.timing_status -ne 'pass' })
$badPolicy = @($rows | Where-Object { $_.policy_status -ne 'pass' })

$decision = if ($badFinite.Count -gt 0 -or $badTiming.Count -gt 0 -or $badPolicy.Count -gt 0) {
  'fmm_e253_jcp_hierarchy_structural_gap'
} elseif ($badCoverage.Count -gt 0) {
  'fmm_e253_jcp_coverage_gap'
} elseif ($badAccuracy.Count -gt 0 -or $q8Status -eq 'fail') {
  'fmm_e253_jcp_channel_accuracy_gap'
} elseif ($badFar.Count -gt 0 -or $badStress.Count -gt 0) {
  if ($Quick) { 'fmm_e253_jcp_leaf_hierarchy_pass_m2m_l2l_pending' } else { 'fmm_e253_jcp_hierarchy_structural_gap' }
} elseif ($Quick) {
  'fmm_e253_jcp_leaf_hierarchy_pass_m2m_l2l_pending'
} else {
  'fmm_e253_jcp_full_hierarchy_synthetic_pass'
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
  q8_confirmation_status = $q8Status
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.3 JCP Full Hierarchy Report',
  '',
  "Decision: $decision",
  '',
  'E25.3 tests a synthetic eight-field K1-K4 hierarchy on deterministic panel clouds only.',
  'It remains diagnostic-only: no real panels, proteins, GMRES, preconditioner, manuscript, or default path is changed.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}: status={1}, order={2}, rms={3}, max_channel={4}:{5}, far_fraction={6}, p2m={7}, m2m={8}, m2l={9}, local={10}, l2l={11}, l2p_targets={12}, accuracy={13}, coverage={14}, structural={15}/{16}' -f
    $row.row_id,$row.status,$row.order,$row.aggregate_rms,$row.max_channel_name,$row.max_channel_rms,$row.far_fraction,$row.p2m_leaf_solves,$row.m2m_translations,$row.m2l_translations,$row.local_node_solves,$row.l2l_translations,$row.l2p_target_count,$row.accuracy_status,$row.coverage_status,$row.far_status,$row.stress_status)
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.3 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e253_jcp_full_hierarchy_synthetic_pass','fmm_e253_jcp_leaf_hierarchy_pass_m2m_l2l_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E25.3 gate did not pass: $decision"
}
