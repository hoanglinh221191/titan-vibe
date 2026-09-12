param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e253b_jcp_accuracy_closure'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_jcp_accuracy_smoke.exe'

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
  throw 'No C++ compiler found for FMM-E25.3b JCP accuracy smoke.'
}

function Parse-KeyValueLine([string]$Line, [string]$Prefix) {
  if (-not $Line.StartsWith($Prefix)) {
    throw "Unexpected line for ${Prefix}: $Line"
  }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring($Prefix.Length).Trim().Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function FiniteDouble([object]$Row, [string]$Name) {
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop) { return $false }
  try {
    $value = [double]$prop.Value
    return (-not [double]::IsNaN($value) -and -not [double]::IsInfinity($value))
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
  $rowId = ('{0}_q{1}_theta{2}_upcheck{3}_localequiv{4}_localcheck{5}' -f
    $CaseId,
    $Order,
    ($Theta.ToString('0.##').Replace('.','p')),
    ($UpwardCheckScale.ToString('0.##').Replace('.','p')),
    ($LocalEquivScale.ToString('0.##').Replace('.','p')),
    ($LocalCheckScale.ToString('0.##').Replace('.','p')))
  $output = & $ExePath @args 2>&1
  $resultLine = @($output | Where-Object { $_ -like 'E253B_RESULT *' } | Select-Object -Last 1)
  if ($LASTEXITCODE -ne 0 -or $resultLine.Count -ne 1) {
    $row = [pscustomobject][ordered]@{
      case = $CaseId
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'execution'
      aggregate_gate = '1.0e-6'
      channel_gate = '5.0e-6'
      order = $Order
      theta = $Theta
      upward_check_scale = $UpwardCheckScale
      local_equiv_scale = $LocalEquivScale
      local_check_scale = $LocalCheckScale
      condition_policy = 'unknown'
      source_box_policy = 'unknown'
      aggregate_rms = 'NaN'
      k1_rms = 'NaN'
      k2_rms = 'NaN'
      k3_rms = 'NaN'
      k4_rms = 'NaN'
      max_channel_rms = 'NaN'
      max_channel_name = 'none'
      k4_abs_rms = 'NaN'
      k4_reference_rms = 'NaN'
      flat_status = 'not_evaluated'
      flat_aggregate_rms = 'NaN'
      flat_k1_rms = 'NaN'
      flat_k2_rms = 'NaN'
      flat_k3_rms = 'NaN'
      flat_k4_rms = 'NaN'
      flat_max_channel_rms = 'NaN'
      flat_max_channel_name = 'none'
      flat_k4_abs_rms = 'NaN'
      hier_vs_flat_aggregate_rms = 'NaN'
      hier_vs_flat_k1_rms = 'NaN'
      hier_vs_flat_k2_rms = 'NaN'
      hier_vs_flat_k3_rms = 'NaN'
      hier_vs_flat_k4_rms = 'NaN'
      hier_vs_flat_max_channel_rms = 'NaN'
      hier_vs_flat_max_channel_name = 'none'
      far_fraction = 'NaN'
      far_pairs = '0'
      direct_pairs = '0'
      coverage = '0'
      source_count = '0'
      target_count = '0'
      m2m_translations = '0'
      m2l_translations = '0'
      l2l_translations = '0'
      l2p_target_count = '0'
      cache_entries = '0'
      total_seconds = 'NaN'
      row_id = $rowId
      coverage_status = 'fail'
      strict_accuracy_status = 'fail'
      flat_comparator_applicable = ($CaseId -like '*separated*')
      flat_accuracy_status = if ($CaseId -like '*separated*') { 'fail' } else { 'not_applicable' }
      hier_vs_flat_status = if ($CaseId -like '*separated*') { 'fail' } else { 'not_applicable' }
      structural_status = 'fail'
      finite_status = 'fail'
      execution_gate_status = 'fail'
      process_exit_code = $LASTEXITCODE
      process_output = ($output -join ' | ')
      optional_order_sweep = (($CaseId -like '*check4*') -or ($CaseId -like '*q10*'))
    }
    return [pscustomobject]@{
      Row = $row
      Fields = @()
      K4 = @()
    }
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E253B_RESULT '
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $LASTEXITCODE
  $row | Add-Member -NotePropertyName process_output -NotePropertyValue ''
  $row | Add-Member -NotePropertyName optional_order_sweep -NotePropertyValue (($CaseId -like '*check4*') -or ($CaseId -like '*q10*'))

  if ($row.execution_status -ne 'pass' -or $row.status -ne 'ok') {
    $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName strict_accuracy_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName flat_comparator_applicable -NotePropertyValue ($CaseId -like '*separated*')
    $row | Add-Member -NotePropertyName flat_accuracy_status -NotePropertyValue $(if ($CaseId -like '*separated*') { 'fail' } else { 'not_applicable' }) -Force
    $row | Add-Member -NotePropertyName hier_vs_flat_status -NotePropertyValue $(if ($CaseId -like '*separated*') { 'fail' } else { 'not_applicable' })
    $row | Add-Member -NotePropertyName structural_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName finite_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName execution_gate_status -NotePropertyValue 'fail'
    return [pscustomobject]@{
      Row = $row
      Fields = @()
      K4 = @()
    }
  }

  $expectedCoverage = [uint64]$row.source_count * [uint64]$row.target_count
  $coverageStatus = if ([uint64]$row.coverage -eq $expectedCoverage) { 'pass' } else { 'fail' }
  $executionStatus = if ($row.execution_status -eq 'pass' -and $row.status -eq 'ok') { 'pass' } else { 'fail' }
  $strictAccuracyStatus = if ([double]$row.aggregate_rms -le [double]$row.aggregate_gate -and [double]$row.max_channel_rms -le [double]$row.channel_gate) { 'pass' } else { 'fail' }
  $flatApplicable = ($CaseId -like '*separated*')
  $flatExecutionStatus = if ($null -ne $row.PSObject.Properties['flat_execution_status']) {
    $row.flat_execution_status
  } elseif ($null -ne $row.PSObject.Properties['flat_status_legacy']) {
    $row.flat_status_legacy
  } elseif ($null -ne $row.PSObject.Properties['flat_status']) {
    $row.flat_status
  } else {
    'not_evaluated'
  }
  $flatAccuracyStatus = if (-not $flatApplicable) {
    'not_applicable'
  } elseif ($flatExecutionStatus -eq 'pass' -and [double]$row.flat_aggregate_rms -le [double]$row.aggregate_gate -and [double]$row.flat_max_channel_rms -le [double]$row.channel_gate) {
    'pass'
  } else {
    'fail'
  }
  $hierFlatStatus = if (-not $flatApplicable) {
    'not_applicable'
  } elseif ([double]$row.hier_vs_flat_aggregate_rms -le [double]$row.aggregate_gate -and [double]$row.hier_vs_flat_max_channel_rms -le [double]$row.channel_gate) {
    'pass'
  } else {
    'fail'
  }
  $structuralStatus = if ($CaseId -like '*all_near*') {
    if ([uint64]$row.far_pairs -eq 0 -and [int]$row.m2l_translations -eq 0) { 'pass' } else { 'fail' }
  } elseif ($CaseId -like '*mixed*') {
    if ([uint64]$row.far_pairs -gt 0 -and [uint64]$row.direct_pairs -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  } elseif ($CaseId -like '*multi_m2l*') {
    if ([int]$row.m2m_translations -gt 0 -and [int]$row.l2l_translations -gt 0 -and [int]$row.m2l_translations -ge 16) { 'pass' } else { 'fail' }
  } elseif ($CaseId -like '*stress*') {
    if ([int]$row.m2m_translations -gt 0 -and [int]$row.l2l_translations -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  } else {
    if ([uint64]$row.far_pairs -gt 0 -and [int]$row.m2l_translations -gt 0) { 'pass' } else { 'fail' }
  }
  $finiteStatus = 'pass'
  foreach ($field in @('aggregate_rms','k1_rms','k2_rms','k3_rms','k4_rms','max_channel_rms','flat_aggregate_rms','flat_max_channel_rms','hier_vs_flat_aggregate_rms','hier_vs_flat_max_channel_rms','total_seconds')) {
    if (-not (FiniteDouble $row $field)) { $finiteStatus = 'fail'; break }
  }
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $row | Add-Member -NotePropertyName strict_accuracy_status -NotePropertyValue $strictAccuracyStatus
  $row | Add-Member -NotePropertyName flat_comparator_applicable -NotePropertyValue $flatApplicable
  $row | Add-Member -NotePropertyName flat_accuracy_status -NotePropertyValue $flatAccuracyStatus -Force
  $row | Add-Member -NotePropertyName hier_vs_flat_status -NotePropertyValue $hierFlatStatus
  $row | Add-Member -NotePropertyName structural_status -NotePropertyValue $structuralStatus
  $row | Add-Member -NotePropertyName finite_status -NotePropertyValue $finiteStatus
  $row | Add-Member -NotePropertyName execution_gate_status -NotePropertyValue $executionStatus

  $fields = @($output | Where-Object { $_ -like 'E253B_FIELD *' } | ForEach-Object {
      $f = Parse-KeyValueLine ([string]$_) 'E253B_FIELD '
      $f | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
      $f
    })
  $k4 = @($output | Where-Object { $_ -like 'E253B_K4 *' } | ForEach-Object {
      $r = Parse-KeyValueLine ([string]$_) 'E253B_K4 '
      $r | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
      $r
    })

  [pscustomobject]@{
    Row = $row
    Fields = $fields
    K4 = $k4
  }
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
  'test\titan_fmm_kifmm_hierarchy_jcp_accuracy_smoke.cpp',
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
    throw 'failed to build titan_fmm_kifmm_hierarchy_jcp_accuracy_smoke'
  }
} finally {
  Pop-Location
}

$cases = @(
  @{ Case='jcp_all_near_q6'; Order=6; Theta=1e-12; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_laplace_yukawa_separated_q6_check3_local3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_laplace_yukawa_separated_q8_check2p4_local2p4'; Order=8; Theta=0.6; Up=2.4; Eq=1.05; LocalEq=2.4; LocalCheck=1.05 },
  @{ Case='jcp_mixed_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 }
)
if (-not $Quick) {
  $cases += @(
    @{ Case='jcp_laplace_yukawa_separated_q8_check3_local3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_laplace_yukawa_separated_q8_check4_local4'; Order=8; Theta=0.6; Up=4.0; Eq=1.05; LocalEq=4.0; LocalCheck=1.05 },
    @{ Case='jcp_laplace_yukawa_separated_q10_check3_local3'; Order=10; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_stress_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_multi_m2l_stress_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 }
  )
}

$rows = New-Object System.Collections.Generic.List[object]
$fieldRows = New-Object System.Collections.Generic.List[object]
$k4Rows = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
  $run = Run-Case $case.Case $case.Order $case.Theta $case.Up $case.Eq $case.LocalEq $case.LocalCheck
  $rows.Add($run.Row)
  foreach ($field in $run.Fields) { $fieldRows.Add($field) }
  foreach ($k4 in $run.K4) { $k4Rows.Add($k4) }
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e253b_matrix.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e253b_channel.csv'
$FieldCsv = Join-Path $WorkRoot 'fmm_e253b_field_matrix.csv'
$K4Csv = Join-Path $WorkRoot 'fmm_e253b_k4_decomposition.csv'
$FlatCsv = Join-Path $WorkRoot 'fmm_e253b_flat_vs_hierarchy.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e253b_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e253b_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E253B_JCP_ACCURACY_CLOSURE_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows | Select-Object row_id,aggregate_gate,channel_gate,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,max_channel_rms,max_channel_name,accuracy_status,strict_accuracy_status,accuracy_failure | Export-Csv -LiteralPath $ChannelCsv -NoTypeInformation
$fieldRows | Export-Csv -LiteralPath $FieldCsv -NoTypeInformation
$k4Rows | Export-Csv -LiteralPath $K4Csv -NoTypeInformation
$rows | Select-Object row_id,flat_status,flat_aggregate_rms,flat_k1_rms,flat_k2_rms,flat_k3_rms,flat_k4_rms,flat_max_channel_rms,flat_max_channel_name,hier_vs_flat_aggregate_rms,hier_vs_flat_k1_rms,hier_vs_flat_k2_rms,hier_vs_flat_k3_rms,hier_vs_flat_k4_rms,hier_vs_flat_max_channel_rms,hier_vs_flat_max_channel_name,flat_accuracy_status,hier_vs_flat_status | Export-Csv -LiteralPath $FlatCsv -NoTypeInformation

$badExecution = @($rows | Where-Object { $_.execution_gate_status -ne 'pass' -or $_.finite_status -ne 'pass' })
$badRequiredExecution = @($badExecution | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badOptionalExecution = @($badExecution | Where-Object { [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badCoverage = @($rows | Where-Object { $_.coverage_status -ne 'pass' })
$badRequiredCoverage = @($badCoverage | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badStructure = @($rows | Where-Object { $_.structural_status -ne 'pass' })
$badRequiredStructure = @($badStructure | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badAccuracy = @($rows | Where-Object { $_.strict_accuracy_status -ne 'pass' })
$badRequiredAccuracy = @($badAccuracy | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badFlat = @($rows | Where-Object { $_.flat_accuracy_status -eq 'fail' })
$badHierFlat = @($rows | Where-Object { $_.hier_vs_flat_status -eq 'fail' })

$decision = if ($badRequiredExecution.Count -gt 0) {
  'fmm_e253b_execution_gap'
} elseif ($badRequiredCoverage.Count -gt 0) {
  'fmm_e253b_jcp_coverage_gap'
} elseif ($badRequiredStructure.Count -gt 0) {
  'fmm_e253b_jcp_hierarchy_structural_gap'
} elseif ($badRequiredAccuracy.Count -gt 0) {
  if ($badFlat.Count -gt 0 -and $badHierFlat.Count -gt 0) {
    if ($badOptionalExecution.Count -gt 0) { 'fmm_e253b_jcp_flat_and_hierarchy_accuracy_gap_with_order_sweep_execution_gap' } else { 'fmm_e253b_jcp_flat_and_hierarchy_accuracy_gap' }
  } elseif ($badFlat.Count -eq 0 -and $badHierFlat.Count -gt 0) {
    if ($badOptionalExecution.Count -gt 0) { 'fmm_e253b_jcp_hierarchy_translation_accuracy_gap_with_order_sweep_execution_gap' } else { 'fmm_e253b_jcp_hierarchy_translation_accuracy_gap' }
  } elseif ($badFlat.Count -gt 0) {
    if ($badOptionalExecution.Count -gt 0) { 'fmm_e253b_jcp_flat_source_equivalent_accuracy_gap_with_order_sweep_execution_gap' } else { 'fmm_e253b_jcp_flat_source_equivalent_accuracy_gap' }
  } else {
    if ($badOptionalExecution.Count -gt 0) { 'fmm_e253b_jcp_k1k4_accuracy_gap_with_order_sweep_execution_gap' } else { 'fmm_e253b_jcp_k1k4_accuracy_gap' }
  }
} elseif ($badOptionalExecution.Count -gt 0) {
  'fmm_e253b_order_sweep_execution_gap'
} elseif ($Quick) {
  'fmm_e253b_jcp_strict_quick_pass_full_pending'
} else {
  'fmm_e253b_jcp_strict_accuracy_pass'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  execution_gap_rows = $badExecution.Count
  required_execution_gap_rows = $badRequiredExecution.Count
  optional_execution_gap_rows = $badOptionalExecution.Count
  coverage_gap_rows = $badCoverage.Count
  structural_gap_rows = $badStructure.Count
  strict_accuracy_gap_rows = $badAccuracy.Count
  required_strict_accuracy_gap_rows = $badRequiredAccuracy.Count
  flat_accuracy_gap_rows = $badFlat.Count
  hier_vs_flat_gap_rows = $badHierFlat.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "field_matrix=$FieldCsv",
  "k4_decomposition=$K4Csv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.3b JCP Accuracy Closure Report',
  '',
  "Decision: $decision",
  '',
  'E25.3b is a strict synthetic K1-K4 accuracy closure gate for the full hierarchy.',
  'It keeps E25.3 as the structural hierarchy checkpoint and applies the stricter gates: aggregate RMS <= 1e-6 and max channel RMS <= 5e-6.',
  'Flat source-equivalent rows are interpreted only for separated far-only rows where sampled targets are outside the source check surface; all-near and mixed rows keep flat diagnostics as non-gating context.',
  ''
)
foreach ($row in $rows) {
  $report += ('  {0}: exec={1}, strict={2}, rms={3}, max={4}:{5}, flat={6}/{7}, hier_vs_flat={8}/{9}, m2m={10}, m2l={11}, l2l={12}, structural={13}' -f
    $row.row_id,$row.execution_gate_status,$row.strict_accuracy_status,$row.aggregate_rms,$row.max_channel_name,$row.max_channel_rms,$row.flat_aggregate_rms,$row.flat_max_channel_rms,$row.hier_vs_flat_aggregate_rms,$row.hier_vs_flat_max_channel_rms,$row.m2m_translations,$row.m2l_translations,$row.l2l_translations,$row.structural_status)
}
if ($badOptionalExecution.Count -gt 0) {
  $report += ''
  $report += 'Optional order/scale sweep execution gaps:'
  foreach ($row in $badOptionalExecution) {
    $report += ('  {0}: exit_code={1}, status={2}' -f $row.row_id,$row.process_exit_code,$row.status)
  }
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.3b decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e253b_jcp_strict_accuracy_pass','fmm_e253b_jcp_strict_quick_pass_full_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E25.3b gate did not pass: $decision"
}
