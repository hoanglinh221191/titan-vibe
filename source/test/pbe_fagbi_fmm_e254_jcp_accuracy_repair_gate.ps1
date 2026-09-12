param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e254_jcp_accuracy_repair'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_jcp_accuracy_smoke_e254.exe'

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
  throw 'No C++ compiler found for FMM-E25.4 JCP repair smoke.'
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

function Row-Id(
  [string]$CaseId,
  [int]$Order,
  [double]$Theta,
  [double]$UpwardCheckScale,
  [double]$LocalEquivScale,
  [double]$LocalCheckScale
) {
  '{0}_q{1}_theta{2}_upcheck{3}_localequiv{4}_localcheck{5}' -f
    $CaseId,
    $Order,
    ($Theta.ToString('0.##').Replace('.','p')),
    ($UpwardCheckScale.ToString('0.##').Replace('.','p')),
    ($LocalEquivScale.ToString('0.##').Replace('.','p')),
    ($LocalCheckScale.ToString('0.##').Replace('.','p'))
}

function Failed-Run(
  [string]$CaseId,
  [int]$Order,
  [double]$Theta,
  [double]$UpwardCheckScale,
  [double]$LocalEquivScale,
  [double]$LocalCheckScale,
  [int]$ExitCode,
  [string[]]$Output
) {
  $rowId = Row-Id $CaseId $Order $Theta $UpwardCheckScale $LocalEquivScale $LocalCheckScale
  $row = [pscustomobject][ordered]@{
    row_id = $rowId
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
    aggregate_rms = 'NaN'
    max_channel_rms = 'NaN'
    max_channel_name = 'none'
    k1_rms = 'NaN'
    k2_rms = 'NaN'
    k3_rms = 'NaN'
    k4_rms = 'NaN'
    k4_abs_rms = 'NaN'
    k4_reference_rms = 'NaN'
    flat_status_legacy = 'not_evaluated'
    flat_execution_status = 'not_evaluated'
    flat_accuracy_status = 'not_evaluated'
    flat_accuracy_failure = 'not_evaluated'
    flat_aggregate_rms = 'NaN'
    flat_max_channel_rms = 'NaN'
    hier_vs_flat_aggregate_rms = 'NaN'
    hier_vs_flat_max_channel_rms = 'NaN'
    far_pairs = '0'
    direct_pairs = '0'
    coverage = '0'
    source_count = '0'
    target_count = '0'
    m2m_translations = '0'
    m2l_translations = '0'
    l2l_translations = '0'
    l2p_target_count = '0'
    process_exit_code = $ExitCode
    process_output = ($Output -join ' | ')
    optional_order_sweep = (($CaseId -like '*check4*') -or ($CaseId -like '*q10*'))
    execution_gate_status = 'fail'
    strict_accuracy_status = 'fail'
    coverage_status = 'fail'
    structural_status = 'fail'
    stage_attribution_status = 'not_evaluated'
  }
  [pscustomobject]@{ Row = $row; Stages = @(); K4 = @() }
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
    '--e254-profile',
    "--case=$CaseId",
    "--order=$Order",
    "--theta=$Theta",
    "--upward-check-scale=$UpwardCheckScale",
    "--equiv-scale=$EquivScale",
    "--local-equiv-scale=$LocalEquivScale",
    "--local-check-scale=$LocalCheckScale"
  )
  $output = & $ExePath @args 2>&1
  $resultLine = @($output | Where-Object { $_ -like 'E253B_RESULT *' } | Select-Object -Last 1)
  if ($resultLine.Count -ne 1) {
    return Failed-Run $CaseId $Order $Theta $UpwardCheckScale $LocalEquivScale $LocalCheckScale $LASTEXITCODE @($output)
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E253B_RESULT '
  $rowId = Row-Id $CaseId $Order $Theta $UpwardCheckScale $LocalEquivScale $LocalCheckScale
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $LASTEXITCODE
  $row | Add-Member -NotePropertyName optional_order_sweep -NotePropertyValue (($CaseId -like '*check4*') -or ($CaseId -like '*q10*'))

  $stageRows = @($output | Where-Object { $_ -like 'E254_STAGE *' } | ForEach-Object {
      $stage = Parse-KeyValueLine ([string]$_) 'E254_STAGE '
      $stage | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
      $stage
    })
  $caseK4Rows = @($output | Where-Object { $_ -like 'E253B_K4 *' } | ForEach-Object {
      $k4 = Parse-KeyValueLine ([string]$_) 'E253B_K4 '
      $k4 | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
      $k4
    })
  if ($row.execution_status -ne 'pass' -or $row.status -ne 'ok') {
    $row | Add-Member -NotePropertyName execution_gate_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName strict_accuracy_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName structural_status -NotePropertyValue 'fail'
    $row | Add-Member -NotePropertyName stage_attribution_status -NotePropertyValue 'not_evaluated'
    $row | Add-Member -NotePropertyName finite_status -NotePropertyValue 'pass'
    return [pscustomobject]@{ Row = $row; Stages = $stageRows; K4 = $caseK4Rows }
  }
  if ($stageRows.Count -lt 5) {
    throw "Case $CaseId order $Order emitted $($stageRows.Count) E254_STAGE rows; expected at least 5."
  }

  $expectedCoverage = [uint64]$row.source_count * [uint64]$row.target_count
  $coverageStatus = if ([uint64]$row.coverage -eq $expectedCoverage) { 'pass' } else { 'fail' }
  $executionStatus = if ($row.execution_status -eq 'pass' -and $row.status -eq 'ok') { 'pass' } else { 'fail' }
  $strictAccuracyStatus = if ([double]$row.aggregate_rms -le [double]$row.aggregate_gate -and [double]$row.max_channel_rms -le [double]$row.channel_gate) { 'pass' } else { 'fail' }
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
  $stageStatus = if (@($stageRows | Where-Object { $_.stage_status -eq 'fail' }).Count -eq 0) { 'pass' } else { 'fail' }
  $finiteStatus = 'pass'
  foreach ($field in @('aggregate_rms','k1_rms','k2_rms','k3_rms','k4_rms','max_channel_rms','total_seconds')) {
    if (-not (FiniteDouble $row $field)) { $finiteStatus = 'fail'; break }
  }

  $row | Add-Member -NotePropertyName execution_gate_status -NotePropertyValue $executionStatus
  $row | Add-Member -NotePropertyName strict_accuracy_status -NotePropertyValue $strictAccuracyStatus
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $row | Add-Member -NotePropertyName structural_status -NotePropertyValue $structuralStatus
  $row | Add-Member -NotePropertyName stage_attribution_status -NotePropertyValue $stageStatus
  $row | Add-Member -NotePropertyName finite_status -NotePropertyValue $finiteStatus
  [pscustomobject]@{ Row = $row; Stages = $stageRows; K4 = $caseK4Rows }
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
    throw 'failed to build titan_fmm_kifmm_hierarchy_jcp_accuracy_smoke_e254'
  }
} finally {
  Pop-Location
}

$cases = @(
  @{ Case='jcp_laplace_yukawa_separated_q6_check3_local3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_laplace_yukawa_separated_q8_check3_local3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_stress_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_stress_q8_check3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_multi_m2l_stress_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
  @{ Case='jcp_multi_m2l_stress_q8_check3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 }
)
if (-not $Quick) {
  $cases = @(
    @{ Case='jcp_all_near_q6'; Order=6; Theta=1e-12; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_mixed_q6_check3'; Order=6; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_mixed_q8_check3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_laplace_yukawa_separated_q8_check3_local3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_stress_q8_check3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_multi_m2l_stress_q8_check3'; Order=8; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 },
    @{ Case='jcp_laplace_yukawa_separated_q8_check4_local4'; Order=8; Theta=0.6; Up=4.0; Eq=1.05; LocalEq=4.0; LocalCheck=1.05 },
    @{ Case='jcp_laplace_yukawa_separated_q10_check3_local3'; Order=10; Theta=0.6; Up=3.0; Eq=1.05; LocalEq=3.0; LocalCheck=1.05 }
  )
}

$rows = New-Object System.Collections.Generic.List[object]
$stageRows = New-Object System.Collections.Generic.List[object]
$k4Rows = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
  $run = Run-Case $case.Case $case.Order $case.Theta $case.Up $case.Eq $case.LocalEq $case.LocalCheck
  $rows.Add($run.Row)
  foreach ($stage in $run.Stages) { $stageRows.Add($stage) }
  foreach ($k4 in $run.K4) { $k4Rows.Add($k4) }
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e254_matrix.csv'
$StageCsv = Join-Path $WorkRoot 'fmm_e254_stage_attribution.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e254_channel.csv'
$K4Csv = Join-Path $WorkRoot 'fmm_e254_k4_decomposition.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e254_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e254_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E254_JCP_ACCURACY_REPAIR_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$stageRows | Export-Csv -LiteralPath $StageCsv -NoTypeInformation
$rows | Select-Object row_id,aggregate_gate,channel_gate,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,max_channel_rms,max_channel_name,strict_accuracy_status,accuracy_failure | Export-Csv -LiteralPath $ChannelCsv -NoTypeInformation
$k4Rows | Export-Csv -LiteralPath $K4Csv -NoTypeInformation

$badExecution = @($rows | Where-Object { $_.execution_gate_status -ne 'pass' -or $_.finite_status -ne 'pass' })
$badRequiredExecution = @($badExecution | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badOptionalExecution = @($badExecution | Where-Object { [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badCoverage = @($rows | Where-Object { $_.coverage_status -ne 'pass' })
$badRequiredCoverage = @($badCoverage | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badStructure = @($rows | Where-Object { $_.structural_status -ne 'pass' })
$badRequiredStructure = @($badStructure | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badAccuracy = @($rows | Where-Object { $_.strict_accuracy_status -ne 'pass' })
$badRequiredAccuracy = @($badAccuracy | Where-Object { -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badStage = @($stageRows | Where-Object { $_.stage_status -eq 'fail' -and $_.stage_applicable -ne 'false' })
$q8Rows = @($rows | Where-Object { [int]$_.order -eq 8 -and -not [System.Convert]::ToBoolean($_.optional_order_sweep) })
$badQ8Required = @($q8Rows | Where-Object { $_.strict_accuracy_status -ne 'pass' -or $_.structural_status -ne 'pass' -or $_.execution_gate_status -ne 'pass' })

$dominantStage = 'none'
if ($badStage.Count -gt 0) {
  $dominantStage = (@($badStage | Group-Object stage | Sort-Object Count -Descending | Select-Object -First 1).Name)
}
$dominantChannel = 'none'
if ($badRequiredAccuracy.Count -gt 0) {
  $dominantChannel = (@($badRequiredAccuracy | Group-Object max_channel_name | Sort-Object Count -Descending | Select-Object -First 1).Name)
}

$decision = if ($badRequiredExecution.Count -gt 0) {
  'fmm_e254_high_order_execution_gap'
} elseif ($badRequiredCoverage.Count -gt 0) {
  'fmm_e254_structural_gap'
} elseif ($badRequiredStructure.Count -gt 0) {
  'fmm_e254_structural_gap'
} elseif ($badQ8Required.Count -eq 0 -and $badOptionalExecution.Count -eq 0) {
  'fmm_e254_jcp_q8_strict_accuracy_pass'
} elseif ($badQ8Required.Count -eq 0 -and $badOptionalExecution.Count -gt 0) {
  'fmm_e254_synthetic_k1k4_q8_strict_accuracy_pass_optional_high_order_gap'
} elseif ($dominantStage -eq 'p2m_flat_field') {
  'fmm_e254_jcp_flat_equivalent_gap'
} elseif ($dominantChannel -eq 'K4') {
  'fmm_e254_jcp_k4_mixed_derivative_gap'
} elseif ($badStage.Count -gt 0) {
  'fmm_e254_jcp_m2m_m2l_l2l_translation_gap'
} else {
  'fmm_e254_jcp_m2m_m2l_l2l_translation_gap'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  required_q8_rows = $q8Rows.Count
  q8_gap_rows = $badQ8Required.Count
  execution_gap_rows = $badExecution.Count
  required_execution_gap_rows = $badRequiredExecution.Count
  optional_execution_gap_rows = $badOptionalExecution.Count
  coverage_gap_rows = $badCoverage.Count
  structural_gap_rows = $badStructure.Count
  strict_accuracy_gap_rows = $badAccuracy.Count
  required_strict_accuracy_gap_rows = $badRequiredAccuracy.Count
  stage_gap_rows = $badStage.Count
  dominant_stage = $dominantStage
  dominant_channel = $dominantChannel
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "stage_attribution=$StageCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.4 JCP Accuracy Repair Report',
  '',
  "Decision: $decision",
  '',
  'E25.4 keeps E25.3 as a structural checkpoint and applies the strict synthetic K1-K4 gates before any real-panel adapter work.',
  'The promoted repair candidate is q8/check3/local3. Optional q8/check4 and q10 rows are execution diagnostics only.',
  'E25.4 is a synthetic K1-K4 full-hierarchy accuracy closure; it is not a real-panel, full-panel, or GMRES readiness result.',
  ''
)
foreach ($row in $rows) {
  $report += ('  {0}: exec={1}, strict={2}, rms={3}, max={4}:{5}, m2m={6}, m2l={7}, l2l={8}, structural={9}' -f
    $row.row_id,$row.execution_gate_status,$row.strict_accuracy_status,$row.aggregate_rms,$row.max_channel_name,$row.max_channel_rms,$row.m2m_translations,$row.m2l_translations,$row.l2l_translations,$row.structural_status)
}
if ($badStage.Count -gt 0) {
  $report += ''
  $report += 'Dominant failing stage proxies:'
  foreach ($group in ($badStage | Group-Object stage | Sort-Object Count -Descending)) {
    $report += ('  {0}: {1} rows' -f $group.Name,$group.Count)
  }
}
if ($badOptionalExecution.Count -gt 0) {
  $report += ''
  $report += 'Optional high-order execution gaps:'
  foreach ($row in $badOptionalExecution) {
    $report += ('  {0}: exit_code={1}, status={2}' -f $row.row_id,$row.process_exit_code,$row.status)
  }
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.4 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e254_jcp_q8_strict_accuracy_pass','fmm_e254_synthetic_k1k4_q8_strict_accuracy_pass_optional_high_order_gap') -and -not $AllowFailedDecision) {
  throw "FMM-E25.4 gate did not pass: $decision"
}
