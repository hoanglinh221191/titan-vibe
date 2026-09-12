param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e39_tree_stop_sweep'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e39_matrix.csv'
$TreeCsv = Join-Path $WorkRoot 'fmm_e39_tree.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e39_runtime.csv'
$SweepSummaryCsv = Join-Path $WorkRoot 'fmm_e39_sweep_summary.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e39_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e39_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E39_TREE_STOP_SWEEP_REPORT.txt'

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) { throw "Requested C++ compiler not found: $Requested" }
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
  throw 'No C++ compiler found for FMM-E39 smoke.'
}

function Assert-SmokeSupportsE39 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $smoke = Get-Content -Raw -LiteralPath $SmokeSource
  $headers = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy.h')
  $realHeaders = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.h')
  foreach ($needle in @(
      'tree-stop-policy',
      'tree-panel-radius-factor',
      'tree-occupancy-stop',
      'm2l-eval-policy',
      'tree_stop_policy=',
      'tree_leaf_count=',
      'tree_max_depth_observed=',
      'tree_radius_stop_leaf_count=',
      'tree_occupancy_stop_leaf_count=',
      'm2l_check_points_total=')) {
    if (-not $smoke.Contains($needle) -and -not $headers.Contains($needle) -and -not $realHeaders.Contains($needle)) {
      throw "E38/E39 tree-stop support is not present yet: missing '$needle'"
    }
  }
}

function Parse-KeyValueLine([string]$Line, [string]$Prefix) {
  if (-not $Line.StartsWith($Prefix)) { throw "Unexpected line for ${Prefix}: $Line" }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring($Prefix.Length).Trim().Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function Get-Field([object]$Row, [string]$Name, [string]$Default = 'NA') {
  if ($null -eq $Row) { return $Default }
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop -or $null -eq $prop.Value -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) {
    return $Default
  }
  [string]$prop.Value
}

function As-Double([object]$Row, [string]$Name, [double]$Default = [double]::NaN) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value) -or $value -eq 'nan') { return $Default }
  try { [double]$value } catch { $Default }
}

function As-Int64([object]$Row, [string]$Name, [int64]$Default = [int64]::MinValue) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  try { [int64]$value } catch { $Default }
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) {
    Set-Content -LiteralPath $Path -Value ''
  } else {
    $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE39
  New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
  $Compiler = Resolve-Cxx $Cxx
  $CompileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.cpp',
    'src\tabi\tabi_fmm_kifmm_cartesian_bridge.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp',
    'src\tabi\tabi_fmm_cartesian_sampled_panel_matvec.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) { throw 'failed to build titan_fmm_kifmm_hierarchy_real_panel_smoke' }
  } finally {
    Pop-Location
  }
}

function New-Candidate([string]$Family, [string]$Policy, [double]$Factor, [int]$Occ) {
  $factorText = ([string]$Factor).Replace('.', 'p')
  [pscustomobject]@{
    family = $Family
    policy = $Policy
    factor = $Factor
    occupancy = $Occ
    key = "${Policy}_r${factorText}_occ${Occ}"
  }
}

function New-Case(
  [string]$Id,
  [string]$Fixture,
  [string]$Sample,
  [int]$SampleCount,
  [int]$SurfacePoints,
  [int]$PanelLimit,
  [object]$Candidate
) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    candidate_family = $Candidate.family
    candidate_key = $Candidate.key
    tree_stop_policy = $Candidate.policy
    panel_radius_factor = $Candidate.factor
    occupancy_stop = $Candidate.occupancy
    row_id = "${Id}_$($Candidate.key)"
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.row_id
  Ensure-Field $Row 'case_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'surface_points' $Case.surface_points
  Ensure-Field $Row 'panel_limit' $Case.panel_limit
  Ensure-Field $Row 'candidate_family' $Case.candidate_family
  Ensure-Field $Row 'candidate_key' $Case.candidate_key
  Ensure-Field $Row 'requested_tree_stop_policy' $Case.tree_stop_policy
  Ensure-Field $Row 'requested_panel_radius_factor' $Case.panel_radius_factor
  Ensure-Field $Row 'requested_occupancy_stop' $Case.occupancy_stop
  foreach ($field in @(
      'tree_stop_policy',
      'tree_panel_radius_factor',
      'tree_occupancy_stop',
      'tree_leaf_count',
      'tree_node_count',
      'tree_max_depth_observed',
      'tree_radius_stop_leaf_count',
      'tree_occupancy_stop_leaf_count',
      'tree_leaf_size_stop_count',
      'tree_max_depth_stop_count',
      'tree_unsplittable_stop_count',
      'tree_min_leaf_half_width_over_panel_radius',
      'tree_max_leaf_panel_radius',
      'tree_mean_leaf_source_count',
      'tree_mean_leaf_target_count',
      'm2l_eval_policy',
      'm2l_seconds',
      'm2l_translations',
      'm2l_check_points_total',
      'm2l_equivalent_points',
      'm2l_kernel_evals_estimated',
      'm2l_seconds_per_translation',
      'm2l_seconds_per_check_point',
      'aggregate_rms',
      'max_channel_rms',
      'k1_rms',
      'k2_rms',
      'k3_rms',
      'k4_rms',
      'far_fraction',
      'coverage',
      'full_coverage_expected',
      'production_invalid_pairs',
      'nonproduction_fallback_pairs',
      'far_panel_pairs',
      'direct_panel_pairs',
      'total_seconds',
      'output_hash',
      'policy_hash',
      'worklist_hash')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $Row
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.row_id) (timeout ${CaseTimeoutSeconds}s)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--surface-points=$($Case.surface_points)",
    "--panel-limit=$($Case.panel_limit)",
    '--mode=full-panel-profile',
    '--emit-e36-stage-attribution=1',
    '--hier-debug-mode=dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian',
    '--m2l-eval-policy=target-check-batched',
    "--tree-stop-policy=$($Case.tree_stop_policy)",
    "--tree-panel-radius-factor=$($Case.panel_radius_factor)",
    "--tree-occupancy-stop=$($Case.occupancy_stop)",
    '--order=8',
    '--theta=0.6',
    '--equiv-scale=1.05',
    '--upward-check-scale=3.0',
    '--local-equiv-scale=3.0',
    '--local-check-scale=1.05',
    '--condition-policy=qr-rank-only',
    '--real-panel-guard=strict-local-check',
    '--guard-margin-factor=0',
    "--pqr-1sh1=$Pqr1sh1"
  )
  $stdout = Join-Path $WorkRoot "$($Case.row_id).stdout.txt"
  $stderr = Join-Path $WorkRoot "$($Case.row_id).stderr.txt"
  Remove-Item -LiteralPath $stdout, $stderr -Force -ErrorAction SilentlyContinue
  $proc = Start-Process -FilePath $ExePath -ArgumentList $args -NoNewWindow -PassThru `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr -WorkingDirectory $RepoRoot
  $completed = $proc.WaitForExit([Math]::Max(1, $CaseTimeoutSeconds) * 1000)
  if (-not $completed) {
    Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
    $proc.WaitForExit()
    $row = [pscustomobject][ordered]@{
      row_id = $Case.row_id
      case_id = $Case.id
      status = 'timeout'
      execution_status = 'timeout'
      stdout = $stdout
      stderr = $stderr
    }
    return Complete-Row $row $Case
  }
  $output = @()
  if (Test-Path -LiteralPath $stdout) { $output += Get-Content -LiteralPath $stdout }
  if (Test-Path -LiteralPath $stderr) { $output += Get-Content -LiteralPath $stderr }
  $line = @($output | Where-Object { $_.StartsWith('E255_RESULT ') } | Select-Object -Last 1)
  if ([string]::IsNullOrWhiteSpace($line)) {
    $row = [pscustomObject][ordered]@{
      row_id = $Case.row_id
      case_id = $Case.id
      status = 'process_failed'
      execution_status = 'process_failed'
      stdout = $stdout
      stderr = $stderr
      process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }
    return Complete-Row $row $Case
  }
  $parsed = Parse-KeyValueLine $line 'E255_RESULT '
  $parsed | Add-Member -NotePropertyName stdout -NotePropertyValue $stdout
  $parsed | Add-Member -NotePropertyName stderr -NotePropertyValue $stderr
  Complete-Row $parsed $Case
}

function Is-Accurate([object]$Row) {
  (Get-Field $Row 'execution_status') -eq 'pass' -and
  (As-Double $Row 'aggregate_rms' 1.0) -le 1.0e-6 -and
  (As-Double $Row 'max_channel_rms' 1.0) -le 5.0e-6
}

function Is-Covered([object]$Row) {
  (As-Int64 $Row 'production_invalid_pairs' 1) -eq 0 -and
  (As-Int64 $Row 'nonproduction_fallback_pairs' 1) -eq 0 -and
  (As-Int64 $Row 'coverage' -1) -eq (As-Int64 $Row 'full_coverage_expected' -2)
}

function Is-GatedPass([object]$Row) {
  (Is-Accurate $Row) -and
  (Is-Covered $Row) -and
  (As-Double $Row 'far_fraction' 0.0) -ge 0.10 -and
  (As-Int64 $Row 'far_panel_pairs' 0) -gt 0
}

function Case-Row([array]$Rows, [string]$CaseId, [string]$CandidateKey) {
  return ,@($Rows | Where-Object { $_.case_id -eq $CaseId -and (Get-Field $_ 'candidate_key') -eq $CandidateKey } | Select-Object -First 1)
}

function Summarize-Candidate([array]$Rows, [array]$CaseIds, [object]$Candidate, [string]$Scope) {
  $caseIdList = @()
  if ($null -ne $CaseIds) { $caseIdList = @($CaseIds) }
  $requiredCount = $caseIdList.Count
  $passCount = 0
  $rowCount = 0
  $baselinePassCount = 0
  $candidateSeconds = 0.0
  $baselineSeconds = 0.0
  $candidateM2l = 0.0
  $baselineM2l = 0.0
  $candidateTranslations = 0.0
  $baselineTranslations = 0.0
  $minFarFraction = [double]::PositiveInfinity
  foreach ($caseId in $caseIdList) {
    $candidateRow = Case-Row $Rows $caseId $Candidate.key
    $baselineRow = Case-Row $Rows $caseId 'baseline_r0_occ0'
    if ($candidateRow.Count -ne 1 -or $baselineRow.Count -ne 1) { continue }
    $rowCount += 1
    if (Is-GatedPass $candidateRow[0]) { $passCount += 1 }
    if (Is-GatedPass $baselineRow[0]) { $baselinePassCount += 1 }
    $candidateSeconds += [double](As-Double $candidateRow[0] 'total_seconds' 0.0)
    $baselineSeconds += [double](As-Double $baselineRow[0] 'total_seconds' 0.0)
    $candidateM2l += [double](As-Double $candidateRow[0] 'm2l_check_points_total' 0.0)
    $baselineM2l += [double](As-Double $baselineRow[0] 'm2l_check_points_total' 0.0)
    $candidateTranslations += [double](As-Double $candidateRow[0] 'm2l_translations' 0.0)
    $baselineTranslations += [double](As-Double $baselineRow[0] 'm2l_translations' 0.0)
    $minFarFraction = [Math]::Min($minFarFraction, [double](As-Double $candidateRow[0] 'far_fraction' 0.0))
  }
  [pscustomobject][ordered]@{
    scope = $Scope
    candidate_key = $Candidate.key
    candidate_family = $Candidate.family
    tree_stop_policy = $Candidate.policy
    panel_radius_factor = $Candidate.factor
    occupancy_stop = $Candidate.occupancy
    rows_seen = $rowCount
    rows_required = $requiredCount
    pass_count = $passCount
    baseline_pass_count = $baselinePassCount
    gated_status = if ($rowCount -eq $requiredCount -and $passCount -eq $requiredCount) { 'pass' } else { 'fail' }
    baseline_status = if ($rowCount -eq $requiredCount -and $baselinePassCount -eq $requiredCount) { 'pass' } else { 'fail' }
    total_seconds_sum = $candidateSeconds
    baseline_seconds_sum = $baselineSeconds
    total_speedup_vs_baseline = if ($candidateSeconds -gt 0.0) { $baselineSeconds / $candidateSeconds } else { 0.0 }
    m2l_check_point_reduction = if ($baselineM2l -gt 0.0) { ($baselineM2l - $candidateM2l) / $baselineM2l } else { 0.0 }
    m2l_translation_reduction = if ($baselineTranslations -gt 0.0) { ($baselineTranslations - $candidateTranslations) / $baselineTranslations } else { 0.0 }
    min_far_fraction = if ([double]::IsPositiveInfinity($minFarFraction)) { 0.0 } else { $minFarFraction }
  }
}

function Best-Summary([array]$Summaries, [string]$Family) {
  return ,@($Summaries | Where-Object {
      $_.candidate_family -eq $Family -and $_.gated_status -eq 'pass'
    } | Sort-Object {[double]$_.total_seconds_sum}, candidate_key | Select-Object -First 1)
}

function Decide-E39([array]$Rows, [array]$QuickCaseIds, [array]$FullCaseIds, [array]$SweepSummary, [bool]$IsQuick) {
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e39_runtime_limited'
  }
  $baselineQuick = @($SweepSummary | Where-Object { $_.scope -eq 'quick' -and $_.candidate_family -eq 'baseline' } | Select-Object -First 1)
  if ($baselineQuick.Count -ne 1 -or $baselineQuick[0].gated_status -ne 'pass') {
    return 'fmm_e39_accuracy_gap'
  }
  $quickBest = @($SweepSummary | Where-Object {
      $_.scope -eq 'quick' -and $_.candidate_family -ne 'baseline' -and $_.gated_status -eq 'pass'
    } | Sort-Object {[double]$_.total_seconds_sum}, candidate_key | Select-Object -First 1)
  if ($IsQuick) {
    if ($quickBest.Count -eq 1) { return 'fmm_e39_quick_pass_full_pending' }
    $fallbackOnly = @($SweepSummary | Where-Object {
        $_.scope -eq 'quick' -and $_.candidate_family -ne 'baseline' -and
        [double]$_.min_far_fraction -lt 0.10 -and [int]$_.pass_count -gt 0
      })
    if ($fallbackOnly.Count -gt 0) { return 'fmm_e39_fallback_dominated' }
    return 'fmm_e39_no_speed_gain'
  }
  $baselineFull = @($SweepSummary | Where-Object { $_.scope -eq 'full_all' -and $_.candidate_family -eq 'baseline' } | Select-Object -First 1)
  if ($baselineFull.Count -ne 1 -or $baselineFull[0].gated_status -ne 'pass') {
    return 'fmm_e39_accuracy_gap'
  }
  $fullUniform = @($SweepSummary | Where-Object {
      $_.scope -eq 'full_uniform' -and $_.candidate_family -ne 'baseline' -and $_.gated_status -eq 'pass'
    } | Sort-Object {[double]$_.total_seconds_sum}, candidate_key | Select-Object -First 1)
  if ($fullUniform.Count -eq 1 -and [double]$fullUniform[0].total_speedup_vs_baseline -ge 1.5) {
    return 'fmm_e39_tree_stop_policy_pass'
  }
  $fallbackFull = @($SweepSummary | Where-Object {
      $_.scope -eq 'full_all' -and $_.candidate_family -ne 'baseline' -and
      [double]$_.min_far_fraction -lt 0.10 -and [int]$_.pass_count -gt 0
    })
  if ($fallbackFull.Count -gt 0) { return 'fmm_e39_fallback_dominated' }
  'fmm_e39_no_speed_gain'
}

function Update-DerivedCsvs([array]$Rows, [array]$Candidates, [array]$QuickCaseIds, [array]$FullCaseIds) {
  $treeRows = @()
  $runtimeRows = @()
  foreach ($row in $Rows) {
    $treeRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      candidate_key = Get-Field $row 'candidate_key'
      candidate_family = Get-Field $row 'candidate_family'
      tree_stop_policy = Get-Field $row 'tree_stop_policy'
      tree_node_count = Get-Field $row 'tree_node_count'
      tree_leaf_count = Get-Field $row 'tree_leaf_count'
      tree_max_depth_observed = Get-Field $row 'tree_max_depth_observed'
      tree_radius_stop_leaf_count = Get-Field $row 'tree_radius_stop_leaf_count'
      tree_occupancy_stop_leaf_count = Get-Field $row 'tree_occupancy_stop_leaf_count'
      tree_leaf_size_stop_count = Get-Field $row 'tree_leaf_size_stop_count'
      tree_max_depth_stop_count = Get-Field $row 'tree_max_depth_stop_count'
      tree_min_leaf_half_width_over_panel_radius = Get-Field $row 'tree_min_leaf_half_width_over_panel_radius'
      tree_max_leaf_panel_radius = Get-Field $row 'tree_max_leaf_panel_radius'
      tree_mean_leaf_source_count = Get-Field $row 'tree_mean_leaf_source_count'
      tree_mean_leaf_target_count = Get-Field $row 'tree_mean_leaf_target_count'
      far_fraction = Get-Field $row 'far_fraction'
    }
    $runtimeRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      candidate_key = Get-Field $row 'candidate_key'
      candidate_family = Get-Field $row 'candidate_family'
      tree_stop_policy = Get-Field $row 'tree_stop_policy'
      m2l_eval_policy = Get-Field $row 'm2l_eval_policy'
      m2l_seconds = Get-Field $row 'm2l_seconds'
      m2l_translations = Get-Field $row 'm2l_translations'
      m2l_check_points_total = Get-Field $row 'm2l_check_points_total'
      m2l_kernel_evals_estimated = Get-Field $row 'm2l_kernel_evals_estimated'
      total_seconds = Get-Field $row 'total_seconds'
      far_fraction = Get-Field $row 'far_fraction'
      aggregate_rms = Get-Field $row 'aggregate_rms'
      max_channel_rms = Get-Field $row 'max_channel_rms'
    }
  }
  $summaries = @()
  foreach ($candidate in $Candidates) {
    $summaries += Summarize-Candidate $Rows $QuickCaseIds $candidate 'quick'
    if ($null -ne $FullCaseIds -and @($FullCaseIds).Count -gt 0) {
      $summaries += Summarize-Candidate $Rows $FullCaseIds $candidate 'full_all'
      $summaries += Summarize-Candidate $Rows @('1SH1_s8_full_uniform128') $candidate 'full_uniform'
    }
  }
  Export-Table $Rows $MatrixCsv
  Export-Table $treeRows $TreeCsv
  Export-Table $runtimeRows $RuntimeCsv
  Export-Table $summaries $SweepSummaryCsv
  $summaries
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
if ($Force) {
  Remove-Item -LiteralPath $MatrixCsv, $TreeCsv, $RuntimeCsv, $SweepSummaryCsv, $DecisionCsv, $SummaryDat, $ReportPath -Force -ErrorAction SilentlyContinue
}
Build-Smoke

$baseline = New-Candidate 'baseline' 'baseline' 0.0 0
$radiusOnly = @(1.2, 1.5, 2.0 | ForEach-Object { New-Candidate 'radius_only' 'panel-radius' $_ 0 })
$radiusOccupancy = @()
foreach ($factor in @(1.2, 1.5, 2.0)) {
  foreach ($occ in @(8, 16, 32)) {
    $radiusOccupancy += New-Candidate 'radius_occupancy' 'panel-radius-occupancy' $factor $occ
  }
}
$allCandidates = @($baseline) + $radiusOnly + $radiusOccupancy

$quickSpecs = @(
  @('1SH1_s4_p512_uniform64', '1SH1', 'uniform_hash', 64, 4, 512),
  @('1SH1_s8_p512_uniform64', '1SH1', 'uniform_hash', 64, 8, 512)
)
$fullSpecs = @(
  @('1SH1_s8_full_uniform128', '1SH1', 'uniform_hash', 128, 8, 0),
  @('1SH1_s8_full_closest64', '1SH1', 'closest_opposing_surface', 64, 8, 0),
  @('1SH1_s8_full_high_aspect64', '1SH1', 'high_aspect', 64, 8, 0)
)
$quickCaseIds = @($quickSpecs | ForEach-Object { $_[0] })
$fullCaseIds = if ($Quick) { @() } else { @($fullSpecs | ForEach-Object { $_[0] }) }

$rows = @()
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  $rows = @(Import-Csv -LiteralPath $MatrixCsv)
}

function Invoke-Cases([array]$Cases) {
  foreach ($case in $Cases) {
    $existing = @($rows | Where-Object { $_.row_id -eq $case.row_id } | Select-Object -First 1)
    $shouldRun = $true
    if ($existing.Count -eq 1 -and -not $Force) {
      $shouldRun = $false
      if ($RerunFailed -and (Get-Field $existing[0] 'execution_status') -ne 'pass') {
        $shouldRun = $true
      }
    }
    if ($shouldRun) {
      $newRow = Run-Case $case
      $script:rows = @($script:rows | Where-Object { $_.row_id -ne $case.row_id })
      $script:rows += $newRow
      Export-Table $script:rows $MatrixCsv
    } else {
      Write-Host "resume: $($case.row_id)"
    }
  }
}

$quickCases = @()
foreach ($spec in $quickSpecs) {
  foreach ($candidate in $allCandidates) {
    $quickCases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] $candidate
  }
}
Invoke-Cases $quickCases

$quickSummaries = Update-DerivedCsvs $rows $allCandidates $quickCaseIds @()
$bestRadius = Best-Summary $quickSummaries 'radius_only'
$bestOcc = Best-Summary $quickSummaries 'radius_occupancy'
$selectedFullCandidates = @($baseline)
if ($bestRadius.Count -eq 1) {
  $selectedFullCandidates += @($allCandidates | Where-Object { $_.key -eq $bestRadius[0].candidate_key } | Select-Object -First 1)
}
if ($bestOcc.Count -eq 1) {
  $selectedFullCandidates += @($allCandidates | Where-Object { $_.key -eq $bestOcc[0].candidate_key } | Select-Object -First 1)
}
$selectedFullCandidates = @($selectedFullCandidates | Sort-Object key -Unique)

if (-not $Quick -and $selectedFullCandidates.Count -gt 1) {
  $fullCases = @()
  foreach ($spec in $fullSpecs) {
    foreach ($candidate in $selectedFullCandidates) {
      $fullCases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] $candidate
    }
  }
  Invoke-Cases $fullCases
}

$sweepSummary = Update-DerivedCsvs $rows $allCandidates $quickCaseIds $fullCaseIds
$decision = Decide-E39 $rows $quickCaseIds $fullCaseIds $sweepSummary ([bool]$Quick)

$overallSelection = 'none'
if ($bestRadius.Count -eq 1 -and $bestOcc.Count -eq 1) {
  $radiusTime = [double]$bestRadius[0].total_seconds_sum
  $occTime = [double]$bestOcc[0].total_seconds_sum
  if ($occTime -lt ($radiusTime * 0.95)) {
    $overallSelection = $bestOcc[0].candidate_key
  } else {
    $overallSelection = $bestRadius[0].candidate_key
  }
} elseif ($bestRadius.Count -eq 1) {
  $overallSelection = $bestRadius[0].candidate_key
} elseif ($bestOcc.Count -eq 1) {
  $overallSelection = $bestOcc[0].candidate_key
}

[pscustomobject][ordered]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  best_radius_only = if ($bestRadius.Count -eq 1) { $bestRadius[0].candidate_key } else { 'none' }
  best_radius_occupancy = if ($bestOcc.Count -eq 1) { $bestOcc[0].candidate_key } else { 'none' }
  selected_overall = $overallSelection
  generated_at = (Get-Date).ToString('s')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

@"
FMM-E39 tree-stop sweep full gate
decision=$decision
mode=$(if ($Quick) { 'quick' } else { 'full' })
best_radius_only=$(if ($bestRadius.Count -eq 1) { $bestRadius[0].candidate_key } else { 'none' })
best_radius_occupancy=$(if ($bestOcc.Count -eq 1) { $bestOcc[0].candidate_key } else { 'none' })
selected_overall=$overallSelection
matrix=$MatrixCsv
tree=$TreeCsv
runtime=$RuntimeCsv
sweep_summary=$SweepSummaryCsv
"@ | Set-Content -LiteralPath $SummaryDat

@"
TITAN FMM-E39 TREE-STOP SWEEP REPORT

Decision: $decision

E39 validates the opt-in E38 panel-aware tree-stop policies on broader 1SH1
rows. The KIFMM math, q8 policy, strict-local-check guard, target-check-batched
M2L, sampled dense production-direct audit, production TABI/grid defaults, and
manuscript scope are unchanged.

Candidate sweep:
- baseline: baseline factor=0 occ=0
- radius-only: panel-radius factor=1.2,1.5,2.0
- radius+occupancy: panel-radius-occupancy factor=1.2,1.5,2.0 and occ=8,16,32

Candidate selection:
- Best quick radius-only: $(if ($bestRadius.Count -eq 1) { $bestRadius[0].candidate_key } else { 'none' })
- Best quick radius+occupancy: $(if ($bestOcc.Count -eq 1) { $bestOcc[0].candidate_key } else { 'none' })
- Selected overall: $overallSelection

Outputs:
- Matrix: $MatrixCsv
- Tree telemetry: $TreeCsv
- Runtime telemetry: $RuntimeCsv
- Sweep summary: $SweepSummaryCsv
- Decision: $DecisionCsv

Decision labels:
- fmm_e39_tree_stop_policy_pass: full rows pass and the selected tree-stop
  policy improves total time by at least 1.5x on the full uniform row.
- fmm_e39_quick_pass_full_pending: quick sweep passes, full validation not run.
- fmm_e39_accuracy_gap, fmm_e39_fallback_dominated,
  fmm_e39_no_speed_gain, fmm_e39_runtime_limited: fail-closed
  classifications.
"@ | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E39 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Sweep summary: $SweepSummaryCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."
if ($decision -notin @('fmm_e39_tree_stop_policy_pass', 'fmm_e39_quick_pass_full_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E39 gate did not pass: $decision"
}
