param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e38_tree_worklist_compression'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e38_matrix.csv'
$TreeCsv = Join-Path $WorkRoot 'fmm_e38_tree.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e38_runtime.csv'
$ComparisonCsv = Join-Path $WorkRoot 'fmm_e38_comparison.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e38_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e38_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E38_TREE_WORKLIST_COMPRESSION_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E38 smoke.'
}

function Assert-SmokeSupportsE38 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $smoke = Get-Content -Raw -LiteralPath $SmokeSource
  $headers = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy.h')
  $realHeaders = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.h')
  foreach ($needle in @(
      'tree-stop-policy',
      'tree_stop_policy=',
      'tree_leaf_count=',
      'tree_max_depth_observed=',
      'tree_radius_stop_leaf_count=',
      'tree_occupancy_stop_leaf_count=')) {
    if (-not $smoke.Contains($needle) -and -not $headers.Contains($needle) -and -not $realHeaders.Contains($needle)) {
      throw "E38 tree/worklist compression support is not present yet: missing '$needle'"
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
  Assert-SmokeSupportsE38
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

function New-Case(
  [string]$Id,
  [string]$Fixture,
  [string]$Sample,
  [int]$SampleCount,
  [int]$SurfacePoints,
  [int]$PanelLimit,
  [string]$TreeStopPolicy,
  [double]$PanelRadiusFactor,
  [int]$OccupancyStop
) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    tree_stop_policy = $TreeStopPolicy
    panel_radius_factor = $PanelRadiusFactor
    occupancy_stop = $OccupancyStop
    row_id = "${Id}_${TreeStopPolicy}_r${PanelRadiusFactor}_occ${OccupancyStop}".Replace('.', 'p')
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
  Ensure-Field $Row 'requested_tree_stop_policy' $Case.tree_stop_policy
  Ensure-Field $Row 'requested_panel_radius_factor' $Case.panel_radius_factor
  Ensure-Field $Row 'requested_occupancy_stop' $Case.occupancy_stop
  foreach ($field in @(
      'tree_stop_policy',
      'tree_leaf_count',
      'tree_node_count',
      'tree_max_depth_observed',
      'tree_radius_stop_leaf_count',
      'tree_occupancy_stop_leaf_count',
      'tree_leaf_size_stop_count',
      'tree_max_depth_stop_count',
      'tree_min_leaf_half_width_over_panel_radius',
      'tree_max_leaf_panel_radius',
      'tree_mean_leaf_source_count',
      'tree_mean_leaf_target_count',
      'm2l_eval_policy',
      'm2l_seconds',
      'm2l_translations',
      'm2l_check_points_total',
      'm2l_kernel_evals_estimated',
      'aggregate_rms',
      'max_channel_rms',
      'k1_rms',
      'k2_rms',
      'k3_rms',
      'k4_rms',
      'far_fraction',
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
    "--mode=full-panel-profile",
    "--emit-e36-stage-attribution=1",
    "--hier-debug-mode=dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian",
    "--m2l-eval-policy=target-check-batched",
    "--tree-stop-policy=$($Case.tree_stop_policy)",
    "--tree-panel-radius-factor=$($Case.panel_radius_factor)",
    "--tree-occupancy-stop=$($Case.occupancy_stop)",
    "--order=8",
    "--theta=0.6",
    "--equiv-scale=1.05",
    "--upward-check-scale=3.0",
    "--local-equiv-scale=3.0",
    "--local-check-scale=1.05",
    "--condition-policy=qr-rank-only",
    "--real-panel-guard=strict-local-check",
    "--guard-margin-factor=0"
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
  $line = ''
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

function Decide([array]$Rows, [array]$Comparisons) {
  $required = @($Rows | Where-Object { $_.case_id -eq '1SH1_s4_p512_uniform64' })
  if (@($required | Where-Object { (Get-Field $_ 'execution_status') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e38_runtime_limited'
  }
  if (@($Rows | Where-Object { -not (Is-Accurate $_) -or -not (Is-Covered $_) }).Count -gt 0) {
    return 'fmm_e38_accuracy_or_coverage_gap'
  }
  $best = @($Comparisons | Where-Object { $_.case_id -eq '1SH1_s4_p512_uniform64' -and $_.candidate_policy -ne 'baseline' } | Sort-Object {[double]$_.m2l_check_point_reduction} -Descending | Select-Object -First 1)
  if ($best.Count -eq 1) {
    $m2lReduction = [double]$best[0].m2l_check_point_reduction
    $timeSpeedup = [double]$best[0].total_speedup_vs_baseline
    $farFraction = [double]$best[0].candidate_far_fraction
    if ($m2lReduction -ge 0.25 -and $timeSpeedup -ge 1.15 -and $farFraction -ge 0.10) {
      return 'fmm_e38_tree_worklist_compression_pass'
    }
    if ($m2lReduction -ge 0.25 -and $farFraction -lt 0.10) {
      return 'fmm_e38_fallback_dominated'
    }
  }
  'fmm_e38_accuracy_pass_no_compression_gain'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$caseSpecs = @(
  @('1SH1_s4_p128_uniform32', '1SH1', 'uniform_hash', 32, 4, 128),
  @('1SH1_s4_p256_uniform32', '1SH1', 'uniform_hash', 32, 4, 256),
  @('1SH1_s4_p512_uniform64', '1SH1', 'uniform_hash', 64, 4, 512)
)
if (-not $Quick) {
  $caseSpecs += @(
    @('1SH1_s8_full_uniform128', '1SH1', 'uniform_hash', 128, 8, 0),
    @('1SH1_s8_full_closest64', '1SH1', 'closest_opposing_surface', 64, 8, 0),
    @('1SH1_s8_full_high_aspect64', '1SH1', 'high_aspect', 64, 8, 0)
  )
}

$cases = @()
foreach ($spec in $caseSpecs) {
  $cases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] 'baseline' 0.0 0
  $cases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] 'panel-radius' 1.5 0
  $cases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] 'panel-radius-occupancy' 1.5 16
}

$rows = @()
foreach ($case in $cases) {
  $existing = @()
  if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
    $existing = @(Import-Csv -LiteralPath $MatrixCsv | Where-Object { $_.row_id -eq $case.row_id } | Select-Object -First 1)
  }
  if ($existing.Count -eq 1 -and -not ($RerunFailed -and (Get-Field $existing[0] 'execution_status') -ne 'pass')) {
    Write-Host "resume: $($case.row_id)"
    $rows += $existing[0]
  } else {
    $rows += Run-Case $case
    Export-Table $rows $MatrixCsv
  }
}

$treeRows = @()
$runtimeRows = @()
foreach ($row in $rows) {
  $treeRows += [pscustomobject][ordered]@{
    row_id = Get-Field $row 'row_id'
    case_id = Get-Field $row 'case_id'
    tree_stop_policy = Get-Field $row 'tree_stop_policy'
    tree_node_count = Get-Field $row 'tree_node_count'
    tree_leaf_count = Get-Field $row 'tree_leaf_count'
    tree_max_depth_observed = Get-Field $row 'tree_max_depth_observed'
    tree_radius_stop_leaf_count = Get-Field $row 'tree_radius_stop_leaf_count'
    tree_occupancy_stop_leaf_count = Get-Field $row 'tree_occupancy_stop_leaf_count'
    tree_leaf_size_stop_count = Get-Field $row 'tree_leaf_size_stop_count'
    tree_max_depth_stop_count = Get-Field $row 'tree_max_depth_stop_count'
    tree_min_leaf_half_width_over_panel_radius = Get-Field $row 'tree_min_leaf_half_width_over_panel_radius'
    tree_mean_leaf_source_count = Get-Field $row 'tree_mean_leaf_source_count'
    tree_mean_leaf_target_count = Get-Field $row 'tree_mean_leaf_target_count'
    far_fraction = Get-Field $row 'far_fraction'
  }
  $runtimeRows += [pscustomobject][ordered]@{
    row_id = Get-Field $row 'row_id'
    case_id = Get-Field $row 'case_id'
    tree_stop_policy = Get-Field $row 'tree_stop_policy'
    m2l_seconds = Get-Field $row 'm2l_seconds'
    m2l_translations = Get-Field $row 'm2l_translations'
    m2l_check_points_total = Get-Field $row 'm2l_check_points_total'
    m2l_kernel_evals_estimated = Get-Field $row 'm2l_kernel_evals_estimated'
    total_seconds = Get-Field $row 'total_seconds'
    aggregate_rms = Get-Field $row 'aggregate_rms'
    max_channel_rms = Get-Field $row 'max_channel_rms'
  }
}

$comparisons = @()
foreach ($caseId in @($rows | Select-Object -ExpandProperty case_id -Unique)) {
  $base = @($rows | Where-Object { $_.case_id -eq $caseId -and (Get-Field $_ 'tree_stop_policy') -eq 'baseline' } | Select-Object -First 1)
  if ($base.Count -ne 1) { continue }
  foreach ($candidate in @($rows | Where-Object { $_.case_id -eq $caseId -and (Get-Field $_ 'tree_stop_policy') -ne 'baseline' })) {
    $baseCheck = [double](As-Double $base[0] 'm2l_check_points_total' 0.0)
    $candCheck = [double](As-Double $candidate 'm2l_check_points_total' 0.0)
    $baseTime = [double](As-Double $base[0] 'total_seconds' 0.0)
    $candTime = [double](As-Double $candidate 'total_seconds' 0.0)
    $comparisons += [pscustomobject][ordered]@{
      case_id = $caseId
      baseline_row = Get-Field $base[0] 'row_id'
      candidate_row = Get-Field $candidate 'row_id'
      candidate_policy = Get-Field $candidate 'tree_stop_policy'
      m2l_check_point_reduction = if ($baseCheck -gt 0.0) { ($baseCheck - $candCheck) / $baseCheck } else { 0.0 }
      m2l_translation_reduction = if ((As-Double $base[0] 'm2l_translations' 0.0) -gt 0.0) { ((As-Double $base[0] 'm2l_translations' 0.0) - (As-Double $candidate 'm2l_translations' 0.0)) / (As-Double $base[0] 'm2l_translations' 1.0) } else { 0.0 }
      total_speedup_vs_baseline = if ($candTime -gt 0.0) { $baseTime / $candTime } else { 0.0 }
      candidate_far_fraction = As-Double $candidate 'far_fraction' 0.0
      candidate_accuracy_status = if (Is-Accurate $candidate) { 'pass' } else { 'fail' }
      candidate_coverage_status = if (Is-Covered $candidate) { 'pass' } else { 'fail' }
    }
  }
}

Export-Table $rows $MatrixCsv
Export-Table $treeRows $TreeCsv
Export-Table $runtimeRows $RuntimeCsv
Export-Table $comparisons $ComparisonCsv
$decision = Decide $rows $comparisons
[pscustomobject][ordered]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  generated_at = (Get-Date).ToString('s')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

@"
FMM-E38 tree/worklist compression gate
decision=$decision
matrix=$MatrixCsv
tree=$TreeCsv
runtime=$RuntimeCsv
comparison=$ComparisonCsv
"@ | Set-Content -LiteralPath $SummaryDat

@"
TITAN FMM-E38 TREE/WORKLIST COMPRESSION REPORT

Decision: $decision

This diagnostic gate compares the E37 target-check-batched M2L baseline
against opt-in panel-aware leaf-stop policies. The gate is intentionally scoped
to tree/worklist compression: KIFMM math, q8 policy, strict-local-check guard,
sampled dense production-direct audit, and production defaults are unchanged.

Outputs:
- Matrix: $MatrixCsv
- Tree telemetry: $TreeCsv
- Runtime telemetry: $RuntimeCsv
- Baseline/candidate comparison: $ComparisonCsv
- Decision: $DecisionCsv
"@ | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E38 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
if ($decision -notin @('fmm_e38_tree_worklist_compression_pass', 'fmm_e38_accuracy_pass_no_compression_gain') -and -not $AllowFailedDecision) {
  throw "FMM-E38 gate did not pass: $decision"
}
