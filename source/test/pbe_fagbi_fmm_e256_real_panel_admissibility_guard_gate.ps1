param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e256_real_panel_admissibility'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e256_matrix.csv'
$GuardCsv = Join-Path $WorkRoot 'fmm_e256_guard.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e256_channel.csv'
$CoverageCsv = Join-Path $WorkRoot 'fmm_e256_coverage.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e256_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e256_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E256_REAL_PANEL_ADMISSIBILITY_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E25.6 real-panel admissibility smoke.'
}

function Assert-SmokeSupportsE256 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) {
    throw "missing smoke source: $SmokeSource"
  }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'real-panel-guard',
      'guard-margin-factor',
      'negative_target_check_excess_pairs')) {
    if (-not $text.Contains($needle)) {
      throw "E25.6 smoke support is not present yet: missing '$needle' in $SmokeSource"
    }
  }
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

function Get-Field([object]$Row, [string]$Name, [string]$Default = 'NA') {
  if ($null -eq $Row) { return $Default }
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop -or $null -eq $prop.Value -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) {
    return $Default
  }
  return [string]$prop.Value
}

function As-Double([object]$Row, [string]$Name, [double]$Default = [double]::NaN) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  return [double]$value
}

function As-Int64([object]$Row, [string]$Name, [int64]$Default = [int64]::MinValue) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  return [int64]$value
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE256
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
    'src\tabi\tabi_fmm_kifmm_cartesian_bridge.cpp',
    'src\tabi\tabi_fmm_cartesian_sampled_panel_matvec.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) {
      throw 'failed to build titan_fmm_kifmm_hierarchy_real_panel_smoke'
    }
  } finally {
    Pop-Location
  }
}

function New-Case(
  [string]$Id,
  [string]$Fixture,
  [string]$Sample,
  [int]$SampleCount,
  [string]$Guard,
  [double]$GuardMarginFactor,
  [string]$Family
) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    guard = $Guard
    guard_margin_factor = $GuardMarginFactor
    family = $Family
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'order' 8
  Ensure-Field $Row 'theta' '0.6'
  Ensure-Field $Row 'equiv_scale' '1.05'
  Ensure-Field $Row 'upward_check_scale' '3.0'
  Ensure-Field $Row 'local_equiv_scale' '3.0'
  Ensure-Field $Row 'local_check_scale' '1.05'
  Ensure-Field $Row 'condition_policy' 'qr-rank-only'
  Ensure-Field $Row 'real_panel_guard' $Case.guard
  Ensure-Field $Row 'guard_margin_factor' $Case.guard_margin_factor
  Ensure-Field $Row 'family' $Case.family
  foreach ($field in @(
      'guard_margin_max',
      'real_guard_rejected_pairs',
      'real_guard_split_pairs',
      'real_guard_direct_pairs',
      'negative_target_check_excess_pairs',
      'negative_target_check_excess_panel_pairs',
      'min_leaf_half_width_over_panel_radius',
      'panels_crossing_node_boundary')) {
    Ensure-Field $Row $field 'MISSING'
  }
  return $Row
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.id)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--order=8",
    "--theta=0.6",
    "--equiv-scale=1.05",
    "--upward-check-scale=3.0",
    "--local-equiv-scale=3.0",
    "--local-check-scale=1.05",
    "--condition-policy=qr-rank-only",
    "--real-panel-guard=$($Case.guard)",
    "--guard-margin-factor=$($Case.guard_margin_factor)",
    "--mode=profile",
    "--pqr-1sh1=$Pqr1sh1",
    "--pqr-1a63=$Pqr1a63"
  )
  $stdoutPath = Join-Path $WorkRoot "$($Case.id).stdout.txt"
  $stderrPath = Join-Path $WorkRoot "$($Case.id).stderr.txt"
  Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue
  $proc = Start-Process -FilePath $ExePath -ArgumentList $args -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
  $output = @()
  if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath }
  if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath }
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($proc.ExitCode -ne 0 -or $resultLine.Count -ne 1) {
    $row = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      family = $Case.family
      order = 8
      theta = '0.6'
      equiv_scale = '1.05'
      upward_check_scale = '3.0'
      local_equiv_scale = '3.0'
      local_check_scale = '1.05'
      condition_policy = 'qr-rank-only'
      real_panel_guard = $Case.guard
      guard_margin_factor = $Case.guard_margin_factor
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'execution'
      aggregate_rms = 'NaN'
      max_channel_rms = 'NaN'
      max_channel_name = 'none'
      k4_rms = 'NaN'
      target_check_excess_min = 'NaN'
      negative_target_check_excess_pairs = 'NaN'
      production_invalid_pairs = 'NaN'
      nonproduction_fallback_pairs = 'NaN'
      m2l_translations = '0'
      far_fraction = 'NaN'
      process_exit_code = $proc.ExitCode
      process_output = ($output -join ' | ')
    }
    return Complete-Row $row $Case
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E255_RESULT '
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $Case.id -Force
  $row | Add-Member -NotePropertyName sample_count -NotePropertyValue $Case.sample_count -Force
  $row | Add-Member -NotePropertyName family -NotePropertyValue $Case.family -Force
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $proc.ExitCode -Force
  $row | Add-Member -NotePropertyName process_output -NotePropertyValue '' -Force
  Complete-Row $row $Case
}

if ((Test-Path -LiteralPath $WorkRoot) -and $Force) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$cases = @(
  (New-Case '1SH1_uniform_128_baseline_off' '1SH1' 'uniform_hash' 128 'off' 0.0 'baseline'),
  (New-Case '1SH1_uniform_128_guard_margin0' '1SH1' 'uniform_hash' 128 'strict' 0.0 'guard'),
  (New-Case '1SH1_uniform_128_guard_margin0p25diag' '1SH1' 'uniform_hash' 128 'strict' 0.25 'margin_sweep'),
  (New-Case '1SH1_uniform_128_guard_margin0p5diag' '1SH1' 'uniform_hash' 128 'strict' 0.5 'margin_sweep'),
  (New-Case '1SH1_closest_128_guard_margin0' '1SH1' 'closest_opposing_surface' 128 'strict' 0.0 'guard'),
  (New-Case '1A63_high_aspect_128_guard_margin0' '1A63' 'high_aspect' 128 'strict' 0.0 'guard')
)
if (-not $Quick) {
  $cases += @(
    (New-Case '1SH1_high_aspect_128_guard_margin0' '1SH1' 'high_aspect' 128 'strict' 0.0 'guard_full'),
    (New-Case '1SH1_groove_128_guard_margin0' '1SH1' 'inward_normal_groove' 128 'strict' 0.0 'guard_full'),
    (New-Case '1A63_uniform_128_guard_margin0' '1A63' 'uniform_hash' 128 'strict' 0.0 'guard_full'),
    (New-Case '1A63_closest_128_guard_margin0' '1A63' 'closest_opposing_surface' 128 'strict' 0.0 'guard_full'),
    (New-Case '1SH1_uniform_512_guard_margin0' '1SH1' 'uniform_hash' 512 'strict' 0.0 'guard_full'),
    (New-Case '1A63_uniform_512_guard_margin0' '1A63' 'uniform_hash' 512 'strict' 0.0 'guard_full')
  )
}

$rows = New-Object System.Collections.Generic.List[object]
$existingRows = @{}
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in (Import-Csv -LiteralPath $MatrixCsv)) {
    if ($row.PSObject.Properties['row_id'] -and -not [string]::IsNullOrWhiteSpace($row.row_id)) {
      $existingRows[$row.row_id] = $row
    }
  }
}
foreach ($case in $cases) {
  if ($existingRows.ContainsKey($case.id) -and
      (-not $RerunFailed -or ($existingRows[$case.id].execution_status -eq 'pass' -and $existingRows[$case.id].accuracy_status -eq 'pass'))) {
    Write-Host "resume: $($case.id)"
    $run = $existingRows[$case.id]
  } else {
    $run = Run-Case $case
  }
  $rows.Add($run)
  $rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows | Select-Object row_id,real_panel_guard,guard_margin_factor,guard_margin_max,real_guard_rejected_pairs,real_guard_split_pairs,real_guard_direct_pairs,negative_target_check_excess_pairs,negative_target_check_excess_panel_pairs,min_leaf_half_width_over_panel_radius,panels_crossing_node_boundary,target_check_excess_min |
  Export-Csv -LiteralPath $GuardCsv -NoTypeInformation
$rows | Select-Object row_id,aggregate_gate,channel_gate,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,max_channel_rms,max_channel_name,k4_abs_rms,k4_reference_rms,accuracy_status,accuracy_failure |
  Export-Csv -LiteralPath $ChannelCsv -NoTypeInformation
$rows | Select-Object row_id,source_panels,sampled_rows,sampled_dense_valid_pairs,production_fallback_pairs,production_invalid_pairs,nonproduction_fallback_pairs,far_panel_pairs,direct_panel_pairs,coverage,expected_coverage,far_fraction,m2l_translations |
  Export-Csv -LiteralPath $CoverageCsv -NoTypeInformation

$requiredRows = @($rows | Where-Object { $_.real_panel_guard -eq 'strict' })
$requiredExecutionFailures = @($requiredRows | Where-Object { $_.execution_status -ne 'pass' })
$passRows = @($requiredRows | Where-Object {
    $_.execution_status -eq 'pass' -and
    $_.accuracy_status -eq 'pass' -and
    (As-Double $_ 'aggregate_rms') -le 1.0e-6 -and
    (As-Double $_ 'max_channel_rms') -le 5.0e-6 -and
    (As-Double $_ 'k4_rms') -le 5.0e-6 -and
    (As-Int64 $_ 'negative_target_check_excess_pairs') -eq 0 -and
    (As-Double $_ 'target_check_excess_min') -ge 0.0 -and
    (As-Int64 $_ 'production_invalid_pairs') -eq 0 -and
    (As-Int64 $_ 'nonproduction_fallback_pairs') -eq 0 -and
    (As-Int64 $_ 'm2l_translations') -gt 0
  })
$strictNegative = @($requiredRows | Where-Object {
    $_.execution_status -eq 'pass' -and
    ((As-Int64 $_ 'negative_target_check_excess_pairs') -gt 0 -or
     (As-Double $_ 'target_check_excess_min') -lt 0.0)
  })
$accuracyFailsNoNegative = @($requiredRows | Where-Object {
    $_.execution_status -eq 'pass' -and
    $_.accuracy_status -ne 'pass' -and
    (As-Int64 $_ 'negative_target_check_excess_pairs') -eq 0 -and
    (As-Double $_ 'target_check_excess_min') -ge 0.0
  })
$fallbackHeavyPass = @($passRows | Where-Object { (As-Double $_ 'far_fraction') -lt 0.10 })

$decision = if ($requiredExecutionFailures.Count -gt 0) {
  'fmm_e256_execution_gap'
} elseif ($strictNegative.Count -gt 0) {
  'fmm_e256_negative_check_excess_persists'
} elseif ($accuracyFailsNoNegative.Count -gt 0) {
  'fmm_e256_k4_persists_after_guard'
} elseif ($passRows.Count -gt 0 -and $fallbackHeavyPass.Count -gt 0) {
  'fmm_e256_accuracy_pass_fallback_heavy'
} elseif ($passRows.Count -gt 0) {
  'fmm_e256_real_panel_admissibility_guard_pass'
} else {
  'fmm_e256_k4_persists_after_guard'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  required_rows = $requiredRows.Count
  pass_rows = $passRows.Count
  fallback_heavy_pass_rows = $fallbackHeavyPass.Count
  execution_gap_rows = $requiredExecutionFailures.Count
  negative_check_excess_rows = $strictNegative.Count
  accuracy_fail_no_negative_rows = $accuracyFailsNoNegative.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "guard=$GuardCsv",
  "channel=$ChannelCsv",
  "coverage=$CoverageCsv",
  "decision_csv=$DecisionCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.6 Real-Panel Admissibility Report',
  '',
  "Decision: $decision",
  '',
  'E25.6 evaluates strict real-panel admissibility guards for sampled KIFMM hierarchy rows.',
  'Pass requires clean execution, q8 channel accuracy, no negative target-check excess, no production invalid/fallback leakage, and active M2L work.',
  ''
)
foreach ($row in $rows) {
  $report += ('  {0}: guard={1}, margin={2}, exec={3}, acc={4}, rms={5}, max={6}:{7}, k4={8}, neg_pairs={9}, check_min={10}, far_fraction={11}, m2l={12}, direct_guard={13}, total_s={14}' -f
    (Get-Field $row 'row_id'),
    (Get-Field $row 'real_panel_guard'),
    (Get-Field $row 'guard_margin_factor'),
    (Get-Field $row 'execution_status'),
    (Get-Field $row 'accuracy_status'),
    (Get-Field $row 'aggregate_rms'),
    (Get-Field $row 'max_channel_name'),
    (Get-Field $row 'max_channel_rms'),
    (Get-Field $row 'k4_rms'),
    (Get-Field $row 'negative_target_check_excess_pairs'),
    (Get-Field $row 'target_check_excess_min'),
    (Get-Field $row 'far_fraction'),
    (Get-Field $row 'm2l_translations'),
    (Get-Field $row 'real_guard_direct_pairs'),
    (Get-Field $row 'total_seconds'))
}
$report += ''
$report += 'Interpretation: the guard passes only when strict rows remove negative check-surface excess without losing K4/channel accuracy. A fallback-heavy pass means the geometry guard is correct but too conservative for production use.'
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.6 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Guard: $GuardCsv"
Write-Host "Channel: $ChannelCsv"
Write-Host "Coverage: $CoverageCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e256_real_panel_admissibility_guard_pass','fmm_e256_accuracy_pass_fallback_heavy') -and -not $AllowFailedDecision) {
  throw "FMM-E25.6 gate did not pass: $decision"
}
