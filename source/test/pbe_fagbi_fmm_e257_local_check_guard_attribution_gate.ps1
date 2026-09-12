param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e257_local_check_guard_attribution'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e257_matrix.csv'
$GuardCsv = Join-Path $WorkRoot 'fmm_e257_local_check_guard.csv'
$AttributionCsv = Join-Path $WorkRoot 'fmm_e257_flat_vs_hierarchy.csv'
$ContainmentCsv = Join-Path $WorkRoot 'fmm_e257_l2l_containment.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e257_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e257_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E257_LOCAL_CHECK_ATTRIBUTION_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E25.7 smoke.'
}

function Assert-SmokeSupportsE257 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'strict-local-check',
      'target_local_check_excess_min',
      'flat_vs_dense_aggregate',
      'hier_vs_flat_k4',
      'child_local_check_inside_parent_local_check_min_margin')) {
    if (-not $text.Contains($needle)) {
      throw "E25.7 smoke support is not present yet: missing '$needle' in $SmokeSource"
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
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  [double]$value
}

function As-Int64([object]$Row, [string]$Name, [int64]$Default = [int64]::MinValue) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  [int64]$value
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE257
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
  [string]$Mode,
  [string]$Guard,
  [double]$GuardMarginFactor
) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    mode = $Mode
    guard = $Guard
    guard_margin_factor = $GuardMarginFactor
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'mode' $Case.mode
  Ensure-Field $Row 'order' 8
  Ensure-Field $Row 'theta' '0.6'
  Ensure-Field $Row 'upward_check_scale' '3.0'
  Ensure-Field $Row 'local_equiv_scale' '3.0'
  Ensure-Field $Row 'local_check_scale' '1.05'
  Ensure-Field $Row 'condition_policy' 'qr-rank-only'
  Ensure-Field $Row 'real_panel_guard' $Case.guard
  Ensure-Field $Row 'guard_margin_factor' $Case.guard_margin_factor
  foreach ($field in @(
      'target_panel_envelope_excess_min',
      'target_local_check_excess_min',
      'negative_target_local_check_excess_pairs',
      'negative_target_local_check_excess_panel_pairs',
      'flat_vs_dense_aggregate',
      'flat_vs_dense_k4',
      'flat_vs_dense_max_channel',
      'hier_vs_dense_aggregate',
      'hier_vs_dense_k4',
      'hier_vs_flat_aggregate',
      'hier_vs_flat_k4',
      'child_local_check_inside_parent_local_check_min_margin',
      'negative_l2l_local_containment_pairs')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $Row
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
      mode = $Case.mode
      order = 8
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
      process_exit_code = $proc.ExitCode
      process_output = ($output -join ' | ')
    }
    return Complete-Row $row $Case
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E255_RESULT '
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $Case.id -Force
  $row | Add-Member -NotePropertyName sample_count -NotePropertyValue $Case.sample_count -Force
  $row | Add-Member -NotePropertyName mode -NotePropertyValue $Case.mode -Force
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $proc.ExitCode -Force
  $row | Add-Member -NotePropertyName process_output -NotePropertyValue '' -Force
  Complete-Row $row $Case
}

if ((Test-Path -LiteralPath $WorkRoot) -and $Force) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$baseCases = @(
  @{ Fixture='1SH1'; Sample='uniform_hash'; Count=128; Prefix='1SH1_uniform_128' },
  @{ Fixture='1SH1'; Sample='closest_opposing_surface'; Count=128; Prefix='1SH1_closest_128' },
  @{ Fixture='1A63'; Sample='high_aspect'; Count=128; Prefix='1A63_high_aspect_128' }
)
$cases = @()
foreach ($base in $baseCases) {
  $cases += New-Case "$($base.Prefix)_panel_envelope_guard_current" $base.Fixture $base.Sample $base.Count 'panel_envelope_guard_current' 'strict' 0.0
  $cases += New-Case "$($base.Prefix)_local_check_guard_margin0" $base.Fixture $base.Sample $base.Count 'local_check_guard_margin0' 'strict-local-check' 0.0
  $cases += New-Case "$($base.Prefix)_local_check_guard_margin0p25diag" $base.Fixture $base.Sample $base.Count 'local_check_guard_margin0p25diag' 'strict-local-check' 0.25
  $cases += New-Case "$($base.Prefix)_local_check_guard_margin0p5diag" $base.Fixture $base.Sample $base.Count 'local_check_guard_margin0p5diag' 'strict-local-check' 0.5
}

$existingRows = @{}
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in (Import-Csv -LiteralPath $MatrixCsv)) {
    if ($row.PSObject.Properties['row_id'] -and -not [string]::IsNullOrWhiteSpace($row.row_id)) {
      $existingRows[$row.row_id] = $row
    }
  }
}
$rows = New-Object System.Collections.Generic.List[object]
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
$rows | Select-Object row_id,mode,real_panel_guard,guard_margin_factor,target_panel_envelope_excess_min,target_local_check_excess_min,negative_target_check_excess_pairs,negative_target_local_check_excess_pairs,negative_target_local_check_excess_panel_pairs,real_guard_rejected_pairs,real_guard_split_pairs,real_guard_direct_pairs |
  Export-Csv -LiteralPath $GuardCsv -NoTypeInformation
$rows | Select-Object row_id,mode,flat_vs_dense_aggregate,flat_vs_dense_k4,flat_vs_dense_max_channel,hier_vs_dense_aggregate,hier_vs_dense_k4,hier_vs_flat_aggregate,hier_vs_flat_k4,aggregate_rms,max_channel_rms,k4_rms,far_fraction,m2l_translations |
  Export-Csv -LiteralPath $AttributionCsv -NoTypeInformation
$rows | Select-Object row_id,mode,child_local_check_inside_parent_local_check_min_margin,negative_l2l_local_containment_pairs |
  Export-Csv -LiteralPath $ContainmentCsv -NoTypeInformation

$requiredRows = @($rows | Where-Object { $_.real_panel_guard -eq 'strict-local-check' })
$executionGap = @($requiredRows | Where-Object { $_.execution_status -ne 'pass' })
$coverageGap = @($requiredRows | Where-Object {
    (As-Int64 $_ 'sampled_dense_valid_pairs') -ne (As-Int64 $_ 'expected_coverage') -or
    (As-Int64 $_ 'coverage') -ne (As-Int64 $_ 'expected_coverage') -or
    (As-Int64 $_ 'production_invalid_pairs') -ne 0 -or
    (As-Int64 $_ 'nonproduction_fallback_pairs') -ne 0
  })
$localCheckGap = @($requiredRows | Where-Object {
    (As-Int64 $_ 'negative_target_local_check_excess_pairs') -ne 0 -or
    (As-Double $_ 'target_local_check_excess_min') -lt 0.0
  })
$flatPassRows = @($requiredRows | Where-Object {
    (As-Double $_ 'flat_vs_dense_aggregate') -le 1.0e-6 -and
    (As-Double $_ 'flat_vs_dense_max_channel') -le 5.0e-6 -and
    (As-Double $_ 'flat_vs_dense_k4') -le 5.0e-6
  })
$hierPassRows = @($requiredRows | Where-Object {
    $_.accuracy_status -eq 'pass' -and
    (As-Double $_ 'aggregate_rms') -le 1.0e-6 -and
    (As-Double $_ 'max_channel_rms') -le 5.0e-6 -and
    (As-Double $_ 'k4_rms') -le 5.0e-6 -and
    (As-Int64 $_ 'm2l_translations') -gt 0
  })
$flatPassHierarchyFail = @($requiredRows | Where-Object {
    (As-Double $_ 'flat_vs_dense_aggregate') -le 1.0e-6 -and
    (As-Double $_ 'flat_vs_dense_max_channel') -le 5.0e-6 -and
    (As-Double $_ 'flat_vs_dense_k4') -le 5.0e-6 -and
    -not ($_.accuracy_status -eq 'pass' -and (As-Double $_ 'k4_rms') -le 5.0e-6)
  })
$fallbackHeavyPass = @($hierPassRows | Where-Object { (As-Double $_ 'far_fraction') -lt 0.10 })

$decision = if ($executionGap.Count -gt 0) {
  'fmm_e257_execution_gap'
} elseif ($coverageGap.Count -gt 0) {
  'fmm_e257_coverage_gap'
} elseif ($localCheckGap.Count -gt 0) {
  'fmm_e257_local_check_guard_internal_gap'
} elseif ($hierPassRows.Count -eq $requiredRows.Count -and $fallbackHeavyPass.Count -gt 0) {
  'fmm_e257_fallback_heavy_accuracy_pass'
} elseif ($hierPassRows.Count -eq $requiredRows.Count) {
  'fmm_e257_local_check_guard_pass'
} elseif ($flatPassHierarchyFail.Count -gt 0) {
  'fmm_e257_flat_pass_hierarchy_fail'
} elseif ($flatPassRows.Count -eq 0) {
  'fmm_e257_worklist_still_too_aggressive'
} else {
  'fmm_e257_k4_persists_after_local_check_guard'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  required_rows = $requiredRows.Count
  execution_gap_rows = $executionGap.Count
  coverage_gap_rows = $coverageGap.Count
  local_check_gap_rows = $localCheckGap.Count
  flat_pass_rows = $flatPassRows.Count
  hierarchy_pass_rows = $hierPassRows.Count
  fallback_heavy_pass_rows = $fallbackHeavyPass.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "guard=$GuardCsv",
  "attribution=$AttributionCsv",
  "containment=$ContainmentCsv",
  "decision_csv=$DecisionCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.7 Local-Check Attribution Report',
  '',
  "Decision: $decision",
  '',
  'E25.7 tests local-check-aware M2L admissibility and compares flat-on-same-worklist against the full hierarchy.',
  ''
)
foreach ($row in $rows) {
  $report += ('  {0}: mode={1}, exec={2}, acc={3}, agg={4}, k4={5}, flat_agg={6}, flat_k4={7}, hier_flat_k4={8}, local_check_min={9}, neg_local={10}, far_fraction={11}, m2l={12}' -f
    (Get-Field $row 'row_id'),
    (Get-Field $row 'mode'),
    (Get-Field $row 'execution_status'),
    (Get-Field $row 'accuracy_status'),
    (Get-Field $row 'aggregate_rms'),
    (Get-Field $row 'k4_rms'),
    (Get-Field $row 'flat_vs_dense_aggregate'),
    (Get-Field $row 'flat_vs_dense_k4'),
    (Get-Field $row 'hier_vs_flat_k4'),
    (Get-Field $row 'target_local_check_excess_min'),
    (Get-Field $row 'negative_target_local_check_excess_pairs'),
    (Get-Field $row 'far_fraction'),
    (Get-Field $row 'm2l_translations'))
}
$report += ''
$report += 'Interpretation: flat-pass/hierarchy-fail points to propagation/local-solve issues; flat-fail/hierarchy-fail points to admissibility/worklist aggressiveness or real-panel field construction.'
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.7 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Guard: $GuardCsv"
Write-Host "Attribution: $AttributionCsv"
Write-Host "Containment: $ContainmentCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e257_local_check_guard_pass','fmm_e257_flat_pass_hierarchy_fail','fmm_e257_worklist_still_too_aggressive','fmm_e257_k4_persists_after_local_check_guard','fmm_e257_fallback_heavy_accuracy_pass') -and -not $AllowFailedDecision) {
  throw "FMM-E25.7 gate did not produce an interpretable diagnostic decision: $decision"
}
