param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 300,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e34_full_panel_dual_worklist'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e34_matrix.csv'
$AuditCsv = Join-Path $WorkRoot 'fmm_e34_sampled_audit.csv'
$CoverageCsv = Join-Path $WorkRoot 'fmm_e34_coverage.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e34_runtime.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e34_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e34_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E34_FULL_PANEL_DUAL_WORKLIST_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E34 smoke.'
}

function Assert-SmokeSupportsE34 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'full-panel-profile',
      'run_kifmm_hierarchy_real_panel_full_apply_sampled_audit',
      'target_panels=',
      'full_coverage_expected=',
      'sampled_coverage_expected=',
      'build_direct_reference')) {
    if (-not $text.Contains($needle)) {
      throw "E34 full-panel sampled-audit support is not present yet: missing '$needle' in $SmokeSource"
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
  Assert-SmokeSupportsE34
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

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'surface_points' $Case.surface_points
  Ensure-Field $Row 'panel_limit' $Case.panel_limit
  Ensure-Field $Row 'mode' $Case.mode
  Ensure-Field $Row 'timing_scope' 'MISSING'
  Ensure-Field $Row 'target_panels' 'MISSING'
  Ensure-Field $Row 'full_coverage_expected' 'MISSING'
  Ensure-Field $Row 'sampled_coverage_expected' 'MISSING'
  Ensure-Field $Row 'sampled_dense_valid_pairs' 'MISSING'
  Ensure-Field $Row 'cartesian_k4_accuracy_status' 'MISSING'
  Ensure-Field $Row 'dual_worklist_status' 'MISSING'
  Ensure-Field $Row 'channel_recombination_status' 'MISSING'
  Ensure-Field $Row 'worklist_independence_status' 'MISSING'
  $Row
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.id) (timeout ${CaseTimeoutSeconds}s)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--surface-points=$($Case.surface_points)",
    "--panel-limit=$($Case.panel_limit)",
    "--order=8",
    "--theta=0.6",
    "--equiv-scale=1.05",
    "--upward-check-scale=3.0",
    "--local-equiv-scale=3.0",
    "--local-check-scale=1.05",
    "--condition-policy=qr-rank-only",
    "--real-panel-guard=strict-local-check",
    "--guard-margin-factor=0",
    "--hier-debug-mode=dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian",
    "--mode=full-panel-profile",
    "--pqr-1sh1=$Pqr1sh1",
    "--pqr-1a63=$Pqr1a63"
  )
  $stdoutPath = Join-Path $WorkRoot "$($Case.id).stdout.txt"
  $stderrPath = Join-Path $WorkRoot "$($Case.id).stderr.txt"
  Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue
  $proc = Start-Process -FilePath $ExePath -ArgumentList $args -NoNewWindow -PassThru `
    -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
  $completed = $proc.WaitForExit([Math]::Max(1, $CaseTimeoutSeconds) * 1000)
  if (-not $completed) {
    Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
    $proc.WaitForExit()
  }
  $output = @()
  if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath }
  if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath }
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($resultLine.Count -ne 1) {
    $failedRow = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      mode = $Case.mode
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = $(if ($completed) { 'process_failed' } else { 'timeout' })
      status_detail = $(if ($completed) { 'process_failed' } else { 'timeout' })
      process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }
    return (Complete-Row $failedRow $Case)
  }
  $row = Complete-Row (Parse-KeyValueLine ([string]$resultLine[-1]) 'E255_RESULT') $Case
  Ensure-Field $row 'process_exit_code' $proc.ExitCode
  $row
}

function New-Case([string]$Id, [string]$Fixture, [string]$Sample, [int]$SampleCount, [int]$SurfacePoints, [int]$PanelLimit, [string]$Mode) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    mode = $Mode
  }
}

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) {
    Set-Content -LiteralPath $Path -Value ''
  } else {
    $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path
  }
}

function Row-Passes([object]$Row) {
  if ((Get-Field $Row 'execution_status') -ne 'pass') { return $false }
  if ((Get-Field $Row 'accuracy_status') -ne 'pass') { return $false }
  if ((Get-Field $Row 'timing_scope') -ne 'full_panel_apply_sampled_audit') { return $false }
  if ((Get-Field $Row 'coverage') -ne (Get-Field $Row 'full_coverage_expected')) { return $false }
  if ((Get-Field $Row 'sampled_dense_valid_pairs') -ne (Get-Field $Row 'sampled_coverage_expected')) { return $false }
  if ((As-Double $Row 'aggregate_rms') -gt 1.0e-6) { return $false }
  if ((As-Double $Row 'max_channel_rms') -gt 5.0e-6) { return $false }
  if ((As-Double $Row 'k4_rms') -gt 5.0e-6) { return $false }
  if ((As-Int64 $Row 'production_invalid_pairs' 1) -ne 0) { return $false }
  if ((As-Int64 $Row 'nonproduction_fallback_pairs' 1) -ne 0) { return $false }
  if ((As-Double $Row 'far_fraction' 0.0) -lt 0.10) { return $false }
  foreach ($field in @('worklist_hash','policy_hash','panel_geometry_hash','unknown_vector_hash','sampled_target_hash','output_hash')) {
    if ((As-Int64 $Row $field 0) -eq 0) { return $false }
  }
  return $true
}

function Decide([array]$Rows, [int]$RequiredRows) {
  if ($Rows.Count -lt $RequiredRows) { return 'fmm_e34_runtime_limited' }
  if (@($Rows | Where-Object { (Get-Field $_ 'status_detail') -eq 'timeout' -or (Get-Field $_ 'accuracy_failure') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e34_runtime_limited'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e34_execution_gap'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'coverage') -ne (Get-Field $_ 'full_coverage_expected') -or (Get-Field $_ 'sampled_dense_valid_pairs') -ne (Get-Field $_ 'sampled_coverage_expected') }).Count -gt 0) {
    return 'fmm_e34_full_coverage_gap'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'cartesian_k4_accuracy_status') -ne 'pass' -or (As-Double $_ 'k4_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e34_cartesian_k4_gap'
  }
  if (@($Rows | Where-Object { (As-Double $_ 'k1_rms') -gt 5.0e-6 -or (As-Double $_ 'k2_rms') -gt 5.0e-6 -or (As-Double $_ 'k3_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e34_k123_hierarchy_gap'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'accuracy_status') -eq 'pass' -and (As-Double $_ 'far_fraction' 0.0) -lt 0.10 }).Count -gt 0) {
    return 'fmm_e34_fallback_heavy_accuracy_pass'
  }
  if (@($Rows | Where-Object { -not (Row-Passes $_) }).Count -eq 0) {
    return 'fmm_e34_full_panel_dual_worklist_pass'
  }
  return 'fmm_e34_sampled_audit_gap'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$cases = @(
  (New-Case '1SH1_uniform_hash_s4_p128_32_full_panel_dual' '1SH1' 'uniform_hash' 32 4 128 'full_panel_dual_worklist_channel_hybrid')
)
if (-not $Quick) {
  $cases += @(
    (New-Case '1SH1_uniform_hash_s8_128_full_panel_dual' '1SH1' 'uniform_hash' 128 8 0 'full_panel_dual_worklist_channel_hybrid'),
    (New-Case '1SH1_closest_s8_64_full_panel_dual' '1SH1' 'closest_opposing_surface' 64 8 0 'full_panel_dual_worklist_channel_hybrid'),
    (New-Case '1SH1_high_aspect_s8_64_full_panel_dual' '1SH1' 'high_aspect' 64 8 0 'full_panel_dual_worklist_channel_hybrid')
  )
}

$existing = @{}
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in @(Import-Csv -LiteralPath $MatrixCsv)) {
    $existing[(Get-Field $row 'row_id')] = $row
  }
}

$matrix = @()
foreach ($case in $cases) {
  $old = $existing[$case.id]
  $needsRun = $Force -or $null -eq $old
  if ($RerunFailed -and $null -ne $old -and (Get-Field $old 'execution_status') -ne 'pass') {
    $needsRun = $true
  }
  if (-not $needsRun) {
    Write-Host "resume: $($case.id)"
    $matrix += (Complete-Row $old $case)
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

Export-Table ($matrix | Select-Object row_id,fixture,sample,surface_points,panel_limit,timing_scope,status,execution_status,accuracy_status,accuracy_failure,aggregate_rms,max_channel_rms,max_channel_name,k1_rms,k2_rms,k3_rms,k4_rms,k4_abs_rms,k4_reference_rms,cartesian_k4_accuracy_status) $AuditCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,surface_points,panel_limit,source_panels,target_panels,sampled_rows,coverage,full_coverage_expected,sampled_dense_valid_pairs,sampled_coverage_expected,production_invalid_pairs,nonproduction_fallback_pairs,far_panel_pairs,direct_panel_pairs,far_fraction,worklist_hash,policy_hash,panel_geometry_hash,unknown_vector_hash,sampled_target_hash,output_hash) $CoverageCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,hierarchy_seconds,production_fallback_seconds,dense_oracle_seconds,cartesian_k4_total_seconds,total_seconds) $RuntimeCsv

$decision = Decide $matrix $cases.Count
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'full' })
  rows = $matrix.Count
  required_rows = $cases.Count
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$MatrixCsv"
)

$report = @(
  '# TITAN FMM Engine E34 Full-Panel Dual-Worklist Channel-Hybrid Report',
  '',
  "Decision: $decision",
  '',
  'E34 evaluates all target panels with the dual-worklist channel-hybrid q8 operator, then verifies only deterministic sampled target rows with dense production direct rows.',
  'K1/K2/K3 are sourced from full hierarchical KiFMM q8 strict-local-check; K4 is sourced from guarded Cartesian p5 geometry0p4 using its own independent worklist. No GMRES, full dense protein oracle, production backend, or default change is included.'
)
if ($decision -eq 'fmm_e34_fallback_heavy_accuracy_pass') {
  $report += 'Interpretation: the full-panel sampled-audit contract passes on the quick real-panel subset, but the far fraction is below the 0.10 compression gate. This is a correctness scaffold, not a useful KIFMM-compression result.'
} elseif ($decision -eq 'fmm_e34_full_panel_dual_worklist_pass') {
  $report += 'Interpretation: full-panel sampled-audit accuracy, coverage, and far-fraction gates pass on the required rows.'
} elseif ($decision -eq 'fmm_e34_runtime_limited') {
  $report += 'Interpretation: at least one required full-panel row timed out before a complete sampled dense audit record was produced.'
}
foreach ($row in $matrix) {
  $report += ('  {0}: exec={1}, acc={2}, scope={3}, agg={4}, max={5}/{6}, K1={7}, K2={8}, K3={9}, K4={10}, coverage={11}/{12}, sampled={13}/{14}, total={15}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'execution_status'), (Get-Field $row 'accuracy_status'), `
      (Get-Field $row 'timing_scope'), (Get-Field $row 'aggregate_rms'), `
      (Get-Field $row 'max_channel_name'), (Get-Field $row 'max_channel_rms'), `
      (Get-Field $row 'k1_rms'), (Get-Field $row 'k2_rms'), (Get-Field $row 'k3_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'coverage'), (Get-Field $row 'full_coverage_expected'), `
      (Get-Field $row 'sampled_dense_valid_pairs'), (Get-Field $row 'sampled_coverage_expected'), `
      (Get-Field $row 'total_seconds'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E34 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Sampled audit: $AuditCsv"
Write-Host "Coverage: $CoverageCsv"
Write-Host "Runtime: $RuntimeCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e34_full_panel_dual_worklist_pass','fmm_e34_k123_hierarchy_gap')) {
  throw "FMM-E34 decision failed: $decision"
}
