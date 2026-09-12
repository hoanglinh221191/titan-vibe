param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e28_split_hybrid_local_bridge'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e28_matrix.csv'
$SplitCsv = Join-Path $WorkRoot 'fmm_e28_k4_split.csv'
$ComparisonCsv = Join-Path $WorkRoot 'fmm_e28_comparison.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e28_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e28_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E28_SPLIT_HYBRID_LOCAL_BRIDGE_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E28 smoke.'
}

function Assert-SmokeSupportsE28 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'hybrid-local-bridge',
      'hybrid-local-split-bridge',
      'hybrid_mode',
      'k4_decomposition_status',
      'split_k4_laplace',
      'split_k4_yukawa',
      'split_k4_combined',
      'split_vs_dense_k4',
      'split_vs_flat_k4',
      'l2p_gradient_method')) {
    if (-not $text.Contains($needle)) {
      throw "E28 split hybrid support is not present yet: missing '$needle' in $SmokeSource"
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
  Assert-SmokeSupportsE28
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
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
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

function New-Case([string]$Id, [string]$Fixture, [string]$Sample, [int]$SampleCount, [string]$Mode) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    mode = $Mode
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'mode' $Case.mode
  Ensure-Field $Row 'hier_debug_mode' $Case.mode
  foreach ($field in @(
      'hybrid_mode',
      'm2e_residual',
      'local_equiv_residual',
      'e2l_coeff_rms',
      'e2l_coeff_status',
      'hierarchy_inner_status',
      'hierarchy_inner_accuracy_failure',
      'k4_decomposition_status',
      'l2p_gradient_method',
      'hybrid_vs_dense_k4',
      'hybrid_vs_flat_k4',
      'split_k4_laplace',
      'split_k4_yukawa',
      'split_k4_combined',
      'split_vs_dense_k4',
      'split_vs_flat_k4',
      'pure_hierarchy_vs_dense_k4',
      'flat_vs_dense_k4')) {
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
    "--real-panel-guard=strict-local-check",
    "--guard-margin-factor=0",
    "--hier-debug-mode=$($Case.mode)",
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
  if ($resultLine.Count -ne 1) {
    $failedRow = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      mode = $Case.mode
      hier_debug_mode = $Case.mode
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'process_failed'
      process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }
    return (Complete-Row $failedRow $Case)
  }
  $row = Complete-Row (Parse-KeyValueLine ([string]$resultLine[-1]) 'E255_RESULT') $Case
  Ensure-Field $row 'process_exit_code' $proc.ExitCode
  $row
}

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) {
    Set-Content -LiteralPath $Path -Value ''
  } else {
    $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path
  }
}

function Mode-Pass([array]$Rows, [string]$Mode) {
  $selected = @($Rows | Where-Object { $_.mode -eq $Mode })
  if ($selected.Count -eq 0) { return $false }
  foreach ($row in $selected) {
    if ((Get-Field $row 'execution_status') -ne 'pass') { return $false }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0) { return $false }
    if ((As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0) { return $false }
    if ((Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage')) { return $false }
    if ((As-Double $row 'far_fraction' 0.0) -le 0.0) { return $false }
    if ((Get-Field $row 'accuracy_status') -ne 'pass') { return $false }
    if ((As-Double $row 'k4_rms') -gt 5.0e-6) { return $false }
  }
  return $true
}

function Split-Pass([array]$Rows) {
  $selected = @($Rows | Where-Object { $_.mode -eq 'hybrid-local-split-bridge' })
  if ($selected.Count -ne 3) { return $false }
  foreach ($row in $selected) {
    if ((Get-Field $row 'execution_status') -ne 'pass') { return $false }
    if ((Get-Field $row 'accuracy_status') -ne 'pass') { return $false }
    if ((As-Double $row 'aggregate_rms') -gt 1.0e-6) { return $false }
    if ((As-Double $row 'max_channel_rms') -gt 5.0e-6) { return $false }
    if ((As-Double $row 'k4_rms') -gt 5.0e-6) { return $false }
    if ((Get-Field $row 'k4_decomposition_status') -ne 'laplace_yukawa_split') { return $false }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0) { return $false }
    if ((As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0) { return $false }
    if ((Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage')) { return $false }
  }
  return $true
}

function Decide([array]$Rows) {
  if ($Rows.Count -lt 9) { return 'fmm_e28_runtime_limited' }
  if (Split-Pass $Rows) { return 'fmm_e28_split_hybrid_local_bridge_pass' }
  $splitRows = @($Rows | Where-Object { $_.mode -eq 'hybrid-local-split-bridge' })
  if (@($splitRows | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e28_split_bridge_execution_gap'
  }
  if (@($splitRows | Where-Object { (Get-Field $_ 'k4_decomposition_status') -ne 'laplace_yukawa_split' }).Count -gt 0) {
    return 'fmm_e28_split_decomposition_gap'
  }
  if (@($splitRows | Where-Object { (As-Double $_ 'k4_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e28_split_k4_accuracy_gap'
  }
  'fmm_e28_split_aggregate_or_channel_gap'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$modes = @('full', 'hybrid-local-bridge', 'hybrid-local-split-bridge')
$samples = @(
  @{ fixture = '1SH1'; sample = 'uniform_hash'; count = 128; prefix = '1SH1_uniform_128' },
  @{ fixture = '1SH1'; sample = 'closest_opposing_surface'; count = 128; prefix = '1SH1_closest_128' },
  @{ fixture = '1A63'; sample = 'high_aspect'; count = 128; prefix = '1A63_high_aspect_128' }
)
$cases = @()
foreach ($sample in $samples) {
  foreach ($mode in $modes) {
    $cases += New-Case "$($sample.prefix)_$($mode.Replace('-', '_'))" $sample.fixture $sample.sample $sample.count $mode
  }
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

$splitRows = $matrix | Select-Object row_id,fixture,sample,mode,status,execution_status,accuracy_status,aggregate_rms,max_channel_rms,k4_rms,k4_decomposition_status,split_k4_laplace,split_k4_yukawa,split_k4_combined,split_vs_dense_k4,split_vs_flat_k4,flat_vs_dense_k4,hierarchy_inner_status,e2l_coeff_status,l2p_gradient_method
$comparisonRows = $matrix | Select-Object row_id,fixture,sample,mode,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,far_fraction,direct_pairs,total_seconds,accuracy_status
Export-Table $splitRows $SplitCsv
Export-Table $comparisonRows $ComparisonCsv

$decision = Decide $matrix
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'quick' })
  rows = $matrix.Count
  required_rows = 9
  full_pass = (Mode-Pass $matrix 'full')
  hybrid_local_pass = (Mode-Pass $matrix 'hybrid-local-bridge')
  hybrid_split_pass = (Split-Pass $matrix)
  hybrid_local_k4_decomposition_status = ((@($matrix | Where-Object { $_.mode -eq 'hybrid-local-bridge' } | ForEach-Object { Get-Field $_ 'k4_decomposition_status' }) | Select-Object -Unique) -join ';')
  hybrid_split_k4_decomposition_status = ((@($matrix | Where-Object { $_.mode -eq 'hybrid-local-split-bridge' } | ForEach-Object { Get-Field $_ 'k4_decomposition_status' }) | Select-Object -Unique) -join ';')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$MatrixCsv"
)

$report = @(
  '# TITAN FMM Engine E28 Split Hybrid Local Bridge Report',
  '',
  "Decision: $decision",
  '',
  'E28 tests a sampled-row diagnostic mode only. It does not promote KIFMM, hybrid bridge, Cartesian local expansions, or any production TABI default.',
  'The question is whether preserving Laplace/Yukawa K4 decomposition through hybrid E2L/L2P reduces the E27 K4 gap on the same real-panel sampled rows.'
)
if ($decision -eq 'fmm_e28_split_hybrid_local_bridge_pass') {
  $report += 'Interpretation: split K4 decomposition closes the sampled hybrid-local K4 gate on the three E27 rows.'
} elseif ($decision -eq 'fmm_e28_split_k4_accuracy_gap') {
  $report += 'Interpretation: split mode executes and emits laplace_yukawa_split telemetry, but K4 remains above gate. The E27 gap is not explained by decomposition bookkeeping alone.'
} else {
  $report += 'Interpretation: split K4 decomposition is diagnostic evidence only; inspect fmm_e28_k4_split.csv before planning broader real-panel coverage.'
}
$report += 'The split_lap and split_yuk component RMS values are diagnostic only: large component errors can coexist with cancellation, so E28 gates only the recombined sampled K4 row against dense production direct and flat-worklist references.'
foreach ($row in $matrix) {
  $report += ('  {0}: mode={1}, exec={2}, acc={3}, inner={4}, agg={5}, k4={6}, split_vs_flat_k4={7}, flat_k4={8}, split_lap={9}, split_yuk={10}, e2l_status={11}, k4_decomp={12}, l2p={13}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'mode'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'hierarchy_inner_status'), `
      (Get-Field $row 'aggregate_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'split_vs_flat_k4'), (Get-Field $row 'flat_vs_dense_k4'), `
      (Get-Field $row 'split_k4_laplace'), (Get-Field $row 'split_k4_yukawa'), `
      (Get-Field $row 'e2l_coeff_status'), (Get-Field $row 'k4_decomposition_status'), `
      (Get-Field $row 'l2p_gradient_method'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E28 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Split: $SplitCsv"
Write-Host "Comparison: $ComparisonCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @(
    'fmm_e28_split_hybrid_local_bridge_pass')) {
  throw "FMM-E28 decision failed: $decision"
}
