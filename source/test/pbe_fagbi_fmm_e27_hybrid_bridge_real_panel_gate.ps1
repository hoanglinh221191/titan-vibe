param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e27_hybrid_bridge_real_panel'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e27_matrix.csv'
$BridgeCsv = Join-Path $WorkRoot 'fmm_e27_bridge_stage.csv'
$K4Csv = Join-Path $WorkRoot 'fmm_e27_k4_decomposition.csv'
$ComparisonCsv = Join-Path $WorkRoot 'fmm_e27_comparison.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e27_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e27_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E27_HYBRID_BRIDGE_REAL_PANEL_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E27 smoke.'
}

function Assert-SmokeSupportsE27 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'hybrid-local-bridge',
      'hybrid-source-local-bridge',
      'hybrid_mode',
      'hybrid_vs_dense_k4',
      'hybrid_vs_flat_k4',
      'm2e_residual',
      'local_equiv_residual',
      'e2l_coeff_rms',
      'e2l_coeff_status',
      'hierarchy_inner_status',
      'k4_decomposition_status',
      'l2p_gradient_method')) {
    if (-not $text.Contains($needle)) {
      throw "E27 smoke support is not present yet: missing '$needle' in $SmokeSource"
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
  Assert-SmokeSupportsE27
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

function Decide([array]$Rows) {
  if ($Rows.Count -lt 9) { return 'fmm_e27_runtime_limited' }
  $fullAndLocal = @($Rows | Where-Object { $_.mode -eq 'full' -or $_.mode -eq 'hybrid-local-bridge' })
  foreach ($row in $fullAndLocal) {
    if ((Get-Field $row 'execution_status') -ne 'pass') { return 'fmm_e27_runtime_limited' }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0 -or
        (Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage') -or
        (As-Double $row 'far_fraction' 0.0) -le 0.0) {
      return 'fmm_e27_hybrid_structural_gap'
    }
  }
  if (-not (Mode-Pass $Rows 'hybrid-local-bridge')) {
    return 'fmm_e27_hybrid_local_bridge_gap'
  }
  $sourceRows = @($Rows | Where-Object { $_.mode -eq 'hybrid-source-local-bridge' })
  foreach ($row in $sourceRows) {
    if ((Get-Field $row 'execution_status') -ne 'pass') {
      return 'fmm_e27_hybrid_source_bridge_gap'
    }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0 -or
        (Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage') -or
        (As-Double $row 'far_fraction' 0.0) -le 0.0) {
      return 'fmm_e27_hybrid_structural_gap'
    }
  }
  if (-not (Mode-Pass $Rows 'hybrid-source-local-bridge')) {
    return 'fmm_e27_hybrid_source_bridge_gap'
  }
  'fmm_e27_hybrid_source_local_pass'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$modes = @('full', 'hybrid-local-bridge', 'hybrid-source-local-bridge')
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

$bridgeRows = $matrix | Select-Object row_id,fixture,sample,mode,hybrid_mode,hierarchy_inner_status,hierarchy_inner_accuracy_failure,m2e_residual,local_equiv_residual,e2l_coeff_rms,e2l_coeff_status,l2p_gradient_method,aggregate_rms,max_channel_rms,k4_rms,hybrid_vs_dense_k4,hybrid_vs_flat_k4,pure_hierarchy_vs_dense_k4,flat_vs_dense_k4,accuracy_status
$k4Rows = $matrix | Select-Object row_id,fixture,sample,mode,k4_decomposition_status,k4_rms,k4_abs_rms,k4_reference_rms,hybrid_vs_dense_k4,hybrid_vs_flat_k4,pure_hierarchy_vs_dense_k4,flat_vs_dense_k4,max_channel_name,max_channel_rms
$comparisonRows = $matrix | Select-Object row_id,fixture,sample,mode,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,far_fraction,direct_pairs,total_seconds,accuracy_status
Export-Table $bridgeRows $BridgeCsv
Export-Table $k4Rows $K4Csv
Export-Table $comparisonRows $ComparisonCsv

$decision = Decide $matrix
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'quick' })
  rows = $matrix.Count
  required_rows = 9
  full_pass = (Mode-Pass $matrix 'full')
  hybrid_local_pass = (Mode-Pass $matrix 'hybrid-local-bridge')
  hybrid_source_local_pass = (Mode-Pass $matrix 'hybrid-source-local-bridge')
  hybrid_source_inner_status = ((@($matrix | Where-Object { $_.mode -eq 'hybrid-source-local-bridge' } | ForEach-Object { Get-Field $_ 'hierarchy_inner_status' }) | Select-Object -Unique) -join ';')
  hybrid_local_k4_decomposition_status = ((@($matrix | Where-Object { $_.mode -eq 'hybrid-local-bridge' } | ForEach-Object { Get-Field $_ 'k4_decomposition_status' }) | Select-Object -Unique) -join ';')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$MatrixCsv"
)

$report = @(
  '# TITAN FMM Engine E27 Hybrid Bridge Real-Panel Report',
  '',
  "Decision: $decision",
  '',
  'E27 tests whether the E26 Cartesian-KiFMM bridge closes the E25.8 real-panel K4 local solve/L2P gap on sampled rows.',
  'E27.1 closeout wording: this is a sampled-row diagnostic harness only, not a production backend claim.'
)
if ($decision -eq 'fmm_e27_hybrid_local_bridge_gap') {
  $report += 'Interpretation: hybrid-local-bridge executes but fails the sampled K4 gate, so the local bridge does not close the real-panel gap in this integration.'
}
$sourceInner = @($matrix | Where-Object { $_.mode -eq 'hybrid-source-local-bridge' } | Select-Object -ExpandProperty hierarchy_inner_status -Unique)
if ($sourceInner.Count -gt 0) {
  $report += ('Hybrid source-local inner status observed: {0}' -f ($sourceInner -join ','))
}
foreach ($row in $matrix) {
  $report += ('  {0}: mode={1}, exec={2}, acc={3}, inner={4}, agg={5}, k4={6}, hybrid_vs_flat_k4={7}, flat_k4={8}, m2e={9}, local_equiv={10}, e2l_status={11}, k4_decomp={12}, l2p={13}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'mode'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'hierarchy_inner_status'), `
      (Get-Field $row 'aggregate_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'hybrid_vs_flat_k4'), (Get-Field $row 'flat_vs_dense_k4'), `
      (Get-Field $row 'm2e_residual'), (Get-Field $row 'local_equiv_residual'), `
      (Get-Field $row 'e2l_coeff_status'), (Get-Field $row 'k4_decomposition_status'), `
      (Get-Field $row 'l2p_gradient_method'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E27 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Bridge: $BridgeCsv"
Write-Host "K4: $K4Csv"
Write-Host "Comparison: $ComparisonCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @(
    'fmm_e27_hybrid_local_l2p_pass',
    'fmm_e27_hybrid_source_local_pass')) {
  throw "FMM-E27 decision failed: $decision"
}
