param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e258_hierarchy_attribution'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e258_matrix.csv'
$DebugCsv = Join-Path $WorkRoot 'fmm_e258_debug_modes.csv'
$FieldCsv = Join-Path $WorkRoot 'fmm_e258_field_attribution.csv'
$LocalRhsCsv = Join-Path $WorkRoot 'fmm_e258_local_check_rhs.csv'
$L2pCsv = Join-Path $WorkRoot 'fmm_e258_l2p_oracle.csv'
$K4Csv = Join-Path $WorkRoot 'fmm_e258_k4_decomposition.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e258_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e258_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E258_HIERARCHY_ATTRIBUTION_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E25.8 smoke.'
}

function Assert-SmokeSupportsE258 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'hier-debug-mode',
      'm2l-direct-to-targets',
      'local_check_rhs_relative_rms',
      'l2p_oracle_k4',
      'k4_laplace_grad_relative_rms',
      'E258_FIELD')) {
    if (-not $text.Contains($needle)) {
      throw "E25.8 smoke support is not present yet: missing '$needle' in $SmokeSource"
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
  Assert-SmokeSupportsE258
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
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
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
  Ensure-Field $Row 'order' 8
  Ensure-Field $Row 'theta' '0.6'
  Ensure-Field $Row 'real_panel_guard' 'strict-local-check'
  Ensure-Field $Row 'guard_margin_factor' '0'
  foreach ($field in @(
      'local_check_rhs_relative_rms',
      'local_check_rhs_max_abs_error',
      'l2p_oracle_aggregate',
      'l2p_oracle_k4',
      'k4_laplace_grad_relative_rms',
      'k4_yukawa_grad_relative_rms',
      'k4_combined_relative_rms',
      'flat_vs_dense_k4',
      'hier_vs_dense_k4',
      'hier_vs_flat_k4')) {
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
  if ($proc.ExitCode -ne 0 -or $resultLine.Count -ne 1) {
    $row = [pscustomobject][ordered]@{
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
    return [pscustomobject]@{ row = $row; fields = @() }
  }
  $row = Complete-Row (Parse-KeyValueLine ([string]$resultLine[-1]) 'E255_RESULT') $Case
  Ensure-Field $row 'process_exit_code' $proc.ExitCode
  Ensure-Field $row 'process_output' 'captured'
  $fields = @()
  foreach ($line in @($output | Where-Object { $_ -like 'E258_FIELD *' })) {
    $field = Parse-KeyValueLine ([string]$line) 'E258_FIELD'
    Ensure-Field $field 'row_id' $Case.id
    Ensure-Field $field 'fixture' $Case.fixture
    Ensure-Field $field 'sample' "$($Case.sample)_$($Case.sample_count)"
    Ensure-Field $field 'mode' $Case.mode
    $fields += $field
  }
  [pscustomobject]@{ row = $row; fields = $fields }
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
    if ((Get-Field $row 'accuracy_status') -ne 'pass') { return $false }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0) { return $false }
    if ((As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0) { return $false }
  }
  return $true
}

function Decide([array]$Rows) {
  $requiredRows = 18
  if ($Rows.Count -lt $requiredRows) { return 'fmm_e258_runtime_limited' }
  foreach ($row in $Rows) {
    if ((Get-Field $row 'execution_status') -ne 'pass') { return 'fmm_e258_runtime_limited' }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0) {
      return 'fmm_e258_worklist_regression_gap'
    }
    if ((Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage')) {
      return 'fmm_e258_worklist_regression_gap'
    }
  }
  if (-not (Mode-Pass $Rows 'm2l-direct-to-targets')) {
    return 'fmm_e258_m2m_upward_gap'
  }
  if (-not (Mode-Pass $Rows 'local-check-oracle')) {
    return 'fmm_e258_m2l_rhs_gap'
  }
  if (-not (Mode-Pass $Rows 'l2p-oracle')) {
    return 'fmm_e258_local_solve_l2p_gap'
  }
  if ((Mode-Pass $Rows 'no-l2l') -and -not (Mode-Pass $Rows 'full')) {
    return 'fmm_e258_l2l_propagation_gap'
  }
  if ((Mode-Pass $Rows 'no-m2m') -and -not (Mode-Pass $Rows 'full')) {
    return 'fmm_e258_m2m_upward_gap'
  }
  if (-not (Mode-Pass $Rows 'full')) {
    return 'fmm_e258_unresolved_hierarchy_gap'
  }
  'fmm_e258_hierarchy_debug_modes_pass'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$modes = @('full', 'm2l-direct-to-targets', 'no-m2m', 'no-l2l', 'local-check-oracle', 'l2p-oracle')
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
$fieldRows = @()
foreach ($case in $cases) {
  $old = $existing[$case.id]
  $needsRun = $Force -or $null -eq $old
  if ($RerunFailed -and $null -ne $old -and (Get-Field $old 'execution_status') -ne 'pass') {
    $needsRun = $true
  }
  if (-not $needsRun) {
    Write-Host "resume: $($case.id)"
    $matrix += $old
    continue
  }
  $result = Run-Case $case
  $matrix += $result.row
  $fieldRows += $result.fields
  Export-Table $matrix $MatrixCsv
  Export-Table $fieldRows $FieldCsv
}

$debugRows = $matrix | Select-Object row_id,fixture,sample,mode,aggregate_rms,max_channel_rms,k4_rms,flat_vs_dense_k4,hier_vs_dense_k4,hier_vs_flat_k4,m2m_translations,m2l_translations,l2l_translations,local_node_solves,l2p_target_count,accuracy_status
$localRows = @($fieldRows | Where-Object { (Get-Field $_ 'scope') -eq 'local-check-rhs' })
$l2pRows = $matrix | Select-Object row_id,fixture,sample,mode,l2p_oracle_aggregate,l2p_oracle_k4,accuracy_status
$k4Rows = $matrix | Select-Object row_id,fixture,sample,mode,k4_laplace_grad_relative_rms,k4_yukawa_grad_relative_rms,k4_combined_relative_rms,k4_abs_rms,k4_reference_rms
Export-Table $debugRows $DebugCsv
Export-Table $localRows $LocalRhsCsv
Export-Table $l2pRows $L2pCsv
Export-Table $k4Rows $K4Csv

$decision = Decide $matrix
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'quick' })
  rows = $matrix.Count
  required_rows = 18
  full_pass = (Mode-Pass $matrix 'full')
  m2l_direct_pass = (Mode-Pass $matrix 'm2l-direct-to-targets')
  no_m2m_pass = (Mode-Pass $matrix 'no-m2m')
  no_l2l_pass = (Mode-Pass $matrix 'no-l2l')
  local_check_oracle_pass = (Mode-Pass $matrix 'local-check-oracle')
  l2p_oracle_pass = (Mode-Pass $matrix 'l2p-oracle')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$MatrixCsv"
)

$report = @(
  '# TITAN FMM Engine E25.8 Hierarchy Attribution Report',
  '',
  "Decision: $decision",
  '',
  'E25.8 isolates the remaining hierarchy propagation/local-representation K4 gap after E25.7.'
)
foreach ($row in $matrix) {
  $report += ('  {0}: mode={1}, exec={2}, acc={3}, agg={4}, k4={5}, local_rhs={6}, l2p_k4={7}, k4_lap={8}, k4_yuk={9}, m2m={10}, m2l={11}, l2l={12}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'mode'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'aggregate_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'local_check_rhs_relative_rms'), (Get-Field $row 'l2p_oracle_k4'), `
      (Get-Field $row 'k4_laplace_grad_relative_rms'), (Get-Field $row 'k4_yukawa_grad_relative_rms'), `
      (Get-Field $row 'm2m_translations'), (Get-Field $row 'm2l_translations'), (Get-Field $row 'l2l_translations'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E25.8 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Debug modes: $DebugCsv"
Write-Host "Field attribution: $FieldCsv"
Write-Host "Local RHS: $LocalRhsCsv"
Write-Host "L2P oracle: $L2pCsv"
Write-Host "K4: $K4Csv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @(
    'fmm_e258_m2m_upward_gap',
    'fmm_e258_m2l_rhs_gap',
    'fmm_e258_local_solve_l2p_gap',
    'fmm_e258_l2l_propagation_gap',
    'fmm_e258_k4_recomposition_gap',
    'fmm_e258_unresolved_hierarchy_gap',
    'fmm_e258_hierarchy_debug_modes_pass')) {
  throw "FMM-E25.8 decision failed: $decision"
}
