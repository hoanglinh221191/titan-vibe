param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e42_grouped_gemm'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e42_matrix.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e42_runtime.csv'
$GroupedCsv = Join-Path $WorkRoot 'fmm_e42_grouped_gemm.csv'
$ParityCsv = Join-Path $WorkRoot 'fmm_e42_parity.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e42_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e42_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E42_GROUPED_GEMM_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E42 smoke.'
}

function Assert-SmokeSupportsE42 {
  $smoke = Get-Content -Raw -LiteralPath $SmokeSource
  $header = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy.h')
  $jcp = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp')
  foreach ($needle in @('m2l-gemm-batch-policy', 'group-by-cache-key', 'm2l_gemm_batch_policy=')) {
    if (-not ($smoke.Contains($needle) -or $header.Contains($needle) -or $jcp.Contains($needle))) {
      throw "E42 grouped-GEMM support is not present yet: missing '$needle'"
    }
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE42
  New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
  $Compiler = Resolve-Cxx $Cxx
  $CompileArgs = @(
    '-std=c++23', '-O2', '-mavx2', '-mfma', '-Isrc', '-Isrc\eigen',
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
    '-o', $ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) { throw 'failed to build titan_fmm_kifmm_hierarchy_real_panel_smoke' }
  } finally {
    Pop-Location
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
  if ($null -eq $prop -or $null -eq $prop.Value -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) { return $Default }
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
  if ($Rows.Count -eq 0) { Set-Content -LiteralPath $Path -Value '' }
  else { $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path }
}

function New-Case(
  [string]$Id,
  [string]$Sample,
  [int]$SampleCount,
  [int]$SurfacePoints,
  [int]$PanelLimit,
  [string]$GemmPolicy
) {
  [pscustomobject]@{
    id = $Id
    fixture = '1SH1'
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    m2l_gemm_batch_policy = $GemmPolicy
    row_id = "${Id}_relative-cache_${GemmPolicy}_panel-radius-occupancy_r1p5_occ32".Replace('.', 'p')
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
  Ensure-Field $Row 'requested_m2l_gemm_batch_policy' $Case.m2l_gemm_batch_policy
  foreach ($field in @(
      'm2l_gemm_batch_policy', 'm2l_gemm_group_count',
      'm2l_gemm_grouped_pair_count', 'm2l_gemm_grouped_rhs_columns',
      'm2l_gemm_max_group_size', 'm2l_matrix_cache_policy',
      'm2l_matrix_cache_hit_rate', 'm2l_matrix_cache_entries',
      'm2l_matrix_cache_hits', 'm2l_matrix_cache_misses',
      'm2l_matrix_gemm_seconds', 'm2l_matrix_cache_build_seconds',
      'tree_stop_policy', 'tree_panel_radius_factor', 'tree_occupancy_stop',
      'm2l_eval_policy', 'm2l_seconds', 'm2l_translations',
      'aggregate_rms', 'max_channel_rms', 'k1_rms', 'k2_rms', 'k3_rms',
      'k4_rms', 'coverage', 'full_coverage_expected',
      'production_invalid_pairs', 'nonproduction_fallback_pairs',
      'far_fraction', 'far_panel_pairs', 'total_seconds', 'output_hash',
      'policy_hash', 'worklist_hash')) {
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
    '--m2l-matrix-cache-policy=relative-geometry-cache',
    "--m2l-gemm-batch-policy=$($Case.m2l_gemm_batch_policy)",
    '--tree-stop-policy=panel-radius-occupancy',
    '--tree-panel-radius-factor=1.5',
    '--tree-occupancy-stop=32',
    '--order=8', '--theta=0.6',
    '--equiv-scale=1.05', '--upward-check-scale=3.0',
    '--local-equiv-scale=3.0', '--local-check-scale=1.05',
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
    return Complete-Row ([pscustomobject][ordered]@{
      row_id = $Case.row_id
      case_id = $Case.id
      status = 'timeout'
      execution_status = 'timeout'
      stdout = $stdout
      stderr = $stderr
    }) $Case
  }
  $output = @()
  if (Test-Path -LiteralPath $stdout) { $output += Get-Content -LiteralPath $stdout }
  if (Test-Path -LiteralPath $stderr) { $output += Get-Content -LiteralPath $stderr }
  $line = @($output | Where-Object { $_.StartsWith('E255_RESULT ') } | Select-Object -Last 1)
  if ([string]::IsNullOrWhiteSpace($line)) {
    return Complete-Row ([pscustomobject][ordered]@{
      row_id = $Case.row_id
      case_id = $Case.id
      status = 'process_failed'
      execution_status = 'process_failed'
      stdout = $stdout
      stderr = $stderr
      process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }) $Case
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

function Case-Row([array]$Rows, [string]$CaseId, [string]$Policy) {
  return ,@($Rows | Where-Object { $_.case_id -eq $CaseId -and (Get-Field $_ 'requested_m2l_gemm_batch_policy') -eq $Policy } | Select-Object -First 1)
}

function Build-ParityRows([array]$Rows, [array]$CaseIds) {
  $parityRows = @()
  foreach ($caseId in $CaseIds) {
    $base = Case-Row $Rows $caseId 'per-pair'
    $grouped = Case-Row $Rows $caseId 'group-by-cache-key'
    if ($base.Count -ne 1 -or $grouped.Count -ne 1) { continue }
    $maxDelta = 0.0
    foreach ($field in @('aggregate_rms', 'k1_rms', 'k2_rms', 'k3_rms', 'k4_rms', 'max_channel_rms')) {
      $maxDelta = [Math]::Max($maxDelta, [Math]::Abs((As-Double $base[0] $field 0.0) - (As-Double $grouped[0] $field 0.0)))
    }
    $groupedM2l = As-Double $grouped[0] 'm2l_seconds' 0.0
    $groupedTotal = As-Double $grouped[0] 'total_seconds' 0.0
    $parityRows += [pscustomobject][ordered]@{
      case_id = $caseId
      baseline_row = Get-Field $base[0] 'row_id'
      grouped_row = Get-Field $grouped[0] 'row_id'
      max_channel_rms_delta = $maxDelta
      parity_status = if ($maxDelta -le 1.0e-12) { 'pass' } else { 'fail' }
      grouped_group_count = As-Int64 $grouped[0] 'm2l_gemm_group_count' 0
      grouped_pair_count = As-Int64 $grouped[0] 'm2l_gemm_grouped_pair_count' 0
      grouped_rhs_columns = As-Int64 $grouped[0] 'm2l_gemm_grouped_rhs_columns' 0
      grouped_max_group_size = As-Int64 $grouped[0] 'm2l_gemm_max_group_size' 0
      baseline_m2l_seconds = As-Double $base[0] 'm2l_seconds' 0.0
      grouped_m2l_seconds = $groupedM2l
      m2l_speedup = if ($groupedM2l -gt 0.0) { (As-Double $base[0] 'm2l_seconds' 0.0) / $groupedM2l } else { 0.0 }
      baseline_total_seconds = As-Double $base[0] 'total_seconds' 0.0
      grouped_total_seconds = $groupedTotal
      total_speedup = if ($groupedTotal -gt 0.0) { (As-Double $base[0] 'total_seconds' 0.0) / $groupedTotal } else { 0.0 }
    }
  }
  $parityRows
}

function Decide-E42([array]$Rows, [array]$ParityRows, [bool]$IsQuick) {
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -eq 'timeout' }).Count -gt 0) { return 'fmm_e42_runtime_limited' }
  if (@($Rows | Where-Object { -not (Is-Accurate $_) -or -not (Is-Covered $_) }).Count -gt 0) { return 'fmm_e42_accuracy_gap' }
  if (@($ParityRows | Where-Object { $_.parity_status -ne 'pass' }).Count -gt 0) { return 'fmm_e42_grouped_gemm_parity_gap' }
  if (@($ParityRows | Where-Object { [int64]$_.grouped_group_count -le 0 -or [int64]$_.grouped_max_group_size -lt 2 }).Count -gt 0) { return 'fmm_e42_grouping_coverage_gap' }
  if (@($ParityRows | Where-Object { [double]$_.m2l_speedup -lt 1.10 }).Count -gt 0) { return 'fmm_e42_accuracy_pass_partial_speedup' }
  if ($IsQuick) { return 'fmm_e42_grouped_gemm_quick_pass_full_pending' }
  $fullUniform = @($ParityRows | Where-Object { $_.case_id -eq '1SH1_s8_full_uniform128' } | Select-Object -First 1)
  if ($fullUniform.Count -eq 1 -and [double]$fullUniform[0].m2l_speedup -ge 1.20 -and [double]$fullUniform[0].total_speedup -ge 1.05) {
    return 'fmm_e42_grouped_gemm_pass'
  }
  'fmm_e42_accuracy_pass_partial_speedup'
}

function Update-DerivedCsvs([array]$Rows, [array]$CaseIds) {
  $runtimeRows = @()
  $groupRows = @()
  foreach ($row in $Rows) {
    $runtimeRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      m2l_gemm_batch_policy = Get-Field $row 'm2l_gemm_batch_policy'
      m2l_seconds = Get-Field $row 'm2l_seconds'
      m2l_matrix_gemm_seconds = Get-Field $row 'm2l_matrix_gemm_seconds'
      total_seconds = Get-Field $row 'total_seconds'
      aggregate_rms = Get-Field $row 'aggregate_rms'
      max_channel_rms = Get-Field $row 'max_channel_rms'
      far_fraction = Get-Field $row 'far_fraction'
    }
    $groupRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      m2l_gemm_batch_policy = Get-Field $row 'm2l_gemm_batch_policy'
      m2l_gemm_group_count = Get-Field $row 'm2l_gemm_group_count'
      m2l_gemm_grouped_pair_count = Get-Field $row 'm2l_gemm_grouped_pair_count'
      m2l_gemm_grouped_rhs_columns = Get-Field $row 'm2l_gemm_grouped_rhs_columns'
      m2l_gemm_max_group_size = Get-Field $row 'm2l_gemm_max_group_size'
      m2l_matrix_cache_hit_rate = Get-Field $row 'm2l_matrix_cache_hit_rate'
      m2l_matrix_cache_entries = Get-Field $row 'm2l_matrix_cache_entries'
      m2l_matrix_cache_hits = Get-Field $row 'm2l_matrix_cache_hits'
      m2l_matrix_cache_misses = Get-Field $row 'm2l_matrix_cache_misses'
    }
  }
  $parityRows = Build-ParityRows $Rows $CaseIds
  Export-Table $Rows $MatrixCsv
  Export-Table $runtimeRows $RuntimeCsv
  Export-Table $groupRows $GroupedCsv
  Export-Table $parityRows $ParityCsv
  $parityRows
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
if ($Force) {
  Remove-Item -LiteralPath $MatrixCsv, $RuntimeCsv, $GroupedCsv, $ParityCsv, $DecisionCsv, $SummaryDat, $ReportPath -Force -ErrorAction SilentlyContinue
}
Build-Smoke

$quickSpecs = @(
  @('1SH1_s4_p512_uniform64', 'uniform_hash', 64, 4, 512),
  @('1SH1_s8_p512_uniform64', 'uniform_hash', 64, 8, 512)
)
$fullSpecs = @(
  @('1SH1_s8_full_uniform128', 'uniform_hash', 128, 8, 0),
  @('1SH1_s8_full_closest64', 'closest_opposing_surface', 64, 8, 0),
  @('1SH1_s8_full_high_aspect64', 'high_aspect', 64, 8, 0)
)
$caseSpecs = @($quickSpecs)
if (-not $Quick) { $caseSpecs += $fullSpecs }
$caseIds = @($caseSpecs | ForEach-Object { $_[0] })

$rows = @()
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  $rows = @(Import-Csv -LiteralPath $MatrixCsv)
}

foreach ($spec in $caseSpecs) {
  foreach ($policy in @('per-pair', 'group-by-cache-key')) {
    $case = New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $policy
    $existing = @($rows | Where-Object { $_.row_id -eq $case.row_id } | Select-Object -First 1)
    $shouldRun = $true
    if ($existing.Count -eq 1 -and -not $Force) {
      $shouldRun = $false
      if ($RerunFailed -and (Get-Field $existing[0] 'execution_status') -ne 'pass') { $shouldRun = $true }
    }
    if ($shouldRun) {
      $newRow = Run-Case $case
      $rows = @($rows | Where-Object { $_.row_id -ne $case.row_id })
      $rows += $newRow
      Export-Table $rows $MatrixCsv
    } else {
      Write-Host "resume: $($case.row_id)"
    }
  }
}

$parityRows = Update-DerivedCsvs $rows $caseIds
$decision = Decide-E42 $rows $parityRows ([bool]$Quick)

[pscustomobject][ordered]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  generated_at = (Get-Date).ToString('s')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

@"
FMM-E42 grouped GEMM gate
decision=$decision
mode=$(if ($Quick) { 'quick' } else { 'full' })
matrix=$MatrixCsv
runtime=$RuntimeCsv
grouped=$GroupedCsv
parity=$ParityCsv
"@ | Set-Content -LiteralPath $SummaryDat

@"
TITAN FMM-E42 GROUP-BY-CACHE-KEY GEMM REPORT

Decision: $decision

E42 tests an opt-in grouped GEMM path for target-check-batched KIFMM M2L.
It keeps the E39 tree-stop candidate, relative-geometry matrix cache, q8,
strict-local-check, sampled dense production-direct audit, and all production
defaults unchanged.

Outputs:
- Matrix: $MatrixCsv
- Runtime: $RuntimeCsv
- Group telemetry: $GroupedCsv
- Grouped-vs-per-pair parity: $ParityCsv
- Decision: $DecisionCsv
"@ | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E42 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Grouped: $GroupedCsv"
Write-Host "Parity: $ParityCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."
if ($decision -notin @('fmm_e42_grouped_gemm_pass', 'fmm_e42_grouped_gemm_quick_pass_full_pending', 'fmm_e42_accuracy_pass_partial_speedup') -and -not $AllowFailedDecision) {
  throw "FMM-E42 gate did not pass: $decision"
}
