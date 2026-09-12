param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e37_m2l_batching'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e37_matrix.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e37_m2l_runtime.csv'
$ParityCsv = Join-Path $WorkRoot 'fmm_e37_parity.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e37_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e37_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E37_M2L_BATCHING_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E37 smoke.'
}

function Assert-SmokeSupportsE37 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $smoke = Get-Content -Raw -LiteralPath $SmokeSource
  $headers = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.h')
  foreach ($needle in @(
      'm2l-eval-policy',
      'm2l_eval_policy=',
      'm2l_target_nodes=',
      'm2l_kernel_evals_estimated=')) {
    if (-not $smoke.Contains($needle) -and -not $headers.Contains($needle)) {
      throw "E37 M2L batching support is not present yet: missing '$needle'"
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
  Assert-SmokeSupportsE37
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

function New-Case([string]$Id, [string]$Fixture, [string]$Sample, [int]$SampleCount, [int]$SurfacePoints, [int]$PanelLimit, [string]$Policy) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    policy = $Policy
    row_id = "${Id}_${Policy}"
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
  Ensure-Field $Row 'requested_m2l_eval_policy' $Case.policy
  foreach ($field in @(
      'm2l_eval_policy',
      'm2l_seconds',
      'm2l_seconds_per_translation',
      'm2l_seconds_per_check_point',
      'm2l_target_nodes',
      'm2l_source_nodes',
      'm2l_check_points_total',
      'm2l_equivalent_points',
      'm2l_kernel_evals_estimated',
      'aggregate_rms',
      'max_channel_rms',
      'k1_rms',
      'k2_rms',
      'k3_rms',
      'k4_rms',
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
    "--m2l-eval-policy=$($Case.policy)",
    "--order=8",
    "--theta=0.6",
    "--equiv-scale=1.05",
    "--upward-check-scale=3.0",
    "--local-equiv-scale=3.0",
    "--local-check-scale=1.05",
    "--condition-policy=qr-rank-only",
    "--real-panel-guard=strict-local-check",
    "--guard-margin-factor=0",
    "--pqr-1sh1=$Pqr1sh1",
    "--pqr-1a63=$Pqr1a63"
  )
  $stdoutPath = Join-Path $WorkRoot "$($Case.row_id).stdout.txt"
  $stderrPath = Join-Path $WorkRoot "$($Case.row_id).stderr.txt"
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
    $failedRow = [pscustomObject][ordered]@{
      row_id = $Case.row_id
      case_id = $Case.id
      fixture = $Case.fixture
      requested_m2l_eval_policy = $Case.policy
      m2l_eval_policy = $Case.policy
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
  Ensure-Field $row 'status_detail' 'completed'
  $row
}

function Build-ParityRows([array]$MatrixRows) {
  $rows = @()
  foreach ($caseId in @($MatrixRows | ForEach-Object { Get-Field $_ 'case_id' } | Sort-Object -Unique)) {
    $scalar = @($MatrixRows | Where-Object { (Get-Field $_ 'case_id') -eq $caseId -and (Get-Field $_ 'requested_m2l_eval_policy') -eq 'scalar-existing' } | Select-Object -First 1)
    $batched = @($MatrixRows | Where-Object { (Get-Field $_ 'case_id') -eq $caseId -and (Get-Field $_ 'requested_m2l_eval_policy') -eq 'target-check-batched' } | Select-Object -First 1)
    if ($scalar.Count -ne 1 -or $batched.Count -ne 1) { continue }
    $scalarM2l = As-Double $scalar[0] 'm2l_seconds'
    $batchedM2l = As-Double $batched[0] 'm2l_seconds'
    $scalarTotal = As-Double $scalar[0] 'total_seconds'
    $batchedTotal = As-Double $batched[0] 'total_seconds'
    $speedup = if ($batchedM2l -gt 0.0) { $scalarM2l / $batchedM2l } else { [double]::NaN }
    $totalSpeedup = if ($batchedTotal -gt 0.0) { $scalarTotal / $batchedTotal } else { [double]::NaN }
    $channelDiff = [Math]::Max(
      [Math]::Max([Math]::Abs((As-Double $scalar[0] 'k1_rms') - (As-Double $batched[0] 'k1_rms')),
                  [Math]::Abs((As-Double $scalar[0] 'k2_rms') - (As-Double $batched[0] 'k2_rms'))),
      [Math]::Max([Math]::Abs((As-Double $scalar[0] 'k3_rms') - (As-Double $batched[0] 'k3_rms')),
                  [Math]::Abs((As-Double $scalar[0] 'k4_rms') - (As-Double $batched[0] 'k4_rms'))))
    $rows += [pscustomobject][ordered]@{
      case_id = $caseId
      scalar_output_hash = Get-Field $scalar[0] 'output_hash'
      batched_output_hash = Get-Field $batched[0] 'output_hash'
      scalar_policy_hash = Get-Field $scalar[0] 'policy_hash'
      batched_policy_hash = Get-Field $batched[0] 'policy_hash'
      scalar_m2l_seconds = $scalarM2l
      batched_m2l_seconds = $batchedM2l
      m2l_speedup = $speedup
      scalar_total_seconds = $scalarTotal
      batched_total_seconds = $batchedTotal
      total_speedup = $totalSpeedup
      aggregate_rms_diff = [Math]::Abs((As-Double $scalar[0] 'aggregate_rms') - (As-Double $batched[0] 'aggregate_rms'))
      max_channel_rms_diff = [Math]::Abs((As-Double $scalar[0] 'max_channel_rms') - (As-Double $batched[0] 'max_channel_rms'))
      channel_rms_max_diff = $channelDiff
      parity_status = $(if ($channelDiff -le 1.0e-10) { 'pass' } else { 'fail' })
    }
  }
  $rows
}

function Decide([array]$MatrixRows, [array]$ParityRows, [int]$RequiredRows) {
  if ($MatrixRows.Count -lt $RequiredRows) { return 'fmm_e37_runtime_limited' }
  if (@($MatrixRows | Where-Object { (Get-Field $_ 'status_detail') -eq 'timeout' -or (Get-Field $_ 'accuracy_failure') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e37_runtime_limited'
  }
  if (@($MatrixRows | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e37_m2l_accuracy_gap'
  }
  if (@($MatrixRows | Where-Object {
        (Get-Field $_ 'accuracy_status') -ne 'pass' -or
        (As-Double $_ 'aggregate_rms') -gt 1.0e-6 -or
        (As-Double $_ 'max_channel_rms') -gt 5.0e-6 -or
        (Get-Field $_ 'coverage') -ne (Get-Field $_ 'full_coverage_expected') -or
        (Get-Field $_ 'sampled_dense_valid_pairs') -ne (Get-Field $_ 'sampled_coverage_expected') -or
        (As-Int64 $_ 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $_ 'nonproduction_fallback_pairs' 1) -ne 0
      }).Count -gt 0) {
    return 'fmm_e37_m2l_accuracy_gap'
  }
  if ($ParityRows.Count -eq 0 -or @($ParityRows | Where-Object { (Get-Field $_ 'parity_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e37_m2l_parity_gap'
  }
  $driver = @($ParityRows | Where-Object { (Get-Field $_ 'case_id') -eq '1SH1_s4_p512_uniform64' } | Select-Object -First 1)
  if ($driver.Count -eq 1 -and (As-Double $driver[0] 'm2l_speedup' 0.0) -ge 1.5 -and (As-Double $driver[0] 'total_speedup' 0.0) -ge 1.25) {
    return 'fmm_e37_m2l_batched_pass'
  }
  return 'fmm_e37_accuracy_pass_partial_speedup'
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
  foreach ($policy in @('scalar-existing', 'target-check-batched')) {
    $cases += New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $spec[5] $policy
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
  $old = $existing[$case.row_id]
  $needsRun = $Force -or $null -eq $old
  if ($RerunFailed -and $null -ne $old -and (Get-Field $old 'execution_status') -ne 'pass') {
    $needsRun = $true
  }
  if (-not $needsRun) {
    Write-Host "resume: $($case.row_id)"
    $matrix += (Complete-Row $old $case)
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

$parity = Build-ParityRows $matrix
$decision = Decide $matrix $parity $cases.Count

Export-Table $matrix $MatrixCsv
Export-Table ($matrix | Select-Object row_id,case_id,fixture,sample,m2l_eval_policy,m2l_seconds,m2l_seconds_per_translation,m2l_seconds_per_check_point,m2l_target_nodes,m2l_source_nodes,m2l_translations,m2l_check_points_total,m2l_equivalent_points,m2l_kernel_evals_estimated,total_seconds) $RuntimeCsv
Export-Table $parity $ParityCsv

[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'full' })
  rows = $matrix.Count
  required_rows = $cases.Count
  parity_rows = $parity.Count
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($matrix.Count)",
  "required_rows=$($cases.Count)",
  "matrix=$MatrixCsv",
  "runtime=$RuntimeCsv",
  "parity=$ParityCsv"
)

$report = @(
  '# TITAN FMM Engine E37 M2L Batching Report',
  '',
  "Decision: $decision",
  "Mode: $(if ($Quick) { 'quick' } else { 'full' })",
  "Rows: $($matrix.Count)/$($cases.Count)",
  '',
  'E37 compares the existing scalar M2L RHS accumulation path with an opt-in target-check-batched double-precision path. It does not change q-order, theta, tree policy, GMRES, preconditioners, production TABI, manuscript, or defaults.',
  ''
)
foreach ($row in $matrix) {
  $report += ('  {0}: policy={1}, exec={2}, acc={3}, m2l={4}s, per_translation={5}, per_check={6}, total={7}s, aggregate={8}, max_channel={9}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'm2l_eval_policy'), `
      (Get-Field $row 'execution_status'), (Get-Field $row 'accuracy_status'), `
      (Get-Field $row 'm2l_seconds'), (Get-Field $row 'm2l_seconds_per_translation'), `
      (Get-Field $row 'm2l_seconds_per_check_point'), (Get-Field $row 'total_seconds'), `
      (Get-Field $row 'aggregate_rms'), (Get-Field $row 'max_channel_rms'))
}
foreach ($row in $parity) {
  $report += ('  parity {0}: status={1}, m2l_speedup={2}, total_speedup={3}, channel_diff={4}' -f `
      (Get-Field $row 'case_id'), (Get-Field $row 'parity_status'), `
      (Get-Field $row 'm2l_speedup'), (Get-Field $row 'total_speedup'), `
      (Get-Field $row 'channel_rms_max_diff'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E37 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Runtime: $RuntimeCsv"
Write-Host "Parity: $ParityCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e37_m2l_batched_pass','fmm_e37_accuracy_pass_partial_speedup')) {
  throw "FMM-E37 decision failed: $decision"
}
