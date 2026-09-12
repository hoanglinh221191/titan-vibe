param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 600,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e36_hierarchy_runtime'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e36_matrix.csv'
$StageCsv = Join-Path $WorkRoot 'fmm_e36_stage_timing.csv'
$CompressionCsv = Join-Path $WorkRoot 'fmm_e36_compression.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e36_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e36_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E36_HIERARCHY_RUNTIME_REPORT.txt'
$StageRows = @()

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
  throw 'No C++ compiler found for FMM-E36 smoke.'
}

function Assert-SmokeSupportsE36 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'emit-e36-stage-attribution',
      'p2m_seconds=',
      'm2m_seconds=',
      'm2l_seconds=',
      'local_solve_seconds=',
      'l2l_seconds=',
      'l2p_seconds=',
      'E36_STAGE')) {
    if (-not $text.Contains($needle)) {
      throw "E36 hierarchy runtime attribution support is not present yet: missing '$needle' in $SmokeSource"
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
  try {
    [double]$value
  } catch {
    return $Default
  }
}

function As-Int64([object]$Row, [string]$Name, [int64]$Default = [int64]::MinValue) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  try {
    [int64]$value
  } catch {
    return $Default
  }
}

function Field-IsFinite([object]$Row, [string]$Name) {
  $value = As-Double $Row $Name
  return -not [double]::IsNaN($value) -and -not [double]::IsInfinity($value)
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
  Assert-SmokeSupportsE36
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
  foreach ($field in @(
      'timing_scope',
      'execution_status',
      'accuracy_status',
      'compression_status',
      'full_apply_stage_status',
      'kifmm_compression_ratio',
      'cartesian_k4_compression_ratio',
      'p2m_seconds',
      'm2m_seconds',
      'm2l_seconds',
      'local_solve_seconds',
      'l2l_seconds',
      'l2p_seconds',
      'stage_timing_sum',
      'hierarchy_bottleneck_stage',
      'hierarchy_bottleneck_fraction')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $Row
}

function New-Case([string]$Id, [string]$Fixture, [string]$Sample, [int]$SampleCount, [int]$SurfacePoints, [int]$PanelLimit) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
  }
}

function Read-StageRowsFromStdout([string]$Path, [string]$RowId) {
  if (-not (Test-Path -LiteralPath $Path)) { return }
  foreach ($stageLine in @(Get-Content -LiteralPath $Path | Where-Object { $_ -like 'E36_STAGE *' })) {
    $stage = Parse-KeyValueLine ([string]$stageLine) 'E36_STAGE '
    Ensure-Field $stage 'row_id' $RowId
    $script:StageRows += $stage
  }
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.id) (timeout ${CaseTimeoutSeconds}s)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--surface-points=$($Case.surface_points)",
    "--panel-limit=$($Case.panel_limit)",
    "--mode=full-panel-profile",
    "--emit-e36-stage-attribution=1",
    "--hier-debug-mode=dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian",
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
  Read-StageRowsFromStdout $stdoutPath $Case.id
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($resultLine.Count -ne 1) {
    $failedRow = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      surface_points = $Case.surface_points
      panel_limit = $Case.panel_limit
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
  Ensure-Field $row 'status_detail' 'completed'
  $row
}

function Largest-Row([array]$Rows) {
  @($Rows | Sort-Object @{ Expression = { As-Int64 $_ 'target_panels' 0 }; Ascending = $true },
                         @{ Expression = { As-Int64 $_ 'sampled_rows' 0 }; Ascending = $true }) | Select-Object -Last 1
}

function Decide([array]$Rows, [int]$RequiredRows, [int]$StageRowCount) {
  if ($Rows.Count -lt $RequiredRows) { return 'fmm_e36_runtime_limited' }
  if (@($Rows | Where-Object { (Get-Field $_ 'status_detail') -eq 'timeout' -or (Get-Field $_ 'accuracy_failure') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e36_runtime_limited'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e36_execution_gap'
  }
  if (@($Rows | Where-Object {
        (Get-Field $_ 'coverage') -ne (Get-Field $_ 'full_coverage_expected') -or
        (Get-Field $_ 'sampled_dense_valid_pairs') -ne (Get-Field $_ 'sampled_coverage_expected') -or
        (As-Int64 $_ 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $_ 'nonproduction_fallback_pairs' 1) -ne 0
      }).Count -gt 0) {
    return 'fmm_e36_coverage_gap'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'accuracy_status') -ne 'pass' -or (As-Double $_ 'aggregate_rms') -gt 1.0e-6 -or (As-Double $_ 'max_channel_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e36_accuracy_gap'
  }
  if ($StageRowCount -lt $Rows.Count) {
    return 'fmm_e36_stage_attribution_gap'
  }
  if (@($Rows | Where-Object {
        (Get-Field $_ 'timing_scope') -ne 'full_panel_apply_sampled_audit' -or
        (Get-Field $_ 'full_apply_stage_status') -ne 'pass' -or
        -not (Field-IsFinite $_ 'p2m_seconds') -or
        -not (Field-IsFinite $_ 'm2m_seconds') -or
        -not (Field-IsFinite $_ 'm2l_seconds') -or
        -not (Field-IsFinite $_ 'local_solve_seconds') -or
        -not (Field-IsFinite $_ 'l2l_seconds') -or
        -not (Field-IsFinite $_ 'l2p_seconds') -or
        -not (Field-IsFinite $_ 'stage_timing_sum') -or
        (Get-Field $_ 'hierarchy_bottleneck_stage') -eq 'MISSING'
      }).Count -gt 0) {
    return 'fmm_e36_stage_attribution_gap'
  }
  $driver = Largest-Row $Rows
  if ((As-Double $driver 'far_fraction' 0.0) -lt 0.10) {
    return 'fmm_e36_fallback_dominated'
  }
  $stage = Get-Field $driver 'hierarchy_bottleneck_stage'
  $fraction = As-Double $driver 'hierarchy_bottleneck_fraction' 0.0
  if ($fraction -ge 0.50) {
    if ($stage -eq 'l2p') { return 'fmm_e36_l2p_dominated' }
    if ($stage -eq 'm2l') { return 'fmm_e36_m2l_dominated' }
    if ($stage -eq 'local_solve') { return 'fmm_e36_local_solve_dominated' }
    if ($stage -eq 'p2m') { return 'fmm_e36_p2m_dominated' }
    if ($stage -eq 'hierarchy_direct_fallback') { return 'fmm_e36_fallback_dominated' }
  }
  return 'fmm_e36_stage_attribution_pass'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$cases = @(
  (New-Case '1SH1_s4_p128_uniform32' '1SH1' 'uniform_hash' 32 4 128),
  (New-Case '1SH1_s4_p256_uniform32' '1SH1' 'uniform_hash' 32 4 256),
  (New-Case '1SH1_s4_p512_uniform64' '1SH1' 'uniform_hash' 64 4 512)
)
if (-not $Quick) {
  $cases += @(
    (New-Case '1SH1_s8_full_uniform128' '1SH1' 'uniform_hash' 128 8 0),
    (New-Case '1SH1_s8_full_closest64' '1SH1' 'closest_opposing_surface' 64 8 0),
    (New-Case '1SH1_s8_full_high_aspect64' '1SH1' 'high_aspect' 64 8 0)
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
    Read-StageRowsFromStdout (Join-Path $WorkRoot "$($case.id).stdout.txt") $case.id
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

Export-Table $matrix $MatrixCsv
Export-Table $StageRows $StageCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,surface_points,panel_limit,source_panels,target_panels,far_panel_pairs,direct_panel_pairs,far_fraction,kifmm_compression_ratio,cartesian_k4_compression_ratio,compression_status,cartesian_k4_far_pairs,cartesian_k4_fallback_pairs) $CompressionCsv

$decision = Decide $matrix $cases.Count $StageRows.Count
$driver = Largest-Row $matrix
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'full' })
  rows = $matrix.Count
  required_rows = $cases.Count
  driver_row = Get-Field $driver 'row_id'
  driver_bottleneck_stage = Get-Field $driver 'hierarchy_bottleneck_stage'
  driver_bottleneck_fraction = Get-Field $driver 'hierarchy_bottleneck_fraction'
  driver_far_fraction = Get-Field $driver 'far_fraction'
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($matrix.Count)",
  "required_rows=$($cases.Count)",
  "driver_row=$(Get-Field $driver 'row_id')",
  "driver_bottleneck_stage=$(Get-Field $driver 'hierarchy_bottleneck_stage')",
  "driver_bottleneck_fraction=$(Get-Field $driver 'hierarchy_bottleneck_fraction')",
  "matrix=$MatrixCsv",
  "stage=$StageCsv",
  "compression=$CompressionCsv"
)

$report = @(
  '# TITAN FMM Engine E36 Hierarchy Runtime Attribution Report',
  '',
  "Decision: $decision",
  "Mode: $(if ($Quick) { 'quick' } else { 'full' })",
  "Rows: $($matrix.Count)/$($cases.Count)",
  '',
  'E36 profiles the existing full-panel dual-worklist/channel-hybrid q8 path and attributes hierarchy runtime by stage. It does not change math, GMRES, preconditioners, production TABI, manuscript, or defaults.',
  ('Driver row: {0}, bottleneck={1}, fraction={2}, far_fraction={3}' -f `
      (Get-Field $driver 'row_id'), (Get-Field $driver 'hierarchy_bottleneck_stage'), `
      (Get-Field $driver 'hierarchy_bottleneck_fraction'), (Get-Field $driver 'far_fraction'))
)
foreach ($row in $matrix) {
  $report += ('  {0}: exec={1}, acc={2}, panels={3}, far={4}, bottleneck={5}/{6}, hierarchy={7}, p2m={8}, m2m={9}, m2l={10}, local={11}, l2l={12}, l2p={13}, total={14}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'execution_status'), (Get-Field $row 'accuracy_status'), `
      (Get-Field $row 'target_panels'), (Get-Field $row 'far_fraction'), `
      (Get-Field $row 'hierarchy_bottleneck_stage'), (Get-Field $row 'hierarchy_bottleneck_fraction'), `
      (Get-Field $row 'hierarchy_seconds'), (Get-Field $row 'p2m_seconds'), `
      (Get-Field $row 'm2m_seconds'), (Get-Field $row 'm2l_seconds'), `
      (Get-Field $row 'local_solve_seconds'), (Get-Field $row 'l2l_seconds'), `
      (Get-Field $row 'l2p_seconds'), (Get-Field $row 'total_seconds'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E36 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Stage timing: $StageCsv"
Write-Host "Compression: $CompressionCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e36_stage_attribution_pass','fmm_e36_l2p_dominated','fmm_e36_m2l_dominated','fmm_e36_local_solve_dominated','fmm_e36_p2m_dominated','fmm_e36_fallback_dominated')) {
  throw "FMM-E36 decision failed: $decision"
}
