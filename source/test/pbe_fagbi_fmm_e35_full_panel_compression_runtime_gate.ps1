param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 300,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e35_compression_runtime'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e35_matrix.csv'
$CompressionCsv = Join-Path $WorkRoot 'fmm_e35_compression.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e35_runtime.csv'
$AttributionCsv = Join-Path $WorkRoot 'fmm_e35_attribution.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e35_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e35_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E35_COMPRESSION_RUNTIME_REPORT.txt'
$AttributionRows = @()

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
  throw 'No C++ compiler found for FMM-E35 smoke.'
}

function Assert-SmokeSupportsE35 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'emit-e35-attribution',
      'compression_status=',
      'kifmm_compression_ratio=',
      'cartesian_k4_compression_ratio=',
      'full_panel_apply_sampled_audit')) {
    if (-not $text.Contains($needle)) {
      throw "E35 compression/runtime attribution support is not present yet: missing '$needle' in $SmokeSource"
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
  [int64]$value
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE35
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
  Ensure-Field $Row 'timing_scope' 'MISSING'
  Ensure-Field $Row 'compression_status' 'MISSING'
  Ensure-Field $Row 'full_apply_stage_status' 'MISSING'
  Ensure-Field $Row 'kifmm_compression_ratio' 'MISSING'
  Ensure-Field $Row 'cartesian_k4_compression_ratio' 'MISSING'
  Ensure-Field $Row 'cartesian_k4_accuracy_status' 'MISSING'
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
    "--mode=full-panel-profile",
    "--emit-e35-attribution=1",
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
  foreach ($attrLine in @($output | Where-Object { $_ -like 'E35_ATTR *' })) {
    $attr = Parse-KeyValueLine ([string]$attrLine) 'E35_ATTR '
    Ensure-Field $attr 'row_id' $Case.id
    $script:AttributionRows += $attr
  }
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

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) {
    Set-Content -LiteralPath $Path -Value ''
  } else {
    $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path
  }
}

function Field-IsFinite([object]$Row, [string]$Name) {
  $value = As-Double $Row $Name
  return -not [double]::IsNaN($value) -and -not [double]::IsInfinity($value)
}

function Decide([array]$Rows, [int]$RequiredRows) {
  if ($Rows.Count -lt $RequiredRows) { return 'fmm_e35_runtime_limited' }
  if (@($Rows | Where-Object { (Get-Field $_ 'status_detail') -eq 'timeout' -or (Get-Field $_ 'accuracy_failure') -eq 'timeout' }).Count -gt 0) {
    return 'fmm_e35_runtime_limited'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e35_execution_gap'
  }
  if (@($Rows | Where-Object {
        (Get-Field $_ 'timing_scope') -ne 'full_panel_apply_sampled_audit' -or
        (Get-Field $_ 'full_apply_stage_status') -ne 'pass' -or
        (Get-Field $_ 'compression_status') -notin @('pass','fallback_dominated') -or
        -not (Field-IsFinite $_ 'kifmm_compression_ratio') -or
        -not (Field-IsFinite $_ 'cartesian_k4_compression_ratio')
      }).Count -gt 0) {
    return 'fmm_e35_execution_gap'
  }
  if (@($Rows | Where-Object {
        (Get-Field $_ 'coverage') -ne (Get-Field $_ 'full_coverage_expected') -or
        (Get-Field $_ 'sampled_dense_valid_pairs') -ne (Get-Field $_ 'sampled_coverage_expected') -or
        (As-Int64 $_ 'production_invalid_pairs' 1) -ne 0 -or
        (As-Int64 $_ 'nonproduction_fallback_pairs' 1) -ne 0
      }).Count -gt 0) {
    return 'fmm_e35_coverage_gap'
  }
  if (@($Rows | Where-Object { (As-Double $_ 'k4_rms') -gt 5.0e-6 -or (Get-Field $_ 'cartesian_k4_accuracy_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e35_cartesian_k4_gap'
  }
  if (@($Rows | Where-Object { (As-Double $_ 'k1_rms') -gt 5.0e-6 -or (As-Double $_ 'k2_rms') -gt 5.0e-6 -or (As-Double $_ 'k3_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e35_k123_accuracy_gap'
  }
  if (@($Rows | Where-Object { (Get-Field $_ 'accuracy_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e35_k123_accuracy_gap'
  }
  if (@($Rows | Where-Object { (As-Double $_ 'far_fraction' 0.0) -lt 0.10 }).Count -gt 0) {
    return 'fmm_e35_fallback_dominated'
  }
  return 'fmm_e35_compression_runtime_pass'
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
    $attrPath = Join-Path $WorkRoot "$($case.id).stdout.txt"
    if (Test-Path -LiteralPath $attrPath) {
      foreach ($attrLine in @(Get-Content -LiteralPath $attrPath | Where-Object { $_ -like 'E35_ATTR *' })) {
        $attr = Parse-KeyValueLine ([string]$attrLine) 'E35_ATTR '
        Ensure-Field $attr 'row_id' $case.id
        $script:AttributionRows += $attr
      }
    }
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

Export-Table $matrix $MatrixCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,surface_points,panel_limit,compression_status,full_apply_stage_status,far_fraction,kifmm_compression_ratio,cartesian_k4_compression_ratio,far_panel_pairs,direct_panel_pairs,cartesian_k4_far_pairs,cartesian_k4_fallback_pairs) $CompressionCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,hierarchy_seconds,cartesian_k4_total_seconds,production_fallback_seconds,dense_oracle_seconds,total_seconds,kifmm_compression_ratio,cartesian_k4_compression_ratio) $RuntimeCsv
Export-Table $AttributionRows $AttributionCsv

$decision = Decide $matrix $cases.Count
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'full' })
  rows = $matrix.Count
  required_rows = $cases.Count
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($matrix.Count)",
  "required_rows=$($cases.Count)",
  "matrix=$MatrixCsv",
  "compression=$CompressionCsv",
  "runtime=$RuntimeCsv",
  "attribution=$AttributionCsv"
)

$report = @(
  '# TITAN FMM Engine E35 Full-Panel Compression/Runtime Attribution Report',
  '',
  "Decision: $decision",
  "Mode: $(if ($Quick) { 'quick' } else { 'full' })",
  "Rows: $($matrix.Count)/$($cases.Count)",
  '',
  'E35 runs the full-panel dual-worklist channel-hybrid sampled audit and exports compression/runtime attribution for KiFMM K1-K3 plus guarded Cartesian K4 under q8 strict-local-check policy.'
)
if ($decision -eq 'fmm_e35_fallback_dominated') {
  $report += 'Interpretation: sampled audit passes but operator remains fallback-dominated under current q8 strict-local-check policy.'
} elseif ($decision -eq 'fmm_e35_compression_runtime_pass') {
  $report += 'Interpretation: sampled audit accuracy, coverage, compression, and runtime attribution gates pass on the required rows.'
} elseif ($decision -eq 'fmm_e35_runtime_limited') {
  $report += 'Interpretation: at least one required row timed out or did not produce a complete sampled audit record.'
}
foreach ($row in $matrix) {
  $report += ('  {0}: exec={1}, acc={2}, comp={3}, full_apply={4}, far_fraction={5}, kifmm_ratio={6}, cart_k4_ratio={7}, K1={8}, K2={9}, K3={10}, K4={11}, coverage={12}/{13}, sampled={14}/{15}, total={16}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'execution_status'), (Get-Field $row 'accuracy_status'), `
      (Get-Field $row 'compression_status'), (Get-Field $row 'full_apply_stage_status'), `
      (Get-Field $row 'far_fraction'), (Get-Field $row 'kifmm_compression_ratio'), `
      (Get-Field $row 'cartesian_k4_compression_ratio'), (Get-Field $row 'k1_rms'), `
      (Get-Field $row 'k2_rms'), (Get-Field $row 'k3_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'coverage'), (Get-Field $row 'full_coverage_expected'), `
      (Get-Field $row 'sampled_dense_valid_pairs'), (Get-Field $row 'sampled_coverage_expected'), `
      (Get-Field $row 'total_seconds'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E35 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Compression: $CompressionCsv"
Write-Host "Runtime: $RuntimeCsv"
Write-Host "Attribution: $AttributionCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -ne 'fmm_e35_compression_runtime_pass') {
  throw "FMM-E35 decision failed: $decision"
}
