param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e20_linearity_repeatability'),
  [switch]$Quick,
  [switch]$Full,
  [string[]]$RerunRows = @(),
  [switch]$RerunFailed,
  [switch]$Force,
  [switch]$AllowFailedDecision,
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) { throw "missing PQR fixture: $pqr" }
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) { throw "requested C++ compiler not found: $Requested" }
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
  throw 'No C++ compiler found for KIFMM E20.1 gate.'
}

function Build-GateExe([string]$Exe) {
  $compiler = Resolve-Cxx $Cxx
  $args = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_full_panel_matvec_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $Exe
  )
  Push-Location $RepoRoot
  try {
    & $compiler @args
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E20.1 KIFMM frozen-apply driver' }
  } finally {
    Pop-Location
  }
}

function Parse-KeyValueLine([string]$Line) {
  $out = [ordered]@{}
  foreach ($part in ($Line -split '\s+')) {
    $idx = $part.IndexOf('=')
    if ($idx -le 0) { continue }
    $out[$part.Substring(0, $idx)] = $part.Substring($idx + 1)
  }
  return $out
}

function D([hashtable]$Map, [string]$Key, [double]$Default = [double]::NaN) {
  if (-not $Map.Contains($Key)) { return $Default }
  $v = 0.0
  if ([double]::TryParse($Map[$Key], [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$v)) { return $v }
  return $Default
}

function U64([hashtable]$Map, [string]$Key) {
  if (-not $Map.Contains($Key)) { return [uint64]0 }
  $v = [uint64]0
  if ([uint64]::TryParse($Map[$Key], [ref]$v)) { return $v }
  return [uint64]0
}

function S([hashtable]$Map, [string]$Key, [string]$Default = '') {
  if ($Map.Contains($Key)) { return [string]$Map[$Key] }
  return $Default
}

function Is-Finite([double]$Value) {
  return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Assert-E201CanonicalSources {
  $fullHeader = Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.h'
  $fullSource = Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp'
  $smoke = Join-Path $RepoRoot 'test\titan_fmm_kifmm_full_panel_matvec_smoke.cpp'
  foreach ($path in @($fullHeader, $fullSource, $smoke)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "missing E20.1 input: $path" }
  }
  $headerText = Get-Content -Raw -LiteralPath $fullHeader
  foreach ($needle in @(
      'KifmmFullPanelPreparedState',
      'frozen_worklist',
      'frozen_equivalent_system_cache',
      'prepared_equivalent_system_cache_entries')) {
    if (-not $headerText.Contains($needle)) { throw "expected '$needle' in $fullHeader" }
  }
  $sourceText = Get-Content -Raw -LiteralPath $fullSource
  foreach ($needle in @(
      'prebuild_equivalent_system_cache',
      'state.equivalent_system_cache',
      'state.pairs',
      'state.surface_rule')) {
    if (-not $sourceText.Contains($needle)) { throw "expected '$needle' in $fullSource" }
  }
  $stale = @(Get-ChildItem -Recurse -File (Join-Path $RepoRoot 'src\tabi'),(Join-Path $RepoRoot 'test') |
    Where-Object { $_.Name -match '(?i)(copy|stage|stale|old|bak|tmp)' -and $_.Name -match '(?i)kifmm' })
  if ($stale.Count -gt 0) {
    throw "stale KIFMM source/header candidates found: $($stale.FullName -join '; ')"
  }
}

$quickRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_128_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=128; order=6 },
  [pscustomobject]@{ row_id='1A63_high_aspect_128_q6'; fixture='1A63'; sample='high_aspect'; sample_count=128; order=6 }
)
$fullRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_512_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1SH1_closest_512_q6'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1A63_uniform_512_q6'; fixture='1A63'; sample='uniform_hash'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1A63_high_aspect_512_q6'; fixture='1A63'; sample='high_aspect'; sample_count=512; order=6 }
)

$mode = if ($Full -and -not $Quick) { 'full' } else { 'quick' }
$rowsToRun = if ($mode -eq 'full') { $fullRows } else { $quickRows }
$matrixPath = Join-Path $WorkRoot 'fmm_e201_matrix.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e201_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e201_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E201_FROZEN_PREPARED_APPLY_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e201_frozen_prepared_apply.exe'

Assert-E201CanonicalSources
Build-GateExe -Exe $exe

function Invoke-E201Row($row) {
  $log = Join-Path $WorkRoot ($row.row_id + '.log')
  $line = @()
  if ((Test-Path -LiteralPath $log) -and -not $Force -and -not ($RerunRows -contains $row.row_id)) {
    $line = @(Get-Content -LiteralPath $log | Where-Object { [string]$_ -match '^KIFMM_E20_CASE ' } | Select-Object -First 1)
  }
  if ($line.Count -gt 0) {
    $m0 = Parse-KeyValueLine ([string]$line[0])
    $status0 = S $m0 'status'
    $reuse0 = U64 $m0 'prepared_context_reused'
    if (-not ($RerunFailed -and ($status0 -ne 'pass' -or $reuse0 -ne 1))) {
      Write-Host "resume: $($row.row_id)"
      return Convert-E201MapToRow $row $m0 0 $log
    }
  }

  Write-Host "run: $($row.row_id)"
  $args = @(
    $Pqr1sh1,
    $Pqr1a63,
    '--e20-profile',
    "--fixture=$($row.fixture)",
    "--sample=$($row.sample)",
    "--sample-count=$($row.sample_count)",
    "--order=$($row.order)"
  )
  Push-Location $RepoRoot
  try {
    $output = & $exe @args 2>&1
    $code = $LASTEXITCODE
  } finally {
    Pop-Location
  }
  $output | Set-Content -LiteralPath $log -Encoding ascii
  $line = @($output | Where-Object { [string]$_ -match '^KIFMM_E20_CASE ' } | Select-Object -First 1)
  if ($line.Count -eq 0) {
    return [pscustomobject]@{
      row_id=$row.row_id; mode=$mode; fixture=$row.fixture; sample_set=$row.sample; sample_count=$row.sample_count; order=$row.order
      status='run_failed'; process_exit_code=$code; linearity_relative_error=[double]::NaN
      prepared_context_reused=0; frozen_worklist=0; frozen_equivalent_system_cache=0
      prepared_worklist_pairs=0; prepared_far_source_nodes=0; prepared_equivalent_system_cache_entries=0
      repeat_output_hash_stable=0; repeat_worklist_hash_stable=0; repeat_policy_hash_stable=0
      full_pair_coverage=0; sampled_direct_valid_pairs=0; sampled_invalid_pairs=0; nonproduction_fallback_pairs=0; far_fraction=0.0
      tree_build_seconds=[double]::NaN; worklist_build_seconds=[double]::NaN
      equivalent_factor_cache_misses=[uint64]0; equivalent_factor_cache_hits=[uint64]0
      setup_seconds=[double]::NaN; apply_seconds=[double]::NaN; total_seconds=[double]::NaN
      worklist_hash=''; policy_hash=''; panel_geometry_hash=''; sampled_target_hash=''; log=$log
    }
  }
  $m = Parse-KeyValueLine ([string]$line[0])
  return Convert-E201MapToRow $row $m $code $log
}

function Convert-E201MapToRow($row, [hashtable]$m, [int]$code, [string]$log) {
  [pscustomobject]@{
    row_id=$row.row_id
    mode=$mode
    fixture=$row.fixture
    sample_set=$row.sample
    sample_count=(U64 $m 'sampled_rows')
    order=(U64 $m 'order')
    status=(S $m 'status')
    process_exit_code=$code
    operator_status_x=(S $m 'operator_status_x')
    operator_status_y=(S $m 'operator_status_y')
    operator_status_z=(S $m 'operator_status_z')
    linearity_relative_error=(D $m 'linearity_relative_error')
    prepare_build_count=(U64 $m 'prepare_build_count')
    apply_count=(U64 $m 'apply_count')
    prepared_context_reused=(U64 $m 'prepared_context_reused')
    frozen_worklist=(U64 $m 'frozen_worklist')
    frozen_equivalent_system_cache=(U64 $m 'frozen_equivalent_system_cache')
    prepared_worklist_pairs=(U64 $m 'prepared_worklist_pairs')
    prepared_far_source_nodes=(U64 $m 'prepared_far_source_nodes')
    prepared_equivalent_system_cache_entries=(U64 $m 'prepared_equivalent_system_cache_entries')
    repeat_output_hash_stable=(U64 $m 'repeat_output_hash_stable')
    repeat_worklist_hash_stable=(U64 $m 'repeat_worklist_hash_stable')
    repeat_policy_hash_stable=(U64 $m 'repeat_policy_hash_stable')
    full_pair_coverage=(U64 $m 'full_pair_coverage')
    sampled_direct_valid_pairs=(U64 $m 'sampled_direct_valid_pairs')
    sampled_invalid_pairs=(U64 $m 'sampled_invalid_pairs')
    nonproduction_fallback_pairs=(U64 $m 'nonproduction_fallback_pairs')
    far_fraction=(D $m 'far_fraction')
    tree_build_seconds=(D $m 'tree_build_seconds' 0.0)
    worklist_build_seconds=(D $m 'worklist_build_seconds' 0.0)
    equivalent_factor_cache_misses=(U64 $m 'equivalent_factor_cache_misses')
    equivalent_factor_cache_hits=(U64 $m 'equivalent_factor_cache_hits')
    setup_seconds=(D $m 'setup_seconds')
    apply_seconds=(D $m 'apply_seconds')
    total_seconds=(D $m 'total_seconds')
    worklist_hash=(S $m 'worklist_hash')
    policy_hash=(S $m 'policy_hash')
    panel_geometry_hash=(S $m 'panel_geometry_hash')
    sampled_target_hash=(S $m 'sampled_target_hash')
    log=$log
  }
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($row in $rowsToRun) { $matrix.Add((Invoke-E201Row $row)) }
$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii

$required = @($matrix.ToArray())
$runtimeLimited = @($required | Where-Object { $_.status -eq 'run_failed' }).Count -gt 0
$linearityGap = @($required | Where-Object { -not (Is-Finite ([double]$_.linearity_relative_error)) -or [double]$_.linearity_relative_error -gt 1.0e-10 }).Count -gt 0
$hashGap = @($required | Where-Object { [uint64]$_.repeat_output_hash_stable -ne 1 -or [uint64]$_.repeat_worklist_hash_stable -ne 1 -or [uint64]$_.repeat_policy_hash_stable -ne 1 }).Count -gt 0
$coverageGap = @($required | Where-Object { [uint64]$_.sampled_invalid_pairs -ne 0 -or [uint64]$_.nonproduction_fallback_pairs -ne 0 -or [uint64]$_.full_pair_coverage -eq 0 -or [uint64]$_.sampled_direct_valid_pairs -eq 0 }).Count -gt 0
$fallbackDominated = @($required | Where-Object { [double]$_.far_fraction -lt 0.10 }).Count -gt 0
$lifecycleGap = @($required | Where-Object { [uint64]$_.prepare_build_count -ne 1 -or [uint64]$_.apply_count -lt 6 -or [uint64]$_.prepared_context_reused -ne 1 -or [uint64]$_.frozen_worklist -ne 1 -or [uint64]$_.frozen_equivalent_system_cache -ne 1 -or [uint64]$_.prepared_worklist_pairs -eq 0 -or [uint64]$_.prepared_far_source_nodes -eq 0 -or [uint64]$_.prepared_equivalent_system_cache_entries -eq 0 -or [uint64]$_.equivalent_factor_cache_misses -ne 0 }).Count -gt 0

$decision = 'fmm_e201_frozen_prepared_apply_pass'
if ($runtimeLimited) {
  $decision = 'fmm_e201_runtime_limited'
} elseif ($linearityGap) {
  $decision = 'fmm_e201_linearity_gap'
} elseif ($hashGap) {
  $decision = 'fmm_e201_hash_instability'
} elseif ($coverageGap) {
  $decision = 'fmm_e201_coverage_gap'
} elseif ($fallbackDominated) {
  $decision = 'fmm_e201_fallback_dominated'
} elseif ($lifecycleGap) {
  $decision = 'fmm_e201_cache_lifecycle_gap'
}

[pscustomobject]@{
  decision=$decision
  mode=$mode
  rows=$matrix.Count
  max_linearity_error=(@($matrix | Measure-Object -Property linearity_relative_error -Maximum).Maximum)
  min_far_fraction=(@($matrix | Measure-Object -Property far_fraction -Minimum).Minimum)
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$mode",
  "rows=$($matrix.Count)",
  "matrix=$matrixPath",
  "decision_csv=$decisionPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$report = @(
  'TITAN FMM Engine E20.1 Frozen Prepared Apply Report',
  '',
  "Decision: $decision",
  "Mode: $mode",
  "Rows: $($matrix.Count)",
  '',
  'E20.1 verifies that the full-panel KIFMM-lite diagnostic operator has a frozen geometry/worklist/equivalent-system lifecycle before any GMRES wiring.',
  'The apply path consumes only a new unknown vector, recomputes source-node RHS/charges, and reuses the prepared worklist and equivalent-system cache.',
  'This is still diagnostic infrastructure only: no production TABI matvec, GMRES path, preconditioner, manuscript, or default behavior is changed.',
  '',
  'Rows:',
  '  row_id,status,linearity_relative_error,prepared_context_reused,cache_misses,far_fraction,apply_seconds,total_seconds'
)
foreach ($row in $matrix) {
  $report += ("  {0},{1},{2},{3},{4},{5},{6},{7}" -f $row.row_id,$row.status,$row.linearity_relative_error,$row.prepared_context_reused,$row.equivalent_factor_cache_misses,$row.far_fraction,$row.apply_seconds,$row.total_seconds)
}
$report += @(
  '',
  'Artifacts:',
  "  Matrix: $matrixPath",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$report | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E20.1 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -ne 'fmm_e201_frozen_prepared_apply_pass' -and -not $AllowFailedDecision) {
  throw "E20.1 did not pass: $decision"
}
