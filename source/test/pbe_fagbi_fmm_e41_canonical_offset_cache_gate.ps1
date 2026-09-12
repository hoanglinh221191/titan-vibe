param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e41_canonical_offset_cache'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e41_matrix.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e41_runtime.csv'
$CacheCsv = Join-Path $WorkRoot 'fmm_e41_cache.csv'
$ParityCsv = Join-Path $WorkRoot 'fmm_e41_parity.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e41_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e41_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E41_CANONICAL_OFFSET_CACHE_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E41 smoke.'
}

function Assert-SmokeSupportsE41 {
  $smoke = Get-Content -Raw -LiteralPath $SmokeSource
  $headers = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy.cpp')
  $jcp = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp')
  foreach ($needle in @('canonical-offset-cache', 'm2l-matrix-cache-policy', 'm2l_matrix_cache_hit_rate=')) {
    if (-not ($smoke.Contains($needle) -or $headers.Contains($needle) -or $jcp.Contains($needle))) {
      throw "E41 canonical-offset cache support is not present yet: missing '$needle'"
    }
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE41
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
  if ($Rows.Count -eq 0) { Set-Content -LiteralPath $Path -Value '' }
  else { $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path }
}

function Resolve-TreeStopCandidate {
  [pscustomobject]@{
    key = 'panel-radius-occupancy_r1p5_occ32'
    policy = 'panel-radius-occupancy'
    factor = 1.5
    occupancy = 32
    source = 'e39_quick_best_fixed_for_e41'
  }
}

function New-Case(
  [string]$Id,
  [string]$Sample,
  [int]$SampleCount,
  [int]$SurfacePoints,
  [int]$PanelLimit,
  [string]$MatrixCachePolicy,
  [object]$TreeCandidate
) {
  [pscustomobject]@{
    id = $Id
    fixture = '1SH1'
    sample = $Sample
    sample_count = $SampleCount
    surface_points = $SurfacePoints
    panel_limit = $PanelLimit
    m2l_matrix_cache_policy = $MatrixCachePolicy
    tree_stop_policy = $TreeCandidate.policy
    panel_radius_factor = $TreeCandidate.factor
    occupancy_stop = $TreeCandidate.occupancy
    tree_candidate_key = $TreeCandidate.key
    row_id = "${Id}_${MatrixCachePolicy}_$($TreeCandidate.key)".Replace('.', 'p')
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
  Ensure-Field $Row 'requested_m2l_matrix_cache_policy' $Case.m2l_matrix_cache_policy
  Ensure-Field $Row 'requested_tree_candidate_key' $Case.tree_candidate_key
  foreach ($field in @(
      'm2l_matrix_cache_policy', 'm2l_matrix_cache_entries',
      'm2l_matrix_cache_hits', 'm2l_matrix_cache_misses',
      'm2l_matrix_cache_hit_rate', 'm2l_matrix_cache_estimated_bytes',
      'm2l_matrix_cache_build_seconds', 'm2l_matrix_gemm_seconds',
      'm2l_matrix_cache_key_status', 'tree_stop_policy',
      'tree_panel_radius_factor', 'tree_occupancy_stop', 'm2l_eval_policy',
      'm2l_seconds', 'm2l_translations', 'aggregate_rms',
      'max_channel_rms', 'k1_rms', 'k2_rms', 'k3_rms', 'k4_rms',
      'coverage', 'full_coverage_expected', 'production_invalid_pairs',
      'nonproduction_fallback_pairs', 'far_fraction', 'far_panel_pairs',
      'total_seconds', 'output_hash', 'policy_hash', 'worklist_hash')) {
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
    "--m2l-matrix-cache-policy=$($Case.m2l_matrix_cache_policy)",
    "--tree-stop-policy=$($Case.tree_stop_policy)",
    "--tree-panel-radius-factor=$($Case.panel_radius_factor)",
    "--tree-occupancy-stop=$($Case.occupancy_stop)",
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
  return ,@($Rows | Where-Object { $_.case_id -eq $CaseId -and (Get-Field $_ 'requested_m2l_matrix_cache_policy') -eq $Policy } | Select-Object -First 1)
}

function Build-ParityRows([array]$Rows, [array]$CaseIds) {
  $parityRows = @()
  foreach ($caseId in $CaseIds) {
    $off = Case-Row $Rows $caseId 'off'
    $relative = Case-Row $Rows $caseId 'relative-geometry-cache'
    $canonical = Case-Row $Rows $caseId 'canonical-offset-cache'
    if ($off.Count -ne 1 -or $canonical.Count -ne 1) { continue }
    $maxDelta = 0.0
    foreach ($field in @('aggregate_rms', 'k1_rms', 'k2_rms', 'k3_rms', 'k4_rms', 'max_channel_rms')) {
      $maxDelta = [Math]::Max($maxDelta, [Math]::Abs((As-Double $off[0] $field 0.0) - (As-Double $canonical[0] $field 0.0)))
    }
    $relativeHitRate = if ($relative.Count -eq 1) { As-Double $relative[0] 'm2l_matrix_cache_hit_rate' 0.0 } else { 0.0 }
    $canonicalM2l = As-Double $canonical[0] 'm2l_seconds' 0.0
    $canonicalTotal = As-Double $canonical[0] 'total_seconds' 0.0
    $parityRows += [pscustomobject][ordered]@{
      case_id = $caseId
      uncached_row = Get-Field $off[0] 'row_id'
      canonical_row = Get-Field $canonical[0] 'row_id'
      output_hash_match = if ((Get-Field $off[0] 'output_hash') -eq (Get-Field $canonical[0] 'output_hash')) { 'pass' } else { 'fail' }
      max_channel_rms_delta = $maxDelta
      parity_status = if ($maxDelta -le 1.0e-12) { 'pass' } else { 'fail' }
      relative_hit_rate = $relativeHitRate
      canonical_hit_rate = As-Double $canonical[0] 'm2l_matrix_cache_hit_rate' 0.0
      hit_rate_gain_vs_relative = (As-Double $canonical[0] 'm2l_matrix_cache_hit_rate' 0.0) - $relativeHitRate
      uncached_m2l_seconds = As-Double $off[0] 'm2l_seconds' 0.0
      canonical_m2l_seconds = $canonicalM2l
      m2l_speedup = if ($canonicalM2l -gt 0.0) { (As-Double $off[0] 'm2l_seconds' 0.0) / $canonicalM2l } else { 0.0 }
      uncached_total_seconds = As-Double $off[0] 'total_seconds' 0.0
      canonical_total_seconds = $canonicalTotal
      total_speedup = if ($canonicalTotal -gt 0.0) { (As-Double $off[0] 'total_seconds' 0.0) / $canonicalTotal } else { 0.0 }
    }
  }
  $parityRows
}

function Decide-E41([array]$Rows, [array]$ParityRows, [bool]$IsQuick) {
  if (@($Rows | Where-Object { (Get-Field $_ 'execution_status') -eq 'timeout' }).Count -gt 0) { return 'fmm_e41_runtime_limited' }
  if (@($Rows | Where-Object { -not (Is-Accurate $_) -or -not (Is-Covered $_) }).Count -gt 0) { return 'fmm_e41_accuracy_gap' }
  if (@($ParityRows | Where-Object { $_.parity_status -ne 'pass' }).Count -gt 0) { return 'fmm_e41_canonical_parity_gap' }
  if (@($ParityRows | Where-Object { [double]$_.canonical_hit_rate -lt 0.70 -or [double]$_.hit_rate_gain_vs_relative -lt 0.20 }).Count -gt 0) { return 'fmm_e41_canonical_cache_key_gap' }
  if ($IsQuick) { return 'fmm_e41_quick_pass_full_pending' }
  $fullUniform = @($ParityRows | Where-Object { $_.case_id -eq '1SH1_s8_full_uniform128' } | Select-Object -First 1)
  if ($fullUniform.Count -eq 1 -and [double]$fullUniform[0].m2l_speedup -ge 1.5 -and [double]$fullUniform[0].total_speedup -ge 1.15) {
    return 'fmm_e41_canonical_offset_cache_pass'
  }
  'fmm_e41_accuracy_pass_partial_speedup'
}

function Update-DerivedCsvs([array]$Rows, [array]$CaseIds) {
  $runtimeRows = @()
  $cacheRows = @()
  foreach ($row in $Rows) {
    $runtimeRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      m2l_matrix_cache_policy = Get-Field $row 'm2l_matrix_cache_policy'
      m2l_seconds = Get-Field $row 'm2l_seconds'
      m2l_translations = Get-Field $row 'm2l_translations'
      total_seconds = Get-Field $row 'total_seconds'
      aggregate_rms = Get-Field $row 'aggregate_rms'
      max_channel_rms = Get-Field $row 'max_channel_rms'
      far_fraction = Get-Field $row 'far_fraction'
    }
    $cacheRows += [pscustomobject][ordered]@{
      row_id = Get-Field $row 'row_id'
      case_id = Get-Field $row 'case_id'
      m2l_matrix_cache_policy = Get-Field $row 'm2l_matrix_cache_policy'
      m2l_matrix_cache_entries = Get-Field $row 'm2l_matrix_cache_entries'
      m2l_matrix_cache_hits = Get-Field $row 'm2l_matrix_cache_hits'
      m2l_matrix_cache_misses = Get-Field $row 'm2l_matrix_cache_misses'
      m2l_matrix_cache_hit_rate = Get-Field $row 'm2l_matrix_cache_hit_rate'
      m2l_matrix_cache_estimated_bytes = Get-Field $row 'm2l_matrix_cache_estimated_bytes'
      m2l_matrix_cache_build_seconds = Get-Field $row 'm2l_matrix_cache_build_seconds'
      m2l_matrix_gemm_seconds = Get-Field $row 'm2l_matrix_gemm_seconds'
      m2l_matrix_cache_key_status = Get-Field $row 'm2l_matrix_cache_key_status'
    }
  }
  $parityRows = Build-ParityRows $Rows $CaseIds
  Export-Table $Rows $MatrixCsv
  Export-Table $runtimeRows $RuntimeCsv
  Export-Table $cacheRows $CacheCsv
  Export-Table $parityRows $ParityCsv
  $parityRows
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
if ($Force) {
  Remove-Item -LiteralPath $MatrixCsv, $RuntimeCsv, $ParityCsv, $CacheCsv, $DecisionCsv, $SummaryDat, $ReportPath -Force -ErrorAction SilentlyContinue
}
Build-Smoke

$treeCandidate = Resolve-TreeStopCandidate
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
  foreach ($policy in @('off', 'relative-geometry-cache', 'canonical-offset-cache')) {
    $case = New-Case $spec[0] $spec[1] $spec[2] $spec[3] $spec[4] $policy $treeCandidate
    $existing = @($rows | Where-Object { $_.row_id -eq $case.row_id } | Select-Object -First 1)
    $shouldRun = $true
    if ($existing.Count -eq 1 -and -not $Force) {
      $shouldRun = $false
      if ($RerunFailed -and (Get-Field $existing[0] 'execution_status') -ne 'pass') {
        $shouldRun = $true
      }
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
$decision = Decide-E41 $rows $parityRows ([bool]$Quick)

[pscustomobject][ordered]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  tree_candidate = $treeCandidate.key
  generated_at = (Get-Date).ToString('s')
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

@"
FMM-E41 canonical-offset M2L matrix-cache gate
decision=$decision
mode=$(if ($Quick) { 'quick' } else { 'full' })
tree_candidate=$($treeCandidate.key)
matrix=$MatrixCsv
runtime=$RuntimeCsv
cache=$CacheCsv
parity=$ParityCsv
"@ | Set-Content -LiteralPath $SummaryDat

@"
TITAN FMM-E41 CANONICAL-OFFSET M2L MATRIX CACHE REPORT

Decision: $decision

E41 tests an opt-in canonical-offset M2L matrix cache. It keeps the E39
tree-stop candidate, q8, strict-local-check, target-check-batched M2L,
sampled dense production-direct audit, and all production defaults unchanged.

Outputs:
- Matrix: $MatrixCsv
- Runtime: $RuntimeCsv
- Cache telemetry: $CacheCsv
- Canonical-vs-uncached parity: $ParityCsv
- Decision: $DecisionCsv
"@ | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E41 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Cache: $CacheCsv"
Write-Host "Parity: $ParityCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."
if ($decision -notin @('fmm_e41_canonical_offset_cache_pass', 'fmm_e41_quick_pass_full_pending', 'fmm_e41_accuracy_pass_partial_speedup') -and -not $AllowFailedDecision) {
  throw "FMM-E41 gate did not pass: $decision"
}
