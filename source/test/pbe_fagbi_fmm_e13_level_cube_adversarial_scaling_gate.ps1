param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e13_level_cube_scaling'),
  [string]$Cxx = '',
  [switch]$Quick,
  [switch]$Full,
  [string[]]$RerunRows = @(),
  [switch]$RerunFailed,
  [switch]$Force,
  [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E13_LEVEL_CUBE_SCALING_REPORT.txt'
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) { throw "missing PQR fixture: $pqr" }
}
if (-not $Quick -and -not $Full) { $Full = $true }

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$BuildRoot = Join-Path $WorkRoot 'build'
$RunRoot = Join-Path $WorkRoot 'runs'
New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null

$MatrixPath = Join-Path $WorkRoot 'fmm_e13_matrix.csv'
$CacheGeometryPath = Join-Path $WorkRoot 'fmm_e13_cache_geometry.csv'
$ThetaOrderPath = Join-Path $WorkRoot 'fmm_e13_theta_order.csv'
$AdversarialPath = Join-Path $WorkRoot 'fmm_e13_adversarial.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e13_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e13_summary.dat'

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
  throw 'No C++ compiler found for KIFMM-lite E13 gate.'
}

function Build-Driver {
  $compiler = Resolve-Cxx $Cxx
  $exe = Join-Path $BuildRoot 'titan_fmm_kifmm_sampled_panel_matvec_smoke.exe'
  $args = @(
    '-std=c++23',
    '-O2',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_sampled_panel_matvec_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $exe
  )
  Push-Location $RepoRoot
  try {
    & $compiler @args
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E13 sampled profile driver' }
  } finally {
    Pop-Location
  }
  return $exe
}

function Format-ThetaLabel([double]$Theta) {
  return ('theta{0}' -f (($Theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture)) -replace '\.', 'p'))
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [int]$Order, [double]$Theta, [string]$Policy,
                 [string]$Role) {
  $sampleLabel = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  $policyLabel = $Policy -replace '[^A-Za-z0-9]+', '_'
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_label = $sampleLabel
    sample_count = $SampleCount
    order = $Order
    theta = $Theta
    condition_policy = 'qr-rank-only'
    source_box_policy = $Policy
    role = $Role
    row_id = ('{0}_{1}_{2}_q{3}_{4}_{5}' -f $Fixture, $sampleLabel, $SampleCount, $Order, (Format-ThetaLabel $Theta), $policyLabel)
  }
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @(
      'fixture','sample','status','source_panels','sampled_rows','order','theta',
      'condition_policy','source_box_policy','aggregate_rms','k1','k2','k3','k4',
      'max_equiv_residual','far_panel_pairs','fallback_panel_pairs',
      'production_direct_reference_valid_pairs','invalid_pairs',
      'nonproduction_fallback_pairs','equivalent_build_solve_seconds',
      'total_seconds','cache_entries','equivalent_factor_cache_hits',
      'equivalent_factor_cache_misses','cache_hit_rate','unique_half_widths',
      'unique_kappa_h_values','laplace_systems','yukawa_systems',
      'unique_far_source_nodes','cached_source_nodes','source_tree_levels',
      'unique_source_node_half_widths','level_cube_geometry','worklist_hash',
      'policy_hash','panel_geometry_hash','unknown_vector_hash',
      'sampled_target_hash')) {
    if (-not $fields.Contains($required)) { throw "malformed E13 profile row, missing $required`: $Line" }
  }
  [pscustomobject]$fields
}

function Parse-ResultFromLog([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path)) { return $null }
  $line = Get-Content -LiteralPath $Path | Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } | Select-Object -Last 1
  if ($null -eq $line) { return $null }
  try { Parse-CaseLine ([string]$line) } catch { return $null }
}

function As-Double($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return [double]::NaN }
  [double]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}

function As-Long($Value) {
  if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return [int64]0 }
  [int64]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}

function Is-AcceptableRow($Parsed) {
  if ($null -eq $Parsed) { return $false }
  return [string]$Parsed.status -eq 'pass'
}

function Run-Row($Exe, $Row) {
  $rowDir = Join-Path $RunRoot $Row.row_id
  New-Item -ItemType Directory -Path $rowDir -Force | Out-Null
  $rawLog = Join-Path $rowDir 'raw.log'
  $parsed = Parse-ResultFromLog $rawLog
  $shouldRun = $Force -or $null -eq $parsed
  if ($RerunFailed -and -not (Is-AcceptableRow $parsed)) { $shouldRun = $true }
  if ($RerunRows.Count -gt 0) { $shouldRun = $true }
  if (-not $shouldRun) {
    Write-Host "resume: $($Row.row_id)"
  } else {
    Write-Host "run: $($Row.row_id)"
    $runArgs = @(
      $Pqr1sh1,
      $Pqr1a63,
      '--e10-factor-profile',
      "--fixture=$($Row.fixture)",
      "--sample=$($Row.sample)",
      "--sample-count=$($Row.sample_count)",
      "--order=$($Row.order)",
      "--theta=$($Row.theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture))",
      '--min-far-panel-fraction=0.10',
      "--condition-policy=$($Row.condition_policy)",
      "--source-box-policy=$($Row.source_box_policy)"
    )
    Push-Location $RepoRoot
    try {
      $output = & $Exe @runArgs 2>&1
      $code = $LASTEXITCODE
    } finally {
      Pop-Location
    }
    $output | Set-Content -LiteralPath $rawLog -Encoding ascii
    $parsed = Parse-ResultFromLog $rawLog
    if ($code -ne 0 -and $null -eq $parsed) {
      $parsed = [pscustomobject]@{
        fixture = $Row.fixture
        sample = $Row.sample
        status = 'runtime_failed'
        source_panels = 0
        sampled_rows = $Row.sample_count
        order = $Row.order
        theta = $Row.theta
      }
    }
  }
  if ($null -eq $parsed) { throw "no parseable result for $($Row.row_id)" }

  $farPairs = As-Long $parsed.far_panel_pairs
  $fallbackPairs = As-Long $parsed.fallback_panel_pairs
  $totalPairs = $farPairs + $fallbackPairs
  $uniqueFar = As-Long $parsed.unique_far_source_nodes
  $cacheMisses = As-Long $parsed.equivalent_factor_cache_misses
  $equivSeconds = As-Double $parsed.equivalent_build_solve_seconds
  $laplaceEntries = As-Long $parsed.laplace_systems
  $yukawaEntries = As-Long $parsed.yukawa_systems
  [pscustomobject]@{
    row_id = $Row.row_id
    fixture = $Row.fixture
    sample_set = $Row.sample
    sample_count = $Row.sample_count
    role = $Row.role
    order = As-Long $parsed.order
    theta = As-Double $parsed.theta
    condition_policy = [string]$parsed.condition_policy
    source_box_policy = [string]$parsed.source_box_policy
    status = [string]$parsed.status
    source_panels = As-Long $parsed.source_panels
    sampled_rows = As-Long $parsed.sampled_rows
    aggregate_rms = As-Double $parsed.aggregate_rms
    k1_rms = As-Double $parsed.k1
    k2_rms = As-Double $parsed.k2
    k3_rms = As-Double $parsed.k3
    k4_rms = As-Double $parsed.k4
    max_channel_rms = [Math]::Max([Math]::Max((As-Double $parsed.k1), (As-Double $parsed.k2)), [Math]::Max((As-Double $parsed.k3), (As-Double $parsed.k4)))
    max_equiv_residual = As-Double $parsed.max_equiv_residual
    max_condition = As-Double $parsed.max_condition
    invalid_pairs = As-Long $parsed.invalid_pairs
    nonproduction_fallback_pairs = As-Long $parsed.nonproduction_fallback_pairs
    production_direct_reference_valid_pairs = As-Long $parsed.production_direct_reference_valid_pairs
    far_panel_pairs = $farPairs
    fallback_panel_pairs = $fallbackPairs
    far_panel_fraction = if ($totalPairs -eq 0) { 0.0 } else { [double]$farPairs / [double]$totalPairs }
    fallback_panel_fraction = if ($totalPairs -eq 0) { 0.0 } else { [double]$fallbackPairs / [double]$totalPairs }
    equivalent_build_solve_seconds = $equivSeconds
    equivalent_seconds_per_unique_far_source_node = if ($uniqueFar -eq 0) { [double]::NaN } else { $equivSeconds / [double]$uniqueFar }
    equivalent_seconds_per_cache_miss = if ($cacheMisses -eq 0) { [double]::NaN } else { $equivSeconds / [double]$cacheMisses }
    total_seconds = As-Double $parsed.total_seconds
    direct_reference_seconds = As-Double $parsed.direct_reference_seconds
    far_eval_seconds = As-Double $parsed.far_eval_seconds
    fallback_eval_seconds = As-Double $parsed.fallback_eval_seconds
    equivalent_system_cache_entries = As-Long $parsed.cache_entries
    equivalent_factor_cache_hits = As-Long $parsed.equivalent_factor_cache_hits
    equivalent_factor_cache_misses = $cacheMisses
    equivalent_cache_hit_rate = As-Double $parsed.cache_hit_rate
    unique_equivalent_half_widths = As-Long $parsed.unique_half_widths
    unique_equivalent_kappa_h_values = As-Long $parsed.unique_kappa_h_values
    laplace_cache_entries = $laplaceEntries
    yukawa_cache_entries = $yukawaEntries
    cache_entries_per_kernel = if (($laplaceEntries + $yukawaEntries) -eq 0) { 0.0 } else { [double](As-Long $parsed.cache_entries) / 2.0 }
    unique_far_source_nodes = $uniqueFar
    cached_source_nodes = As-Long $parsed.cached_source_nodes
    source_tree_levels = As-Long $parsed.source_tree_levels
    unique_source_node_half_widths = As-Long $parsed.unique_source_node_half_widths
    level_cube_geometry = [string]$parsed.level_cube_geometry
    worklist_hash = [string]$parsed.worklist_hash
    policy_hash = [string]$parsed.policy_hash
    panel_geometry_hash = [string]$parsed.panel_geometry_hash
    unknown_vector_hash = [string]$parsed.unknown_vector_hash
    sampled_target_hash = [string]$parsed.sampled_target_hash
  }
}

$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 6 0.6 'tight-bbox' 'tight_baseline'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.6 'level-cube' 'theta_tuning'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.7 'level-cube' 'theta_tuning'),
  (New-Row '1SH1' 'uniform_hash' 128 6 0.8 'level-cube' 'theta_tuning'),
  (New-Row '1SH1' 'uniform_hash' 128 8 0.6 'level-cube' 'q8_confirm')
)
$fullRows = @(
  $quickRows
  (New-Row '1SH1' 'uniform_hash' 512 6 0.6 'level-cube' 'uniform_scaling'),
  (New-Row '1SH1' 'closest_opposing_surface' 128 6 0.6 'level-cube' 'adversarial'),
  (New-Row '1SH1' 'high_aspect' 128 6 0.6 'level-cube' 'adversarial'),
  (New-Row '1SH1' 'inward_normal_groove' 128 6 0.6 'level-cube' 'adversarial'),
  (New-Row '1A63' 'uniform_hash' 128 6 0.6 'level-cube' 'cross_fixture'),
  (New-Row '1A63' 'uniform_hash' 512 6 0.6 'level-cube' 'cross_fixture')
)

$rows = if ($Quick) { $quickRows } else { $fullRows }
$seen = @{}
$rows = @($rows | Where-Object {
    if ($seen.ContainsKey($_.row_id)) { $false } else { $seen[$_.row_id] = $true; $true }
  })
if ($RerunRows.Count -gt 0) {
  $wanted = @{}
  foreach ($id in $RerunRows) { $wanted[$id] = $true }
  $rows = @($rows | Where-Object { $wanted.ContainsKey($_.row_id) })
  if ($rows.Count -eq 0) { throw 'RerunRows did not match any E13 row id' }
}

$exe = Build-Driver
$parsedRows = @()
foreach ($row in $rows) {
  $parsedRows += Run-Row $exe $row
}

$parsedRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation -Encoding ascii
$parsedRows |
  Select-Object row_id,fixture,sample_set,sample_count,role,order,theta,source_box_policy,status,equivalent_system_cache_entries,equivalent_factor_cache_hits,equivalent_factor_cache_misses,equivalent_cache_hit_rate,unique_equivalent_half_widths,unique_equivalent_kappa_h_values,laplace_cache_entries,yukawa_cache_entries,cache_entries_per_kernel,unique_far_source_nodes,cached_source_nodes,source_tree_levels,unique_source_node_half_widths,level_cube_geometry,equivalent_build_solve_seconds,equivalent_seconds_per_unique_far_source_node,equivalent_seconds_per_cache_miss,total_seconds,far_panel_pairs,fallback_panel_pairs,far_panel_fraction,fallback_panel_fraction |
  Export-Csv -LiteralPath $CacheGeometryPath -NoTypeInformation -Encoding ascii
$parsedRows |
  Where-Object { $_.sample_set -eq 'uniform_hash' -and $_.fixture -eq '1SH1' } |
  Sort-Object sample_count,order,theta,source_box_policy |
  Select-Object row_id,fixture,sample_count,role,order,theta,source_box_policy,status,aggregate_rms,max_channel_rms,total_seconds,equivalent_build_solve_seconds,far_panel_fraction,fallback_panel_fraction,equivalent_cache_hit_rate |
  Export-Csv -LiteralPath $ThetaOrderPath -NoTypeInformation -Encoding ascii
$parsedRows |
  Where-Object { $_.sample_set -ne 'uniform_hash' } |
  Select-Object row_id,fixture,sample_set,sample_count,order,theta,status,aggregate_rms,max_channel_rms,far_panel_fraction,fallback_panel_fraction,equivalent_cache_hit_rate,total_seconds |
  Export-Csv -LiteralPath $AdversarialPath -NoTypeInformation -Encoding ascii

function Test-AccuracyGate($row) {
  if ($null -eq $row) { return $false }
  $expected = $row.source_panels * $row.sampled_rows
  return $row.status -eq 'pass' -and
    $row.aggregate_rms -le 1.0e-6 -and
    $row.max_channel_rms -le 5.0e-6 -and
    $row.max_equiv_residual -le 1.0e-10 -and
    $row.invalid_pairs -eq 0 -and
    $row.nonproduction_fallback_pairs -eq 0 -and
    (($row.far_panel_pairs + $row.fallback_panel_pairs) -eq $expected) -and
    $row.production_direct_reference_valid_pairs -eq $expected
}

$tightBaseline = $parsedRows | Where-Object { $_.row_id -eq '1SH1_uniform_128_q6_theta0p6_tight_bbox' } | Select-Object -First 1
$levelQ6 = $parsedRows | Where-Object { $_.row_id -eq '1SH1_uniform_128_q6_theta0p6_level_cube' } | Select-Object -First 1
$levelQ8 = $parsedRows | Where-Object { $_.row_id -eq '1SH1_uniform_128_q8_theta0p6_level_cube' } | Select-Object -First 1
$levelRows = @($parsedRows | Where-Object { $_.source_box_policy -eq 'level-cube' })
$requiredFullRows = @($parsedRows | Where-Object {
    $_.source_box_policy -eq 'level-cube' -and
    ($_.role -in @('theta_tuning','q8_confirm','uniform_scaling','adversarial','cross_fixture'))
  })
$accuracyGap = @($levelRows | Where-Object { -not (Test-AccuracyGate $_) })
$cacheGap = @($levelRows | Where-Object { $_.equivalent_cache_hit_rate -lt 0.90 -or $_.level_cube_geometry -ne 'true' })
$fallbackDominated = @($levelRows | Where-Object { $_.far_panel_fraction -lt 0.10 })
$q8Stable = $false
if ($null -ne $levelQ6 -and $null -ne $levelQ8 -and (Test-AccuracyGate $levelQ8)) {
  $q8Stable = ($levelQ8.aggregate_rms -le ($levelQ6.aggregate_rms * 1.10)) -and
    ($levelQ8.max_channel_rms -le ($levelQ6.max_channel_rms * 1.10))
}
$timeGate = $null -ne $tightBaseline -and $null -ne $levelQ6 -and $levelQ6.total_seconds -lt $tightBaseline.total_seconds

$thetaRows = @($parsedRows | Where-Object {
    $_.fixture -eq '1SH1' -and $_.sample_set -eq 'uniform_hash' -and
    $_.sample_count -eq 128 -and $_.order -eq 6 -and
    $_.source_box_policy -eq 'level-cube'
  })
$thetaPassCandidates = @($thetaRows | Where-Object {
    (Test-AccuracyGate $_) -and
    $_.equivalent_cache_hit_rate -ge 0.90 -and
    $_.far_panel_fraction -ge 0.10
  })
$bestThetaCandidate = $thetaPassCandidates | Sort-Object total_seconds | Select-Object -First 1
$bestThetaIsAbove06 = $false
if ($null -ne $bestThetaCandidate) {
  $bestThetaIsAbove06 = $bestThetaCandidate.theta -gt (0.6 + 1.0e-9)
}

$decision = 'fmm_e13_runtime_limited'
if ($accuracyGap.Count -gt 0) {
  $decision = 'fmm_e13_level_cube_accuracy_gap'
} elseif ($cacheGap.Count -gt 0) {
  $decision = 'fmm_e13_cache_geometry_gap'
} elseif ($fallbackDominated.Count -gt 0) {
  $decision = 'fmm_e13_fallback_dominated'
} elseif (-not $timeGate -or -not $q8Stable) {
  $decision = 'fmm_e13_level_cube_accuracy_gap'
} elseif ($bestThetaIsAbove06) {
  $decision = 'fmm_e13_theta07_or_08_candidate'
} elseif ($Quick) {
  $decision = 'fmm_e13_level_cube_quick_pass_full_gate_pending'
} elseif ($requiredFullRows.Count -ge 10) {
  $decision = 'fmm_e13_level_cube_1sh1_adversarial_1a63_uniform_pass'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $parsedRows.Count
  level_rows = $levelRows.Count
  level_accuracy_gaps = $accuracyGap.Count
  cache_geometry_gaps = $cacheGap.Count
  fallback_dominated_rows = $fallbackDominated.Count
  tight_total_seconds = if ($null -eq $tightBaseline) { 0.0 } else { $tightBaseline.total_seconds }
  level_q6_total_seconds = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.total_seconds }
  level_q6_cache_hit_rate = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.equivalent_cache_hit_rate }
  level_q6_far_fraction = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.far_panel_fraction }
  q8_stable_or_improved = $q8Stable
  theta_candidate_count = $thetaPassCandidates.Count
  theta_candidate = if ($null -eq $bestThetaCandidate) { '' } else { $bestThetaCandidate.theta.ToString('0.0', [Globalization.CultureInfo]::InvariantCulture) }
  theta_candidate_row_id = if ($null -eq $bestThetaCandidate) { '' } else { $bestThetaCandidate.row_id }
  theta_candidate_rule = 'lowest_total_seconds_among_passing_level_cube_q6_uniform128_theta_rows'
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($parsedRows.Count)",
  "matrix=$MatrixPath",
  "cache_geometry=$CacheGeometryPath",
  "theta_order=$ThetaOrderPath",
  "adversarial=$AdversarialPath"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

$tightSeconds = if ($null -eq $tightBaseline) { 0.0 } else { $tightBaseline.total_seconds }
$levelSeconds = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.total_seconds }
$levelCacheHit = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.equivalent_cache_hit_rate }
$levelFarFraction = if ($null -eq $levelQ6) { 0.0 } else { $levelQ6.far_panel_fraction }
$thetaNote = if ($null -eq $bestThetaCandidate) {
  'No theta row passed the full quick theta-candidate gate.'
} else {
  "Best theta candidate: theta=$($bestThetaCandidate.theta), row=$($bestThetaCandidate.row_id), far_fraction=$($bestThetaCandidate.far_panel_fraction), total_seconds=$($bestThetaCandidate.total_seconds)."
}

@"
TITAN FMM-E13 Level-Stable KIFMM Cache Policy Adversarial And Scaling Gate

Decision: $decision
Mode: $(if ($Quick) { 'quick' } else { 'full' })

E13 extends the E12 level-stable source-box cache policy from quick evidence to
adversarial sampled rows, sample-count scaling, q-order confirmation, and theta
tuning. This remains diagnostic-only; no KIFMM production TABI matvec, GMRES,
preconditioner, parser default, or manuscript claim is changed.

Baseline comparison:
- tight-bbox 1SH1 uniform-128 q6 total seconds: $tightSeconds
- level-cube 1SH1 uniform-128 q6 total seconds: $levelSeconds
- level-cube 1SH1 uniform-128 q6 cache hit rate: $levelCacheHit
- level-cube 1SH1 uniform-128 q6 far fraction: $levelFarFraction
- q8 stable or improved versus q6: $q8Stable
- $thetaNote

Scope boundary:
The adversarial rows in E13 are 1SH1 closest-opposing-surface, high-aspect,
and inward-normal-groove samples. The 1A63 coverage in this scope is uniform
128/512 cross-fixture sampling, not 1A63 adversarial sampling.

Artifacts:
- Matrix: $MatrixPath
- Cache geometry: $CacheGeometryPath
- Theta/order rows: $ThetaOrderPath
- Adversarial rows: $AdversarialPath
- Decision: $DecisionPath
- Summary: $SummaryPath

Terminology:
The source-box policy is level-stable half-width quantization for equivalent
system reuse. It is not a full spatial octree KIFMM or a production local
upward/downward hierarchy.

Gate notes:
The level-stable policy must keep sampled dense-row accuracy below aggregate
RMS 1e-6 and max-channel RMS 5e-6, keep direct-reference coverage exact, keep
nonproduction fallback at zero, maintain cache hit rate at least 0.90, retain
far panel fraction at least 0.10, and beat the same-sample tight-bbox baseline.
Under condition_policy=qr-rank-only, max_condition=-1 is a sentinel meaning the
condition estimate was intentionally not computed; rank/finite checks and the
equivalent residual gate remain active.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E13 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Cache geometry: $CacheGeometryPath"
Write-Host "Theta/order: $ThetaOrderPath"
Write-Host "Adversarial: $AdversarialPath"
if ($decision -notin @('fmm_e13_level_cube_quick_pass_full_gate_pending','fmm_e13_level_cube_1sh1_adversarial_1a63_uniform_pass','fmm_e13_theta07_or_08_candidate') -and -not $AllowFailedDecision) {
  throw "FMM E13 gate did not pass: $decision"
}
