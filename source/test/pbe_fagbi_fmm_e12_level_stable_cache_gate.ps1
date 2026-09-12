param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e12_level_stable_boxes'),
  [string]$Cxx = '',
  [switch]$Quick,
  [switch]$Full,
  [switch]$Force,
  [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E12_LEVEL_STABLE_BOX_REPORT.txt'
if ([string]::IsNullOrWhiteSpace($Pqr1sh1)) {
  $Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
}
if ([string]::IsNullOrWhiteSpace($Pqr1a63)) {
  $Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
}
foreach ($pqr in @($Pqr1sh1, $Pqr1a63)) {
  if (-not (Test-Path -LiteralPath $pqr)) { throw "missing PQR fixture: $pqr" }
}
if (-not $Quick -and -not $Full) { $Quick = $true }

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$BuildRoot = Join-Path $WorkRoot 'build'
$RunRoot = Join-Path $WorkRoot 'runs'
New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null

$MatrixPath = Join-Path $WorkRoot 'fmm_e12_matrix.csv'
$CachePath = Join-Path $WorkRoot 'fmm_e12_cache_geometry.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e12_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e12_summary.dat'

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
  throw 'No C++ compiler found for KIFMM-lite E12 gate.'
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E12 profile driver' }
  } finally {
    Pop-Location
  }
  return $exe
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [string]$Policy, [string]$Role) {
  $sampleLabel = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    order = 6
    theta = 0.6
    condition_policy = 'qr-rank-only'
    source_box_policy = $Policy
    role = $Role
    row_id = ('{0}_{1}_{2}_q6_{3}' -f $Fixture, $sampleLabel, $SampleCount, ($Policy -replace '[^A-Za-z0-9]+','_'))
  }
}

function Parse-CaseLine([string]$Line) {
  $fields = [ordered]@{}
  foreach ($token in ($Line -split '\s+')) {
    if ($token -eq 'KIFMM_SAMPLED_PANEL_CASE') { continue }
    $parts = $token -split '=', 2
    if ($parts.Count -eq 2) { $fields[$parts[0]] = $parts[1] }
  }
  foreach ($required in @('fixture','sample','status','sampled_rows','order','theta','condition_policy','source_box_policy','condition_computed','aggregate_rms','k1','k2','k3','k4','max_equiv_residual','far_panel_pairs','fallback_panel_pairs','invalid_pairs','nonproduction_fallback_pairs','equivalent_build_solve_seconds','total_seconds','cache_entries','equivalent_factor_cache_hits','equivalent_factor_cache_misses','cache_hit_rate','unique_half_widths','unique_kappa_h_values','source_tree_levels','unique_source_node_half_widths','source_root_half_width','level_cube_geometry','worklist_hash','policy_hash')) {
    if (-not $fields.Contains($required)) { throw "malformed E12 profile row, missing $required`: $Line" }
  }
  [pscustomobject]$fields
}

function As-Double($Value) {
  [double]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}
function As-Long($Value) {
  [int64]::Parse([string]$Value, [Globalization.CultureInfo]::InvariantCulture)
}

function Parse-ResultFromLog([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path)) { return $null }
  $line = Get-Content -LiteralPath $Path | Where-Object { $_ -like 'KIFMM_SAMPLED_PANEL_CASE*' } | Select-Object -Last 1
  if ($null -eq $line) { return $null }
  try { Parse-CaseLine ([string]$line) } catch { return $null }
}

function Run-Row($Exe, $Row) {
  $rowDir = Join-Path $RunRoot $Row.row_id
  New-Item -ItemType Directory -Path $rowDir -Force | Out-Null
  $rawLog = Join-Path $rowDir 'raw.log'
  $parsed = if ($Force) { $null } else { Parse-ResultFromLog $rawLog }
  if ($null -ne $parsed) {
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
      "--theta=$($Row.theta)",
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
      throw "E12 row failed before producing parseable output: $($Row.row_id)`n$($output -join [Environment]::NewLine)"
    }
  }
  if ($null -eq $parsed) { throw "no parseable result for $($Row.row_id)" }
  [pscustomobject]@{
    row_id = $Row.row_id
    fixture = $Row.fixture
    sample_set = $Row.sample
    sample_count = $Row.sample_count
    role = $Row.role
    order = $Row.order
    theta = $Row.theta
    source_box_policy = [string]$parsed.source_box_policy
    status = [string]$parsed.status
    aggregate_rms = As-Double $parsed.aggregate_rms
    max_channel_rms = [Math]::Max([Math]::Max((As-Double $parsed.k1), (As-Double $parsed.k2)), [Math]::Max((As-Double $parsed.k3), (As-Double $parsed.k4)))
    condition_policy = [string]$parsed.condition_policy
    condition_computed = [string]$parsed.condition_computed
    invalid_pairs = As-Long $parsed.invalid_pairs
    nonproduction_fallback_pairs = As-Long $parsed.nonproduction_fallback_pairs
    unique_far_source_nodes = As-Long $parsed.unique_far_source_nodes
    cached_source_nodes = As-Long $parsed.cached_source_nodes
    equivalent_system_cache_entries = As-Long $parsed.cache_entries
    equivalent_factor_cache_hits = As-Long $parsed.equivalent_factor_cache_hits
    equivalent_factor_cache_misses = As-Long $parsed.equivalent_factor_cache_misses
    equivalent_cache_hit_rate = As-Double $parsed.cache_hit_rate
    unique_equivalent_half_widths = As-Long $parsed.unique_half_widths
    unique_equivalent_kappa_h_values = As-Long $parsed.unique_kappa_h_values
    source_tree_levels = As-Long $parsed.source_tree_levels
    unique_source_node_half_widths = As-Long $parsed.unique_source_node_half_widths
    source_root_half_width = As-Double $parsed.source_root_half_width
    level_cube_geometry = [string]$parsed.level_cube_geometry
    equivalent_build_solve_seconds = As-Double $parsed.equivalent_build_solve_seconds
    total_seconds = As-Double $parsed.total_seconds
    far_panel_pairs = As-Long $parsed.far_panel_pairs
    fallback_panel_pairs = As-Long $parsed.fallback_panel_pairs
    far_panel_fraction = if ((As-Long $parsed.far_panel_pairs) + (As-Long $parsed.fallback_panel_pairs) -eq 0) { 0.0 } else { [double](As-Long $parsed.far_panel_pairs) / [double]((As-Long $parsed.far_panel_pairs) + (As-Long $parsed.fallback_panel_pairs)) }
    fallback_panel_fraction = if ((As-Long $parsed.far_panel_pairs) + (As-Long $parsed.fallback_panel_pairs) -eq 0) { 0.0 } else { [double](As-Long $parsed.fallback_panel_pairs) / [double]((As-Long $parsed.far_panel_pairs) + (As-Long $parsed.fallback_panel_pairs)) }
    max_equiv_residual = As-Double $parsed.max_equiv_residual
    worklist_hash = [string]$parsed.worklist_hash
    policy_hash = [string]$parsed.policy_hash
  }
}

$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 'tight-bbox' 'quick_baseline'),
  (New-Row '1SH1' 'uniform_hash' 128 'level-cube' 'quick_candidate')
)
$fullRows = @(
  $quickRows
  (New-Row '1SH1' 'uniform_hash' 512 'level-cube' 'uniform_512'),
  (New-Row '1SH1' 'inward_normal_groove' 128 'level-cube' 'groove'),
  (New-Row '1SH1' 'closest_opposing_surface' 128 'level-cube' 'opposing'),
  (New-Row '1A63' 'uniform_hash' 128 'level-cube' 'uniform_1a63'),
  (New-Row '1A63' 'high_aspect' 128 'level-cube' 'aspect_1a63')
)
$rows = if ($Full) { $fullRows } else { $quickRows }
$seen = @{}
$rows = @($rows | Where-Object {
    if ($seen.ContainsKey($_.row_id)) { $false } else { $seen[$_.row_id] = $true; $true }
  })

$exe = Build-Driver
$parsedRows = @()
foreach ($row in $rows) {
  $parsedRows += Run-Row $exe $row
}

$parsedRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
$parsedRows | Select-Object row_id,source_box_policy,source_tree_levels,unique_source_node_half_widths,source_root_half_width,equivalent_system_cache_entries,equivalent_factor_cache_hits,equivalent_factor_cache_misses,equivalent_cache_hit_rate,unique_equivalent_half_widths,unique_equivalent_kappa_h_values,equivalent_build_solve_seconds,total_seconds,far_panel_pairs,fallback_panel_pairs,far_panel_fraction,fallback_panel_fraction | Export-Csv -LiteralPath $CachePath -NoTypeInformation

$tight = $parsedRows | Where-Object { $_.role -eq 'quick_baseline' } | Select-Object -First 1
$level = $parsedRows | Where-Object { $_.role -eq 'quick_candidate' } | Select-Object -First 1
$requiredRows = @($parsedRows | Where-Object { $_.source_box_policy -eq 'level-cube' })
$accuracyPass = ($requiredRows.Count -gt 0) -and (@($requiredRows | Where-Object {
      $_.status -ne 'pass' -or $_.aggregate_rms -gt 1.0e-6 -or $_.max_channel_rms -gt 5.0e-6 -or $_.invalid_pairs -ne 0 -or $_.nonproduction_fallback_pairs -ne 0 -or $_.max_equiv_residual -gt 1.0e-10
    }).Count -eq 0)
$levelGeometry = $null -ne $level -and $level.level_cube_geometry -eq 'true'
$halfWidthGate = $null -ne $tight -and $null -ne $level -and $level.unique_source_node_half_widths -le [Math]::Max(1, [Math]::Floor($tight.unique_source_node_half_widths / 4))
$cacheHitGate = $null -ne $level -and $level.equivalent_cache_hit_rate -ge 0.25
$timeGate = $null -ne $tight -and $null -ne $level -and $level.equivalent_build_solve_seconds -lt $tight.equivalent_build_solve_seconds

$decision = if (-not $accuracyPass) {
  'fmm_e12_accuracy_gap'
} elseif (-not $levelGeometry) {
  'fmm_e12_worklist_gap'
} elseif ($halfWidthGate -and $cacheHitGate -and $timeGate) {
  if ($Full) { 'fmm_e12_level_cube_cache_pass' } else { 'fmm_e12_level_cube_cache_quick_pass_full_gate_pending' }
} elseif ($halfWidthGate -and -not $cacheHitGate) {
  'fmm_e12_cache_key_gap'
} else {
  'fmm_e12_cache_geometry_partial'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Full) { 'full' } else { 'quick' }
  rows = $parsedRows.Count
  tight_unique_source_node_half_widths = if ($null -eq $tight) { 0 } else { $tight.unique_source_node_half_widths }
  level_unique_source_node_half_widths = if ($null -eq $level) { 0 } else { $level.unique_source_node_half_widths }
  level_cache_hit_rate = if ($null -eq $level) { 0.0 } else { $level.equivalent_cache_hit_rate }
  level_equivalent_seconds = if ($null -eq $level) { 0.0 } else { $level.equivalent_build_solve_seconds }
  tight_far_panel_fraction = if ($null -eq $tight) { 0.0 } else { $tight.far_panel_fraction }
  level_far_panel_fraction = if ($null -eq $level) { 0.0 } else { $level.far_panel_fraction }
  tight_total_seconds = if ($null -eq $tight) { 0.0 } else { $tight.total_seconds }
  level_total_seconds = if ($null -eq $level) { 0.0 } else { $level.total_seconds }
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation

@(
  "decision $decision",
  "mode $(if ($Full) { 'full' } else { 'quick' })",
  "rows $($parsedRows.Count)",
  "tight_unique_source_node_half_widths $(if ($null -eq $tight) { 0 } else { $tight.unique_source_node_half_widths })",
  "level_unique_source_node_half_widths $(if ($null -eq $level) { 0 } else { $level.unique_source_node_half_widths })",
  "level_cache_hit_rate $(if ($null -eq $level) { 0.0 } else { $level.equivalent_cache_hit_rate })",
  "level_equivalent_seconds $(if ($null -eq $level) { 0.0 } else { $level.equivalent_build_solve_seconds })",
  "tight_far_panel_fraction $(if ($null -eq $tight) { 0.0 } else { $tight.far_panel_fraction })",
  "level_far_panel_fraction $(if ($null -eq $level) { 0.0 } else { $level.far_panel_fraction })"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

$tightCacheEntries = if ($null -eq $tight) { 0 } else { $tight.equivalent_system_cache_entries }
$levelCacheEntries = if ($null -eq $level) { 0 } else { $level.equivalent_system_cache_entries }
$tightCacheHitRate = if ($null -eq $tight) { 0.0 } else { $tight.equivalent_cache_hit_rate }
$levelCacheHitRate = if ($null -eq $level) { 0.0 } else { $level.equivalent_cache_hit_rate }
$tightFarPairs = if ($null -eq $tight) { 0 } else { $tight.far_panel_pairs }
$levelFarPairs = if ($null -eq $level) { 0 } else { $level.far_panel_pairs }
$tightFallbackPairs = if ($null -eq $tight) { 0 } else { $tight.fallback_panel_pairs }
$levelFallbackPairs = if ($null -eq $level) { 0 } else { $level.fallback_panel_pairs }
$tightTotalSeconds = if ($null -eq $tight) { 0.0 } else { $tight.total_seconds }
$levelTotalSeconds = if ($null -eq $level) { 0.0 } else { $level.total_seconds }

@"
TITAN FMM-E12 Level-Stable KIFMM Box Cache Gate

Decision: $decision
Mode: $(if ($Full) { 'full' } else { 'quick' })

E12 compares the historical tight-bounding-box source tree against a
level-stable enclosing-box policy for KIFMM-lite sampled production-panel rows.
The level-cube policy preserves the existing deterministic binary split order
but rounds each source box half-width up to a root-relative power-of-two value,
so equivalent-system factorization can be reused across source nodes with the
same half-width.

Quick measured comparison:
- tight unique source-node half-widths: $(if ($null -eq $tight) { 0 } else { $tight.unique_source_node_half_widths })
- level-cube unique source-node half-widths: $(if ($null -eq $level) { 0 } else { $level.unique_source_node_half_widths })
- level-cube cache hit rate: $(if ($null -eq $level) { 0.0 } else { $level.equivalent_cache_hit_rate })
- tight equivalent build/solve seconds: $(if ($null -eq $tight) { 0.0 } else { $tight.equivalent_build_solve_seconds })
- level-cube equivalent build/solve seconds: $(if ($null -eq $level) { 0.0 } else { $level.equivalent_build_solve_seconds })
- level-cube aggregate RMS: $(if ($null -eq $level) { 0.0 } else { $level.aggregate_rms })
- level-cube max channel RMS: $(if ($null -eq $level) { 0.0 } else { $level.max_channel_rms })

Artifacts:
- Matrix: $MatrixPath
- Cache geometry: $CachePath
- Decision: $DecisionPath
- Summary: $SummaryPath

Scope:
E12 remains diagnostic-only. It changes KIFMM-lite sampled tree geometry behind
an explicit source-box policy and does not wire KIFMM into production TABI
matvec, GMRES, preconditioners, parser defaults, or manuscript claims.

Terminology:
The E12 "level-cube" policy is a level-stable source-box cache policy, not a
full octree KIFMM or a KIFMM-local upward/downward hierarchy. It preserves the
existing deterministic binary split order and tight node centers, but rounds
source box half-widths up to root-relative power-of-two enclosing values. This
is sufficient for equivalent-matrix reuse because the canonical equivalent
matrix depends on order, kernel, scale, condition policy, and half-width, not on
absolute center.

Work decomposition:
The quick row reduces equivalent-system cache entries from $tightCacheEntries
to $levelCacheEntries and improves cache hit rate from
$tightCacheHitRate to $levelCacheHitRate. It
also changes work decomposition: accepted far panel pairs move from
$tightFarPairs to $levelFarPairs, while fallback pairs move
from $tightFallbackPairs to $levelFallbackPairs. The
measured total-time change, $tightTotalSeconds seconds to
$levelTotalSeconds seconds, is therefore a combined effect of cache reuse
and a more conservative far/fallback split, not a cache-only speedup claim.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E12 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Cache geometry: $CachePath"
if (-not $AllowFailedDecision -and $decision -notin @('fmm_e12_level_cube_cache_pass','fmm_e12_level_cube_cache_quick_pass_full_gate_pending')) {
  throw "FMM E12 gate did not pass: $decision"
}
