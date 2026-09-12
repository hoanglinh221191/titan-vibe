param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e19_full_panel_matvec'),
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
  throw 'No C++ compiler found for KIFMM E19 gate.'
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E19 full-panel KIFMM driver' }
  } finally {
    Pop-Location
  }
}

function Assert-CanonicalPair([string]$Prefix) {
  $tabiDir = Join-Path $RepoRoot 'src\tabi'
  $allowed = @("$Prefix.h", "$Prefix.cpp")
  $matches = @(Get-ChildItem -LiteralPath $tabiDir -File -Filter "$Prefix*" | Select-Object -ExpandProperty Name)
  $extra = @($matches | Where-Object { $allowed -notcontains $_ })
  foreach ($required in $allowed) {
    if ($matches -notcontains $required) { throw "missing canonical source/header for ${Prefix}: $required" }
  }
  if ($extra.Count -gt 0) { throw "stale duplicate source/header for ${Prefix}: $($extra -join ', ')" }
}

function Assert-E19CanonicalSources {
  foreach ($prefix in @(
      'tabi_fmm_kifmm_lite',
      'tabi_fmm_kifmm_real_panel_cloud',
      'tabi_fmm_kifmm_full_panel_matvec')) {
    Assert-CanonicalPair $prefix
  }
  $fullHeader = Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.h'
  $text = Get-Content -Raw -LiteralPath $fullHeader
  foreach ($needle in @(
      'KifmmFullPanelMatvecResult',
      'full_pair_coverage',
      'sampled_direct_reference_valid_pairs',
      'full_output_hash_phi',
      'full_output_hash_phin',
      'full_channel_output_hash',
      'setup_seconds',
      'apply_seconds',
      'first_apply_seconds',
      'run_kifmm_lite_full_panel_matvec')) {
    if (-not $text.Contains($needle)) { throw "expected '$needle' in $fullHeader" }
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

$quickRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_128_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=128; order=6 },
  [pscustomobject]@{ row_id='1A63_high_aspect_128_q6'; fixture='1A63'; sample='high_aspect'; sample_count=128; order=6 }
)
$fullRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_512_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1SH1_closest_opposing_surface_512_q6'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1SH1_high_aspect_512_q6'; fixture='1SH1'; sample='high_aspect'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1SH1_inward_normal_groove_512_q6'; fixture='1SH1'; sample='inward_normal_groove'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1A63_uniform_512_q6'; fixture='1A63'; sample='uniform_hash'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1A63_closest_opposing_surface_512_q6'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=512; order=6 },
  [pscustomobject]@{ row_id='1A63_high_aspect_512_q6'; fixture='1A63'; sample='high_aspect'; sample_count=512; order=6 }
)

$mode = if ($Full -and -not $Quick) { 'full' } else { 'quick' }
$rowsToRun = if ($mode -eq 'full') { $fullRows } else { $quickRows }
$matrixPath = Join-Path $WorkRoot 'fmm_e19_matrix.csv'
$auditPath = Join-Path $WorkRoot 'fmm_e19_sampled_audit.csv'
$coveragePath = Join-Path $WorkRoot 'fmm_e19_coverage.csv'
$runtimePath = Join-Path $WorkRoot 'fmm_e19_runtime_cache.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e19_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e19_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E19_FULL_PANEL_MATVEC_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e19_full_panel_matvec.exe'

Assert-E19CanonicalSources
Build-GateExe -Exe $exe

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($old in @(Import-Csv -LiteralPath $matrixPath)) { $existing[$old.row_id] = $old }
}

function Invoke-E19Row($row) {
  if ($existing.ContainsKey($row.row_id)) {
    $old = $existing[$row.row_id]
    $rerun = $Force -or ($RerunRows -contains $row.row_id) -or
      (-not ($old.PSObject.Properties.Name -contains 'full_pair_coverage')) -or
      (-not ($old.PSObject.Properties.Name -contains 'full_output_hash_phi')) -or
      (-not ($old.PSObject.Properties.Name -contains 'full_output_hash_phin')) -or
      (-not ($old.PSObject.Properties.Name -contains 'full_channel_output_hash')) -or
      (-not ($old.PSObject.Properties.Name -contains 'setup_seconds')) -or
      (-not ($old.PSObject.Properties.Name -contains 'apply_seconds')) -or
      (-not ($old.PSObject.Properties.Name -contains 'first_apply_seconds')) -or
      ($RerunFailed -and ($old.status -ne 'pass'))
    if (-not $rerun) {
      Write-Host "resume: $($row.row_id)"
      return $old
    }
  }

  Write-Host "run: $($row.row_id)"
  $log = Join-Path $WorkRoot ($row.row_id + '.log')
  $args = @(
    $Pqr1sh1,
    $Pqr1a63,
    '--e19-profile',
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
  $line = @($output | Where-Object { [string]$_ -match '^KIFMM_FULL_PANEL_CASE ' } | Select-Object -First 1)
  if ($code -ne 0 -or $line.Count -eq 0) {
    return [pscustomobject]@{
      row_id=$row.row_id; mode=$mode; fixture=$row.fixture; sample_set=$row.sample; sample_count=$row.sample_count; order=$row.order
      status='run_failed'; source_panels=0; target_panels=0; full_pair_coverage=0; sampled_direct_valid_pairs=0; sampled_invalid_pairs=0; nonproduction_fallback_pairs=0; production_fallback_valid_pairs=0
      far_panel_pairs=0; fallback_panel_pairs=0; far_fraction=0.0; aggregate_rms=[double]::NaN; max_channel_rms=[double]::NaN; max_channel_name='none'; k1=[double]::NaN; k2=[double]::NaN; k3=[double]::NaN; k4=[double]::NaN
      worklist_hash=''; policy_hash=''; panel_geometry_hash=''; unknown_vector_hash=''; sampled_target_hash=''; full_output_hash_phi=''; full_output_hash_phin=''; full_channel_output_hash=''; setup_seconds=[double]::NaN; apply_seconds=[double]::NaN; first_apply_seconds=[double]::NaN; repeat_apply_seconds=[double]::NaN; repeat_apply_min_seconds=[double]::NaN; repeat_apply_max_seconds=[double]::NaN; total_seconds=[double]::NaN; log=$log
    }
  }
  $m = Parse-KeyValueLine ([string]$line[0])
  [pscustomobject]@{
    row_id=$row.row_id
    mode=$mode
    fixture=$row.fixture
    sample_set=$row.sample
    sample_count=(U64 $m 'sampled_rows')
    order=(U64 $m 'order')
    status=(S $m 'status')
    source_panels=(U64 $m 'source_panels')
    target_panels=(U64 $m 'target_panels')
    full_pair_coverage=(U64 $m 'full_pair_coverage')
    sampled_direct_valid_pairs=(U64 $m 'sampled_direct_valid_pairs')
    sampled_invalid_pairs=(U64 $m 'sampled_invalid_pairs')
    nonproduction_fallback_pairs=(U64 $m 'nonproduction_fallback_pairs')
    production_fallback_valid_pairs=(U64 $m 'production_fallback_valid_pairs')
    far_panel_pairs=(U64 $m 'far_panel_pairs')
    fallback_panel_pairs=(U64 $m 'fallback_panel_pairs')
    far_fraction=(D $m 'far_fraction')
    aggregate_rms=(D $m 'aggregate_rms')
    max_channel_rms=(D $m 'max_channel_rms')
    max_channel_name=(S $m 'max_channel_name')
    k1=(D $m 'k1')
    k2=(D $m 'k2')
    k3=(D $m 'k3')
    k4=(D $m 'k4')
    max_equiv_residual=(D $m 'max_equiv_residual')
    max_condition=(D $m 'max_condition')
    cache_entries=(U64 $m 'cache_entries')
    cache_hit_rate=(D $m 'cache_hit_rate')
    unique_far_source_nodes=(U64 $m 'unique_far_source_nodes')
    unique_half_widths=(U64 $m 'unique_half_widths')
    unique_kappa_h_values=(U64 $m 'unique_kappa_h_values')
    tree_build_seconds=(D $m 'tree_build_seconds')
    worklist_build_seconds=(D $m 'worklist_build_seconds')
    equivalent_build_solve_seconds=(D $m 'equivalent_build_solve_seconds')
    full_far_eval_seconds=(D $m 'full_far_eval_seconds')
    full_fallback_eval_seconds=(D $m 'full_fallback_eval_seconds')
    sampled_dense_audit_seconds=(D $m 'sampled_dense_audit_seconds')
    setup_seconds=(D $m 'setup_seconds')
    apply_seconds=(D $m 'apply_seconds')
    first_apply_seconds=(D $m 'first_apply_seconds')
    repeat_apply_seconds=(D $m 'repeat_apply_seconds')
    repeat_apply_min_seconds=(D $m 'repeat_apply_min_seconds')
    repeat_apply_max_seconds=(D $m 'repeat_apply_max_seconds')
    total_seconds=(D $m 'total_seconds')
    panel_geometry_hash=(S $m 'panel_geometry_hash')
    unknown_vector_hash=(S $m 'unknown_vector_hash')
    sampled_target_hash=(S $m 'sampled_target_hash')
    full_output_hash_phi=(S $m 'full_output_hash_phi')
    full_output_hash_phin=(S $m 'full_output_hash_phin')
    full_channel_output_hash=(S $m 'full_channel_output_hash')
    worklist_hash=(S $m 'worklist_hash')
    policy_hash=(S $m 'policy_hash')
    log=$log
  }
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($row in $rowsToRun) { $matrix.Add((Invoke-E19Row $row)) }

if ($mode -eq 'full') {
  $q6Passable = @($matrix | Where-Object { $_.status -eq 'pass' })
  if ($q6Passable.Count -gt 0) {
    $worst = @($q6Passable | Sort-Object {[double]$_.max_channel_rms} -Descending | Select-Object -First 1)[0]
    $q8Row = [pscustomobject]@{
      row_id=("{0}_{1}_{2}_q8_confirm" -f $worst.fixture,$worst.sample_set,$worst.sample_count)
      fixture=$worst.fixture
      sample=$worst.sample_set
      sample_count=[int]$worst.sample_count
      order=8
    }
    $matrix.Add((Invoke-E19Row $q8Row))
  }
}

$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,sample_count,order,status,aggregate_rms,max_channel_rms,max_channel_name,k1,k2,k3,k4,max_equiv_residual | Export-Csv -LiteralPath $auditPath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,source_panels,target_panels,full_pair_coverage,sampled_direct_valid_pairs,sampled_invalid_pairs,nonproduction_fallback_pairs,production_fallback_valid_pairs,far_panel_pairs,fallback_panel_pairs,far_fraction,worklist_hash,policy_hash,panel_geometry_hash,unknown_vector_hash,sampled_target_hash,full_output_hash_phi,full_output_hash_phin,full_channel_output_hash | Export-Csv -LiteralPath $coveragePath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,total_seconds,setup_seconds,apply_seconds,first_apply_seconds,repeat_apply_seconds,repeat_apply_min_seconds,repeat_apply_max_seconds,tree_build_seconds,worklist_build_seconds,equivalent_build_solve_seconds,full_far_eval_seconds,full_fallback_eval_seconds,sampled_dense_audit_seconds,cache_entries,cache_hit_rate,unique_far_source_nodes,unique_half_widths,unique_kappa_h_values | Export-Csv -LiteralPath $runtimePath -NoTypeInformation -Encoding ascii

$required = @($matrix | Where-Object { [int]$_.order -eq 6 })
$runFailed = @($required | Where-Object { $_.status -eq 'run_failed' }).Count -gt 0
$coverageGap = @($required | Where-Object {
    [uint64]$_.full_pair_coverage -ne ([uint64]$_.source_panels * [uint64]$_.target_panels) -or
    [uint64]$_.sampled_direct_valid_pairs -ne ([uint64]$_.source_panels * [uint64]$_.sample_count) -or
    [uint64]$_.sampled_invalid_pairs -ne 0 -or
    [uint64]$_.nonproduction_fallback_pairs -ne 0 -or
    [uint64]$_.production_fallback_valid_pairs -ne [uint64]$_.fallback_panel_pairs
  }).Count -gt 0
$auditGap = @($required | Where-Object { $_.status -ne 'pass' -or [double]$_.aggregate_rms -gt 1.0e-6 -or [double]$_.max_channel_rms -gt 5.0e-6 }).Count -gt 0
$fallbackDominated = @($required | Where-Object { [double]$_.far_fraction -lt 0.10 }).Count -gt 0
$hashGap = @($required | Where-Object {
    [string]::IsNullOrWhiteSpace($_.worklist_hash) -or
    [string]::IsNullOrWhiteSpace($_.policy_hash) -or
    [string]::IsNullOrWhiteSpace($_.panel_geometry_hash) -or
    [string]::IsNullOrWhiteSpace($_.unknown_vector_hash) -or
    [string]::IsNullOrWhiteSpace($_.sampled_target_hash) -or
    [string]::IsNullOrWhiteSpace($_.full_output_hash_phi) -or
    [string]::IsNullOrWhiteSpace($_.full_output_hash_phin) -or
    [string]::IsNullOrWhiteSpace($_.full_channel_output_hash)
  }).Count -gt 0
$timingGap = @($required | Where-Object {
    -not (Is-Finite ([double]$_.setup_seconds)) -or
    -not (Is-Finite ([double]$_.apply_seconds)) -or
    -not (Is-Finite ([double]$_.first_apply_seconds)) -or
    [double]$_.apply_seconds -le 0.0 -or
    [double]$_.first_apply_seconds -le 0.0
  }).Count -gt 0
$q8Gap = $false
if ($mode -eq 'full') {
  $q8Rows = @($matrix | Where-Object { [int]$_.order -eq 8 })
  foreach ($q8 in $q8Rows) {
    $q6 = @($matrix | Where-Object { $_.fixture -eq $q8.fixture -and $_.sample_set -eq $q8.sample_set -and [int]$_.sample_count -eq [int]$q8.sample_count -and [int]$_.order -eq 6 } | Select-Object -First 1)
    if ($q6.Count -eq 0 -or $q8.status -ne 'pass' -or [double]$q8.aggregate_rms -gt 1.25 * [double]$q6[0].aggregate_rms -or [double]$q8.max_channel_rms -gt 1.25 * [double]$q6[0].max_channel_rms) {
      $q8Gap = $true
    }
  }
}

$decision = if ($mode -eq 'quick') { 'fmm_e19_full_panel_quick_pass_linearity_pending' } else { 'fmm_e19_full_panel_matvec_pass_linearity_pending' }
if ($runFailed) {
  $decision = 'fmm_e19_runtime_limited'
} elseif ($coverageGap) {
  $decision = 'fmm_e19_full_coverage_gap'
} elseif ($auditGap) {
  $decision = 'fmm_e19_sampled_audit_gap'
} elseif ($fallbackDominated) {
  $decision = 'fmm_e19_fallback_dominated'
} elseif ($hashGap) {
  $decision = 'fmm_e19_hash_instability'
} elseif ($timingGap) {
  $decision = 'fmm_e19_timing_schema_gap'
} elseif ($q8Gap) {
  $decision = 'fmm_e19_sampled_audit_gap'
}

[pscustomobject]@{
  decision=$decision
  mode=$mode
  rows=$matrix.Count
  required_rows=$required.Count
  max_aggregate_rms=(@($required | Measure-Object -Property aggregate_rms -Maximum).Maximum)
  max_channel_rms=(@($required | Measure-Object -Property max_channel_rms -Maximum).Maximum)
  min_far_fraction=(@($required | Measure-Object -Property far_fraction -Minimum).Minimum)
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$mode",
  "rows=$($matrix.Count)",
  "required_rows=$($required.Count)",
  "matrix=$matrixPath",
  "sampled_audit=$auditPath",
  "coverage=$coveragePath",
  "runtime_cache=$runtimePath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$report = @(
  'TITAN FMM Engine E19 Full-Panel Matvec Report',
  '',
  "Decision: $decision",
  "Mode: $mode",
  "Rows: $($matrix.Count)",
  '',
  'E19 evaluates KIFMM-lite on all target panels and verifies the full-panel result with sampled dense production-direct rows.',
  'This remains diagnostic infrastructure only: no production TABI matvec, GMRES path, preconditioner, manuscript, or default behavior is changed.',
  'KIFMM rows use q6 level-cube and qr-rank-only by default. Condition sentinel fields are not interpreted when condition diagnostics are disabled.',
  'E19.1 bounds the result as a quick/full-panel matvec gate with linearity and repeatability still pending before GMRES.',
  'The measured bottleneck is full far evaluation rather than equivalent-system setup; target-side batched far evaluation is the next performance direction.',
  '',
  'Rows:',
  '  row_id,status,aggregate_rms,max_channel_rms,far_fraction,setup_seconds,apply_seconds,repeat_apply_seconds,far_eval_fraction,total_seconds'
)
foreach ($row in $matrix) {
  $farEvalFraction = if ([double]$row.total_seconds -gt 0.0) {
    [double]$row.full_far_eval_seconds / [double]$row.total_seconds
  } else {
    [double]::NaN
  }
  $report += ("  {0},{1},{2},{3},{4},{5},{6},{7},{8},{9}" -f $row.row_id,$row.status,$row.aggregate_rms,$row.max_channel_rms,$row.far_fraction,$row.setup_seconds,$row.apply_seconds,$row.repeat_apply_seconds,$farEvalFraction,$row.total_seconds)
}
$report += @(
  '',
  'Artifacts:',
  "  Matrix: $matrixPath",
  "  Sampled audit: $auditPath",
  "  Coverage: $coveragePath",
  "  Runtime/cache: $runtimePath",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$report | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E19 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -notin @('fmm_e19_full_panel_quick_pass_linearity_pending','fmm_e19_full_panel_matvec_pass_linearity_pending') -and -not $AllowFailedDecision) {
  throw "E19 did not pass an accepted full-panel decision: $decision"
}
