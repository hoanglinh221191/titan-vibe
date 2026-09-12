param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e17_cartesian_k4_gap'),
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
  throw 'No C++ compiler found for KIFMM E17 gate.'
}

function Build-GateExe([string]$Exe) {
  $compiler = Resolve-Cxx $Cxx
  $args = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_ENABLE_E15_CARTESIAN_COMPARATOR=1',
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
    $Exe
  )
  Push-Location $RepoRoot
  try {
    & $compiler @args
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E17 Cartesian K4 diagnostic driver' }
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

function RequiredD([hashtable]$Map, [string]$Key) {
  if (-not $Map.Contains($Key)) { throw "missing required numeric field: $Key" }
  return D $Map $Key
}

function RequiredS([hashtable]$Map, [string]$Key) {
  if (-not $Map.Contains($Key)) { throw "missing required string field: $Key" }
  return S $Map $Key
}

function Max-Channel([hashtable]$Map) {
  return (@('k1','k2','k3','k4') | ForEach-Object { D $Map $_ 0.0 } | Measure-Object -Maximum).Maximum
}

$quickSamples = @(
  [pscustomobject]@{ row_key='1A63_closest_opposing_surface_128'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ row_key='1SH1_uniform_128'; fixture='1SH1'; sample='uniform_hash'; sample_count=128 }
)
$fullSamples = @(
  [pscustomobject]@{ row_key='1SH1_uniform_128'; fixture='1SH1'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ row_key='1SH1_closest_opposing_surface_128'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ row_key='1SH1_high_aspect_128'; fixture='1SH1'; sample='high_aspect'; sample_count=128 },
  [pscustomobject]@{ row_key='1SH1_inward_normal_groove_128'; fixture='1SH1'; sample='inward_normal_groove'; sample_count=128 },
  [pscustomobject]@{ row_key='1A63_uniform_128'; fixture='1A63'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ row_key='1A63_uniform_512'; fixture='1A63'; sample='uniform_hash'; sample_count=512 },
  [pscustomobject]@{ row_key='1A63_high_aspect_128'; fixture='1A63'; sample='high_aspect'; sample_count=128 },
  [pscustomobject]@{ row_key='1A63_closest_opposing_surface_128'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 }
)
$quickVariants = @(
  [pscustomobject]@{ variant='baseline_p4_theta0p6'; order=4; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p4_theta0p5'; order=4; theta=0.5; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p4_theta0p4'; order=4; theta=0.4; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p5_theta0p6'; order=5; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p5_geometry0p4'; order=5; theta=0.6; tail=0.002; geom=0.4; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p4_geometry0p4'; order=4; theta=0.6; tail=0.002; geom=0.4; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p4_tail0p001'; order=4; theta=0.6; tail=0.001; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ variant='p4_green_jet_tiled'; order=4; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-tiled' }
)

$matrixPath = Join-Path $WorkRoot 'fmm_e17_matrix.csv'
$attrPath = Join-Path $WorkRoot 'fmm_e17_k4_attribution.csv'
$policyPath = Join-Path $WorkRoot 'fmm_e17_policy_sweep.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e17_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e17_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E17_CARTESIAN_K4_GAP_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e17_cartesian_k4_gap.exe'

Build-GateExe -Exe $exe

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($old in @(Import-Csv -LiteralPath $matrixPath)) { $existing[$old.row_id] = $old }
}

$samples = if ($Full -and -not $Quick) { $fullSamples } else { $quickSamples }
$variants = $quickVariants
if ($Full -and -not $Quick) {
  $quickRows = @()
  if (Test-Path -LiteralPath $matrixPath) {
    $quickRows = @(Import-Csv -LiteralPath $matrixPath | Where-Object {
        $_.row_id -match '^(1A63_closest_opposing_surface_128|1SH1_uniform_128)__'
      })
  }
  $bestVariantName = 'baseline_p4_theta0p6'
  $candidateRows = @($quickRows | Where-Object { $_.variant -ne 'baseline_p4_theta0p6' -and $_.cartesian_execution_status -eq 'pass' } |
    Sort-Object @{Expression={[double]$_.cartesian_k4}; Ascending=$true})
  if ($candidateRows.Count -gt 0) { $bestVariantName = $candidateRows[0].variant }
  $variants = @($quickVariants | Where-Object { $_.variant -eq 'baseline_p4_theta0p6' -or $_.variant -eq $bestVariantName })
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($sample in $samples) {
  foreach ($variant in $variants) {
    $rowId = "$($sample.row_key)__$($variant.variant)"
    if ($existing.ContainsKey($rowId)) {
      $old = $existing[$rowId]
      $rerun = $Force -or ($RerunRows -contains $rowId) -or
        (-not ($old.PSObject.Properties.Name -contains 'k4_abs_rms')) -or
        (-not ($old.PSObject.Properties.Name -contains 'cartesian_max_channel_rms')) -or
        (-not ($old.PSObject.Properties.Name -contains 'cartesian_max_channel_name')) -or
        ($RerunFailed -and ($old.cartesian_execution_status -ne 'pass'))
      if (-not $rerun) {
        Write-Host "resume: $rowId"
        $matrix.Add($old)
        continue
      }
    }
    Write-Host "run: $rowId"
    $log = Join-Path $WorkRoot ($rowId + '.log')
    $args = @(
      $Pqr1sh1,
      $Pqr1a63,
      '--e17-profile',
      "--fixture=$($sample.fixture)",
      "--sample=$($sample.sample)",
      "--sample-count=$($sample.sample_count)",
      "--cartesian-order=$($variant.order)",
      "--cartesian-theta=$($variant.theta)",
      "--cartesian-tail-tol=$($variant.tail)",
      "--cartesian-geometry-theta=$($variant.geom)",
      "--cartesian-cache-mode=$($variant.cache)"
    )
    Push-Location $RepoRoot
    try {
      $output = & $exe @args 2>&1
      $code = $LASTEXITCODE
    } finally {
      Pop-Location
    }
    $output | Set-Content -LiteralPath $log -Encoding ascii
    $kifmmLine = @($output | Where-Object { [string]$_ -match '^KIFMM_SAMPLED_PANEL_CASE ' } | Select-Object -First 1)
    $cartLine = @($output | Where-Object { [string]$_ -match '^CARTESIAN_SAMPLE_RESULT ' } | Select-Object -First 1)
    if ($code -ne 0 -or $cartLine.Count -eq 0 -or $kifmmLine.Count -eq 0) {
      $matrix.Add([pscustomobject]@{
        row_id=$rowId; fixture=$sample.fixture; sample_set=$sample.sample; sample_count=$sample.sample_count
        variant=$variant.variant; cartesian_order=$variant.order; theta=$variant.theta; tail_tol=$variant.tail; geometry_theta=$variant.geom; cache_mode=$variant.cache
        kifmm_status='run_failed'; cartesian_status='run_failed'; cartesian_execution_status='run_failed'; cartesian_accuracy_status='not_evaluated'; cartesian_accuracy_failure='not_evaluated'
        cartesian_aggregate_rms=[double]::NaN; cartesian_max_channel_rms=[double]::NaN; cartesian_max_channel_name='none'; cartesian_k4=[double]::NaN; k4_abs_rms=[double]::NaN; k4_reference_rms=[double]::NaN
        cartesian_k1=[double]::NaN; cartesian_k2=[double]::NaN; cartesian_k3=[double]::NaN
        far_panel_pairs=0; fallback_panel_pairs=0; expected_pairs=0; valid_pairs=0; invalid_pairs=0; nonproduction_fallback_pairs=0
        linf_max_abs_tail=[double]::NaN; linf_max_relative_tail=[double]::NaN; geometry_guard_rejected_cells=0; geometry_guard_direct_fallback_pairs=0; geometry_ratio_max=[double]::NaN
        fallback_refined_pairs=0; fallback_base_near_pairs=0; total_seconds=[double]::NaN; log=$log
      })
      continue
    }
    $k = Parse-KeyValueLine ([string]$kifmmLine[0])
    $c = Parse-KeyValueLine ([string]$cartLine[0])
    $sourcePanels = [uint64](U64 $c 'source_panels')
    $sampledRows = [uint64](U64 $c 'sampled_rows')
    $expectedPairs = $sourcePanels * $sampledRows
    $far = [uint64](U64 $c 'sampled_far_panel_pairs')
    $fallback = [uint64](U64 $c 'sampled_fallback_panel_pairs')
    $matrix.Add([pscustomobject]@{
      row_id=$rowId
      fixture=$sample.fixture
      sample_set=$sample.sample
      sample_count=$sampledRows
      variant=$variant.variant
      cartesian_order=(D $c 'order')
      theta=(D $c 'theta')
      tail_tol=(D $c 'tail_tol')
      geometry_theta=(D $c 'geometry_theta')
      cache_mode=(S $c 'cache_mode')
      kifmm_status=(S $k 'status')
      cartesian_status=(S $c 'status')
      cartesian_execution_status=(S $c 'execution_status')
      cartesian_accuracy_status=(S $c 'accuracy_status')
      cartesian_accuracy_failure=(S $c 'accuracy_failure')
      cartesian_aggregate_rms=(D $c 'aggregate_rms')
      cartesian_max_channel_rms=(RequiredD $c 'max_channel_rms')
      cartesian_max_channel_name=(RequiredS $c 'max_channel_name')
      cartesian_k1=(D $c 'k1')
      cartesian_k2=(D $c 'k2')
      cartesian_k3=(D $c 'k3')
      cartesian_k4=(D $c 'k4')
      k4_abs_rms=(D $c 'k4_abs_rms')
      k4_reference_rms=(D $c 'k4_reference_rms')
      far_panel_pairs=$far
      fallback_panel_pairs=$fallback
      expected_pairs=$expectedPairs
      valid_pairs=(U64 $c 'valid_pairs')
      invalid_pairs=(U64 $c 'invalid_pairs')
      nonproduction_fallback_pairs=(U64 $c 'nonproduction_fallback_pairs')
      linf_max_abs_tail=(D $c 'linf_max_abs_tail')
      linf_max_relative_tail=(D $c 'linf_max_relative_tail')
      geometry_guard_rejected_cells=(U64 $c 'geometry_guard_rejected_cells')
      geometry_guard_direct_fallback_pairs=(U64 $c 'geometry_guard_direct_fallback_pairs')
      geometry_ratio_max=(D $c 'geometry_ratio_max')
      fallback_refined_pairs=(U64 $c 'fallback_refined_pairs')
      fallback_base_near_pairs=(U64 $c 'fallback_base_near_pairs')
      total_seconds=(D $c 'total_seconds')
      log=$log
    })
  }
}

$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii

$attr = foreach ($m in $matrix) {
  $total = [double]([uint64]$m.expected_pairs)
  $farFraction = if ($total -gt 0.0) { [double]([uint64]$m.far_panel_pairs) / $total } else { [double]::NaN }
  $fallbackFraction = if ($total -gt 0.0) { [double]([uint64]$m.fallback_panel_pairs) / $total } else { [double]::NaN }
  [pscustomobject]@{
    row_id=$m.row_id; fixture=$m.fixture; sample_set=$m.sample_set; variant=$m.variant
    source_depth='mixed'; mac_bin='sampled_total'; attribution_scope='global_sampled_row_error'; pair_count=$m.expected_pairs
    far_fraction=$farFraction; fallback_fraction=$fallbackFraction
    k4_abs_rms=$m.k4_abs_rms; k4_reference_rms=$m.k4_reference_rms; k4_relative_rms=$m.cartesian_k4
    geometry_ratio_max=$m.geometry_ratio_max; linf_max_relative_tail=$m.linf_max_relative_tail
  }
  [pscustomobject]@{
    row_id=$m.row_id; fixture=$m.fixture; sample_set=$m.sample_set; variant=$m.variant
    source_depth='mixed'; mac_bin='accepted_far_work'; attribution_scope='work_bucket_pair_coverage'; pair_count=$m.far_panel_pairs
    far_fraction=$farFraction; fallback_fraction=0.0
    k4_abs_rms=$m.k4_abs_rms; k4_reference_rms=$m.k4_reference_rms; k4_relative_rms=$m.cartesian_k4
    geometry_ratio_max=$m.geometry_ratio_max; linf_max_relative_tail=$m.linf_max_relative_tail
  }
  [pscustomobject]@{
    row_id=$m.row_id; fixture=$m.fixture; sample_set=$m.sample_set; variant=$m.variant
    source_depth='mixed'; mac_bin='fallback_direct_work'; attribution_scope='work_bucket_pair_coverage'; pair_count=$m.fallback_panel_pairs
    far_fraction=0.0; fallback_fraction=$fallbackFraction
    k4_abs_rms=0.0; k4_reference_rms=$m.k4_reference_rms; k4_relative_rms=0.0
    geometry_ratio_max=$m.geometry_ratio_max; linf_max_relative_tail=$m.linf_max_relative_tail
  }
}
$attr | Export-Csv -LiteralPath $attrPath -NoTypeInformation -Encoding ascii

$policy = $matrix | Group-Object variant | ForEach-Object {
  $rows = @($_.Group)
  $baselineRows = @($matrix | Where-Object { $_.variant -eq 'baseline_p4_theta0p6' })
  [pscustomobject]@{
    variant=$_.Name
    rows=$rows.Count
    pass_rows=@($rows | Where-Object { $_.cartesian_accuracy_status -eq 'pass' }).Count
    mean_k4=($rows | Measure-Object -Property cartesian_k4 -Average).Average
    max_k4=($rows | Measure-Object -Property cartesian_k4 -Maximum).Maximum
    mean_max_channel_rms=($rows | Measure-Object -Property cartesian_max_channel_rms -Average).Average
    mean_total_seconds=($rows | Measure-Object -Property total_seconds -Average).Average
    mean_far_pairs=($rows | Measure-Object -Property far_panel_pairs -Average).Average
    mean_fallback_pairs=($rows | Measure-Object -Property fallback_panel_pairs -Average).Average
    baseline_mean_k4=($baselineRows | Measure-Object -Property cartesian_k4 -Average).Average
  }
}
$policy | Export-Csv -LiteralPath $policyPath -NoTypeInformation -Encoding ascii

$decision = 'fmm_e17_unresolved_cartesian_gap'
$baseline = @($matrix | Where-Object { $_.variant -eq 'baseline_p4_theta0p6' -and $_.cartesian_execution_status -eq 'pass' })
$baselineMeanK4 = ($baseline | Measure-Object -Property cartesian_k4 -Average).Average
$passVariant = @($policy | Where-Object { [int]$_.pass_rows -eq [int]$_.rows -and [int]$_.rows -gt 0 } | Select-Object -First 1)
if (@($matrix | Where-Object { $_.cartesian_execution_status -eq 'run_failed' }).Count -gt 0) {
  $decision = 'fmm_e17_runtime_limited'
} elseif ($passVariant.Count -gt 0) {
  $decision = 'fmm_e17_cartesian_policy_recovers_accuracy'
} elseif ($baselineMeanK4 -gt 0) {
  $p5Geom = @($policy | Where-Object { $_.variant -eq 'p5_geometry0p4' } | Select-Object -First 1)
  $p5 = @($policy | Where-Object { $_.variant -eq 'p5_theta0p6' } | Select-Object -First 1)
  $geom = @($policy | Where-Object { $_.variant -eq 'p4_geometry0p4' } | Select-Object -First 1)
  $tail = @($policy | Where-Object { $_.variant -eq 'p4_tail0p001' } | Select-Object -First 1)
  $theta = @($policy | Where-Object { $_.variant -in @('p4_theta0p5','p4_theta0p4') } | Sort-Object {[double]$_.mean_k4} | Select-Object -First 1)
  $tiled = @($policy | Where-Object { $_.variant -eq 'p4_green_jet_tiled' } | Select-Object -First 1)
  if ($p5Geom.Count -gt 0 -and [double]$p5Geom[0].mean_k4 -le $baselineMeanK4 / 4.0) {
    $decision = 'fmm_e17_k4_derivative_order_gap'
  } elseif ($p5.Count -gt 0 -and [double]$p5[0].mean_k4 -le $baselineMeanK4 / 4.0) {
    $decision = 'fmm_e17_k4_derivative_order_gap'
  } elseif ($geom.Count -gt 0 -and [double]$geom[0].mean_k4 -le $baselineMeanK4 / 4.0) {
    $decision = 'fmm_e17_geometry_guard_gap'
  } elseif ($tail.Count -gt 0 -and [double]$tail[0].mean_k4 -le $baselineMeanK4 / 4.0) {
    $decision = 'fmm_e17_tail_estimator_gap'
  } elseif ($theta.Count -gt 0 -and [double]$theta[0].mean_k4 -le $baselineMeanK4 / 4.0) {
    $decision = 'fmm_e17_k4_high_mac_far_gap'
  } elseif ($tiled.Count -gt 0 -and [math]::Abs([double]$tiled[0].mean_k4 - $baselineMeanK4) -le [math]::Max(1.0e-12, 1.0e-6 * $baselineMeanK4)) {
    $decision = 'fmm_e17_cache_layout_neutral'
  }
}

[pscustomobject]@{
  decision=$decision
  mode=($(if ($Full -and -not $Quick) { 'full' } else { 'quick' }))
  rows=$matrix.Count
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$matrixPath",
  "k4_attribution=$attrPath",
  "policy_sweep=$policyPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$worst = @($matrix | Sort-Object { [double]$_.cartesian_k4 } -Descending | Select-Object -First 1)
$recovery = @($policy | Where-Object { [int]$_.pass_rows -eq [int]$_.rows -and [int]$_.rows -gt 0 } | Select-Object -First 1)
$report = @(
  'TITAN FMM Engine E17 Cartesian K4 Gap Diagnostic Report',
  '',
  "Decision: $decision",
  "Mode: $(if ($Full -and -not $Quick) { 'full' } else { 'quick' })",
  "Rows: $($matrix.Count)",
  '',
  'E17 diagnoses the guarded Cartesian sampled comparator only. KIFMM-lite remains the E16 passing subgate and no production TABI or GMRES path is changed.',
  'The comparison is KIFMM q6 level-cube versus guarded Cartesian p4/p5 sampled comparator, not a same-order comparison.',
  "Worst Cartesian K4 row: $($worst[0].row_id), K4_RMS=$($worst[0].cartesian_k4), max_channel=$($worst[0].cartesian_max_channel_name), max_channel_RMS=$($worst[0].cartesian_max_channel_rms), status=$($worst[0].cartesian_status)",
  '',
  'Policy sweep:',
  '  variant,rows,pass_rows,mean_k4,max_k4,mean_max_channel_rms,mean_total_seconds,mean_far_pairs,mean_fallback_pairs'
)
foreach ($p in $policy) {
  $report += ("  {0},{1},{2},{3},{4},{5},{6},{7},{8}" -f $p.variant,$p.rows,$p.pass_rows,$p.mean_k4,$p.max_k4,$p.mean_max_channel_rms,$p.mean_total_seconds,$p.mean_far_pairs,$p.mean_fallback_pairs)
}
if ($recovery.Count -gt 0) {
  $report += @(
    '',
    ("Recovery candidate: {0} passed {1}/{2} rows with max_k4={3}, mean_total_seconds={4}, mean_far_pairs={5}, mean_fallback_pairs={6}." -f $recovery[0].variant,$recovery[0].pass_rows,$recovery[0].rows,$recovery[0].max_k4,$recovery[0].mean_total_seconds,$recovery[0].mean_far_pairs,$recovery[0].mean_fallback_pairs),
    'This is a sampled-row diagnostic result only; E18 should compare the recovered Cartesian policy against KIFMM q6/q8 before any production wiring.'
  )
}
$report += @(
  '',
  'Artifacts:',
  "  Matrix: $matrixPath",
  "  K4 attribution: $attrPath",
  "  Policy sweep: $policyPath",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$report | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E17 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -notin @('fmm_e17_cartesian_policy_recovers_accuracy','fmm_e17_k4_high_mac_far_gap','fmm_e17_k4_derivative_order_gap','fmm_e17_geometry_guard_gap','fmm_e17_tail_estimator_gap','fmm_e17_cache_layout_neutral') -and -not $AllowFailedDecision) {
  throw "E17 did not pass/classify accepted diagnostic: $decision"
}
