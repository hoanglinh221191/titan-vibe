param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e18_k4_closure'),
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
  throw 'No C++ compiler found for KIFMM E18 gate.'
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E18 K4 closure comparator driver' }
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
  if ($extra.Count -gt 0) {
    throw "stale duplicate source/header for ${Prefix}: $($extra -join ', ')"
  }
}

function Assert-FileContains([string]$Path, [string]$Needle) {
  $text = Get-Content -Raw -LiteralPath $Path
  if (-not $text.Contains($Needle)) { throw "expected '$Needle' in $Path" }
}

function Assert-FileNotContains([string]$Path, [string]$Needle) {
  $text = Get-Content -Raw -LiteralPath $Path
  if ($text.Contains($Needle)) { throw "stale pattern '$Needle' found in $Path" }
}

function Assert-E18CanonicalSources {
  foreach ($prefix in @(
      'tabi_fmm_kifmm_lite',
      'tabi_fmm_kifmm_sampled_panel_matvec',
      'tabi_fmm_kifmm_real_panel_cloud',
      'tabi_fmm_cartesian_sampled_panel_matvec')) {
    Assert-CanonicalPair $prefix
  }
  $tabiDir = Join-Path $RepoRoot 'src\tabi'
  $kifmmHeader = Join-Path $tabiDir 'tabi_fmm_kifmm_lite.h'
  $cartesianHeader = Join-Path $tabiDir 'tabi_fmm_cartesian_sampled_panel_matvec.h'
  Assert-FileContains $kifmmHeader 'condition_policy{"qr-rank-only"}'
  Assert-FileContains $kifmmHeader 'source_box_policy'
  Assert-FileContains $kifmmHeader 'min_far_panel_fraction'
  Assert-FileNotContains $kifmmHeader 'condition_policy{"svd-full-debug"}'
  Assert-FileContains $cartesianHeader 'eps_ratio{80.0}'
  Assert-FileContains $cartesianHeader 'max_aggregate_relative_rms'
  Assert-FileContains $cartesianHeader 'max_channel_relative_rms'
  Assert-FileNotContains $cartesianHeader 'eps_ratio{1.0 / 80.0}'
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

function Max-Channel([hashtable]$Map) {
  return (@('k1','k2','k3','k4') | ForEach-Object { D $Map $_ 0.0 } | Measure-Object -Maximum).Maximum
}

function Max-Channel-Name([hashtable]$Map) {
  $bestName = 'K1'
  $bestValue = D $Map 'k1' 0.0
  foreach ($entry in @(@('K2','k2'),@('K3','k3'),@('K4','k4'))) {
    $value = D $Map $entry[1] 0.0
    if ($value -gt $bestValue) {
      $bestValue = $value
      $bestName = $entry[0]
    }
  }
  return $bestName
}

function Is-Finite([double]$Value) {
  return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Median([object[]]$Values) {
  $items = @($Values | Where-Object { $_ -ne $null -and (Is-Finite ([double]$_)) } | ForEach-Object { [double]$_ } | Sort-Object)
  if ($items.Count -eq 0) { return [double]::NaN }
  $mid = [int]($items.Count / 2)
  if (($items.Count % 2) -eq 1) { return $items[$mid] }
  return 0.5 * ($items[$mid - 1] + $items[$mid])
}

$quickSamples = @(
  [pscustomobject]@{ key='1A63_closest_opposing_surface_128'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ key='1SH1_uniform_128'; fixture='1SH1'; sample='uniform_hash'; sample_count=128 }
)
$fullSamples = @(
  [pscustomobject]@{ key='1SH1_uniform_128'; fixture='1SH1'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ key='1SH1_closest_opposing_surface_128'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ key='1SH1_high_aspect_128'; fixture='1SH1'; sample='high_aspect'; sample_count=128 },
  [pscustomobject]@{ key='1SH1_inward_normal_groove_128'; fixture='1SH1'; sample='inward_normal_groove'; sample_count=128 },
  [pscustomobject]@{ key='1A63_uniform_128'; fixture='1A63'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ key='1A63_uniform_512'; fixture='1A63'; sample='uniform_hash'; sample_count=512 },
  [pscustomobject]@{ key='1A63_high_aspect_128'; fixture='1A63'; sample='high_aspect'; sample_count=128 },
  [pscustomobject]@{ key='1A63_closest_opposing_surface_128'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 }
)
$quickVariants = @(
  [pscustomobject]@{ name='p4_baseline'; order=4; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_theta0p6'; order=5; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_geometry0p4'; order=5; theta=0.6; tail=0.002; geom=0.4; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_tail0p001'; order=5; theta=0.6; tail=0.001; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_geometry0p4_tail0p001'; order=5; theta=0.6; tail=0.001; geom=0.4; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_theta0p4'; order=5; theta=0.4; tail=0.002; geom=0.5; cache='green-jet-compressed' }
)
$fullVariants = @(
  [pscustomobject]@{ name='p4_baseline'; order=4; theta=0.6; tail=0.002; geom=0.5; cache='green-jet-compressed' },
  [pscustomobject]@{ name='p5_geometry0p4'; order=5; theta=0.6; tail=0.002; geom=0.4; cache='green-jet-compressed' }
)
$q8Keys = @('1A63_closest_opposing_surface_128','1A63_uniform_512')

$mode = if ($Full -and -not $Quick) { 'full' } else { 'quick' }
$samples = if ($mode -eq 'full') { $fullSamples } else { $quickSamples }
$variants = if ($mode -eq 'full') { $fullVariants } else { $quickVariants }

$matrixPath = Join-Path $WorkRoot 'fmm_e18_matrix.csv'
$policyPath = Join-Path $WorkRoot 'fmm_e18_policy_comparison.csv'
$channelPath = Join-Path $WorkRoot 'fmm_e18_channel_summary.csv'
$runtimePath = Join-Path $WorkRoot 'fmm_e18_runtime_fallback.csv'
$q8Path = Join-Path $WorkRoot 'fmm_e18_q8_convergence.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e18_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e18_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E18_K4_CLOSURE_COMPARATOR_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e18_k4_closure.exe'

Assert-E18CanonicalSources
Build-GateExe -Exe $exe

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($old in @(Import-Csv -LiteralPath $matrixPath)) { $existing[$old.row_id] = $old }
}

$jobs = New-Object System.Collections.Generic.List[object]
foreach ($sample in $samples) {
  foreach ($variant in $variants) {
    $jobs.Add([pscustomobject]@{ sample=$sample; variant=$variant; kifmm_order=6 })
  }
  if ($q8Keys -contains $sample.key) {
    $jobs.Add([pscustomobject]@{ sample=$sample; variant=($quickVariants | Where-Object { $_.name -eq 'p5_geometry0p4' } | Select-Object -First 1); kifmm_order=8 })
  }
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($job in $jobs) {
  $sample = $job.sample
  $variant = $job.variant
  $rowId = "$($sample.key)__q$($job.kifmm_order)__$($variant.name)"
  if ($existing.ContainsKey($rowId)) {
    $old = $existing[$rowId]
    $rerun = $Force -or ($RerunRows -contains $rowId) -or
      (-not ($old.PSObject.Properties.Name -contains 'cartesian_max_channel_name')) -or
      (-not ($old.PSObject.Properties.Name -contains 'kifmm_condition_policy')) -or
      ($RerunFailed -and ($old.kifmm_status -ne 'pass' -or $old.cartesian_execution_status -ne 'pass'))
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
    '--e18-profile',
    "--fixture=$($sample.fixture)",
    "--sample=$($sample.sample)",
    "--sample-count=$($sample.sample_count)",
    "--kifmm-order=$($job.kifmm_order)",
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
  if ($code -ne 0 -or $kifmmLine.Count -eq 0 -or $cartLine.Count -eq 0) {
    $matrix.Add([pscustomobject]@{
      row_id=$rowId; mode=$mode; fixture=$sample.fixture; sample_set=$sample.sample; sample_count=$sample.sample_count
      kifmm_order=$job.kifmm_order; cartesian_variant=$variant.name; cartesian_order=$variant.order; cartesian_theta=$variant.theta; cartesian_tail_tol=$variant.tail; cartesian_geometry_theta=$variant.geom; cartesian_cache_mode=$variant.cache
      kifmm_status='run_failed'; cartesian_status='run_failed'; cartesian_execution_status='run_failed'; cartesian_accuracy_status='not_evaluated'; cartesian_accuracy_failure='not_evaluated'
      kifmm_condition_policy=''; kifmm_source_box_policy=''; kifmm_eps_ratio=[double]::NaN; cartesian_eps_ratio=[double]::NaN; kifmm_kappa=[double]::NaN; cartesian_kappa=[double]::NaN; kifmm_near_integration=''; cartesian_near_integration=''; kifmm_quad_buffer=''; cartesian_quad_buffer=''
      kifmm_aggregate_rms=[double]::NaN; kifmm_max_channel_rms=[double]::NaN; kifmm_max_channel_name='none'; kifmm_k4=[double]::NaN
      cartesian_aggregate_rms=[double]::NaN; cartesian_max_channel_rms=[double]::NaN; cartesian_max_channel_name='none'; cartesian_k4=[double]::NaN
      kifmm_total_seconds=[double]::NaN; cartesian_total_seconds=[double]::NaN; kifmm_fallback_pairs=0; cartesian_fallback_pairs=0; kifmm_far_pairs=0; cartesian_far_pairs=0
      expected_pairs=0; kifmm_valid_pairs=0; cartesian_valid_pairs=0; invalid_pairs=0; nonproduction_fallback_pairs=0; coverage_exact='false'; parameter_match='false'; sampled_target_hash_match='false'; log=$log
    })
    continue
  }

  $k = Parse-KeyValueLine ([string]$kifmmLine[0])
  $c = Parse-KeyValueLine ([string]$cartLine[0])
  $sourcePanels = [uint64](U64 $c 'source_panels')
  $sampledRows = [uint64](U64 $c 'sampled_rows')
  $expectedPairs = $sourcePanels * $sampledRows
  $kifmmValid = [uint64](U64 $k 'production_direct_reference_valid_pairs')
  $cartValid = [uint64](U64 $c 'valid_pairs')
  $invalidPairs = [uint64](U64 $k 'invalid_pairs') + [uint64](U64 $c 'invalid_pairs')
  $nonProductionPairs = [uint64](U64 $k 'nonproduction_fallback_pairs') + [uint64](U64 $c 'nonproduction_fallback_pairs')
  $coverageExact = ($kifmmValid -eq $expectedPairs -and $cartValid -eq $expectedPairs)
  $parameterMatch = (
    [math]::Abs((D $k 'eps_ratio') - (D $c 'eps_ratio')) -le 1.0e-12 -and
    [math]::Abs((D $k 'kappa') - (D $c 'kappa')) -le 1.0e-12 -and
    (S $k 'near_integration') -eq (S $c 'near_integration') -and
    (S $k 'quad_buffer') -eq (S $c 'quad_buffer')
  )
  $sampledTargetHashMatch = ((S $k 'sampled_target_hash') -eq (S $c 'sampled_target_hash'))
  $matrix.Add([pscustomobject]@{
    row_id=$rowId
    mode=$mode
    fixture=$sample.fixture
    sample_set=$sample.sample
    sample_count=$sampledRows
    kifmm_order=$job.kifmm_order
    cartesian_variant=$variant.name
    cartesian_order=(D $c 'order')
    cartesian_theta=(D $c 'theta')
    cartesian_tail_tol=(D $c 'tail_tol')
    cartesian_geometry_theta=(D $c 'geometry_theta')
    cartesian_cache_mode=(S $c 'cache_mode')
    kifmm_status=(S $k 'status')
    cartesian_status=(S $c 'status')
    cartesian_execution_status=(S $c 'execution_status')
    cartesian_accuracy_status=(S $c 'accuracy_status')
    cartesian_accuracy_failure=(S $c 'accuracy_failure')
    kifmm_condition_policy=(S $k 'condition_policy')
    kifmm_source_box_policy=(S $k 'source_box_policy')
    kifmm_eps_ratio=(D $k 'eps_ratio')
    cartesian_eps_ratio=(D $c 'eps_ratio')
    kifmm_kappa=(D $k 'kappa')
    cartesian_kappa=(D $c 'kappa')
    kifmm_near_integration=(S $k 'near_integration')
    cartesian_near_integration=(S $c 'near_integration')
    kifmm_quad_buffer=(S $k 'quad_buffer')
    cartesian_quad_buffer=(S $c 'quad_buffer')
    kifmm_aggregate_rms=(D $k 'aggregate_rms')
    kifmm_max_channel_rms=(Max-Channel $k)
    kifmm_max_channel_name=(Max-Channel-Name $k)
    kifmm_k1=(D $k 'k1')
    kifmm_k2=(D $k 'k2')
    kifmm_k3=(D $k 'k3')
    kifmm_k4=(D $k 'k4')
    cartesian_aggregate_rms=(D $c 'aggregate_rms')
    cartesian_max_channel_rms=(D $c 'max_channel_rms')
    cartesian_max_channel_name=(S $c 'max_channel_name')
    cartesian_k1=(D $c 'k1')
    cartesian_k2=(D $c 'k2')
    cartesian_k3=(D $c 'k3')
    cartesian_k4=(D $c 'k4')
    cartesian_k4_abs_rms=(D $c 'k4_abs_rms')
    cartesian_k4_reference_rms=(D $c 'k4_reference_rms')
    kifmm_total_seconds=(D $k 'total_seconds')
    cartesian_total_seconds=(D $c 'total_seconds')
    kifmm_far_pairs=(U64 $k 'far_panel_pairs')
    kifmm_fallback_pairs=(U64 $k 'fallback_panel_pairs')
    cartesian_far_pairs=(U64 $c 'sampled_far_panel_pairs')
    cartesian_fallback_pairs=(U64 $c 'sampled_fallback_panel_pairs')
    cartesian_geometry_guard_direct_fallback_pairs=(U64 $c 'geometry_guard_direct_fallback_pairs')
    cartesian_linf_max_abs_tail=(D $c 'linf_max_abs_tail')
    cartesian_linf_max_relative_tail=(D $c 'linf_max_relative_tail')
    expected_pairs=$expectedPairs
    kifmm_valid_pairs=$kifmmValid
    cartesian_valid_pairs=$cartValid
    invalid_pairs=$invalidPairs
    nonproduction_fallback_pairs=$nonProductionPairs
    coverage_exact=([string]$coverageExact).ToLowerInvariant()
    parameter_match=([string]$parameterMatch).ToLowerInvariant()
    sampled_target_hash_match=([string]$sampledTargetHashMatch).ToLowerInvariant()
    log=$log
  })
}

$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii

$policy = $matrix | Where-Object { [int]$_.kifmm_order -eq 6 } | Group-Object cartesian_variant | ForEach-Object {
  $rows = @($_.Group)
  [pscustomobject]@{
    cartesian_variant=$_.Name
    rows=$rows.Count
    cartesian_pass_rows=@($rows | Where-Object { $_.cartesian_accuracy_status -eq 'pass' }).Count
    kifmm_pass_rows=@($rows | Where-Object { $_.kifmm_status -eq 'pass' }).Count
    cartesian_max_channel_max=($rows | Measure-Object -Property cartesian_max_channel_rms -Maximum).Maximum
    cartesian_k4_max=($rows | Measure-Object -Property cartesian_k4 -Maximum).Maximum
    cartesian_median_seconds=(Median @($rows | ForEach-Object { $_.cartesian_total_seconds }))
    kifmm_median_seconds=(Median @($rows | ForEach-Object { $_.kifmm_total_seconds }))
    cartesian_median_fallback_pairs=(Median @($rows | ForEach-Object { $_.cartesian_fallback_pairs }))
    kifmm_median_fallback_pairs=(Median @($rows | ForEach-Object { $_.kifmm_fallback_pairs }))
  }
}
$policy | Export-Csv -LiteralPath $policyPath -NoTypeInformation -Encoding ascii

$channels = foreach ($row in $matrix) {
  foreach ($channel in @('k1','k2','k3','k4')) {
    [pscustomobject]@{
      row_id=$row.row_id; fixture=$row.fixture; sample_set=$row.sample_set; kifmm_order=$row.kifmm_order; cartesian_variant=$row.cartesian_variant; channel=$channel.ToUpperInvariant()
      kifmm_rms=$row.("kifmm_$channel")
      cartesian_rms=$row.("cartesian_$channel")
    }
  }
}
$channels | Export-Csv -LiteralPath $channelPath -NoTypeInformation -Encoding ascii

$runtime = $matrix | Select-Object row_id,fixture,sample_set,kifmm_order,cartesian_variant,kifmm_total_seconds,cartesian_total_seconds,kifmm_far_pairs,kifmm_fallback_pairs,cartesian_far_pairs,cartesian_fallback_pairs,cartesian_geometry_guard_direct_fallback_pairs,cartesian_linf_max_abs_tail,cartesian_linf_max_relative_tail
$runtime | Export-Csv -LiteralPath $runtimePath -NoTypeInformation -Encoding ascii

$q8Rows = @($matrix | Where-Object { [int]$_.kifmm_order -eq 8 })
$q8 = foreach ($row in $q8Rows) {
  $q6 = @($matrix | Where-Object {
      $_.fixture -eq $row.fixture -and $_.sample_set -eq $row.sample_set -and [int]$_.sample_count -eq [int]$row.sample_count -and [int]$_.kifmm_order -eq 6 -and $_.cartesian_variant -eq 'p5_geometry0p4'
    } | Select-Object -First 1)
  $aggRatio = [double]::NaN
  $channelRatio = [double]::NaN
  $stable = 'false'
  if ($q6.Count -gt 0) {
    if ([double]$q6[0].kifmm_aggregate_rms -gt 0.0) { $aggRatio = [double]$row.kifmm_aggregate_rms / [double]$q6[0].kifmm_aggregate_rms }
    if ([double]$q6[0].kifmm_max_channel_rms -gt 0.0) { $channelRatio = [double]$row.kifmm_max_channel_rms / [double]$q6[0].kifmm_max_channel_rms }
    $stable = ([string]($row.kifmm_status -eq 'pass' -and $aggRatio -le 1.25 -and $channelRatio -le 1.25)).ToLowerInvariant()
  }
  [pscustomobject]@{
    row_id=$row.row_id
    q6_row_id=$(if ($q6.Count -gt 0) { $q6[0].row_id } else { '' })
    fixture=$row.fixture
    sample_set=$row.sample_set
    sample_count=$row.sample_count
    q6_aggregate_rms=$(if ($q6.Count -gt 0) { $q6[0].kifmm_aggregate_rms } else { [double]::NaN })
    q8_aggregate_rms=$row.kifmm_aggregate_rms
    q6_max_channel_rms=$(if ($q6.Count -gt 0) { $q6[0].kifmm_max_channel_rms } else { [double]::NaN })
    q8_max_channel_rms=$row.kifmm_max_channel_rms
    aggregate_ratio=$aggRatio
    max_channel_ratio=$channelRatio
    q8_stable=$stable
  }
}
$q8 | Export-Csv -LiteralPath $q8Path -NoTypeInformation -Encoding ascii

$candidateRows = @($matrix | Where-Object { [int]$_.kifmm_order -eq 6 -and $_.cartesian_variant -eq 'p5_geometry0p4' })
$runFailed = @($matrix | Where-Object { $_.kifmm_status -eq 'run_failed' -or $_.cartesian_execution_status -eq 'run_failed' }).Count -gt 0
$coverageGap = @($matrix | Where-Object { $_.coverage_exact -ne 'true' -or $_.parameter_match -ne 'true' -or $_.sampled_target_hash_match -ne 'true' -or [uint64]$_.invalid_pairs -ne 0 -or [uint64]$_.nonproduction_fallback_pairs -ne 0 }).Count -gt 0
$kifmmGap = @($candidateRows | Where-Object { $_.kifmm_status -ne 'pass' -or [double]$_.kifmm_aggregate_rms -gt 1.0e-6 -or [double]$_.kifmm_max_channel_rms -gt 5.0e-6 -or [double]$_.kifmm_k4 -gt 5.0e-6 }).Count -gt 0
$cartGap = @($candidateRows | Where-Object { $_.cartesian_accuracy_status -ne 'pass' -or [double]$_.cartesian_aggregate_rms -gt 1.0e-6 -or [double]$_.cartesian_max_channel_rms -gt 5.0e-6 -or [double]$_.cartesian_k4 -gt 5.0e-6 }).Count -gt 0
$q8Gap = @($q8 | Where-Object { $_.q8_stable -ne 'true' }).Count -gt 0
$cartMedian = Median @($candidateRows | ForEach-Object { $_.cartesian_total_seconds })
$kifmmMedian = Median @($candidateRows | ForEach-Object { $_.kifmm_total_seconds })
$cartFallbackMedian = Median @($candidateRows | ForEach-Object { $_.cartesian_fallback_pairs })
$kifmmFallbackMedian = Median @($candidateRows | ForEach-Object { $_.kifmm_fallback_pairs })
$runtimeRatio = if ((Is-Finite $kifmmMedian) -and $kifmmMedian -gt 0.0) { $cartMedian / $kifmmMedian } else { [double]::NaN }
$fallbackRatio = if ((Is-Finite $kifmmFallbackMedian) -and $kifmmFallbackMedian -gt 0.0) { $cartFallbackMedian / $kifmmFallbackMedian } else { [double]::NaN }

$decision = 'fmm_e18_k4_closure_kifmm_preferred'
if ($runFailed) {
  $decision = 'fmm_e18_runtime_limited'
} elseif ($coverageGap) {
  $decision = 'fmm_e18_runtime_limited'
} elseif ($kifmmGap) {
  $decision = 'fmm_e18_kifmm_accuracy_gap'
} elseif ($cartGap) {
  $decision = 'fmm_e18_cartesian_accuracy_gap'
} elseif ($q8Gap) {
  $decision = 'fmm_e18_q8_convergence_gap'
} elseif ($runtimeRatio -le 1.25 -and $fallbackRatio -le 2.0) {
  $decision = 'fmm_e18_k4_closure_cartesian_competitive'
}

[pscustomobject]@{
  decision=$decision
  mode=$mode
  rows=$matrix.Count
  candidate_rows=$candidateRows.Count
  cartesian_median_seconds=$cartMedian
  kifmm_median_seconds=$kifmmMedian
  runtime_ratio=$runtimeRatio
  cartesian_median_fallback_pairs=$cartFallbackMedian
  kifmm_median_fallback_pairs=$kifmmFallbackMedian
  fallback_ratio=$fallbackRatio
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "mode=$mode",
  "rows=$($matrix.Count)",
  "candidate_rows=$($candidateRows.Count)",
  "runtime_ratio=$runtimeRatio",
  "fallback_ratio=$fallbackRatio",
  "matrix=$matrixPath",
  "policy=$policyPath",
  "channels=$channelPath",
  "runtime=$runtimePath",
  "q8=$q8Path"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$report = @(
  'TITAN FMM Engine E18 K4 Closure Comparator Report',
  '',
  "Decision: $decision",
  "Mode: $mode",
  "Rows: $($matrix.Count)",
  '',
  'E18 compares KIFMM q6/q8 level-cube against true sampled-row guarded Cartesian policies on identical sampled target rows and dense production-direct references.',
  'The preferred-engine classification is limited to the sampled-row diagnostic harness; it is not a preferred production backend claim.',
  'No production TABI matvec, GMRES path, preconditioner, manuscript, or default behavior is changed.',
  'KIFMM rows use qr-rank-only condition policy in the speed path. Condition diagnostics are disabled in that mode; sigma and condition sentinel fields are not interpreted as singular equivalent systems.',
  'Cartesian p5_geometry0p4 is treated as the Cartesian closure baseline for this diagnostic branch only and does not update the manuscript protocol.',
  '',
  "Candidate p5_geometry0p4 rows: $($candidateRows.Count)",
  "Median Cartesian seconds: $cartMedian",
  "Median KIFMM q6 seconds: $kifmmMedian",
  "Runtime ratio Cartesian/KIFMM: $runtimeRatio",
  "Fallback ratio Cartesian/KIFMM: $fallbackRatio",
  '',
  'Policy comparison:',
  '  cartesian_variant,rows,cartesian_pass_rows,kifmm_pass_rows,cartesian_max_channel_max,cartesian_k4_max,cartesian_median_seconds,kifmm_median_seconds,cartesian_median_fallback_pairs,kifmm_median_fallback_pairs'
)
foreach ($p in $policy) {
  $report += ("  {0},{1},{2},{3},{4},{5},{6},{7},{8},{9}" -f $p.cartesian_variant,$p.rows,$p.cartesian_pass_rows,$p.kifmm_pass_rows,$p.cartesian_max_channel_max,$p.cartesian_k4_max,$p.cartesian_median_seconds,$p.kifmm_median_seconds,$p.cartesian_median_fallback_pairs,$p.kifmm_median_fallback_pairs)
}
$report += @(
  '',
  'Artifacts:',
  "  Matrix: $matrixPath",
  "  Policy comparison: $policyPath",
  "  Channel summary: $channelPath",
  "  Runtime/fallback: $runtimePath",
  "  Q8 convergence: $q8Path",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$report | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E18 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -notin @('fmm_e18_k4_closure_cartesian_competitive','fmm_e18_k4_closure_kifmm_preferred') -and -not $AllowFailedDecision) {
  throw "E18 did not pass an accepted comparator decision: $decision"
}
