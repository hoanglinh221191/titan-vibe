param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e16_cartesian_true_sampled_comparator'),
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
  throw 'No C++ compiler found for KIFMM E16 gate.'
}

function Build-GateExe([string]$Exe) {
  $compiler = Resolve-Cxx $Cxx
  $args = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_ENABLE_E16_CARTESIAN_COMPARATOR=1',
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E16 true sampled Cartesian comparator gate driver' }
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

function Max-Channel([hashtable]$Map) {
  return (@('k1','k2','k3','k4') | ForEach-Object { D $Map $_ 0.0 } | Measure-Object -Maximum).Maximum
}

function Accuracy-Status([string]$Status, [double]$AggregateRms, [double]$MaxChannelRms) {
  if ($Status -ne 'pass') { return 'unavailable' }
  $aggregateGap = $AggregateRms -gt 1.0e-6
  $channelGap = $MaxChannelRms -gt 5.0e-6
  if ($aggregateGap -and $channelGap) { return 'aggregate_and_channel_gap' }
  if ($aggregateGap) { return 'aggregate_gap' }
  if ($channelGap) { return 'channel_gap' }
  return 'pass'
}

function Nearly-Equal([double]$A, [double]$B, [double]$Tol = 1.0e-12) {
  return [math]::Abs($A - $B) -le $Tol
}

$quickRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_128_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ row_id='1A63_high_aspect_128_q6'; fixture='1A63'; sample='high_aspect'; sample_count=128 }
)
$fullRows = @(
  $quickRows[0],
  [pscustomobject]@{ row_id='1SH1_closest_opposing_surface_128_q6'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ row_id='1SH1_high_aspect_128_q6'; fixture='1SH1'; sample='high_aspect'; sample_count=128 },
  [pscustomobject]@{ row_id='1SH1_inward_normal_groove_128_q6'; fixture='1SH1'; sample='inward_normal_groove'; sample_count=128 },
  [pscustomobject]@{ row_id='1A63_uniform_128_q6'; fixture='1A63'; sample='uniform_hash'; sample_count=128 },
  [pscustomobject]@{ row_id='1A63_uniform_512_q6'; fixture='1A63'; sample='uniform_hash'; sample_count=512 },
  $quickRows[1],
  [pscustomobject]@{ row_id='1A63_closest_opposing_surface_128_q6'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 }
)
$rows = if ($Quick -and -not $Full) { $quickRows } else { $fullRows }

$matrixPath = Join-Path $WorkRoot 'fmm_e16_matrix.csv'
$coveragePath = Join-Path $WorkRoot 'fmm_e16_coverage.csv'
$runtimePath = Join-Path $WorkRoot 'fmm_e16_runtime.csv'
$channelPath = Join-Path $WorkRoot 'fmm_e16_channel_summary.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e16_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e16_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E16_CARTESIAN_TRUE_SAMPLED_COMPARATOR_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e16_cartesian_true_sampled_comparator_gate.exe'

Build-GateExe -Exe $exe

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($old in @(Import-Csv -LiteralPath $matrixPath)) { $existing[$old.row_id] = $old }
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
  if ($existing.ContainsKey($row.row_id)) {
    $old = $existing[$row.row_id]
    $rerun = $Force -or ($RerunRows -contains $row.row_id) -or
      (-not ($old.PSObject.Properties.Name -contains 'cartesian_timing_scope')) -or
      (-not ($old.PSObject.Properties.Name -contains 'cartesian_execution_status')) -or
      ($RerunFailed -and ($old.kifmm_status -ne 'pass' -or $old.cartesian_status -ne 'pass'))
    if (-not $rerun) {
      Write-Host "resume: $($row.row_id)"
      $matrix.Add($old)
      continue
    }
  }

  Write-Host "run: $($row.row_id)"
  $log = Join-Path $WorkRoot ($row.row_id + '.log')
  $args = @(
    $Pqr1sh1,
    $Pqr1a63,
    '--e16-profile',
    "--fixture=$($row.fixture)",
    "--sample=$($row.sample)",
    "--sample-count=$($row.sample_count)"
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
      row_id=$row.row_id; fixture=$row.fixture; sample_set=$row.sample; sample_count=$row.sample_count; order=6
      kifmm_status='run_failed'; cartesian_status='run_failed'; cartesian_execution_status='run_failed'; cartesian_accuracy_status='not_evaluated'; cartesian_accuracy_failure='not_evaluated'
      kifmm_aggregate_rms=[double]::NaN; cartesian_aggregate_rms=[double]::NaN
      kifmm_max_channel_rms=[double]::NaN; cartesian_max_channel_rms=[double]::NaN
      kifmm_total_seconds=[double]::NaN; cartesian_total_seconds=[double]::NaN
      cartesian_sampled_apply_seconds=[double]::NaN; cartesian_direct_reference_seconds=[double]::NaN
      kifmm_far_panel_pairs=0; kifmm_fallback_panel_pairs=0
      cartesian_sampled_far_panel_pairs=0; cartesian_sampled_fallback_panel_pairs=0
      expected_pairs=0; cartesian_valid_pairs=0; kifmm_valid_pairs=0
      coverage_exact='false'; invalid_pairs=0; nonproduction_fallback_pairs=0
      parameter_match='false'; sampled_target_hash_match='false'
      cartesian_timing_scope='missing'; cartesian_pair_count_scope='missing'
      kifmm_eps_ratio=[double]::NaN; cartesian_eps_ratio=[double]::NaN
      kifmm_kappa=[double]::NaN; cartesian_kappa=[double]::NaN
      kifmm_theta=[double]::NaN; cartesian_theta=[double]::NaN
      kifmm_near_integration=''; cartesian_near_integration=''
      kifmm_quad_buffer=''; cartesian_quad_buffer=''
      kifmm_fallback_kernel=''; cartesian_fallback_kernel=''
      cartesian_fallback_kernel_effective=''
      cartesian_k1=[double]::NaN; cartesian_k2=[double]::NaN; cartesian_k3=[double]::NaN; cartesian_k4=[double]::NaN
      kifmm_k1=[double]::NaN; kifmm_k2=[double]::NaN; kifmm_k3=[double]::NaN; kifmm_k4=[double]::NaN
      log=$log
    })
    continue
  }

  $k = Parse-KeyValueLine ([string]$kifmmLine[0])
  $c = Parse-KeyValueLine ([string]$cartLine[0])
  $sourcePanels = [uint64](U64 $c 'source_panels')
  $sampledRows = [uint64](U64 $c 'sampled_rows')
  $expectedPairs = $sourcePanels * $sampledRows
  $cartValid = [uint64](U64 $c 'valid_pairs')
  $kifmmValid = [uint64](U64 $k 'production_direct_reference_valid_pairs')
  $cartFarPairs = [uint64](U64 $c 'sampled_far_panel_pairs')
  $cartFallbackPairs = [uint64](U64 $c 'sampled_fallback_panel_pairs')
  $kifmmFarPairs = [uint64](U64 $k 'far_panel_pairs')
  $kifmmFallbackPairs = [uint64](U64 $k 'fallback_panel_pairs')
  $invalidPairs = [uint64](U64 $c 'invalid_pairs') + [uint64](U64 $k 'invalid_pairs')
  $nonProduction = [uint64](U64 $c 'nonproduction_fallback_pairs') + [uint64](U64 $k 'nonproduction_fallback_pairs')
  $hashMatch = (U64 $c 'sampled_target_hash') -eq (U64 $k 'sampled_target_hash')
  $coverageExact = ($cartValid -eq $expectedPairs) -and ($kifmmValid -eq $expectedPairs) -and
    (($cartFarPairs + $cartFallbackPairs) -eq $expectedPairs) -and
    (($kifmmFarPairs + $kifmmFallbackPairs) -eq $expectedPairs) -and
    ((S $c 'timing_scope') -eq 'sampled_rows_only') -and ((S $c 'pair_count_scope') -eq 'sampled_rows_only')
  $kifmmAggregate = D $k 'aggregate_rms'
  $cartAggregate = D $c 'aggregate_rms'
  $kifmmMaxChannel = Max-Channel $k
  $cartMaxChannel = Max-Channel $c
  $cartExecution = S $c 'execution_status' (S $c 'status')
  $cartAccuracy = S $c 'accuracy_status'
  if ([string]::IsNullOrWhiteSpace($cartAccuracy)) {
    $cartAccuracy = Accuracy-Status $cartExecution $cartAggregate $cartMaxChannel
    if ($cartAccuracy -ne 'pass') { $cartAccuracy = 'fail' }
  }
  $cartAccuracyFailure = S $c 'accuracy_failure'
  if ([string]::IsNullOrWhiteSpace($cartAccuracyFailure)) {
    $cartAccuracyFailure = Accuracy-Status $cartExecution $cartAggregate $cartMaxChannel
  }
  $epsMatch = Nearly-Equal (D $k 'eps_ratio') (D $c 'eps_ratio')
  $kappaMatch = Nearly-Equal (D $k 'kappa') (D $c 'kappa')
  $thetaMatch = Nearly-Equal (D $k 'theta') (D $c 'theta')
  $nearMatch = ((S $k 'near_integration') -eq 'analytic_singular') -and ((S $c 'near_integration') -eq 'analytic_singular')
  $quadMatch = ((S $k 'quad_buffer') -eq 'off') -and ((S $c 'quad_buffer') -eq 'off')
  $parameterMatch = $epsMatch -and $kappaMatch -and $thetaMatch -and $nearMatch -and $quadMatch

  $matrix.Add([pscustomobject]@{
    row_id=$row.row_id
    fixture=$row.fixture
    sample_set=$row.sample
    sample_count=$sampledRows
    order=6
    kifmm_status=(S $k 'status')
    cartesian_status=(S $c 'status')
    cartesian_execution_status=$cartExecution
    cartesian_accuracy_status=$cartAccuracy
    cartesian_accuracy_failure=$cartAccuracyFailure
    kifmm_aggregate_rms=$kifmmAggregate
    cartesian_aggregate_rms=$cartAggregate
    kifmm_max_channel_rms=$kifmmMaxChannel
    cartesian_max_channel_rms=$cartMaxChannel
    kifmm_total_seconds=(D $k 'total_seconds')
    cartesian_total_seconds=(D $c 'total_seconds')
    cartesian_sampled_apply_seconds=(D $c 'sampled_apply_seconds')
    cartesian_direct_reference_seconds=(D $c 'direct_reference_seconds')
    kifmm_far_panel_pairs=$kifmmFarPairs
    kifmm_fallback_panel_pairs=$kifmmFallbackPairs
    cartesian_sampled_far_panel_pairs=$cartFarPairs
    cartesian_sampled_fallback_panel_pairs=$cartFallbackPairs
    expected_pairs=$expectedPairs
    cartesian_valid_pairs=$cartValid
    kifmm_valid_pairs=$kifmmValid
    coverage_exact=([string]$coverageExact).ToLowerInvariant()
    invalid_pairs=$invalidPairs
    nonproduction_fallback_pairs=$nonProduction
    parameter_match=([string]$parameterMatch).ToLowerInvariant()
    sampled_target_hash_match=([string]$hashMatch).ToLowerInvariant()
    cartesian_timing_scope=(S $c 'timing_scope')
    cartesian_pair_count_scope=(S $c 'pair_count_scope')
    kifmm_eps_ratio=(D $k 'eps_ratio')
    cartesian_eps_ratio=(D $c 'eps_ratio')
    kifmm_kappa=(D $k 'kappa')
    cartesian_kappa=(D $c 'kappa')
    kifmm_theta=(D $k 'theta')
    cartesian_theta=(D $c 'theta')
    kifmm_near_integration=(S $k 'near_integration')
    cartesian_near_integration=(S $c 'near_integration')
    kifmm_quad_buffer=(S $k 'quad_buffer')
    cartesian_quad_buffer=(S $c 'quad_buffer')
    kifmm_fallback_kernel=(S $k 'fallback_kernel')
    cartesian_fallback_kernel=(S $c 'fallback_kernel')
    cartesian_fallback_kernel_effective=(S $c 'fallback_kernel_effective')
    cartesian_k1=(D $c 'k1')
    cartesian_k2=(D $c 'k2')
    cartesian_k3=(D $c 'k3')
    cartesian_k4=(D $c 'k4')
    kifmm_k1=(D $k 'k1')
    kifmm_k2=(D $k 'k2')
    kifmm_k3=(D $k 'k3')
    kifmm_k4=(D $k 'k4')
    log=$log
  })
}

$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,sample_count,expected_pairs,cartesian_valid_pairs,kifmm_valid_pairs,coverage_exact,invalid_pairs,nonproduction_fallback_pairs,sampled_target_hash_match,cartesian_sampled_far_panel_pairs,cartesian_sampled_fallback_panel_pairs,kifmm_far_panel_pairs,kifmm_fallback_panel_pairs |
  Export-Csv -LiteralPath $coveragePath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,kifmm_total_seconds,cartesian_total_seconds,cartesian_sampled_apply_seconds,cartesian_direct_reference_seconds,cartesian_timing_scope |
  Export-Csv -LiteralPath $runtimePath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,kifmm_k1,kifmm_k2,kifmm_k3,kifmm_k4,kifmm_max_channel_rms,cartesian_k1,cartesian_k2,cartesian_k3,cartesian_k4,cartesian_max_channel_rms |
  Export-Csv -LiteralPath $channelPath -NoTypeInformation -Encoding ascii

$kifmmSubgate = 'fmm_e16_kifmm_level_cube_sampled_pass'
if (@($matrix | Where-Object { $_.kifmm_status -ne 'pass' -or [double]$_.kifmm_aggregate_rms -gt 1.0e-6 -or [double]$_.kifmm_max_channel_rms -gt 5.0e-6 }).Count -gt 0) {
  $kifmmSubgate = 'fmm_e16_kifmm_accuracy_gap'
}
$cartesianSubgate = 'fmm_e16_cartesian_sampled_pass'
if (@($matrix | Where-Object { $_.cartesian_execution_status -ne 'pass' }).Count -gt 0) {
  $cartesianSubgate = 'fmm_e16_cartesian_execution_gap'
} elseif (@($matrix | Where-Object { $_.cartesian_accuracy_status -ne 'pass' }).Count -gt 0) {
  $cartesianSubgate = 'fmm_e16_cartesian_sampled_accuracy_gap'
}

$decision = 'fmm_e16_true_sampled_comparator_pass'
if (@($matrix | Where-Object { $_.parameter_match -ne 'true' }).Count -gt 0) {
  $decision = 'fmm_e16_parameter_mismatch'
} elseif (@($matrix | Where-Object { $_.sampled_target_hash_match -ne 'true' -or $_.coverage_exact -ne 'true' -or [uint64]$_.invalid_pairs -ne 0 -or [uint64]$_.nonproduction_fallback_pairs -ne 0 }).Count -gt 0) {
  $decision = 'fmm_e16_pair_coverage_gap'
} elseif ($kifmmSubgate -ne 'fmm_e16_kifmm_level_cube_sampled_pass') {
  $decision = 'fmm_e16_kifmm_accuracy_gap'
} elseif (@($matrix | Where-Object { $_.cartesian_execution_status -ne 'pass' }).Count -gt 0) {
  $decision = 'fmm_e16_sampled_adapter_unavailable'
} elseif ($cartesianSubgate -ne 'fmm_e16_cartesian_sampled_pass') {
  $decision = 'fmm_e16_cartesian_accuracy_gap'
}

[pscustomobject]@{
  decision=$decision
  kifmm_subgate=$kifmmSubgate
  cartesian_subgate=$cartesianSubgate
  mode=($(if ($Quick -and -not $Full) { 'quick' } else { 'full' }))
  rows=$matrix.Count
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "kifmm_subgate=$kifmmSubgate",
  "cartesian_subgate=$cartesianSubgate",
  "rows=$($matrix.Count)",
  "matrix=$matrixPath",
  "coverage=$coveragePath",
  "runtime=$runtimePath",
  "channels=$channelPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$worstK4 = @($matrix | Sort-Object { [double]$_.cartesian_k4 } -Descending | Select-Object -First 1)
$worstK4Line = if ($worstK4.Count -gt 0) {
  "Cartesian K4 maximum: row=$($worstK4[0].row_id), K4_RMS=$($worstK4[0].cartesian_k4)"
} else {
  'Cartesian K4 maximum: unavailable'
}

$reportLines = @(
  'TITAN FMM Engine E16 True Sampled Cartesian Comparator Report',
  '',
  "Decision: $decision",
  "KIFMM subgate: $kifmmSubgate",
  "Cartesian comparator: $cartesianSubgate",
  "Mode: $(if ($Quick -and -not $Full) { 'quick' } else { 'full' })",
  "Rows: $($matrix.Count)",
  '',
  'E16 compares KIFMM-lite and guarded Cartesian p4 on identical sampled target rows, full source panels, matched physics parameters, and the same production dense-direct row oracle.',
  'Cartesian timing and pair counts in this report are sampled-row-only. No production TABI matvec, GMRES path, preconditioner, or default policy is changed.',
  'The Cartesian comparator gap is primarily channel-level, especially in K4, rather than a missing-pair or parameter-mismatch issue.',
  $worstK4Line,
  '',
  'Row summary:',
  '  row_id,kifmm_rms,kifmm_max_channel,cartesian_rms,cartesian_max_channel,cartesian_execution_status,cartesian_accuracy_status,cartesian_accuracy_failure,kifmm_seconds,cartesian_seconds,coverage_exact,parameter_match'
)
foreach ($m in $matrix) {
  $reportLines += ("  {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}" -f
    $m.row_id,
    $m.kifmm_aggregate_rms,
    $m.kifmm_max_channel_rms,
    $m.cartesian_aggregate_rms,
    $m.cartesian_max_channel_rms,
    $m.cartesian_execution_status,
    $m.cartesian_accuracy_status,
    $m.cartesian_accuracy_failure,
    $m.kifmm_total_seconds,
    $m.cartesian_total_seconds,
    $m.coverage_exact,
    $m.parameter_match)
}
$reportLines += @(
  '',
  'Artifact paths:',
  "  Matrix: $matrixPath",
  "  Coverage: $coveragePath",
  "  Runtime: $runtimePath",
  "  Channel summary: $channelPath",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$reportLines | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E16 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -ne 'fmm_e16_true_sampled_comparator_pass' -and -not $AllowFailedDecision) {
  throw "E16 did not pass: $decision"
}
