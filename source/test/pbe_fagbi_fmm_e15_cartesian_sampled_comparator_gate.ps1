param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e15_cartesian_sampled_comparator'),
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
  throw 'No C++ compiler found for KIFMM E15 gate.'
}

function Build-SmokeExe {
  param([string]$Exe)
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E15 Cartesian sampled comparator gate driver' }
  } finally {
    Pop-Location
  }
}

function Parse-KeyValueLine([string]$Line) {
  $out = [ordered]@{}
  foreach ($part in ($Line -split '\s+')) {
    $idx = $part.IndexOf('=')
    if ($idx -le 0) { continue }
    $key = $part.Substring(0, $idx)
    $value = $part.Substring($idx + 1)
    $out[$key] = $value
  }
  return $out
}

function D([hashtable]$Map, [string]$Key, [double]$Default = [double]::NaN) {
  if (-not $Map.Contains($Key)) { return $Default }
  $v = 0.0
  if ([double]::TryParse($Map[$Key], [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$v)) {
    return $v
  }
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
  [pscustomobject]@{ row_id='1A63_closest_opposing_surface_128_q6'; fixture='1A63'; sample='closest_opposing_surface'; sample_count=128 },
  [pscustomobject]@{ row_id='1A63_inward_normal_groove_128_q6'; fixture='1A63'; sample='inward_normal_groove'; sample_count=128 }
)
$rows = if ($Quick -and -not $Full) { $quickRows } else { $fullRows }

$matrixPath = Join-Path $WorkRoot 'fmm_e15_matrix.csv'
$coveragePath = Join-Path $WorkRoot 'fmm_e15_pair_coverage.csv'
$runtimePath = Join-Path $WorkRoot 'fmm_e15_runtime.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e15_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e15_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E15_CARTESIAN_SAMPLED_COMPARATOR_REPORT.txt'
$exe = Join-Path $WorkRoot 'titan_fmm_e15_cartesian_sampled_comparator_gate.exe'

Build-SmokeExe -Exe $exe

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($row in @(Import-Csv -LiteralPath $matrixPath)) {
    $existing[$row.row_id] = $row
  }
}

$matrix = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
  if ($existing.ContainsKey($row.row_id)) {
    $old = $existing[$row.row_id]
    $rerun = $Force -or ($RerunRows -contains $row.row_id) -or
      (-not ($old.PSObject.Properties.Name -contains 'parameter_match')) -or
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
    '--e15-profile',
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
      kifmm_status='run_failed'; cartesian_status='run_failed'
      cartesian_execution_status='run_failed'
      cartesian_accuracy_status='not_evaluated'
      cartesian_accuracy_failure='not_evaluated'
      kifmm_aggregate_rms=[double]::NaN; cartesian_aggregate_rms=[double]::NaN
      kifmm_max_channel_rms=[double]::NaN; cartesian_max_channel_rms=[double]::NaN
      kifmm_total_seconds=[double]::NaN; cartesian_full_apply_context_seconds=[double]::NaN
      kifmm_far_panel_pairs=0; kifmm_fallback_panel_pairs=0
      cartesian_full_apply_context_far_panel_pairs=0
      cartesian_full_apply_context_fallback_panel_pairs=0
      coverage_exact='false'; invalid_pairs=0; nonproduction_fallback_pairs=0
      parameter_match='false'
      kifmm_eps_ratio=[double]::NaN; cartesian_eps_ratio=[double]::NaN
      kifmm_kappa=[double]::NaN; cartesian_kappa=[double]::NaN
      kifmm_theta=[double]::NaN; cartesian_theta=[double]::NaN
      kifmm_near_integration=''; cartesian_near_integration=''
      kifmm_quad_buffer=''; cartesian_quad_buffer=''
      kifmm_fallback_kernel=''; cartesian_fallback_kernel=''
      sampled_target_hash_match='false'; policy_hash_kifmm=0; policy_hash_cartesian=0
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
  $invalidPairs = [uint64](U64 $c 'invalid_pairs') + [uint64](U64 $k 'invalid_pairs')
  $nonProduction = [uint64](U64 $c 'nonproduction_fallback_pairs') + [uint64](U64 $k 'nonproduction_fallback_pairs')
  $hashMatch = (U64 $c 'sampled_target_hash') -eq (U64 $k 'sampled_target_hash')
  $coverageExact = ($cartValid -eq $expectedPairs) -and ($kifmmValid -eq $expectedPairs)
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
  $nearMatch = ((S $k 'near_integration') -eq 'analytic_singular') -and
    ((S $c 'near_integration') -eq 'analytic_singular')
  $quadMatch = ((S $k 'quad_buffer') -eq 'off') -and
    ((S $c 'quad_buffer') -eq 'off')
  $parameterMatch = $epsMatch -and $kappaMatch -and $thetaMatch -and
    $nearMatch -and $quadMatch
  $cartSeconds = D $c 'full_apply_context_seconds'
  if ([double]::IsNaN($cartSeconds) -or [double]::IsInfinity($cartSeconds)) { $cartSeconds = D $c 'total_seconds' }
  $cartFarPairs = U64 $c 'full_apply_context_far_panel_pairs'
  if ($cartFarPairs -eq 0) { $cartFarPairs = U64 $c 'sampled_far_panel_pairs' }
  $cartFallbackPairs = U64 $c 'full_apply_context_fallback_panel_pairs'
  if ($cartFallbackPairs -eq 0) { $cartFallbackPairs = U64 $c 'sampled_fallback_panel_pairs' }
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
    cartesian_full_apply_context_seconds=$cartSeconds
    kifmm_far_panel_pairs=(U64 $k 'far_panel_pairs')
    kifmm_fallback_panel_pairs=(U64 $k 'fallback_panel_pairs')
    cartesian_full_apply_context_far_panel_pairs=$cartFarPairs
    cartesian_full_apply_context_fallback_panel_pairs=$cartFallbackPairs
    coverage_exact=([string]$coverageExact).ToLowerInvariant()
    invalid_pairs=$invalidPairs
    nonproduction_fallback_pairs=$nonProduction
    parameter_match=([string]$parameterMatch).ToLowerInvariant()
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
    sampled_target_hash_match=([string]$hashMatch).ToLowerInvariant()
    policy_hash_kifmm=(U64 $k 'policy_hash')
    policy_hash_cartesian=(U64 $c 'policy_hash')
    log=$log
  })
}

$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,sample_count,coverage_exact,invalid_pairs,nonproduction_fallback_pairs,parameter_match,sampled_target_hash_match,kifmm_far_panel_pairs,kifmm_fallback_panel_pairs,cartesian_full_apply_context_far_panel_pairs,cartesian_full_apply_context_fallback_panel_pairs |
  Export-Csv -LiteralPath $coveragePath -NoTypeInformation -Encoding ascii
$matrix | Select-Object row_id,fixture,sample_set,kifmm_total_seconds,cartesian_full_apply_context_seconds |
  Export-Csv -LiteralPath $runtimePath -NoTypeInformation -Encoding ascii

$decision = 'fmm_e15_cartesian_sampled_comparator_pass'
if (@($matrix | Where-Object { $_.parameter_match -ne 'true' }).Count -gt 0) {
  $decision = 'fmm_e15_parameter_mismatch'
} elseif (@($matrix | Where-Object { $_.sampled_target_hash_match -ne 'true' }).Count -gt 0) {
  $decision = 'fmm_e15_sample_hash_gap'
} elseif (@($matrix | Where-Object { $_.coverage_exact -ne 'true' -or [uint64]$_.invalid_pairs -ne 0 -or [uint64]$_.nonproduction_fallback_pairs -ne 0 }).Count -gt 0) {
  $decision = 'fmm_e15_pair_coverage_gap'
} elseif (@($matrix | Where-Object { $_.kifmm_status -ne 'pass' -or [double]$_.kifmm_aggregate_rms -gt 1.0e-6 -or [double]$_.kifmm_max_channel_rms -gt 5.0e-6 }).Count -gt 0) {
  $decision = 'fmm_e15_kifmm_accuracy_gap'
} elseif (@($matrix | Where-Object { $_.cartesian_execution_status -ne 'pass' }).Count -gt 0) {
  $decision = 'fmm_e15_cartesian_unavailable'
} elseif (@($matrix | Where-Object { $_.cartesian_accuracy_status -ne 'pass' }).Count -gt 0) {
  $decision = 'fmm_e15_cartesian_accuracy_gap'
}

[pscustomobject]@{
  decision=$decision
  mode=($(if ($Quick -and -not $Full) { 'quick' } else { 'full' }))
  rows=$matrix.Count
  timestamp=(Get-Date).ToString('s')
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$matrixPath",
  "coverage=$coveragePath",
  "runtime=$runtimePath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$reportLines = @(
  'TITAN FMM Engine E15 Cartesian Sampled Comparator Report',
  '',
  "Decision: $decision",
  "Mode: $(if ($Quick -and -not $Full) { 'quick' } else { 'full' })",
  "Rows: $($matrix.Count)",
  '',
  'E15 compares sampled rows only. It does not wire either KIFMM-lite or the Cartesian comparator into GMRES, TABI production matvec, or preconditioners. Dense direct production rows remain the reference for each sampled target set.',
  '',
  'Interpretation:',
  '  A cartesian_accuracy_gap decision means KIFMM-lite passed the sampled dense-direct gates, while the guarded Cartesian p4 comparator exceeded the strict sampled-row threshold on at least one identical target set. This is a diagnostic comparator result, not a production-policy change.',
  '  In E16-capable builds the Cartesian comparator reports sampled-row timing and pair counts; older full-apply context columns are retained here only for CSV compatibility.',
  '',
  'Row summary:',
  '  row_id,kifmm_aggregate_rms,kifmm_max_channel_rms,cartesian_aggregate_rms,cartesian_max_channel_rms,cartesian_execution_status,cartesian_accuracy_status,cartesian_accuracy_failure,kifmm_seconds,cartesian_full_apply_context_seconds,parameter_match,coverage_exact,hash_match'
)
foreach ($m in $matrix) {
  $reportLines += ("  {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}" -f
    $m.row_id,
    $m.kifmm_aggregate_rms,
    $m.kifmm_max_channel_rms,
    $m.cartesian_aggregate_rms,
    $m.cartesian_max_channel_rms,
    $m.cartesian_execution_status,
    $m.cartesian_accuracy_status,
    $m.cartesian_accuracy_failure,
    $m.kifmm_total_seconds,
    $m.cartesian_full_apply_context_seconds,
    $m.parameter_match,
    $m.coverage_exact,
    $m.sampled_target_hash_match)
}
$reportLines += @(
  '',
  'Artifact paths:',
  "  Matrix: $matrixPath",
  "  Pair coverage: $coveragePath",
  "  Runtime: $runtimePath",
  "  Decision: $decisionPath",
  "  Summary: $summaryPath"
)
$reportLines | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "E15 decision: $decision"
Write-Host "Matrix: $matrixPath"
if ($decision -ne 'fmm_e15_cartesian_sampled_comparator_pass' -and -not $AllowFailedDecision) {
  throw "E15 did not pass: $decision"
}
