param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [switch]$Include1A63,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e30_shared_worklist_channel_hybrid'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e30_matrix.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e30_channel.csv'
$ProvenanceCsv = Join-Path $WorkRoot 'fmm_e30_provenance.csv'
$CoverageCsv = Join-Path $WorkRoot 'fmm_e30_coverage.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e30_runtime.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e30_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e30_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E30_SHARED_WORKLIST_CHANNEL_HYBRID_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E30 smoke.'
}

function Assert-SmokeSupportsE30 {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @(
      'channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian',
      'channel-hybrid-k123-hierarchy-k4-cartesian',
      'k1_source',
      'k2_source',
      'k3_source',
      'k4_source',
      'cartesian_k4_status',
      'cartesian_k4_accuracy_status',
      'channel_hybrid_double_count_guard',
      'shared_worklist_hash',
      'worklist_consistency_status',
      'direct_avx2_all_channel_pairs')) {
    if (-not $text.Contains($needle)) {
      throw "E30 shared-worklist channel-hybrid support is not present yet: missing '$needle' in $SmokeSource"
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
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  [double]$value
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
  Assert-SmokeSupportsE30
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

function New-Case([string]$Id, [string]$Fixture, [string]$Sample, [int]$SampleCount, [string]$Mode) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    mode = $Mode
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'mode' $Case.mode
  Ensure-Field $Row 'hier_debug_mode' $Case.mode
  foreach ($field in @(
      'k1_source','k2_source','k3_source','k4_source',
      'cartesian_k4_status','cartesian_k4_execution_status',
      'cartesian_k4_accuracy_status','cartesian_k4_accuracy_failure',
      'cartesian_k4_rms','cartesian_k4_abs_rms','cartesian_k4_reference_rms',
      'cartesian_k4_far_pairs','cartesian_k4_fallback_pairs',
      'cartesian_k4_total_seconds','cartesian_k4_policy_hash',
      'channel_hybrid_double_count_guard','channel_hybrid_recombined_hash',
      'hierarchy_k123_max_channel_rms','hierarchy_k123_max_channel_name',
      'hybrid_recombined_aggregate_rms','hybrid_recombined_max_channel_rms',
      'hybrid_recombined_max_channel_name',
      'shared_worklist_hash','shared_theta','shared_guard',
      'shared_far_pairs','shared_direct_pairs','shared_far_fraction',
      'cartesian_k4_shared_far_pairs','cartesian_k4_shared_direct_pairs',
      'direct_avx2_all_channel_pairs','worklist_consistency_status')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $Row
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.id)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--order=8",
    "--theta=0.5",
    "--equiv-scale=1.05",
    "--upward-check-scale=3.0",
    "--local-equiv-scale=3.0",
    "--local-check-scale=1.05",
    "--condition-policy=qr-rank-only",
    "--real-panel-guard=strict-local-check",
    "--guard-margin-factor=0",
    "--hier-debug-mode=$($Case.mode)",
    "--mode=profile",
    "--pqr-1sh1=$Pqr1sh1",
    "--pqr-1a63=$Pqr1a63"
  )
  $stdoutPath = Join-Path $WorkRoot "$($Case.id).stdout.txt"
  $stderrPath = Join-Path $WorkRoot "$($Case.id).stderr.txt"
  Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue
  $proc = Start-Process -FilePath $ExePath -ArgumentList $args -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
  $output = @()
  if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath }
  if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath }
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($resultLine.Count -ne 1) {
    $failedRow = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      mode = $Case.mode
      hier_debug_mode = $Case.mode
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'process_failed'
      process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }
    return (Complete-Row $failedRow $Case)
  }
  $row = Complete-Row (Parse-KeyValueLine ([string]$resultLine[-1]) 'E255_RESULT') $Case
  Ensure-Field $row 'process_exit_code' $proc.ExitCode
  $row
}

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) {
    Set-Content -LiteralPath $Path -Value ''
  } else {
    $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path
  }
}

function Hybrid-Pass([array]$Rows) {
  $selected = @($Rows | Where-Object { $_.mode -eq 'channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian' })
  if ($selected.Count -eq 0) { return $false }
  foreach ($row in $selected) {
    if ((Get-Field $row 'execution_status') -ne 'pass') { return $false }
    if ((Get-Field $row 'accuracy_status') -ne 'pass') { return $false }
    if ((As-Double $row 'aggregate_rms') -gt 1.0e-6) { return $false }
    if ((As-Double $row 'max_channel_rms') -gt 5.0e-6) { return $false }
    if ((As-Double $row 'k4_rms') -gt 5.0e-6) { return $false }
    if ((As-Double $row 'k1_rms') -gt 5.0e-6) { return $false }
    if ((As-Double $row 'k2_rms') -gt 5.0e-6) { return $false }
    if ((As-Double $row 'k3_rms') -gt 5.0e-6) { return $false }
    if ((Get-Field $row 'k1_source') -ne 'hierarchy_q8_shared_worklist') { return $false }
    if ((Get-Field $row 'k2_source') -ne 'hierarchy_q8_shared_worklist') { return $false }
    if ((Get-Field $row 'k3_source') -ne 'hierarchy_q8_shared_worklist') { return $false }
    if ((Get-Field $row 'k4_source') -ne 'cartesian_p5_shared_worklist') { return $false }
    if ((Get-Field $row 'channel_hybrid_double_count_guard') -ne 'pass') { return $false }
    if ((Get-Field $row 'worklist_consistency_status') -ne 'pass') { return $false }
    if ((As-Int64 $row 'production_invalid_pairs' 1) -ne 0) { return $false }
    if ((As-Int64 $row 'nonproduction_fallback_pairs' 1) -ne 0) { return $false }
    if ((Get-Field $row 'coverage') -ne (Get-Field $row 'expected_coverage')) { return $false }
    if ((As-Double $row 'shared_far_fraction' 0.0) -lt 0.10) { return $false }
    if ((As-Int64 $row 'cartesian_k4_shared_far_pairs' 0) -ne (As-Int64 $row 'shared_far_pairs' -1)) { return $false }
    if ((As-Int64 $row 'cartesian_k4_shared_direct_pairs' 0) -ne (As-Int64 $row 'shared_direct_pairs' -1)) { return $false }
    if ((As-Int64 $row 'direct_avx2_all_channel_pairs' 0) -ne (As-Int64 $row 'shared_direct_pairs' -1)) { return $false }
    foreach ($hashField in @('shared_worklist_hash','worklist_hash','policy_hash','panel_geometry_hash','unknown_vector_hash','sampled_target_hash','channel_hybrid_recombined_hash')) {
      if ((As-Int64 $row $hashField 0) -eq 0) { return $false }
    }
  }
  return $true
}

function Decide([array]$Rows, [int]$RequiredRows) {
  if ($Rows.Count -lt $RequiredRows) { return 'fmm_e30_runtime_limited' }
  $hybrid = @($Rows | Where-Object { $_.mode -eq 'channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian' })
  if (@($hybrid | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' -or (Get-Field $_ 'coverage') -ne (Get-Field $_ 'expected_coverage') -or (As-Int64 $_ 'production_invalid_pairs' 1) -ne 0 -or (As-Int64 $_ 'nonproduction_fallback_pairs' 1) -ne 0 -or (Get-Field $_ 'channel_hybrid_double_count_guard') -ne 'pass' -or (Get-Field $_ 'worklist_consistency_status') -ne 'pass' }).Count -gt 0) {
    return 'fmm_e30_worklist_consistency_gap'
  }
  if (@($hybrid | Where-Object { (As-Double $_ 'shared_far_fraction' 0.0) -lt 0.10 }).Count -gt 0) {
    return 'fmm_e30_shared_worklist_too_conservative'
  }
  if (@($hybrid | Where-Object { (Get-Field $_ 'cartesian_k4_accuracy_status') -ne 'pass' -or (As-Double $_ 'cartesian_k4_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e30_k4_cartesian_shared_worklist_gap'
  }
  if (@($hybrid | Where-Object { (As-Double $_ 'k1_rms') -gt 5.0e-6 -or (As-Double $_ 'k2_rms') -gt 5.0e-6 -or (As-Double $_ 'k3_rms') -gt 5.0e-6 }).Count -gt 0) {
    return 'fmm_e30_k123_hierarchy_gap'
  }
  if (Hybrid-Pass $Rows) { return 'fmm_e30_shared_worklist_channel_hybrid_pass' }
  'fmm_e30_channel_recombination_gap'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$modes = @('full', 'channel-hybrid-k123-hierarchy-k4-cartesian', 'channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian')
$sampleSpecs = @(
  @{ fixture = '1SH1'; sample = 'uniform_hash'; count = 128; prefix = '1SH1_uniform_128' },
  @{ fixture = '1SH1'; sample = 'closest_opposing_surface'; count = 128; prefix = '1SH1_closest_128' }
)
if (-not $Quick) {
  $sampleSpecs += @(
    @{ fixture = '1SH1'; sample = 'high_aspect'; count = 128; prefix = '1SH1_high_aspect_128' },
    @{ fixture = '1SH1'; sample = 'inward_normal_groove'; count = 128; prefix = '1SH1_groove_128' }
  )
}
if ($Include1A63) {
  $sampleSpecs += @(
    @{ fixture = '1A63'; sample = 'uniform_hash'; count = 128; prefix = '1A63_uniform_128' },
    @{ fixture = '1A63'; sample = 'closest_opposing_surface'; count = 128; prefix = '1A63_closest_128' },
    @{ fixture = '1A63'; sample = 'high_aspect'; count = 128; prefix = '1A63_high_aspect_128' }
  )
}
$cases = @()
foreach ($sample in $sampleSpecs) {
  foreach ($mode in $modes) {
    $cases += New-Case "$($sample.prefix)_$($mode.Replace('-', '_'))" $sample.fixture $sample.sample $sample.count $mode
  }
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
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

Export-Table ($matrix | Select-Object row_id,fixture,sample,mode,status,execution_status,accuracy_status,accuracy_failure,aggregate_rms,max_channel_rms,max_channel_name,k1_rms,k2_rms,k3_rms,k4_rms,k4_abs_rms,k4_reference_rms) $ChannelCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,mode,k1_source,k2_source,k3_source,k4_source,cartesian_k4_status,cartesian_k4_accuracy_status,cartesian_k4_accuracy_failure,channel_hybrid_double_count_guard,worklist_consistency_status,channel_hybrid_recombined_hash,shared_worklist_hash,shared_theta,shared_guard) $ProvenanceCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,mode,source_panels,sampled_rows,coverage,expected_coverage,production_invalid_pairs,nonproduction_fallback_pairs,far_panel_pairs,direct_panel_pairs,shared_far_pairs,shared_direct_pairs,shared_far_fraction,cartesian_k4_far_pairs,cartesian_k4_fallback_pairs,cartesian_k4_shared_far_pairs,cartesian_k4_shared_direct_pairs,direct_avx2_all_channel_pairs,worklist_hash,shared_worklist_hash,policy_hash,panel_geometry_hash,unknown_vector_hash,sampled_target_hash) $CoverageCsv
Export-Table ($matrix | Select-Object row_id,fixture,sample,mode,hierarchy_seconds,production_fallback_seconds,dense_oracle_seconds,total_seconds,cartesian_k4_total_seconds) $RuntimeCsv

$decision = Decide $matrix $cases.Count
[pscustomobject][ordered]@{
  decision = $decision
  mode = $(if ($Quick) { 'quick' } else { 'full' })
  rows = $matrix.Count
  required_rows = $cases.Count
  hybrid_rows = @($matrix | Where-Object { $_.mode -eq 'channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian' }).Count
  hybrid_pass = (Hybrid-Pass $matrix)
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @(
  "decision=$decision",
  "rows=$($matrix.Count)",
  "matrix=$MatrixCsv"
)

$report = @(
  '# TITAN FMM Engine E30 Shared-Worklist Channel-Hybrid Real-Panel Sampled Report',
  '',
  "Decision: $decision",
  '',
  'E30 tests a sampled-row diagnostic shared-worklist channel hybrid only. It does not promote KIFMM, Cartesian p5, GMRES, production TABI, or any default.',
  'A single strict-local-check worklist at theta=0.5 decides far/direct coverage. K1/K2/K3 are sourced from hierarchy q8 on that shared worklist; K4 is sourced from Cartesian p5 on the same far worklist; direct pairs are counted as all-channel AVX2/direct fallback.',
  'By default this gate uses 1SH1 sampled rows only; 1A63 rows are optional with -Include1A63 because E30 is an architecture gate, not a protein-breadth gate.'
)
if ($decision -eq 'fmm_e30_k123_hierarchy_gap') {
  $report += 'Interpretation: shared-worklist K4 is not the blocker; the remaining sampled-row gap is in hierarchy K1/K2/K3.'
} elseif ($decision -eq 'fmm_e30_k4_cartesian_shared_worklist_gap') {
  $report += 'Interpretation: the shared worklist has exact coverage and enough far work, but Cartesian p5 K4 fails when forced onto the same far/direct partition. E29 K4 closure depended on the independent guarded Cartesian worklist, not merely on the p5 channel formula.'
} elseif ($decision -eq 'fmm_e30_shared_worklist_too_conservative') {
  $report += 'Interpretation: accuracy is not the primary issue, but the shared worklist demotes too much work to direct fallback for a useful KIFMM/Cartesian hybrid.'
} elseif ($decision -eq 'fmm_e30_worklist_consistency_gap') {
  $report += 'Interpretation: the shared-worklist accounting is inconsistent; inspect coverage and provenance before interpreting channel accuracy.'
} elseif ($decision -eq 'fmm_e30_shared_worklist_channel_hybrid_pass') {
  $report += 'Interpretation: shared-worklist channel-hybrid composition closes the sampled K1-K4 gate.'
} else {
  $report += 'Interpretation: inspect the channel, provenance, and coverage CSVs before broadening the shared-worklist channel-hybrid gate.'
}
foreach ($row in $matrix) {
  $report += ('  {0}: mode={1}, exec={2}, acc={3}, agg={4}, max={5}/{6}, K1={7}, K2={8}, K3={9}, K4={10}, src=({11},{12},{13},{14}), consistency={15}, shared=({16},{17},{18}), cartK4={19}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'mode'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'aggregate_rms'), `
      (Get-Field $row 'max_channel_name'), (Get-Field $row 'max_channel_rms'), `
      (Get-Field $row 'k1_rms'), (Get-Field $row 'k2_rms'), (Get-Field $row 'k3_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'k1_source'), (Get-Field $row 'k2_source'), (Get-Field $row 'k3_source'), (Get-Field $row 'k4_source'), `
      (Get-Field $row 'worklist_consistency_status'), (Get-Field $row 'shared_far_pairs'), (Get-Field $row 'shared_direct_pairs'), (Get-Field $row 'shared_far_fraction'), (Get-Field $row 'cartesian_k4_rms'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E30 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Channel: $ChannelCsv"
Write-Host "Provenance: $ProvenanceCsv"
Write-Host "Coverage: $CoverageCsv"
Write-Host "Runtime: $RuntimeCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e30_shared_worklist_channel_hybrid_pass')) {
  throw "FMM-E30 decision failed: $decision"
}
