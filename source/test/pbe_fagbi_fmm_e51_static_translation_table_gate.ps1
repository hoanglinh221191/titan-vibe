param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$Force,
  [int]$Workers = 4,
  [int]$CaseTimeoutSeconds = 900,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e51_static_translation_table'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_hierarchy_real_panel_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e51_matrix.csv'
$CacheCsv = Join-Path $WorkRoot 'fmm_e51_cache.csv'
$ParityCsv = Join-Path $WorkRoot 'fmm_e51_parity.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e51_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e51_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E51_STATIC_TRANSLATION_TABLE_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E51 smoke.'
}

function Assert-Contains([string]$Path, [string[]]$Needles) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Missing file: $Path" }
  $text = Get-Content -Raw -LiteralPath $Path
  foreach ($needle in $Needles) {
    if (-not $text.Contains($needle)) { throw "Missing marker '$needle' in $Path" }
  }
}

function Assert-E51Support {
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_jcp.cpp') @(
    'static-translation-table-cache',
    'static_translation_table_pass',
    'm2l_matrix_cache_key_index'
  )
  Assert-Contains (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy.h') @(
    'static-translation-table-cache'
  )
  Assert-Contains $SmokeSource @(
    'static-translation-table-cache',
    'm2l_matrix_cache_key_status='
  )
}

function Build-Smoke {
  Assert-E51Support
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
    'src\tabi\tabi_fmm_kifmm_real_panel_geometry.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_worklist.cpp',
    'src\tabi\tabi_fmm_kifmm_channel_hybrid_operator.cpp',
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
  if ($null -eq $prop -or $null -eq $prop.Value -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) { return $Default }
  [string]$prop.Value
}

function As-Double([object]$Row, [string]$Name, [double]$Default = [double]::NaN) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value) -or $value -eq 'nan') { return $Default }
  try { [double]$value } catch { $Default }
}

function As-Int64([object]$Row, [string]$Name, [int64]$Default = 0) {
  $value = Get-Field $Row $Name ''
  if ([string]::IsNullOrWhiteSpace($value)) { return $Default }
  try { [int64]$value } catch { $Default }
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.row_id
  Ensure-Field $Row 'case_id' $Case.id
  Ensure-Field $Row 'requested_m2l_matrix_cache_policy' $Case.policy
  foreach ($field in @(
      'status','execution_status','accuracy_status','aggregate_rms',
      'max_channel_rms','k1_rms','k2_rms','k3_rms','k4_rms',
      'coverage','full_coverage_expected','production_invalid_pairs',
      'nonproduction_fallback_pairs','far_fraction','m2l_seconds',
      'm2l_matrix_cache_policy','m2l_matrix_cache_entries',
      'm2l_matrix_cache_peak_entries','m2l_matrix_cache_hits',
      'm2l_matrix_cache_misses','m2l_matrix_cache_hit_rate',
      'm2l_matrix_cache_estimated_bytes',
      'm2l_matrix_cache_peak_estimated_bytes',
      'm2l_matrix_cache_total_built_bytes',
      'm2l_matrix_cache_build_seconds','m2l_matrix_gemm_seconds',
      'm2l_matrix_cache_key_status','hierarchy_seconds','total_seconds',
      'output_hash','policy_hash','worklist_hash')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $Row
}

function New-Case([string]$Policy) {
  [pscustomobject]@{
    id = '1SH1_s4_p512_uniform64'
    row_id = ('1SH1_s4_p512_uniform64_' + $Policy).Replace('-', '_')
    policy = $Policy
  }
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.row_id) (timeout ${CaseTimeoutSeconds}s)"
  $args = @(
    '--fixture=1SH1', '--sample=uniform_hash', '--sample-count=64',
    '--surface-points=4', '--panel-limit=512',
    '--mode=full-panel-profile', '--emit-e36-stage-attribution=1',
    '--hier-debug-mode=dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian',
    '--jcp-channel-policy=k123-only',
    '--m2l-eval-policy=target-check-batched',
    "--m2l-matrix-cache-policy=$($Case.policy)",
    '--m2l-gemm-batch-policy=per-pair',
    '--m2l-parallel-policy=target-node-key-group',
    "--m2l-parallel-workers=$Workers",
    '--tree-stop-policy=panel-radius-occupancy',
    '--tree-panel-radius-factor=1.5', '--tree-occupancy-stop=32',
    '--order=8', '--theta=0.6', '--equiv-scale=1.05',
    '--upward-check-scale=3.0', '--local-equiv-scale=3.0',
    '--local-check-scale=1.05', '--condition-policy=qr-rank-only',
    '--real-panel-guard=strict-local-check', '--guard-margin-factor=0',
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
    return Complete-Row ([pscustomobject][ordered]@{ row_id=$Case.row_id; case_id=$Case.id; status='timeout'; execution_status='timeout'; stdout=$stdout; stderr=$stderr }) $Case
  }
  $output = @()
  if (Test-Path -LiteralPath $stdout) { $output += Get-Content -LiteralPath $stdout }
  if (Test-Path -LiteralPath $stderr) { $output += Get-Content -LiteralPath $stderr }
  $line = @($output | Where-Object { $_.StartsWith('E255_RESULT ') } | Select-Object -Last 1)
  if ([string]::IsNullOrWhiteSpace($line)) {
    return Complete-Row ([pscustomobject][ordered]@{ row_id=$Case.row_id; case_id=$Case.id; status='process_failed'; execution_status='process_failed'; stdout=$stdout; stderr=$stderr; process_exit_code=$proc.ExitCode; process_output=(($output | Select-Object -First 20) -join '|') }) $Case
  }
  $parsed = Parse-KeyValueLine $line 'E255_RESULT '
  $parsed | Add-Member -NotePropertyName stdout -NotePropertyValue $stdout
  $parsed | Add-Member -NotePropertyName stderr -NotePropertyValue $stderr
  Complete-Row $parsed $Case
}

function Is-Row-Pass([object]$Row) {
  (Get-Field $Row 'execution_status') -eq 'pass' -and
  (Get-Field $Row 'accuracy_status') -eq 'pass' -and
  (As-Double $Row 'aggregate_rms' 1.0) -le 1.0e-6 -and
  (As-Double $Row 'max_channel_rms' 1.0) -le 5.0e-6 -and
  (As-Int64 $Row 'coverage' -1) -eq (As-Int64 $Row 'full_coverage_expected' -2) -and
  (As-Int64 $Row 'production_invalid_pairs' 1) -eq 0 -and
  (As-Int64 $Row 'nonproduction_fallback_pairs' 1) -eq 0
}

function Decide([array]$Rows) {
  if ($Rows.Count -lt 2) { return 'fmm_e51_runtime_limited' }
  if (@($Rows | Where-Object { -not (Is-Row-Pass $_) }).Count -gt 0) { return 'fmm_e51_accuracy_gap' }
  $canonical = @($Rows | Where-Object { $_.requested_m2l_matrix_cache_policy -eq 'canonical-offset-cache' } | Select-Object -First 1)
  $static = @($Rows | Where-Object { $_.requested_m2l_matrix_cache_policy -eq 'static-translation-table-cache' } | Select-Object -First 1)
  if ($canonical.Count -ne 1 -or $static.Count -ne 1) { return 'fmm_e51_runtime_limited' }
  $delta = 0.0
  foreach ($field in @('aggregate_rms','k1_rms','k2_rms','k3_rms','k4_rms','max_channel_rms')) {
    $delta = [Math]::Max($delta, [Math]::Abs((As-Double $canonical[0] $field 0.0) - (As-Double $static[0] $field 0.0)))
  }
  if ($delta -gt 1.0e-12) { return 'fmm_e51_static_table_parity_gap' }
  if ((Get-Field $static[0] 'm2l_matrix_cache_key_status') -notmatch 'static_translation_table') { return 'fmm_e51_static_table_key_gap' }
  $canonicalEntries = As-Int64 $canonical[0] 'm2l_matrix_cache_peak_entries' 0
  $staticEntries = As-Int64 $static[0] 'm2l_matrix_cache_peak_entries' ([int64]::MaxValue)
  if ($canonicalEntries -gt 0 -and $staticEntries -gt ($canonicalEntries * 2)) { return 'fmm_e51_static_table_key_gap' }
  return 'fmm_e51_static_translation_table_pass'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
if ($Force) {
  Remove-Item -LiteralPath $MatrixCsv, $CacheCsv, $ParityCsv, $DecisionCsv, $SummaryDat, $ReportPath -Force -ErrorAction SilentlyContinue
}
Build-Smoke
$cases = @(
  (New-Case 'canonical-offset-cache'),
  (New-Case 'static-translation-table-cache')
)
$rows = foreach ($case in $cases) { Run-Case $case }
$rows | Export-Csv -NoTypeInformation -LiteralPath $MatrixCsv
$rows | Select-Object row_id,requested_m2l_matrix_cache_policy,m2l_matrix_cache_policy,m2l_matrix_cache_entries,m2l_matrix_cache_peak_entries,m2l_matrix_cache_hits,m2l_matrix_cache_misses,m2l_matrix_cache_hit_rate,m2l_matrix_cache_estimated_bytes,m2l_matrix_cache_peak_estimated_bytes,m2l_matrix_cache_total_built_bytes,m2l_matrix_cache_build_seconds,m2l_matrix_gemm_seconds,m2l_matrix_cache_key_status | Export-Csv -NoTypeInformation -LiteralPath $CacheCsv
$canonical = @($rows | Where-Object { $_.requested_m2l_matrix_cache_policy -eq 'canonical-offset-cache' } | Select-Object -First 1)
$static = @($rows | Where-Object { $_.requested_m2l_matrix_cache_policy -eq 'static-translation-table-cache' } | Select-Object -First 1)
$delta = 0.0
if ($canonical.Count -eq 1 -and $static.Count -eq 1) {
  foreach ($field in @('aggregate_rms','k1_rms','k2_rms','k3_rms','k4_rms','max_channel_rms')) {
    $delta = [Math]::Max($delta, [Math]::Abs((As-Double $canonical[0] $field 0.0) - (As-Double $static[0] $field 0.0)))
  }
}
@([pscustomobject]@{
  canonical_row = if ($canonical.Count -eq 1) { $canonical[0].row_id } else { 'missing' }
  static_row = if ($static.Count -eq 1) { $static[0].row_id } else { 'missing' }
  max_channel_rms_delta = $delta
  canonical_output_hash = if ($canonical.Count -eq 1) { Get-Field $canonical[0] 'output_hash' } else { 'missing' }
  static_output_hash = if ($static.Count -eq 1) { Get-Field $static[0] 'output_hash' } else { 'missing' }
}) | Export-Csv -NoTypeInformation -LiteralPath $ParityCsv
$decision = Decide $rows
@([pscustomobject]@{ decision=$decision; rows=$rows.Count }) | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv
@(
  "decision=$decision",
  "matrix=$MatrixCsv",
  "cache=$CacheCsv",
  "parity=$ParityCsv"
) | Set-Content -LiteralPath $SummaryDat -Encoding ascii
@(
  'TITAN FMM-E51 static M2L translation-table cache gate',
  '======================================================',
  '',
  "Decision: $decision",
  '',
  'E51 adds an opt-in static-translation-table-cache policy for K123 hierarchy M2L matrices. It uses canonical relative box offsets and a keyed lookup table so repeated translation geometries are not rebuilt in the hot path.',
  "Matrix: $MatrixCsv",
  "Cache: $CacheCsv",
  "Parity: $ParityCsv"
) | Set-Content -LiteralPath $ReportPath -Encoding ascii
Write-Host "FMM-E51 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
if ($decision -notin @('fmm_e51_static_translation_table_pass') -and -not $AllowFailedDecision) {
  throw "FMM-E51 failed decision: $decision"
}
