param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e14_level_cube_policy'),
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
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E14_LEVEL_CUBE_POLICY_REPORT.txt'
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

$MatrixPath = Join-Path $WorkRoot 'fmm_e14_matrix.csv'
$AdversarialPath = Join-Path $WorkRoot 'fmm_e14_adversarial.csv'
$QOrderPath = Join-Path $WorkRoot 'fmm_e14_q_order.csv'
$FarFractionPath = Join-Path $WorkRoot 'fmm_e14_far_fraction.csv'
$WorstRowsPath = Join-Path $WorkRoot 'fmm_e14_worst_rows.csv'
$CartesianContextPath = Join-Path $WorkRoot 'fmm_e14_cartesian_context.csv'
$DecisionPath = Join-Path $WorkRoot 'fmm_e14_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_e14_summary.dat'

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
  throw 'No C++ compiler found for KIFMM-lite E14 gate.'
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
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E14 sampled profile driver' }
  } finally {
    Pop-Location
  }
  return $exe
}

function New-Row([string]$Fixture, [string]$Sample, [int]$SampleCount,
                 [int]$Order, [string]$Role, [bool]$Candidate) {
  $sampleLabel = if ($Sample -eq 'uniform_hash') { 'uniform' } else { $Sample }
  [pscustomobject]@{
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    order = $Order
    theta = 0.6
    role = $Role
    candidate = $Candidate
    row_id = ('{0}_{1}_{2}_q{3}_level_cube' -f $Fixture, $sampleLabel, $SampleCount, $Order)
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
      'fixture','sample','status','source_panels','sampled_rows','order',
      'source_box_policy','min_far_panel_fraction','aggregate_rms','k1','k2',
      'k3','k4','max_equiv_residual','far_panel_pairs','fallback_panel_pairs',
      'production_direct_reference_valid_pairs','invalid_pairs',
      'nonproduction_fallback_pairs','cache_hit_rate','cache_entries',
      'unique_half_widths','total_seconds','worklist_hash','policy_hash')) {
    if (-not $fields.Contains($required)) { throw "malformed E14 row, missing $required`: $Line" }
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

function Run-Row($Exe, $Row) {
  $rowDir = Join-Path $RunRoot $Row.row_id
  New-Item -ItemType Directory -Path $rowDir -Force | Out-Null
  $rawLog = Join-Path $rowDir 'raw.log'
  $parsed = Parse-ResultFromLog $rawLog
  $shouldRun = $Force -or $null -eq $parsed
  if ($RerunFailed -and ($null -eq $parsed -or [string]$parsed.status -notin @('pass','aggregate_rms_gate_failed','channel_rms_gate_failed'))) { $shouldRun = $true }
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
      '--theta=0.6',
      '--min-far-panel-fraction=0.10',
      '--condition-policy=qr-rank-only',
      '--source-box-policy=level-cube'
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
      throw "E14 row failed before producing parseable output: $($Row.row_id)"
    }
  }
  if ($null -eq $parsed) { throw "no parseable result for $($Row.row_id)" }
  $farPairs = As-Long $parsed.far_panel_pairs
  $fallbackPairs = As-Long $parsed.fallback_panel_pairs
  $totalPairs = $farPairs + $fallbackPairs
  [pscustomobject]@{
    row_id = $Row.row_id
    fixture = $Row.fixture
    sample_set = $Row.sample
    sample_count = $Row.sample_count
    role = $Row.role
    candidate = $Row.candidate
    order = As-Long $parsed.order
    theta = As-Double $parsed.theta
    status = [string]$parsed.status
    source_panels = As-Long $parsed.source_panels
    sampled_rows = As-Long $parsed.sampled_rows
    source_box_policy = [string]$parsed.source_box_policy
    min_far_panel_fraction = As-Double $parsed.min_far_panel_fraction
    aggregate_rms = As-Double $parsed.aggregate_rms
    k1_rms = As-Double $parsed.k1
    k2_rms = As-Double $parsed.k2
    k3_rms = As-Double $parsed.k3
    k4_rms = As-Double $parsed.k4
    max_channel_rms = [Math]::Max([Math]::Max((As-Double $parsed.k1), (As-Double $parsed.k2)), [Math]::Max((As-Double $parsed.k3), (As-Double $parsed.k4)))
    max_equiv_residual = As-Double $parsed.max_equiv_residual
    invalid_pairs = As-Long $parsed.invalid_pairs
    nonproduction_fallback_pairs = As-Long $parsed.nonproduction_fallback_pairs
    production_direct_reference_valid_pairs = As-Long $parsed.production_direct_reference_valid_pairs
    far_panel_pairs = $farPairs
    fallback_panel_pairs = $fallbackPairs
    far_panel_fraction = if ($totalPairs -eq 0) { 0.0 } else { [double]$farPairs / [double]$totalPairs }
    fallback_panel_fraction = if ($totalPairs -eq 0) { 0.0 } else { [double]$fallbackPairs / [double]$totalPairs }
    equivalent_cache_hit_rate = As-Double $parsed.cache_hit_rate
    equivalent_system_cache_entries = As-Long $parsed.cache_entries
    unique_equivalent_half_widths = As-Long $parsed.unique_half_widths
    total_seconds = As-Double $parsed.total_seconds
    worklist_hash = [string]$parsed.worklist_hash
    policy_hash = [string]$parsed.policy_hash
  }
}

$quickRows = @(
  (New-Row '1SH1' 'uniform_hash' 128 4 'q_order' $false),
  (New-Row '1SH1' 'uniform_hash' 128 6 'q_order' $true),
  (New-Row '1SH1' 'uniform_hash' 128 8 'q_order' $false),
  (New-Row '1A63' 'high_aspect' 128 6 'quick_1a63_adversarial' $true)
)
$fullRows = @(
  $quickRows
  (New-Row '1SH1' 'closest_opposing_surface' 128 6 '1sh1_adversarial' $true),
  (New-Row '1SH1' 'high_aspect' 128 6 '1sh1_adversarial' $true),
  (New-Row '1SH1' 'inward_normal_groove' 128 6 '1sh1_adversarial' $true),
  (New-Row '1A63' 'uniform_hash' 128 6 '1a63_uniform' $true),
  (New-Row '1A63' 'uniform_hash' 512 6 '1a63_uniform' $true),
  (New-Row '1A63' 'closest_opposing_surface' 128 6 '1a63_adversarial' $true),
  (New-Row '1A63' 'inward_normal_groove' 128 6 '1a63_adversarial' $true)
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
  if ($rows.Count -eq 0) { throw 'RerunRows did not match any E14 row id' }
}

$exe = Build-Driver
$parsedRows = @()
foreach ($row in $rows) {
  $parsedRows += Run-Row $exe $row
}

$parsedRows | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation -Encoding ascii
$parsedRows | Where-Object { $_.sample_set -ne 'uniform_hash' } |
  Export-Csv -LiteralPath $AdversarialPath -NoTypeInformation -Encoding ascii
$parsedRows | Where-Object { $_.fixture -eq '1SH1' -and $_.sample_set -eq 'uniform_hash' -and $_.sample_count -eq 128 } |
  Sort-Object order |
  Export-Csv -LiteralPath $QOrderPath -NoTypeInformation -Encoding ascii
$parsedRows | Select-Object row_id,fixture,sample_set,order,status,min_far_panel_fraction,far_panel_pairs,fallback_panel_pairs,far_panel_fraction,fallback_panel_fraction |
  Export-Csv -LiteralPath $FarFractionPath -NoTypeInformation -Encoding ascii

function Select-WorstRow($Rows, [string]$Metric, [switch]$Ascending) {
  $selected = if ($Ascending) {
    $Rows | Sort-Object -Property $Metric | Select-Object -First 1
  } else {
    $Rows | Sort-Object -Property $Metric -Descending | Select-Object -First 1
  }
  if ($null -eq $selected) { return $null }
  [pscustomobject]@{
    metric = $Metric
    row_id = $selected.row_id
    fixture = $selected.fixture
    sample_set = $selected.sample_set
    order = $selected.order
    status = $selected.status
    value = $selected.$Metric
  }
}

@(
  (Select-WorstRow $parsedRows 'aggregate_rms')
  (Select-WorstRow $parsedRows 'k1_rms')
  (Select-WorstRow $parsedRows 'k2_rms')
  (Select-WorstRow $parsedRows 'k3_rms')
  (Select-WorstRow $parsedRows 'k4_rms')
  (Select-WorstRow $parsedRows 'far_panel_fraction' -Ascending)
  (Select-WorstRow $parsedRows 'total_seconds')
) | Where-Object { $null -ne $_ } |
  Export-Csv -LiteralPath $WorstRowsPath -NoTypeInformation -Encoding ascii

[pscustomobject]@{
  cartesian_context_status = 'not_available_in_current_harness'
  reason = 'No clean sampled-only guarded Cartesian p4 comparator is exposed in the KIFMM E14 worktree; keep for E15.'
  blocks_decision = $false
} | Export-Csv -LiteralPath $CartesianContextPath -NoTypeInformation -Encoding ascii

function Test-CandidatePass($row) {
  if ($null -eq $row) { return $false }
  $expected = $row.source_panels * $row.sampled_rows
  return $row.status -eq 'pass' -and
    $row.aggregate_rms -le 1.0e-6 -and
    $row.max_channel_rms -le 5.0e-6 -and
    $row.max_equiv_residual -le 1.0e-10 -and
    $row.invalid_pairs -eq 0 -and
    $row.nonproduction_fallback_pairs -eq 0 -and
    (($row.far_panel_pairs + $row.fallback_panel_pairs) -eq $expected) -and
    $row.production_direct_reference_valid_pairs -eq $expected -and
    $row.equivalent_cache_hit_rate -ge 0.90 -and
    $row.far_panel_fraction -ge 0.10
}

$candidateRows = @($parsedRows | Where-Object { $_.candidate })
$candidateGaps = @($candidateRows | Where-Object { -not (Test-CandidatePass $_) })
$fallbackDominated = @($candidateRows | Where-Object { $_.status -eq 'fallback_dominated' -or $_.far_panel_fraction -lt 0.10 })
$oneA63Adversarial = @($candidateRows | Where-Object { $_.fixture -eq '1A63' -and $_.role -eq '1a63_adversarial' })
$oneA63HighAspect = @($candidateRows | Where-Object { $_.fixture -eq '1A63' -and $_.role -eq 'quick_1a63_adversarial' })
$oneA63AdversarialLike = @($candidateRows | Where-Object { $_.fixture -eq '1A63' -and ($_.role -eq '1a63_adversarial' -or $_.role -eq 'quick_1a63_adversarial') })
$oneA63AdversarialGaps = @($oneA63Adversarial | Where-Object { -not (Test-CandidatePass $_) })
$oneA63AdversarialLikeGaps = @($oneA63AdversarialLike | Where-Object { -not (Test-CandidatePass $_) })
$q6 = $parsedRows | Where-Object { $_.row_id -eq '1SH1_uniform_128_q6_level_cube' } | Select-Object -First 1
$q8 = $parsedRows | Where-Object { $_.row_id -eq '1SH1_uniform_128_q8_level_cube' } | Select-Object -First 1
$q8Stable = $false
if ($null -ne $q6 -and $null -ne $q8 -and $q8.status -eq 'pass') {
  $q8Stable = ($q8.aggregate_rms -le ($q6.aggregate_rms * 1.10)) -and
    ($q8.max_channel_rms -le ($q6.max_channel_rms * 1.10))
}

$decision = 'fmm_e14_runtime_limited'
if ($fallbackDominated.Count -gt 0) {
  $kifmmSubgateDecision = 'fmm_e14_fallback_dominated'
} elseif ($oneA63AdversarialGaps.Count -gt 0) {
  $kifmmSubgateDecision = 'fmm_e14_1a63_adversarial_gap'
} elseif ($candidateGaps.Count -gt 0) {
  $kifmmSubgateDecision = 'fmm_e14_candidate_gap'
} elseif (-not $q8Stable) {
  $kifmmSubgateDecision = 'fmm_e14_q_order_gap'
} elseif (-not $Quick) {
  $kifmmSubgateDecision = 'fmm_e14_level_cube_1sh1_1a63_adversarial_pass'
} else {
  $kifmmSubgateDecision = 'fmm_e14_level_cube_quick_pass_full_gate_pending'
}

$decision = $kifmmSubgateDecision
if ($kifmmSubgateDecision -eq 'fmm_e14_level_cube_1sh1_1a63_adversarial_pass') {
  $decision = 'fmm_e14_cartesian_context_pending'
} elseif (-not $Quick) {
  $decision = $kifmmSubgateDecision
}

[pscustomobject]@{
  decision = $decision
  kifmm_sampled_row_decision = $kifmmSubgateDecision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $parsedRows.Count
  candidate_rows = $candidateRows.Count
  candidate_gaps = $candidateGaps.Count
  fallback_dominated_rows = $fallbackDominated.Count
  one_a63_adversarial_rows = $oneA63Adversarial.Count
  one_a63_adversarial_gaps = $oneA63AdversarialGaps.Count
  one_a63_high_aspect_rows = $oneA63HighAspect.Count
  one_a63_adversarial_like_rows = $oneA63AdversarialLike.Count
  one_a63_adversarial_like_gaps = $oneA63AdversarialLikeGaps.Count
  q8_stable_or_improved = $q8Stable
  cartesian_context_status = 'not_available_in_current_harness'
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation -Encoding ascii

@(
  "decision=$decision",
  "kifmm_sampled_row_decision=$kifmmSubgateDecision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($parsedRows.Count)",
  "matrix=$MatrixPath",
  "adversarial=$AdversarialPath",
  "q_order=$QOrderPath",
  "far_fraction=$FarFractionPath",
  "worst_rows=$WorstRowsPath",
  "cartesian_context=$CartesianContextPath"
) | Set-Content -LiteralPath $SummaryPath -Encoding ascii

@"
TITAN FMM-E14 Level-Cube Policy Hardening Gate

KIFMM sampled-row decision: $kifmmSubgateDecision
Overall E14 decision: $decision
Mode: $(if ($Quick) { 'quick' } else { 'full' })

E14 hardens the E13 level-stable source-box policy by enforcing a code-level
minimum far-panel fraction of 0.10, extending adversarial sampled rows to 1A63,
and reporting a q4/q6/q8 ladder on 1SH1 uniform-128. KIFMM-lite remains
diagnostic-only and is not wired into TABI production matvec, GMRES, or
preconditioners.

Scientific status:
The KIFMM-lite sampled-row subgate is separated from the overall E14 decision.
When the sampled-row subgate passes, E14 still remains pending if the same-run
guarded Cartesian p4 sampled comparator is unavailable. No KIFMM-vs-Cartesian
runtime claim is made in E14.

Artifacts:
- Matrix: $MatrixPath
- Adversarial rows: $AdversarialPath
- q-order ladder: $QOrderPath
- far-fraction gate: $FarFractionPath
- worst-row summary: $WorstRowsPath
- Cartesian context: $CartesianContextPath
- Decision: $DecisionPath
- Summary: $SummaryPath

Cartesian p4 context:
No clean sampled-only guarded Cartesian p4 comparator is available in the E14
KIFMM worktree. The KIFMM candidate rows are therefore classified independently,
and the Cartesian timing comparison is deferred to E15.
"@ | Set-Content -LiteralPath $ReportPath -Encoding ascii

Write-Host "FMM E14 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Adversarial: $AdversarialPath"
Write-Host "q-order: $QOrderPath"
Write-Host "far-fraction: $FarFractionPath"
if ($decision -notin @('fmm_e14_level_cube_quick_pass_full_gate_pending','fmm_e14_level_cube_cross_fixture_adversarial_pass','fmm_e14_cartesian_context_pending') -and -not $AllowFailedDecision) {
  throw "FMM E14 gate did not pass: $decision"
}
