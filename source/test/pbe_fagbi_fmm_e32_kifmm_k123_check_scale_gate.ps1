param(
  [switch]$AllowFailedDecision,
  [switch]$Force,
  [switch]$RerunFailed,
  [switch]$IncludeClosest,
  [int]$SampleCount = 128,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e32_kifmm_k123_check_scale'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e32_matrix.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e32_channel.csv'
$ScaleCsv = Join-Path $WorkRoot 'fmm_e32_check_scale.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e32_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e32_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E32_KIFMM_K123_CHECK_SCALE_REPORT.txt'

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) { throw "Requested C++ compiler not found: $Requested" }
    return $Requested
  }
  foreach ($candidate in @('C:\msys64\ucrt64\bin\clang++.exe','C:\msys64\mingw64\bin\clang++.exe','C:\msys64\ucrt64\bin\g++.exe','C:\msys64\mingw64\bin\g++.exe')) {
    if (Test-Path -LiteralPath $candidate) { return $candidate }
  }
  $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  $cmd = Get-Command g++ -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  throw 'No C++ compiler found for FMM-E32 smoke.'
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

function Assert-SmokeSupportsE32 {
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @('dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian','cartesian_p5_geometry0p4')) {
    if (-not $text.Contains($needle)) { throw "E32 prerequisite missing '$needle' in $SmokeSource" }
  }
}

function Build-Smoke {
  Assert-SmokeSupportsE32
  New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
  $Compiler = Resolve-Cxx $Cxx
  $CompileArgs = @(
    '-std=c++23','-O2','-mavx2','-mfma','-Isrc','-Isrc\eigen',
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
    '-o',$ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) { throw 'failed to build titan_fmm_kifmm_hierarchy_real_panel_smoke' }
  } finally {
    Pop-Location
  }
}

function Scale-Tag([double]$Value) {
  ('{0:G17}' -f $Value).Replace('.', 'p')
}

function New-Case([string]$Sample, [double]$UpwardCheckScale) {
  $tag = Scale-Tag $UpwardCheckScale
  [pscustomobject]@{
    id = "1SH1_${Sample}_${SampleCount}_q8_upcheck${tag}_dual_worklist"
    fixture = '1SH1'
    sample = $Sample
    sample_count = $SampleCount
    order = 8
    theta = 0.6
    upward_check_scale = $UpwardCheckScale
    mode = 'dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian'
  }
}

function Complete-Row([object]$Row, [object]$Case) {
  Ensure-Field $Row 'row_id' $Case.id
  Ensure-Field $Row 'fixture' $Case.fixture
  Ensure-Field $Row 'sample' "$($Case.sample)_$($Case.sample_count)"
  Ensure-Field $Row 'sample_count' $Case.sample_count
  Ensure-Field $Row 'requested_order' $Case.order
  Ensure-Field $Row 'requested_theta' $Case.theta
  Ensure-Field $Row 'requested_upward_check_scale' $Case.upward_check_scale
  Ensure-Field $Row 'mode' $Case.mode
  foreach ($field in @('status','execution_status','accuracy_status','accuracy_failure','aggregate_rms','max_channel_rms','max_channel_name','k1_rms','k2_rms','k3_rms','k4_rms','k1_source','k2_source','k3_source','k4_source','dual_worklist_status','channel_recombination_status','worklist_independence_status','cartesian_k4_accuracy_status','cartesian_k4_rms','kifmm_far_pairs','kifmm_direct_pairs','kifmm_far_fraction','cartesian_k4_far_pairs','cartesian_k4_fallback_pairs','total_seconds','cartesian_k4_total_seconds')) {
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
    "--order=$($Case.order)",
    "--theta=$($Case.theta)",
    "--equiv-scale=1.05",
    "--upward-check-scale=$($Case.upward_check_scale)",
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
  $proc = Start-Process -FilePath $ExePath -ArgumentList $args -NoNewWindow -Wait -PassThru -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
  $output = @()
  if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath }
  if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath }
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($resultLine.Count -ne 1) {
    return (Complete-Row ([pscustomobject][ordered]@{
      row_id = $Case.id; fixture = $Case.fixture; sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count; requested_order = $Case.order; requested_theta = $Case.theta
      requested_upward_check_scale = $Case.upward_check_scale; mode = $Case.mode
      status = 'process_failed'; execution_status = 'fail'; accuracy_status = 'not_evaluated'
      accuracy_failure = 'process_failed'; process_exit_code = $proc.ExitCode
      process_output = (($output | Select-Object -First 20) -join '|')
    }) $Case)
  }
  $row = Complete-Row (Parse-KeyValueLine ([string]$resultLine[-1]) 'E255_RESULT') $Case
  Ensure-Field $row 'process_exit_code' $proc.ExitCode
  $row
}

function Export-Table([array]$Rows, [string]$Path) {
  if ($Rows.Count -eq 0) { Set-Content -LiteralPath $Path -Value '' } else { $Rows | Export-Csv -NoTypeInformation -LiteralPath $Path }
}

function Row-Provenance-Pass([object]$Row) {
  (Get-Field $Row 'execution_status') -eq 'pass' -and
  (Get-Field $Row 'dual_worklist_status') -eq 'pass' -and
  (Get-Field $Row 'channel_recombination_status') -eq 'pass' -and
  (Get-Field $Row 'worklist_independence_status') -eq 'expected' -and
  (Get-Field $Row 'k1_source') -eq 'kifmm_hierarchy_q8' -and
  (Get-Field $Row 'k2_source') -eq 'kifmm_hierarchy_q8' -and
  (Get-Field $Row 'k3_source') -eq 'kifmm_hierarchy_q8' -and
  (Get-Field $Row 'k4_source') -eq 'cartesian_p5_geometry0p4'
}

function Candidate-Pass([object]$Row) {
  (Row-Provenance-Pass $Row) -and
  (As-Double $Row 'k1_rms') -le 5.0e-6 -and
  (As-Double $Row 'k2_rms') -le 5.0e-6 -and
  (As-Double $Row 'k3_rms') -le 5.0e-6 -and
  (As-Double $Row 'k4_rms') -le 5.0e-6 -and
  (As-Double $Row 'aggregate_rms') -le 1.0e-6 -and
  (Get-Field $Row 'cartesian_k4_accuracy_status') -eq 'pass'
}

function Decide([array]$Rows) {
  if ($Rows.Count -lt 3) { return 'fmm_e32_runtime_limited' }
  $passRows = @($Rows | Where-Object { (Get-Field $_ 'execution_status') -eq 'pass' })
  if ($passRows.Count -eq 0) { return 'fmm_e32_execution_gap' }
  if (@($passRows | Where-Object { -not (Row-Provenance-Pass $_) }).Count -gt 0) { return 'fmm_e32_provenance_gap' }
  if (@($passRows | Where-Object { (Get-Field $_ 'cartesian_k4_accuracy_status') -ne 'pass' -or (As-Double $_ 'k4_rms') -gt 5.0e-6 }).Count -gt 0) { return 'fmm_e32_k4_cartesian_gap' }
  if (@($passRows | Where-Object { Candidate-Pass $_ }).Count -gt 0) { return 'fmm_e32_upward_check_candidate_pass' }
  $baseline = @($passRows | Where-Object { [math]::Abs((As-Double $_ 'requested_upward_check_scale') - 3.0) -lt 1.0e-12 } | Select-Object -First 1)
  if ($baseline.Count -eq 1) {
    $best = ($passRows | Sort-Object @{ Expression = { As-Double $_ 'k2_rms' }; Ascending = $true } | Select-Object -First 1)
    if ((As-Double $best 'k2_rms') -lt (0.75 * (As-Double $baseline[0] 'k2_rms'))) { return 'fmm_e32_upward_check_partial_k2_gain' }
  }
  'fmm_e32_no_upward_check_gain'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$scales = @(3.0, 3.5, 4.0)
$cases = @()
foreach ($scale in $scales) { $cases += New-Case 'uniform_hash' $scale }
if ($IncludeClosest) {
  foreach ($scale in $scales) { $cases += New-Case 'closest_opposing_surface' $scale }
}

$existing = @{}
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in @(Import-Csv -LiteralPath $MatrixCsv)) { $existing[(Get-Field $row 'row_id')] = $row }
}

$matrix = @()
foreach ($case in $cases) {
  $old = $existing[$case.id]
  $needsRun = $Force -or $null -eq $old
  if ($RerunFailed -and $null -ne $old -and (Get-Field $old 'execution_status') -ne 'pass') { $needsRun = $true }
  if (-not $needsRun) {
    Write-Host "resume: $($case.id)"
    $matrix += (Complete-Row $old $case)
    continue
  }
  $row = Run-Case $case
  $matrix += $row
  Export-Table $matrix $MatrixCsv
}

Export-Table ($matrix | Select-Object row_id,fixture,sample,requested_order,requested_theta,requested_upward_check_scale,mode,status,execution_status,accuracy_status,accuracy_failure,aggregate_rms,max_channel_rms,max_channel_name,k1_rms,k2_rms,k3_rms,k4_rms,k1_source,k2_source,k3_source,k4_source,dual_worklist_status,channel_recombination_status,worklist_independence_status,cartesian_k4_accuracy_status,cartesian_k4_rms,kifmm_far_pairs,kifmm_direct_pairs,kifmm_far_fraction,cartesian_k4_far_pairs,cartesian_k4_fallback_pairs,total_seconds,cartesian_k4_total_seconds) $ChannelCsv
Export-Table ($matrix | Select-Object row_id,requested_upward_check_scale,k1_rms,k2_rms,k3_rms,k4_rms,aggregate_rms,max_channel_name,max_channel_rms,kifmm_far_fraction,total_seconds) $ScaleCsv

$decision = Decide $matrix
[pscustomobject][ordered]@{
  decision = $decision
  fixture = '1SH1'
  rows = $matrix.Count
  sample_count = $SampleCount
  theta = 0.6
  order = 8
  include_closest = [bool]$IncludeClosest
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @("decision=$decision","rows=$($matrix.Count)","matrix=$MatrixCsv")

$report = @(
  '# TITAN FMM Engine E32 KIFMM K123 Upward Check-Scale Report',
  '',
  "Decision: $decision",
  '',
  'E32 is a 1SH1-only sampled diagnostic. It keeps KIFMM q8 and theta 0.6 for K1/K2/K3, changes only upward_check_scale, and keeps K4 on the independent Cartesian p5 geometry0p4 branch.',
  'Rows test upward_check_scale = 3.0, 3.5, and 4.0. High-scale execution failures are reported as diagnostics; the decision is based on completed rows unless all rows fail. No production TABI, GMRES, manuscript, or default behavior is changed.'
)
$executionFailures = @($matrix | Where-Object { (Get-Field $_ 'execution_status') -ne 'pass' }).Count
$report += ("Completed rows: {0}/{1}; execution diagnostic failures: {2}" -f ($matrix.Count - $executionFailures), $matrix.Count, $executionFailures)
foreach ($row in $matrix) {
  $report += ('  {0}: upcheck={1}, exec={2}, acc={3}, status={4}, max={5}/{6}, K1={7}, K2={8}, K3={9}, K4={10}, cartK4={11}, far={12}/{13}/{14}, seconds={15}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'requested_upward_check_scale'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'status'), (Get-Field $row 'max_channel_name'), `
      (Get-Field $row 'max_channel_rms'), (Get-Field $row 'k1_rms'), (Get-Field $row 'k2_rms'), `
      (Get-Field $row 'k3_rms'), (Get-Field $row 'k4_rms'), (Get-Field $row 'cartesian_k4_rms'), `
      (Get-Field $row 'kifmm_far_pairs'), (Get-Field $row 'kifmm_direct_pairs'), (Get-Field $row 'kifmm_far_fraction'), (Get-Field $row 'total_seconds'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E32 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Channel: $ChannelCsv"
Write-Host "Scale: $ScaleCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e32_upward_check_candidate_pass')) {
  throw "FMM-E32 decision failed: $decision"
}
