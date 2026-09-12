param(
  [switch]$AllowFailedDecision,
  [switch]$Force,
  [switch]$RerunFailed,
  [int[]]$Orders = @(10),
  [int]$SampleCount = 128,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e33_kifmm_k123_order10'),
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
$MatrixCsv = Join-Path $WorkRoot 'fmm_e33_matrix.csv'
$OrderCsv = Join-Path $WorkRoot 'fmm_e33_order.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e33_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e33_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E33_KIFMM_K123_Q10_DIAGNOSTIC_REPORT.txt'

foreach ($order in $Orders) {
  if ($order -lt 6 -or $order -gt 10) {
    throw "E33 diagnostic accepts only KIFMM hierarchy orders 6..10; got $order"
  }
}

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
  throw 'No C++ compiler found for FMM-E33 smoke.'
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
  if ($value -ieq 'nan') { return [double]::NaN }
  if ($value -ieq 'inf' -or $value -ieq '+inf' -or $value -ieq 'infinity') { return [double]::PositiveInfinity }
  if ($value -ieq '-inf' -or $value -ieq '-infinity') { return [double]::NegativeInfinity }
  [double]$value
}

function Ensure-Field([object]$Row, [string]$Name, [object]$DefaultValue) {
  if ($null -eq $Row.PSObject.Properties[$Name]) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $DefaultValue
  }
}

function Build-Smoke {
  $text = Get-Content -Raw -LiteralPath $SmokeSource
  foreach ($needle in @('dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian','cartesian_p5_geometry0p4','get_arg(args, "order"')) {
    if (-not $text.Contains($needle)) { throw "E33 prerequisite missing '$needle' in $SmokeSource" }
  }
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

function New-Case([int]$Order) {
  [pscustomobject]@{
    id = "1SH1_uniform_hash_${SampleCount}_q${Order}_dual_worklist"
    fixture = '1SH1'
    sample = 'uniform_hash'
    sample_count = $SampleCount
    order = $Order
    theta = 0.6
    upward_check_scale = 3.0
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
  foreach ($field in @('status','execution_status','accuracy_status','accuracy_failure','aggregate_rms','max_channel_rms','max_channel_name','k1_rms','k2_rms','k3_rms','k4_rms','k1_source','k2_source','k3_source','k4_source','dual_worklist_status','channel_recombination_status','worklist_independence_status','cartesian_k4_accuracy_status','cartesian_k4_rms','kifmm_far_pairs','kifmm_direct_pairs','kifmm_far_fraction','total_seconds')) {
    Ensure-Field $Row $field 'MISSING'
  }
  $k2 = As-Double $Row 'k2_rms'
  $Row | Add-Member -Force -NotePropertyName k2_relaxed_gate -NotePropertyValue '1e-5'
  $Row | Add-Member -Force -NotePropertyName k2_relaxed_status -NotePropertyValue $(if ([double]::IsNaN($k2)) { 'not_evaluated' } elseif ($k2 -le 1.0e-5) { 'pass' } else { 'fail' })
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
    "--condition-policy=rank-bypass-diagnostic",
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
  (Get-Field $Row 'k4_source') -eq 'cartesian_p5_geometry0p4'
}

function Strict-Pass([object]$Row) {
  (Row-Provenance-Pass $Row) -and
  (As-Double $Row 'k1_rms') -le 5.0e-6 -and
  (As-Double $Row 'k2_rms') -le 5.0e-6 -and
  (As-Double $Row 'k3_rms') -le 5.0e-6 -and
  (As-Double $Row 'k4_rms') -le 5.0e-6 -and
  (As-Double $Row 'aggregate_rms') -le 1.0e-6 -and
  (Get-Field $Row 'cartesian_k4_accuracy_status') -eq 'pass'
}

function Decide([array]$Rows) {
  if ($Rows.Count -lt $Orders.Count) { return 'fmm_e33_runtime_limited' }
  $passRows = @($Rows | Where-Object { (Get-Field $_ 'execution_status') -eq 'pass' })
  if ($passRows.Count -eq 0) { return 'fmm_e33_execution_gap' }
  if (@($passRows | Where-Object { -not (Row-Provenance-Pass $_) }).Count -gt 0) { return 'fmm_e33_provenance_gap' }
  if (@($passRows | Where-Object { (Get-Field $_ 'cartesian_k4_accuracy_status') -ne 'pass' -or (As-Double $_ 'k4_rms') -gt 5.0e-6 }).Count -gt 0) { return 'fmm_e33_k4_cartesian_gap' }
  if (@($passRows | Where-Object { Strict-Pass $_ }).Count -gt 0) { return 'fmm_e33_q10_strict_k2_pass' }
  if (@($passRows | Where-Object { [int](Get-Field $_ 'requested_order' '0') -eq 10 -and (Get-Field $_ 'k2_relaxed_status') -eq 'pass' }).Count -gt 0) { return 'fmm_e33_q10_relaxed_k2_pass_strict_gap' }
  'fmm_e33_q10_no_k2_gain'
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$cases = @($Orders | ForEach-Object { New-Case $_ })
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

Export-Table ($matrix | Select-Object row_id,fixture,sample,requested_order,requested_theta,requested_upward_check_scale,mode,status,execution_status,accuracy_status,accuracy_failure,aggregate_rms,max_channel_rms,max_channel_name,k1_rms,k2_rms,k2_relaxed_gate,k2_relaxed_status,k3_rms,k4_rms,k1_source,k2_source,k3_source,k4_source,dual_worklist_status,channel_recombination_status,worklist_independence_status,cartesian_k4_accuracy_status,cartesian_k4_rms,kifmm_far_pairs,kifmm_direct_pairs,kifmm_far_fraction,total_seconds) $OrderCsv

$decision = Decide $matrix
[pscustomobject][ordered]@{
  decision = $decision
  fixture = '1SH1'
  rows = $matrix.Count
  sample_count = $SampleCount
  theta = 0.6
  orders = ($Orders -join ';')
  k2_strict_gate = '5e-6'
  k2_relaxed_gate = '1e-5'
} | Export-Csv -NoTypeInformation -LiteralPath $DecisionCsv

Set-Content -LiteralPath $SummaryDat -Value @("decision=$decision","rows=$($matrix.Count)","matrix=$MatrixCsv")

$report = @(
  '# TITAN FMM Engine E33 KIFMM K123 q10 Diagnostic Report',
  '',
  "Decision: $decision",
  '',
  'E33 is a 1SH1-only sampled diagnostic. It keeps KIFMM theta 0.6, upward_check_scale 3.0, strict-local-check, and independent Cartesian p5 geometry0p4 K4.',
  'The harness accepts hierarchy orders 6..10. Production/default policy is unchanged; q6 remains the production-facing baseline outside this diagnostic.',
  'q10 uses condition-policy=rank-bypass-diagnostic so rank-deficient equivalent systems can be measured instead of fail-closed; this is not a production policy.',
  'K2 is reported against the strict 5e-6 channel gate and a diagnostic relaxed 1e-5 K2 gate.'
)
foreach ($row in $matrix) {
  $report += ('  {0}: q={1}, exec={2}, acc={3}, status={4}, max={5}/{6}, K1={7}, K2={8} ({9}), K3={10}, K4={11}, far={12}/{13}/{14}, seconds={15}' -f `
      (Get-Field $row 'row_id'), (Get-Field $row 'requested_order'), (Get-Field $row 'execution_status'), `
      (Get-Field $row 'accuracy_status'), (Get-Field $row 'status'), (Get-Field $row 'max_channel_name'), `
      (Get-Field $row 'max_channel_rms'), (Get-Field $row 'k1_rms'), (Get-Field $row 'k2_rms'), `
      (Get-Field $row 'k2_relaxed_status'), (Get-Field $row 'k3_rms'), (Get-Field $row 'k4_rms'), `
      (Get-Field $row 'kifmm_far_pairs'), (Get-Field $row 'kifmm_direct_pairs'), (Get-Field $row 'kifmm_far_fraction'), (Get-Field $row 'total_seconds'))
}
Set-Content -LiteralPath $ReportPath -Value $report

Write-Host "FMM-E33 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Order: $OrderCsv"
Write-Host "Report: $ReportPath"

if (-not $AllowFailedDecision -and $decision -notin @('fmm_e33_q10_strict_k2_pass','fmm_e33_q10_relaxed_k2_pass_strict_gap')) {
  throw "FMM-E33 decision failed: $decision"
}
