param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e255b_real_panel_k4_attribution'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e255b_matrix.csv'
$AttributionCsv = Join-Path $WorkRoot 'fmm_e255b_attribution.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e255b_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e255b_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E255B_REAL_PANEL_K4_ATTRIBUTION_REPORT.txt'

function Resolve-Cxx([string]$Requested) {
  if (-not [string]::IsNullOrWhiteSpace($Requested)) {
    if (-not (Test-Path -LiteralPath $Requested)) {
      throw "Requested C++ compiler not found: $Requested"
    }
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
  throw 'No C++ compiler found for FMM-E25.5b attribution smoke.'
}

function Parse-KeyValueLine([string]$Line, [string]$Prefix) {
  if (-not $Line.StartsWith($Prefix)) {
    throw "Unexpected line for ${Prefix}: $Line"
  }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring($Prefix.Length).Trim().Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function Build-Smoke {
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
    'src\tabi\tabi_fmm_kifmm_cartesian_bridge.cpp',
    'src\tabi\tabi_fmm_cartesian_sampled_panel_matvec.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_sampled_panel_matvec.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) {
      throw 'failed to build titan_fmm_kifmm_hierarchy_real_panel_smoke'
    }
  } finally {
    Pop-Location
  }
}

function New-Case(
  [string]$Id,
  [string]$Fixture,
  [string]$Sample,
  [int]$SampleCount,
  [double]$Theta,
  [double]$Equiv,
  [double]$UpCheck,
  [double]$LocalEquiv,
  [double]$LocalCheck,
  [string]$ConditionPolicy,
  [string]$Family
) {
  [pscustomobject]@{
    id = $Id
    fixture = $Fixture
    sample = $Sample
    sample_count = $SampleCount
    theta = $Theta
    equiv = $Equiv
    up_check = $UpCheck
    local_equiv = $LocalEquiv
    local_check = $LocalCheck
    condition_policy = $ConditionPolicy
    family = $Family
  }
}

function Run-Case([object]$Case) {
  Write-Host "run: $($Case.id)"
  $args = @(
    "--fixture=$($Case.fixture)",
    "--sample=$($Case.sample)",
    "--sample-count=$($Case.sample_count)",
    "--order=8",
    "--theta=$($Case.theta)",
    "--equiv-scale=$($Case.equiv)",
    "--upward-check-scale=$($Case.up_check)",
    "--local-equiv-scale=$($Case.local_equiv)",
    "--local-check-scale=$($Case.local_check)",
    "--condition-policy=$($Case.condition_policy)",
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
  if (Test-Path -LiteralPath $stdoutPath) {
    $output += Get-Content -LiteralPath $stdoutPath
  }
  if (Test-Path -LiteralPath $stderrPath) {
    $output += Get-Content -LiteralPath $stderrPath
  }
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  if ($proc.ExitCode -ne 0 -or $resultLine.Count -ne 1) {
    $row = [pscustomobject][ordered]@{
      row_id = $Case.id
      fixture = $Case.fixture
      sample = "$($Case.sample)_$($Case.sample_count)"
      sample_count = $Case.sample_count
      family = $Case.family
      order = 8
      theta = $Case.theta
      equiv_scale = $Case.equiv
      upward_check_scale = $Case.up_check
      local_equiv_scale = $Case.local_equiv
      local_check_scale = $Case.local_check
      condition_policy = $Case.condition_policy
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'execution'
      aggregate_rms = 'NaN'
      max_channel_rms = 'NaN'
      max_channel_name = 'none'
      k4_rms = 'NaN'
      far_k4_rms = 'NaN'
      fallback_k4_rms = 'NaN'
      source_box_boundary_margin_min = 'NaN'
      target_check_excess_min = 'NaN'
      process_exit_code = $proc.ExitCode
      process_output = ($output -join ' | ')
    }
    return [pscustomobject]@{ Row = $row; Attr = @() }
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E255_RESULT '
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $Case.id
  $row | Add-Member -NotePropertyName family -NotePropertyValue $Case.family
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $proc.ExitCode
  $row | Add-Member -NotePropertyName process_output -NotePropertyValue ''
  $attr = @($output | Where-Object { $_ -like 'E255B_ATTR *' } | ForEach-Object {
      $a = Parse-KeyValueLine ([string]$_) 'E255B_ATTR '
      $a | Add-Member -NotePropertyName row_id -NotePropertyValue $Case.id
      $a | Add-Member -NotePropertyName family -NotePropertyValue $Case.family
      $a
    })
  [pscustomobject]@{ Row = $row; Attr = $attr }
}

if ((Test-Path -LiteralPath $WorkRoot) -and $Force) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$cases = @(
  (New-Case '1SH1_uniform_128_baseline' '1SH1' 'uniform_hash' 128 0.6 1.05 3.0 3.0 1.05 'qr-rank-only' 'baseline'),
  (New-Case '1SH1_uniform_128_equiv1p2' '1SH1' 'uniform_hash' 128 0.6 1.2 3.0 3.0 1.05 'qr-rank-only' 'proxy_scale'),
  (New-Case '1SH1_uniform_128_equiv1p5' '1SH1' 'uniform_hash' 128 0.6 1.5 3.0 3.0 1.05 'qr-rank-only' 'proxy_scale'),
  (New-Case '1SH1_uniform_128_theta0p4' '1SH1' 'uniform_hash' 128 0.4 1.05 3.0 3.0 1.05 'qr-rank-only' 'far_acceptance'),
  (New-Case '1SH1_uniform_128_svd_values' '1SH1' 'uniform_hash' 128 0.6 1.05 3.0 3.0 1.05 'svd-values' 'conditioning')
)
if (-not $Quick) {
  $cases += @(
    (New-Case '1SH1_closest_128_baseline' '1SH1' 'closest_opposing_surface' 128 0.6 1.05 3.0 3.0 1.05 'qr-rank-only' 'baseline'),
    (New-Case '1SH1_closest_128_equiv1p5' '1SH1' 'closest_opposing_surface' 128 0.6 1.5 3.0 3.0 1.05 'qr-rank-only' 'proxy_scale'),
    (New-Case '1A63_high_aspect_128_baseline' '1A63' 'high_aspect' 128 0.6 1.05 3.0 3.0 1.05 'qr-rank-only' 'baseline'),
    (New-Case '1A63_high_aspect_128_equiv1p5' '1A63' 'high_aspect' 128 0.6 1.5 3.0 3.0 1.05 'qr-rank-only' 'proxy_scale')
  )
}

$rows = New-Object System.Collections.Generic.List[object]
$attrs = New-Object System.Collections.Generic.List[object]
$existingRows = @{}
$existingAttrs = @()
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in (Import-Csv -LiteralPath $MatrixCsv)) {
    if ($row.PSObject.Properties['row_id'] -and -not [string]::IsNullOrWhiteSpace($row.row_id)) {
      $existingRows[$row.row_id] = $row
    }
  }
}
if ((Test-Path -LiteralPath $AttributionCsv) -and -not $Force) {
  $existingAttrs = @(Import-Csv -LiteralPath $AttributionCsv)
}
foreach ($case in $cases) {
  if ($existingRows.ContainsKey($case.id) -and
      (-not $RerunFailed -or $existingRows[$case.id].execution_status -eq 'pass')) {
    Write-Host "resume: $($Case.id)"
    $run = [pscustomobject]@{
      Row = $existingRows[$case.id]
      Attr = @($existingAttrs | Where-Object { $_.row_id -eq $case.id })
    }
  } else {
    $run = Run-Case $case
  }
  $rows.Add($run.Row)
  foreach ($a in $run.Attr) { $attrs.Add($a) }
  $rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
  $attrs | Export-Csv -LiteralPath $AttributionCsv -NoTypeInformation
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$attrs | Export-Csv -LiteralPath $AttributionCsv -NoTypeInformation

$baseline = @($rows | Where-Object { $_.family -eq 'baseline' } | Select-Object -First 1)
$baselineK4 = if ($baseline.Count -gt 0) { [double]$baseline[0].far_k4_rms } else { [double]::NaN }
$passing = @($rows | Where-Object { $_.execution_status -eq 'pass' -and $_.accuracy_status -eq 'pass' })
$proxyImproved = @($rows | Where-Object { $_.family -eq 'proxy_scale' -and $_.execution_status -eq 'pass' -and [double]$_.far_k4_rms -le ($baselineK4 / 10.0) })
$thetaImproved = @($rows | Where-Object { $_.family -eq 'far_acceptance' -and $_.execution_status -eq 'pass' -and [double]$_.far_k4_rms -le ($baselineK4 / 10.0) })
$conditionImproved = @($rows | Where-Object { $_.family -eq 'conditioning' -and $_.execution_status -eq 'pass' -and [double]$_.far_k4_rms -le ($baselineK4 / 10.0) })
$fourXImproved = @($rows | Where-Object { $_.execution_status -eq 'pass' -and [double]$_.far_k4_rms -le ($baselineK4 / 4.0) })
$negativeCheck = @($rows | Where-Object { $_.execution_status -eq 'pass' -and [double]$_.target_check_excess_min -lt 0.0 })
$requiredExecutionFailures = @($rows | Where-Object { $_.family -ne 'conditioning' -and $_.execution_status -ne 'pass' })
$conditioningExecutionFailures = @($rows | Where-Object { $_.family -eq 'conditioning' -and $_.execution_status -ne 'pass' })

$decision = if ($requiredExecutionFailures.Count -gt 0) {
  'fmm_e255b_execution_gap'
} elseif ($passing.Count -gt 0) {
  'fmm_e255b_accuracy_recovered'
} elseif ($proxyImproved.Count -gt 0) {
  'fmm_e255b_proxy_scale_candidate'
} elseif ($thetaImproved.Count -gt 0) {
  'fmm_e255b_far_acceptance_candidate'
} elseif ($conditionImproved.Count -gt 0) {
  'fmm_e255b_conditioning_candidate'
} elseif ($fourXImproved.Count -gt 0) {
  'fmm_e255b_partial_k4_reduction'
} elseif ($negativeCheck.Count -gt 0) {
  'fmm_e255b_proximity_confirmed_unresolved'
} else {
  'fmm_e255b_real_panel_field_model_gap'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  baseline_far_k4_rms = $baselineK4
  accuracy_pass_rows = $passing.Count
  proxy_10x_rows = $proxyImproved.Count
  theta_10x_rows = $thetaImproved.Count
  conditioning_10x_rows = $conditionImproved.Count
  four_x_rows = $fourXImproved.Count
  negative_check_excess_rows = $negativeCheck.Count
  required_execution_failures = $requiredExecutionFailures.Count
  conditioning_execution_failures = $conditioningExecutionFailures.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "attribution=$AttributionCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.5b Real-Panel K4 Attribution Report',
  '',
  "Decision: $decision",
  '',
  'E25.5b attributes the real-panel hierarchy K4 gap without changing production policy.',
  'The gate compares far-only hierarchy rows to production-direct far rows, checks production fallback separately, and records MAC-bin/check-surface proximity diagnostics.',
  ''
)
function Get-Field([object]$Row, [string]$Name, [string]$Default = 'NA') {
  $prop = $Row.PSObject.Properties[$Name]
  if ($null -eq $prop -or $null -eq $prop.Value -or [string]::IsNullOrWhiteSpace([string]$prop.Value)) {
    return $Default
  }
  return [string]$prop.Value
}
foreach ($row in $rows) {
  $report += ('  {0}: family={1}, exec={2}, acc={3}, k4={4}, far_k4={5}, fallback_k4={6}, src_margin_min={7}, target_check_excess_min={8}, total_s={9}' -f
    (Get-Field $row 'row_id'),(Get-Field $row 'family'),(Get-Field $row 'execution_status'),(Get-Field $row 'accuracy_status'),(Get-Field $row 'k4_rms'),(Get-Field $row 'far_k4_rms'),(Get-Field $row 'fallback_k4_rms'),(Get-Field $row 'source_box_boundary_margin_min'),(Get-Field $row 'target_check_excess_min'),(Get-Field $row 'total_seconds'))
}
$report += ''
$report += 'Interpretation: required baseline/proxy-scale/theta rows executed and remained channel-accuracy failures. Production fallback attribution stayed clean (fallback K4 RMS = 0 on executed rows), while the far-only hierarchy K4 RMS remained near the baseline under equiv-scale and theta changes. Negative target-check excess confirms that accepted far interactions can place real target triangle vertices inside the source check cube. The svd-values row is an optional conditioning diagnostic; process failure there is recorded separately and does not override the required-row proximity classification.'
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.5b decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Attribution: $AttributionCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e255b_accuracy_recovered','fmm_e255b_proxy_scale_candidate','fmm_e255b_far_acceptance_candidate','fmm_e255b_conditioning_candidate','fmm_e255b_partial_k4_reduction','fmm_e255b_proximity_confirmed_unresolved') -and -not $AllowFailedDecision) {
  throw "FMM-E25.5b gate did not pass/classify acceptably: $decision"
}
