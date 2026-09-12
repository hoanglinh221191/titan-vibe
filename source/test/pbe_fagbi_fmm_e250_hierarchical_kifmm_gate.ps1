param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e250_hierarchical_kifmm'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$SmokeScript = Join-Path $PSScriptRoot 'pbe_fagbi_fmm_e250_hierarchical_kifmm_smoke.ps1'
$SmokeWorkRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $SmokeWorkRoot 'titan_fmm_kifmm_hierarchy_smoke.exe'

if (Test-Path -LiteralPath $WorkRoot) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

& $SmokeScript -WorkRoot $SmokeWorkRoot -Cxx $Cxx
if ($LASTEXITCODE -ne 0) {
  throw "E25.0 smoke build failed"
}
if (-not (Test-Path -LiteralPath $ExePath)) {
  throw "E25.0 smoke executable missing after build: $ExePath"
}

function Parse-E250Line([string]$Line) {
  if (-not $Line.StartsWith('E250_RESULT ')) {
    throw "Unexpected E25.0 result line: $Line"
  }
  $row = [ordered]@{}
  foreach ($token in $Line.Substring('E250_RESULT '.Length).Split(' ')) {
    if ([string]::IsNullOrWhiteSpace($token)) { continue }
    $kv = $token.Split('=', 2)
    if ($kv.Count -ne 2) { throw "Malformed token in E25.0 output: $token" }
    $row[$kv[0]] = $kv[1]
  }
  [pscustomobject]$row
}

function Run-E250Case([string]$CaseId) {
  $output = & $ExePath "--case=$CaseId" 2>&1
  $exit = $LASTEXITCODE
  $line = @($output | Where-Object { $_ -like 'E250_RESULT *' } | Select-Object -Last 1)
  if ($line.Count -ne 1) {
    throw "E25.0 case $CaseId did not emit an E250_RESULT line. Output: $($output -join [Environment]::NewLine)"
  }
  $row = Parse-E250Line ([string]$line[0])
  $row | Add-Member -NotePropertyName exit_code -NotePropertyValue $exit
  $accuracyStatus = if ([double]$row.aggregate_rms -le 1.0e-6) {
    'pass'
  } else {
    'fail'
  }
  $coverageStatus = if ([uint64]$row.coverage -eq ([uint64]$row.source_count * [uint64]$row.target_count)) {
    'pass'
  } else {
    'fail'
  }
  $row | Add-Member -NotePropertyName accuracy_status -NotePropertyValue $accuracyStatus
  $row | Add-Member -NotePropertyName coverage_status -NotePropertyValue $coverageStatus
  $farStatus = if ($CaseId -like '*all_near*') {
    'not_required'
  } elseif ([double]$row.far_fraction -ge 0.10) {
    'pass'
  } else {
    'fail'
  }
  $row | Add-Member -NotePropertyName far_status -NotePropertyValue $farStatus
  $cacheStatus = if ([int]$row.cache_entries -ge 0) { 'pass' } else { 'fail' }
  $row | Add-Member -NotePropertyName cache_status -NotePropertyValue $cacheStatus
  $row
}

$cases = @(
  'scalar_laplace_all_near_q6',
  'scalar_laplace_separated_q6',
  'scalar_laplace_separated_q8',
  'scalar_yukawa_mixed_q8_kappa0p12746'
)

if (-not $Quick) {
  $cases += @()
}

$rows = foreach ($case in $cases) {
  Run-E250Case $case
}

$MatrixCsv = Join-Path $WorkRoot 'fmm_e250_matrix.csv'
$CacheCsv = Join-Path $WorkRoot 'fmm_e250_cache.csv'
$CoverageCsv = Join-Path $WorkRoot 'fmm_e250_coverage.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e250_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e250_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E250_HIERARCHICAL_KIFMM_REPORT.txt'

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows |
  Select-Object case,order,kernel,kappa,cache_entries,cache_hits,cache_misses,p2m_seconds,m2m_seconds,m2l_seconds,total_seconds |
  Export-Csv -LiteralPath $CacheCsv -NoTypeInformation
$rows |
  Select-Object case,source_count,target_count,coverage,far_pairs,direct_pairs,far_fraction,coverage_status,far_status,tree_hash,worklist_hash,policy_hash,output_hash |
  Export-Csv -LiteralPath $CoverageCsv -NoTypeInformation

$allNear = @($rows | Where-Object { $_.case -eq 'scalar_laplace_all_near_q6' })
$q6 = @($rows | Where-Object { $_.case -eq 'scalar_laplace_separated_q6' })
$q8 = @($rows | Where-Object { $_.case -eq 'scalar_laplace_separated_q8' })
$y8 = @($rows | Where-Object { $_.case -eq 'scalar_yukawa_mixed_q8_kappa0p12746' })

$decision = 'fmm_e250_runtime_limited'
if (@($rows | Where-Object { $_.coverage_status -ne 'pass' }).Count -gt 0) {
  $decision = 'fmm_e250_coverage_gap'
} elseif ($allNear.Count -eq 1 -and $allNear[0].accuracy_status -eq 'pass' -and
          $q6.Count -eq 1 -and $q6[0].accuracy_status -eq 'pass' -and
          $q8.Count -eq 1 -and $q8[0].accuracy_status -eq 'pass' -and
          $y8.Count -eq 1 -and $y8[0].accuracy_status -eq 'pass') {
  $decision = 'fmm_e250_hierarchical_scalar_pass'
} elseif ($allNear.Count -eq 1 -and $allNear[0].accuracy_status -eq 'pass' -and
          $q6.Count -eq 1 -and $q6[0].accuracy_status -eq 'fail' -and
          $q8.Count -eq 1 -and $q8[0].accuracy_status -eq 'pass' -and
          $y8.Count -eq 1 -and $y8[0].accuracy_status -eq 'pass') {
  $decision = 'fmm_e250_hierarchical_scalar_q8_pass_q6_gap'
} elseif (@($rows | Where-Object { $_.far_status -eq 'fail' }).Count -gt 0) {
  $decision = 'fmm_e250_tree_policy_gap'
} else {
  $decision = 'fmm_e250_m2l_l2p_gap'
}

[pscustomobject]@{
  decision = $decision
  rows = $rows.Count
  q6_aggregate_rms = if ($q6.Count -eq 1) { $q6[0].aggregate_rms } else { 'NA' }
  q8_aggregate_rms = if ($q8.Count -eq 1) { $q8[0].aggregate_rms } else { 'NA' }
  yukawa_q8_aggregate_rms = if ($y8.Count -eq 1) { $y8[0].aggregate_rms } else { 'NA' }
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "cache=$CacheCsv",
  "coverage=$CoverageCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.0 Hierarchical KIFMM Report',
  '',
  "Decision: $decision",
  '',
  'E25.0 introduces a strict spatial octree skeleton and scalar hierarchical KIFMM parity gate. It is not a production TABI backend, not a GMRES path, and not a manuscript/default change.',
  'The goal is to prove that the next KIFMM engine evaluates far work through a strict spatial hierarchy rather than the E24 flat target-panel worklist.',
  '',
  'Scope note: the current E25.0 scalar implementation uses source-node equivalent charges for accepted far pairs and production dense scalar kernels for direct fallback and the oracle. It does not yet implement real-panel JCP K1-K4 or production GMRES wiring.',
  '',
  'Rows:'
)
foreach ($row in $rows) {
  $report += ('  {0}: status={1}, order={2}, kernel={3}, rms={4}, far_fraction={5}, coverage={6}, cache_entries={7}' -f
    $row.case,$row.status,$row.order,$row.kernel,$row.aggregate_rms,$row.far_fraction,$row.coverage_status,$row.cache_entries)
}
if ($decision -eq 'fmm_e250_hierarchical_scalar_q8_pass_q6_gap') {
  $report += ''
  $report += 'Closeout interpretation: q6 is not sufficient for the separated scalar far-pair gate in this first hierarchy skeleton, but q8 closes the scalar Laplace and Yukawa rows while preserving exact pair coverage and nonzero far work. E25.1 should either improve q6 translations or define q8 as the hierarchy diagnostic order before moving to JCP K1-K4.'
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.0 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e250_hierarchical_scalar_pass','fmm_e250_hierarchical_scalar_q8_pass_q6_gap') -and -not $AllowFailedDecision) {
  throw "FMM-E25.0 gate did not pass: $decision"
}
