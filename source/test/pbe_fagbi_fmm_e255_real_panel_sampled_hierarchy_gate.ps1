param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [switch]$RerunFailed,
  [switch]$Force,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e255_real_panel_sampled_hierarchy'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Pqr1sh1 = Join-Path $RepoRoot 'test\pbe_inputs_public\1sh1_full_amber_ph7.pqr'
$Pqr1a63 = Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_hierarchy_real_panel_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e255_matrix.csv'
$ChannelCsv = Join-Path $WorkRoot 'fmm_e255_channel.csv'
$CoverageCsv = Join-Path $WorkRoot 'fmm_e255_coverage.csv'
$TreeCsv = Join-Path $WorkRoot 'fmm_e255_tree.csv'
$RuntimeCsv = Join-Path $WorkRoot 'fmm_e255_runtime.csv'
$LiteCsv = Join-Path $WorkRoot 'fmm_e255_comparison_lite.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e255_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e255_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E255_REAL_PANEL_SAMPLED_HIERARCHY_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E25.5 real-panel hierarchy smoke.'
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

function Assert-CanonicalPair([string]$Prefix) {
  $tabiDir = Join-Path $RepoRoot 'src\tabi'
  $allowed = @("$Prefix.h", "$Prefix.cpp")
  $matches = @(Get-ChildItem -LiteralPath $tabiDir -File -Filter "$Prefix*" |
      Select-Object -ExpandProperty Name)
  $extra = @($matches | Where-Object { $allowed -notcontains $_ })
  if ($extra.Count -gt 0) {
    throw "stale source/header variants found for ${Prefix}: $($extra -join ', ')"
  }
}

function Assert-CanonicalSources {
  foreach ($prefix in @(
      'tabi_fmm_kifmm_hierarchy_real_panel',
      'tabi_fmm_kifmm_hierarchy_jcp',
      'tabi_fmm_kifmm_sampled_panel_matvec')) {
    Assert-CanonicalPair $prefix
  }
  $header = Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_hierarchy_real_panel.h'
  if (-not (Test-Path -LiteralPath $header)) {
    throw "missing E25.5 header: $header"
  }
  $text = Get-Content -Raw -LiteralPath $header
  foreach ($needle in @(
      'KifmmHierarchyRealPanelResult',
      'run_kifmm_hierarchy_real_panel_sampled',
      'production_fallback_pairs',
      'nonproduction_fallback_pairs')) {
    if (-not $text.Contains($needle)) {
      throw "expected '$needle' in E25.5 header"
    }
  }
}

function Build-Smoke {
  New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
  Assert-CanonicalSources
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

function Existing-Row([string]$RowId) {
  if (-not (Test-Path -LiteralPath $MatrixCsv)) { return $null }
  $rows = @(Import-Csv -LiteralPath $MatrixCsv)
  $match = @($rows | Where-Object {
      $_.PSObject.Properties['row_id'] -and $_.row_id -eq $RowId
    } | Select-Object -First 1)
  if ($match.Count -eq 0) { return $null }
  return $match[0]
}

function Should-Run([string]$RowId) {
  if ($Force) { return $true }
  $existing = Existing-Row $RowId
  if ($null -eq $existing) { return $true }
  if ($RerunFailed -and ($existing.execution_status -ne 'pass' -or $existing.accuracy_status -ne 'pass')) {
    return $true
  }
  return $false
}

function Run-Profile(
  [string]$Fixture,
  [string]$Sample,
  [int]$SampleCount,
  [int]$Order
) {
  $rowId = ('{0}_{1}_{2}_q{3}' -f $Fixture.ToLowerInvariant(),$Sample,$SampleCount,$Order)
  if (-not (Should-Run $rowId)) {
    Write-Host "resume: $rowId"
    return Existing-Row $rowId
  }
  Write-Host "run: $rowId"
  $args = @(
    "--fixture=$Fixture",
    "--sample=$Sample",
    "--sample-count=$SampleCount",
    "--order=$Order",
    "--mode=profile",
    "--pqr-1sh1=$Pqr1sh1",
    "--pqr-1a63=$Pqr1a63"
  )
  $output = & $ExePath @args 2>&1
  $resultLine = @($output | Where-Object { $_ -like 'E255_RESULT *' } | Select-Object -Last 1)
  $liteLine = @($output | Where-Object { $_ -like 'E255_LITE *' } | Select-Object -Last 1)
  if ($LASTEXITCODE -ne 0 -or $resultLine.Count -ne 1) {
    $row = [pscustomobject][ordered]@{
      row_id = $rowId
      fixture = $Fixture
      sample = $Sample
      sampled_rows = $SampleCount
      order = $Order
      status = 'process_failed'
      execution_status = 'fail'
      accuracy_status = 'not_evaluated'
      accuracy_failure = 'execution'
      aggregate_rms = 'NaN'
      max_channel_rms = 'NaN'
      max_channel_name = 'none'
      k4_rms = 'NaN'
      sampled_dense_valid_pairs = '0'
      production_fallback_pairs = '0'
      production_invalid_pairs = '0'
      nonproduction_fallback_pairs = '0'
      far_panel_pairs = '0'
      direct_panel_pairs = '0'
      coverage = '0'
      expected_coverage = '0'
      m2l_translations = '0'
      l2p_target_count = '0'
      output_hash = '0'
      process_exit_code = $LASTEXITCODE
      process_output = ($output -join ' | ')
    }
    return $row
  }
  $row = Parse-KeyValueLine ([string]$resultLine[0]) 'E255_RESULT '
  $row | Add-Member -NotePropertyName row_id -NotePropertyValue $rowId
  $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $LASTEXITCODE
  $row | Add-Member -NotePropertyName process_output -NotePropertyValue ''
  if ($liteLine.Count -eq 1) {
    $lite = Parse-KeyValueLine ([string]$liteLine[0]) 'E255_LITE '
    foreach ($name in @('status','aggregate_rms','max_channel_rms','total_seconds')) {
      $prop = $lite.PSObject.Properties[$name]
      if ($null -ne $prop) {
        $row | Add-Member -NotePropertyName ("lite_$name") -NotePropertyValue $prop.Value -Force
      }
    }
  }
  $row
}

function Run-InvalidContracts {
  foreach ($mode in @('invalid-duplicates','invalid-range','invalid-nonfinite-unknown','invalid-area','invalid-normal','invalid-triangle')) {
    $output = & $ExePath "--mode=$mode" "--fixture=1SH1" "--sample=uniform_hash" "--sample-count=16" "--order=8" "--pqr-1sh1=$Pqr1sh1" "--pqr-1a63=$Pqr1a63" 2>&1
    $line = @($output | Where-Object { $_ -like 'E255_INVALID *' } | Select-Object -Last 1)
    if ($LASTEXITCODE -ne 0 -or $line.Count -ne 1) {
      throw "invalid contract smoke failed for ${mode}: $($output -join ' | ')"
    }
    $row = Parse-KeyValueLine ([string]$line[0]) 'E255_INVALID '
    if ($row.execution_status -ne 'fail') {
      throw "invalid contract did not fail closed for ${mode}: $($line[0])"
    }
  }
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke
Run-InvalidContracts

$required = @(
  @{ Fixture='1SH1'; Sample='uniform_hash'; Count=128; Order=8 },
  @{ Fixture='1SH1'; Sample='closest_opposing_surface'; Count=128; Order=8 },
  @{ Fixture='1A63'; Sample='high_aspect'; Count=128; Order=8 }
)
if (-not $Quick) {
  $required += @(
    @{ Fixture='1SH1'; Sample='high_aspect'; Count=128; Order=8 },
    @{ Fixture='1SH1'; Sample='inward_normal_groove'; Count=128; Order=8 },
    @{ Fixture='1A63'; Sample='uniform_hash'; Count=128; Order=8 },
    @{ Fixture='1A63'; Sample='closest_opposing_surface'; Count=128; Order=8 },
    @{ Fixture='1SH1'; Sample='uniform_hash'; Count=512; Order=8 },
    @{ Fixture='1A63'; Sample='uniform_hash'; Count=512; Order=8 }
  )
}

$rowsById = @{}
if ((Test-Path -LiteralPath $MatrixCsv) -and -not $Force) {
  foreach ($row in @(Import-Csv -LiteralPath $MatrixCsv)) {
    if ($row.PSObject.Properties['row_id'] -and
        -not [string]::IsNullOrWhiteSpace([string]$row.row_id)) {
      $rowsById[$row.row_id] = $row
    }
  }
}
foreach ($case in $required) {
  $row = Run-Profile $case.Fixture $case.Sample $case.Count $case.Order
  if ($row.PSObject.Properties['row_id'] -and
      -not [string]::IsNullOrWhiteSpace([string]$row.row_id)) {
    $rowsById[$row.row_id] = $row
  }
}

$rows = @($rowsById.GetEnumerator() | ForEach-Object { $_.Value } | Sort-Object row_id)
$requiredRows = @($rows | Where-Object { [int]$_.order -eq 8 })
$requiredPass = $true
foreach ($row in $requiredRows) {
  if ($row.execution_status -ne 'pass' -or $row.accuracy_status -ne 'pass') {
    $requiredPass = $false
  }
}
if ((-not $Quick) -and $requiredPass) {
  foreach ($case in @(
      @{ Fixture='1SH1'; Sample='uniform_hash'; Count=128; Order=6 },
      @{ Fixture='1A63'; Sample='high_aspect'; Count=128; Order=6 })) {
    $row = Run-Profile $case.Fixture $case.Sample $case.Count $case.Order
  $rowsById[$row.row_id] = $row
  }
  $rows = @($rowsById.GetEnumerator() | ForEach-Object { $_.Value } | Sort-Object row_id)
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$rows | Select-Object row_id,aggregate_gate,channel_gate,aggregate_rms,k1_rms,k2_rms,k3_rms,k4_rms,max_channel_rms,max_channel_name,k4_abs_rms,k4_reference_rms,accuracy_status,accuracy_failure | Export-Csv -LiteralPath $ChannelCsv -NoTypeInformation
$rows | Select-Object row_id,source_panels,sampled_rows,sampled_dense_valid_pairs,production_fallback_pairs,production_invalid_pairs,nonproduction_fallback_pairs,far_panel_pairs,direct_panel_pairs,coverage,expected_coverage,far_fraction | Export-Csv -LiteralPath $CoverageCsv -NoTypeInformation
$rows | Select-Object row_id,condition_policy,source_box_policy,m2m_translations,m2l_translations,l2l_translations,l2p_target_count,p2m_leaf_solves,local_node_solves,cache_entries,cache_hits,cache_misses,tree_hash,worklist_hash,policy_hash,panel_geometry_hash,unknown_vector_hash,sampled_target_hash,output_hash | Export-Csv -LiteralPath $TreeCsv -NoTypeInformation
$rows | Select-Object row_id,hierarchy_seconds,production_fallback_seconds,dense_oracle_seconds,total_seconds | Export-Csv -LiteralPath $RuntimeCsv -NoTypeInformation
$rows | Select-Object row_id,lite_status,lite_aggregate_rms,lite_max_channel_rms,lite_total_seconds | Export-Csv -LiteralPath $LiteCsv -NoTypeInformation

$badExecution = @($requiredRows | Where-Object { $_.execution_status -ne 'pass' })
$badCoverage = @($requiredRows | Where-Object {
    [uint64]$_.sampled_dense_valid_pairs -ne [uint64]$_.expected_coverage -or
    [uint64]$_.coverage -ne [uint64]$_.expected_coverage -or
    [uint64]$_.production_invalid_pairs -ne 0 -or
    [uint64]$_.nonproduction_fallback_pairs -ne 0
  })
$badStructure = @($requiredRows | Where-Object { [uint64]$_.far_panel_pairs -eq 0 -or [int]$_.m2l_translations -le 0 -or [uint64]$_.output_hash -eq 0 })
$badAccuracy = @($requiredRows | Where-Object { $_.accuracy_status -ne 'pass' })
$badK4 = @($badAccuracy | Where-Object { $_.max_channel_name -eq 'K4' -or ([double]$_.k4_rms -gt [double]$_.channel_gate) })

$decision = if ($badExecution.Count -gt 0) {
  'fmm_e255_real_panel_execution_gap'
} elseif ($badCoverage.Count -gt 0) {
  'fmm_e255_real_panel_fallback_gap'
} elseif ($badStructure.Count -gt 0) {
  'fmm_e255_real_panel_structural_gap'
} elseif ($badK4.Count -gt 0) {
  'fmm_e255_real_panel_k4_gap'
} elseif ($badAccuracy.Count -gt 0) {
  'fmm_e255_real_panel_aggregate_gap'
} elseif ($Quick) {
  'fmm_e255_real_panel_sampled_hierarchy_quick_pass_full_pending'
} else {
  'fmm_e255_real_panel_sampled_hierarchy_pass'
}

[pscustomobject]@{
  decision = $decision
  mode = if ($Quick) { 'quick' } else { 'full' }
  rows = $rows.Count
  required_rows = $requiredRows.Count
  execution_gap_rows = $badExecution.Count
  coverage_gap_rows = $badCoverage.Count
  structural_gap_rows = $badStructure.Count
  accuracy_gap_rows = $badAccuracy.Count
  k4_gap_rows = $badK4.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
  "decision=$decision",
  "mode=$(if ($Quick) { 'quick' } else { 'full' })",
  "rows=$($rows.Count)",
  "matrix=$MatrixCsv",
  "channel=$ChannelCsv",
  "coverage=$CoverageCsv",
  "tree=$TreeCsv",
  "runtime=$RuntimeCsv",
  "comparison_lite=$LiteCsv"
) | Set-Content -LiteralPath $SummaryDat

$report = @(
  '# TITAN FMM Engine E25.5 Real-Panel Sampled Hierarchy Report',
  '',
  "Decision: $decision",
  '',
  'E25.5 evaluates the strict-octree q8 K1-K4 full hierarchy on sampled rows of real production TABI panels.',
  'The hierarchy contribution is combined only with production analytic_singular direct fallback; the sampled dense oracle also uses the production near-pair kernel with quad_buffer=off.',
  'This is a sampled-row diagnostic gate only. It does not implement full-panel apply, GMRES, preconditioning, production backend policy, or default changes.',
  ''
)
foreach ($row in $rows) {
  $report += ('  {0}: exec={1}, acc={2}, fail={3}, rms={4}, max={5}:{6}, coverage={7}/{8}, m2l={9}, fallback={10}, invalid={11}, total_s={12}' -f
    $row.row_id,$row.execution_status,$row.accuracy_status,$row.accuracy_failure,$row.aggregate_rms,$row.max_channel_name,$row.max_channel_rms,$row.coverage,$row.expected_coverage,$row.m2l_translations,$row.production_fallback_pairs,$row.production_invalid_pairs,$row.total_seconds)
}
$report | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E25.5 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Report: $ReportPath"

if ($decision -notin @('fmm_e255_real_panel_sampled_hierarchy_pass','fmm_e255_real_panel_sampled_hierarchy_quick_pass_full_pending') -and -not $AllowFailedDecision) {
  throw "FMM-E25.5 gate did not pass: $decision"
}
