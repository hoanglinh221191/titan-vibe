param(
  [string]$Pqr1sh1 = '',
  [string]$Pqr1a63 = '',
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e22_batched_far_eval'),
  [switch]$Quick,
  [switch]$Full,
  [string[]]$RerunRows = @(),
  [switch]$RerunFailed,
  [switch]$Force,
  [switch]$AllowFailedDecision,
  [int]$CaseTimeoutSeconds = 900,
  [int]$Workers = 4,
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
  throw 'No C++ compiler found for KIFMM E22 gate.'
}

function Build-GateExe([string]$Exe) {
  $compiler = Resolve-Cxx $Cxx
  $args = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_full_panel_matvec_smoke.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_lite.cpp',
    'src\tabi\tabi_fmm_kifmm_real_panel_cloud.cpp',
    'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp',
    'src\tabi\tabi_fmm_kifmm.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    '-o',
    $Exe
  )
  Push-Location $RepoRoot
  try {
    & $compiler @args
    if ($LASTEXITCODE -ne 0) { throw 'failed to build E22 batched far-eval driver' }
  } finally {
    Pop-Location
  }
}

function Assert-E22CanonicalSources {
  function Assert-E22CanonicalPair {
    param([Parameter(Mandatory = $true)][string]$Prefix)
    $tabiDir = Join-Path $RepoRoot 'src\tabi'
    $allowed = @("$Prefix.h", "$Prefix.cpp")
    $matches = @(Get-ChildItem -LiteralPath $tabiDir -File -Filter "$Prefix*" |
      Select-Object -ExpandProperty Name)
    foreach ($required in $allowed) {
      if ($matches -notcontains $required) {
        throw "missing canonical KIFMM source/header for ${Prefix}: $required"
      }
    }
    $extra = @($matches | Where-Object { $allowed -notcontains $_ })
    if ($extra.Count -gt 0) {
      throw "stale duplicate KIFMM source/header for ${Prefix}: $($extra -join ', ')"
    }
  }
  foreach ($prefix in @(
      'tabi_fmm_kifmm_lite',
      'tabi_fmm_kifmm_real_panel_cloud',
      'tabi_fmm_kifmm_full_panel_matvec',
      'tabi_fmm_kifmm_sampled_panel_matvec')) {
    Assert-E22CanonicalPair -Prefix $prefix
  }
  foreach ($path in @(
      'src\tabi\tabi_fmm_kifmm_lite.h',
      'src\tabi\tabi_fmm_kifmm_real_panel_cloud.h',
      'src\tabi\tabi_fmm_kifmm_full_panel_matvec.h',
      'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp',
      'test\titan_fmm_kifmm_full_panel_matvec_smoke.cpp')) {
    $full = Join-Path $RepoRoot $path
    if (-not (Test-Path -LiteralPath $full)) { throw "missing E22 source: $path" }
  }
  $header = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.h')
  foreach ($needle in @('kifmm_full_apply_policy','kifmm_target_tasks','kifmm_parallel_workers','kifmm_channel_build_seconds')) {
    if (-not $header.Contains($needle)) { throw "expected '$needle' in full-panel header" }
  }
  $source = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot 'src\tabi\tabi_fmm_kifmm_full_panel_matvec.cpp')
  foreach ($needle in @('target-node-parallel','parallel_for','kifmm_channel_build_seconds')) {
    if (-not $source.Contains($needle)) { throw "expected '$needle' in full-panel source" }
  }
  $stale = @(Get-ChildItem -Recurse -File (Join-Path $RepoRoot 'src\tabi'),(Join-Path $RepoRoot 'test') |
    Where-Object {
      $_.Extension -in @('.h','.cpp') -and
      $_.Name -match '(?i)(copy|stage|stale|old|bak|tmp)' -and
      $_.Name -match '(?i)kifmm'
    })
  if ($stale.Count -gt 0) {
    throw "stale KIFMM source/header candidates found: $($stale.FullName -join '; ')"
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

function Invoke-E22Process([string[]]$CmdArgs, [string]$LogPath) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $script:exe
  $psi.WorkingDirectory = $RepoRoot
  $psi.UseShellExecute = $false
  $psi.RedirectStandardOutput = $true
  $psi.RedirectStandardError = $true
  $psi.CreateNoWindow = $true
  $quoted = foreach ($arg in $CmdArgs) {
    '"' + ([string]$arg).Replace('\', '\\').Replace('"', '\"') + '"'
  }
  $psi.Arguments = ($quoted -join ' ')
  $process = New-Object System.Diagnostics.Process
  $process.StartInfo = $psi
  [void]$process.Start()
  $stdoutTask = $process.StandardOutput.ReadToEndAsync()
  $stderrTask = $process.StandardError.ReadToEndAsync()
  if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
    try { $process.Kill() } catch {}
    $stdout = $stdoutTask.Result
    $stderr = $stderrTask.Result
    $output = @()
    if (-not [string]::IsNullOrEmpty($stdout)) { $output += ($stdout -split "`r?`n" | Where-Object { $_ -ne '' }) }
    if (-not [string]::IsNullOrEmpty($stderr)) { $output += ($stderr -split "`r?`n" | Where-Object { $_ -ne '' }) }
    $output | Set-Content -LiteralPath $LogPath -Encoding ascii
    return [pscustomobject]@{ timed_out=$true; exit_code=124; output=$output }
  }
  $stdout = $stdoutTask.Result
  $stderr = $stderrTask.Result
  $output = @()
  if (-not [string]::IsNullOrEmpty($stdout)) { $output += ($stdout -split "`r?`n" | Where-Object { $_ -ne '' }) }
  if (-not [string]::IsNullOrEmpty($stderr)) { $output += ($stderr -split "`r?`n" | Where-Object { $_ -ne '' }) }
  $output | Set-Content -LiteralPath $LogPath -Encoding ascii
  return [pscustomobject]@{ timed_out=$false; exit_code=$process.ExitCode; output=$output }
}

function Convert-CaseToRow($row, [hashtable]$m, [string]$policy, [string]$log, [int]$exitCode, [bool]$timedOut) {
  [pscustomobject]@{
    row_id=('{0}_{1}' -f $row.row_id,$policy)
    pair_id=$row.row_id
    mode=$mode
    fixture=$row.fixture
    sample_set=$row.sample
    sample_count=(U64 $m 'sampled_rows')
    order=(U64 $m 'order')
    full_apply_policy=$policy
    status=$(if ($timedOut) { 'timeout' } else { S $m 'status' 'missing' })
    process_exit_code=$exitCode
    source_panels=(U64 $m 'source_panels')
    target_panels=(U64 $m 'target_panels')
    full_pair_coverage=(U64 $m 'full_pair_coverage')
    sampled_direct_valid_pairs=(U64 $m 'sampled_direct_valid_pairs')
    sampled_invalid_pairs=(U64 $m 'sampled_invalid_pairs')
    nonproduction_fallback_pairs=(U64 $m 'nonproduction_fallback_pairs')
    far_panel_pairs=(U64 $m 'far_panel_pairs')
    fallback_panel_pairs=(U64 $m 'fallback_panel_pairs')
    far_fraction=(D $m 'far_fraction')
    aggregate_rms=(D $m 'aggregate_rms')
    max_channel_rms=(D $m 'max_channel_rms')
    max_channel_name=(S $m 'max_channel_name')
    k1=(D $m 'k1')
    k2=(D $m 'k2')
    k3=(D $m 'k3')
    k4=(D $m 'k4')
    max_equiv_residual=(D $m 'max_equiv_residual')
    full_far_eval_seconds=(D $m 'full_far_eval_seconds')
    full_fallback_eval_seconds=(D $m 'full_fallback_eval_seconds')
    kifmm_channel_build_seconds=(D $m 'kifmm_channel_build_seconds')
    kifmm_far_eval_seconds=(D $m 'kifmm_far_eval_seconds')
    kifmm_fallback_eval_seconds=(D $m 'kifmm_fallback_eval_seconds')
    sampled_dense_audit_seconds=(D $m 'sampled_dense_audit_seconds')
    apply_seconds=(D $m 'apply_seconds')
    total_seconds=(D $m 'total_seconds')
    target_tasks=(U64 $m 'kifmm_target_tasks')
    parallel_workers=(U64 $m 'kifmm_parallel_workers')
    worklist_hash=(S $m 'worklist_hash')
    policy_hash=(S $m 'policy_hash')
    panel_geometry_hash=(S $m 'panel_geometry_hash')
    sampled_target_hash=(S $m 'sampled_target_hash')
    full_output_hash_phi=(S $m 'full_output_hash_phi')
    full_output_hash_phin=(S $m 'full_output_hash_phin')
    full_channel_output_hash=(S $m 'full_channel_output_hash')
    log=$log
  }
}

function Convert-ParityToRow($row, [hashtable]$m, [string]$log, [int]$exitCode, [bool]$timedOut) {
  [pscustomobject]@{
    pair_id=$row.row_id
    mode=$mode
    fixture=$row.fixture
    sample_set=$row.sample
    sample_count=(U64 $m 'sampled_rows')
    order=(U64 $m 'order')
    status=$(if ($timedOut) { 'timeout' } elseif ($exitCode -eq 0) { 'pass' } else { 'failed' })
    process_exit_code=$exitCode
    serial_status=(S $m 'serial_status')
    parallel_status=(S $m 'parallel_status')
    serial_apply_seconds=(D $m 'serial_apply_seconds')
    parallel_apply_seconds=(D $m 'parallel_apply_seconds')
    speedup=(D $m 'speedup')
    row_relative_diff=(D $m 'row_relative_diff')
    channel_relative_diff=(D $m 'channel_relative_diff')
    hash_stable=(U64 $m 'hash_stable')
    serial_worklist_hash=(S $m 'serial_worklist_hash')
    parallel_worklist_hash=(S $m 'parallel_worklist_hash')
    serial_policy_hash=(S $m 'serial_policy_hash')
    parallel_policy_hash=(S $m 'parallel_policy_hash')
    serial_far_panel_pairs=(U64 $m 'serial_far_panel_pairs')
    parallel_far_panel_pairs=(U64 $m 'parallel_far_panel_pairs')
    serial_fallback_panel_pairs=(U64 $m 'serial_fallback_panel_pairs')
    parallel_fallback_panel_pairs=(U64 $m 'parallel_fallback_panel_pairs')
    log=$log
  }
}

$quickRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_128_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=128; order=6; speed_gate=1.5 },
  [pscustomobject]@{ row_id='1A63_high_aspect_128_q6'; fixture='1A63'; sample='high_aspect'; sample_count=128; order=6; speed_gate=2.0 }
)
$fullRows = @(
  [pscustomobject]@{ row_id='1SH1_uniform_512_q6'; fixture='1SH1'; sample='uniform_hash'; sample_count=512; order=6; speed_gate=1.5 },
  [pscustomobject]@{ row_id='1SH1_closest_opposing_512_q6'; fixture='1SH1'; sample='closest_opposing_surface'; sample_count=512; order=6; speed_gate=1.5 },
  [pscustomobject]@{ row_id='1A63_uniform_512_q6'; fixture='1A63'; sample='uniform_hash'; sample_count=512; order=6; speed_gate=2.0 },
  [pscustomobject]@{ row_id='1A63_high_aspect_512_q6'; fixture='1A63'; sample='high_aspect'; sample_count=512; order=6; speed_gate=2.0 }
)

$mode = if ($Full -and -not $Quick) { 'full' } else { 'quick' }
$rowsToRun = if ($mode -eq 'full') { $fullRows } else { $quickRows }
$matrixPath = Join-Path $WorkRoot 'fmm_e22_matrix.csv'
$runtimePath = Join-Path $WorkRoot 'fmm_e22_runtime.csv'
$parityPath = Join-Path $WorkRoot 'fmm_e22_parity.csv'
$coveragePath = Join-Path $WorkRoot 'fmm_e22_coverage.csv'
$decisionPath = Join-Path $WorkRoot 'fmm_e22_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fmm_e22_summary.dat'
$reportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E22_BATCHED_FAR_EVAL_REPORT.txt'
$script:exe = Join-Path $WorkRoot 'titan_fmm_e22_batched_far_eval.exe'

Assert-E22CanonicalSources
Build-GateExe -Exe $script:exe

$existingParity = @{}
if (Test-Path -LiteralPath $parityPath) {
  foreach ($old in @(Import-Csv -LiteralPath $parityPath)) { $existingParity[$old.pair_id] = $old }
}
$existingMatrix = @{}
if (Test-Path -LiteralPath $matrixPath) {
  foreach ($old in @(Import-Csv -LiteralPath $matrixPath)) { $existingMatrix[$old.row_id] = $old }
}

$matrix = New-Object System.Collections.Generic.List[object]
$parity = New-Object System.Collections.Generic.List[object]

foreach ($row in $rowsToRun) {
  $serialId = '{0}_serial' -f $row.row_id
  $parallelId = '{0}_target-node-parallel' -f $row.row_id
  $resumeOk = $existingParity.ContainsKey($row.row_id) -and
    $existingMatrix.ContainsKey($serialId) -and
    $existingMatrix.ContainsKey($parallelId)
  if ($resumeOk -and -not $Force -and -not ($RerunRows -contains $row.row_id) -and
      -not ($RerunRows -contains $serialId) -and -not ($RerunRows -contains $parallelId)) {
    $oldParity = $existingParity[$row.row_id]
    $rerun = $RerunFailed -and ($oldParity.status -ne 'pass')
    if (-not $rerun) {
      Write-Host "resume: $($row.row_id)"
      $matrix.Add($existingMatrix[$serialId])
      $matrix.Add($existingMatrix[$parallelId])
      $parity.Add($oldParity)
      continue
    }
  }

  Write-Host "run: $($row.row_id) serial vs target-node-parallel (timeout ${CaseTimeoutSeconds}s)"
  $log = Join-Path $WorkRoot ($row.row_id + '.log')
  $args = @(
    $Pqr1sh1,
    $Pqr1a63,
    '--e22-profile',
    "--fixture=$($row.fixture)",
    "--sample=$($row.sample)",
    "--sample-count=$($row.sample_count)",
    "--order=$($row.order)",
    "--workers=$Workers"
  )
  $run = Invoke-E22Process -CmdArgs $args -LogPath $log
  $caseLines = @($run.output | Where-Object { [string]$_ -match '^KIFMM_FULL_PANEL_CASE ' })
  $parityLines = @($run.output | Where-Object { [string]$_ -match '^KIFMM_E22_PARITY ' })
  if ($run.timed_out -or $caseLines.Count -lt 2 -or $parityLines.Count -lt 1) {
    $empty = [ordered]@{}
    $matrix.Add((Convert-CaseToRow $row $empty 'serial' $log $run.exit_code $run.timed_out))
    $matrix.Add((Convert-CaseToRow $row $empty 'target-node-parallel' $log $run.exit_code $run.timed_out))
    $parity.Add((Convert-ParityToRow $row $empty $log $run.exit_code $run.timed_out))
    continue
  }
  foreach ($line in $caseLines) {
    $m = Parse-KeyValueLine ([string]$line)
    $policy = S $m 'full_apply_policy' 'unknown'
    $matrix.Add((Convert-CaseToRow $row $m $policy $log $run.exit_code $run.timed_out))
  }
  $parity.Add((Convert-ParityToRow $row (Parse-KeyValueLine ([string]$parityLines[0])) $log $run.exit_code $run.timed_out))
}

$matrix.ToArray() | Export-Csv -LiteralPath $matrixPath -NoTypeInformation
$parity.ToArray() | Export-Csv -LiteralPath $parityPath -NoTypeInformation
$matrix.ToArray() | Select-Object row_id,fixture,full_apply_policy,status,apply_seconds,kifmm_channel_build_seconds,kifmm_far_eval_seconds,kifmm_fallback_eval_seconds,target_tasks,parallel_workers |
  Export-Csv -LiteralPath $runtimePath -NoTypeInformation
$matrix.ToArray() | Select-Object row_id,full_pair_coverage,sampled_direct_valid_pairs,sampled_invalid_pairs,nonproduction_fallback_pairs,far_panel_pairs,fallback_panel_pairs,far_fraction,worklist_hash,policy_hash |
  Export-Csv -LiteralPath $coveragePath -NoTypeInformation

$timeouts = @($parity | Where-Object { $_.status -eq 'timeout' })
$coverageGaps = @($matrix | Where-Object {
    $_.status -ne 'pass' -or
    [uint64]$_.sampled_invalid_pairs -ne 0 -or
    [uint64]$_.nonproduction_fallback_pairs -ne 0 -or
    [uint64]$_.full_pair_coverage -eq 0 -or
    [uint64]$_.sampled_direct_valid_pairs -eq 0
  })
$parityGaps = @($parity | Where-Object {
    $_.status -ne 'pass' -or
    [double]$_.row_relative_diff -gt 1.0e-12 -or
    [double]$_.channel_relative_diff -gt 1.0e-12 -or
    [uint64]$_.hash_stable -ne 1
  })
$speedGaps = @()
foreach ($p in $parity) {
  $gate = @($rowsToRun | Where-Object { $_.row_id -eq $p.pair_id } | Select-Object -First 1)
  if ($gate.Count -eq 0) { continue }
  if ([double]$p.speedup -lt [double]$gate.speed_gate) { $speedGaps += $p }
}

if ($timeouts.Count -gt 0) {
  $decision = 'fmm_e22_runtime_limited'
} elseif ($coverageGaps.Count -gt 0) {
  $decision = 'fmm_e22_coverage_gap'
} elseif ($parityGaps.Count -gt 0) {
  $decision = 'fmm_e22_parallel_parity_gap'
} elseif ($speedGaps.Count -gt 0) {
  $decision = 'fmm_e22_accuracy_pass_partial_speedup'
} else {
  if ($mode -eq 'quick') {
    $decision = 'fmm_e22_target_node_parallel_apply_quick_pass'
  } else {
    $decision = 'fmm_e22_target_node_parallel_apply_pass'
  }
}

[pscustomobject]@{
  decision=$decision
  mode=$mode
  rows=$($rowsToRun.Count)
  workers=$Workers
  timeout_rows=$timeouts.Count
  coverage_gaps=$coverageGaps.Count
  parity_gaps=$parityGaps.Count
  speed_gaps=$speedGaps.Count
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation

@(
  "decision: $decision",
  "mode: $mode",
  "workers: $Workers",
  "rows: $($rowsToRun.Count)",
  "matrix: $matrixPath",
  "parity: $parityPath",
  "coverage: $coveragePath"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$report = @(
  '# FMM-E22: Target-Node Parallel Full-Panel KIFMM Apply Gate',
  '',
  "Decision: $decision",
  '',
  'E22 compares the existing serial full-panel KIFMM-lite apply with a target-node-parallel apply. The source-node channel construction remains serial and deterministic; far/fallback target-row evaluation is parallelized by target node so each worker owns disjoint output rows.',
  '',
  'The worklist_hash must match across serial and target-node-parallel rows. The policy_hash is allowed to differ because full_apply_policy is intentionally part of the policy hash. Panel, unknown, sampled-target and coverage hashes remain matched for the same row pair.',
  '',
  'The measured speedup is sublinear but useful: E22 validates target-node parallel row ownership and deterministic parity, not a SIMD/SoA batched far-evaluation kernel and not a production backend.',
  '',
  'The gate preserves the frozen worklist, panel hash, sampled-target hash, production analytic fallback, and sampled dense-direct audit used by E19-E21.',
  '',
  "Rows: $($rowsToRun.Count)",
  "Workers: $Workers",
  '',
  'Quick/full row summary:',
  $($parity | ForEach-Object {
      "  $($_.pair_id): speedup=$($_.speedup), row_diff=$($_.row_relative_diff), channel_diff=$($_.channel_relative_diff), hash_stable=$($_.hash_stable)"
    }),
  '',
  "Matrix: $matrixPath",
  "Parity: $parityPath"
)
$report | Set-Content -LiteralPath $reportPath -Encoding ascii

Write-Host "FMM-E22 decision: $decision"
Write-Host "Matrix: $matrixPath"
Write-Host "Parity: $parityPath"

$passingDecisions = @(
  'fmm_e22_target_node_parallel_apply_quick_pass',
  'fmm_e22_target_node_parallel_apply_pass'
)
if ($passingDecisions -notcontains $decision -and -not $AllowFailedDecision) {
  throw "FMM-E22 gate did not pass: $decision"
}
