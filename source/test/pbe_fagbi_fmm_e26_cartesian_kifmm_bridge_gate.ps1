param(
  [switch]$Quick,
  [switch]$AllowFailedDecision,
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e26_cartesian_kifmm_bridge'),
  [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$SmokeSource = Join-Path $RepoRoot 'test\titan_fmm_kifmm_cartesian_bridge_smoke.cpp'
$BuildRoot = Join-Path $WorkRoot 'smoke_build'
$ExePath = Join-Path $BuildRoot 'titan_fmm_kifmm_cartesian_bridge_smoke.exe'
$MatrixCsv = Join-Path $WorkRoot 'fmm_e26_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_e26_decision.csv'
$SummaryDat = Join-Path $WorkRoot 'fmm_e26_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_ENGINE_E26_CARTESIAN_KIFMM_BRIDGE_REPORT.txt'

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
  throw 'No C++ compiler found for FMM-E26 smoke.'
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

function Build-Smoke {
  if (-not (Test-Path -LiteralPath $SmokeSource)) { throw "missing smoke source: $SmokeSource" }
  New-Item -ItemType Directory -Path $BuildRoot -Force | Out-Null
  $Compiler = Resolve-Cxx $Cxx
  $CompileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-Isrc',
    '-Isrc\eigen',
    'test\titan_fmm_kifmm_cartesian_bridge_smoke.cpp',
    'src\tabi\tabi_fmm_kifmm_surface.cpp',
    'src\tabi\tabi_fmm_kifmm_kernels.cpp',
    'src\tabi\tabi_fmm_kifmm_solve.cpp',
    'src\tabi\tabi_fmm_kifmm_fields.cpp',
    'src\tabi\tabi_fmm_kifmm_cartesian_bridge.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_near.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
  )
  Push-Location $RepoRoot
  try {
    & $Compiler @CompileArgs
    if ($LASTEXITCODE -ne 0) { throw 'failed to build titan_fmm_kifmm_cartesian_bridge_smoke' }
  } finally {
    Pop-Location
  }
}

function Write-CsvRows([string]$Path, [object[]]$Rows) {
  if ($Rows.Count -eq 0) { throw "no rows for $Path" }
  $Rows | Export-Csv -LiteralPath $Path -NoTypeInformation
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
Build-Smoke

$stdoutPath = Join-Path $WorkRoot 'e26.stdout.txt'
$stderrPath = Join-Path $WorkRoot 'e26.stderr.txt'
Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue
$proc = Start-Process -FilePath $ExePath -NoNewWindow -Wait -PassThru `
  -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
$output = @()
if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath }
if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath }

$rows = @($output | Where-Object { $_ -like 'E26_ROW *' } | ForEach-Object {
    Parse-KeyValueLine $_ 'E26_ROW'
  })
$decisionLine = @($output | Where-Object { $_ -like 'E26_DECISION *' } | Select-Object -Last 1)
if ($decisionLine.Count -ne 1) {
  throw "FMM-E26 smoke did not emit a decision. ExitCode=$($proc.ExitCode). Output=$($output -join '|')"
}
$decision = Parse-KeyValueLine $decisionLine[0] 'E26_DECISION'

Write-CsvRows $MatrixCsv $rows
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
@(
  "decision=$($decision.decision)",
  "rows=$($decision.rows)",
  "required_rows=$($decision.required_rows)",
  "quick=$([bool]$Quick)"
) | Set-Content -LiteralPath $SummaryDat
@(
  'FMM-E26 Cartesian-KiFMM Hybrid Bridge Report',
  "Decision: $($decision.decision)",
  "Rows: $($decision.rows)",
  'Policy: KiFMM cube-surface order 8 with Cartesian Taylor order 5 for bridge feasibility.',
  'Scope: synthetic single-box bridge feasibility only; no real-panel full hierarchy, GMRES, preconditioner, production backend, manuscript, or default change.',
  'Telemetry note: e2l_coeff_rms=-1 means coefficient-level comparison is not computed in this smoke; scalar E2L gradient uses analytic Cartesian-local differentiation.',
  'Rows:',
  (($rows | ForEach-Object {
      "  $($_.row_id): execution=$($_.execution_status) accuracy=$($_.accuracy_status) aggregate=$($_.aggregate_rms) max_channel=$($_.max_channel_name):$($_.max_channel_rms) e2l_coeff_rms=$($_.e2l_coeff_rms) l2p_gradient_method=$($_.l2p_gradient_method)"
    }) -join "`n")
) | Set-Content -LiteralPath $ReportPath

Write-Host "FMM-E26 decision: $($decision.decision)"
Write-Host "Matrix: $MatrixCsv"
if ($proc.ExitCode -ne 0 -and -not $AllowFailedDecision) {
  throw "FMM-E26 smoke failed with exit code $($proc.ExitCode)"
}
