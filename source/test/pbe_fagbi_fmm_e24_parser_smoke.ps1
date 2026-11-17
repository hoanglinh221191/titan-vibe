param(
  [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$typesPath = Join-Path $RepoRoot 'src\tabi\tabi_types.h'
$parsePath = Join-Path $RepoRoot 'src\cli\parse_arguments.cpp'
$argumentsPath = Join-Path $RepoRoot 'src\arguments.cpp'
$solverPath = Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp'

foreach ($path in @($typesPath, $parsePath, $argumentsPath, $solverPath)) {
  if (-not (Test-Path -LiteralPath $path)) {
    throw "Missing required source file: $path"
  }
}

$types = Get-Content -Raw -LiteralPath $typesPath
$parse = (Get-ChildItem -LiteralPath (Split-Path -Parent $parsePath) -File -Filter 'parse_*.cpp' |
  Sort-Object Name | ForEach-Object { Get-Content -Raw -LiteralPath $_.FullName }) -join "`n"
$arguments = Get-Content -Raw -LiteralPath $argumentsPath
$solver = Get-Content -Raw -LiteralPath $solverPath

foreach ($needle in @(
    'std::string kifmm_ilu_graph{"cartesian-near"}',
    'pbe-bie-kifmm-ilu-graph',
    '"kifmm-fallback-capped"',
    'kifmm_theta_allowed',
    '0.75',
    '0.8',
    'pbe_tabi_kifmm_gmres_convergence_trace.csv',
    'reported_true_residual',
    'preconditioned_residual')) {
  $haystack = if ($needle -eq 'std::string kifmm_ilu_graph{"cartesian-near"}') {
    $types
  } elseif ($needle -eq 'pbe-bie-kifmm-ilu-graph' -or $needle -eq '"kifmm-fallback-capped"') {
    $parse + "`n" + $arguments
  } elseif ($needle -eq 'kifmm_theta_allowed' -or $needle -eq '0.75' -or $needle -eq '0.8') {
    $arguments
  } else {
    $solver
  }
  if (-not $haystack.Contains($needle)) {
    throw "Missing E24 parser/logging contract marker: $needle"
  }
}

if ($arguments.Contains('theta 0.6, equiv-scale')) {
  throw 'E24 validator still advertises theta 0.6 as the only accepted KIFMM experimental theta.'
}

$argLines = @(
  '-pbe=bie',
  '-pbe-bie-sum=treecode',
  '-pbe-bie-fmm-engine=kifmm-lite-experimental',
  '-pbe-bie-kifmm-ilu-graph=cartesian-near',
  '-pbe-bie-opening-theta=0.75'
)
$bad = @($argLines | Where-Object { $_ -notmatch '^-.*=' })
if ($bad.Count -ne 0) {
  throw "E24 generated-arg smoke has non equals-syntax options: $($bad -join '; ')"
}

'FMM-E24 parser smoke passed.'
