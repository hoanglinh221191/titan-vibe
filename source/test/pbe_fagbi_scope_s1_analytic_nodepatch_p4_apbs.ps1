param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_s1_analytic_nodepatch_p4_apbs'),
    [string[]]$Fixtures = @('1vii','1ubq','2ci2'),
    [int[]]$Dimes = @(225,257),
    [int]$TitanTimeoutSeconds = 1800,
    [int]$ApbsTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [string]$ApbsExe = '',
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$driver = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4z_apbs_grid_sweep.ps1'
$plotter = Join-Path $PSScriptRoot 'pbe_fagbi_scope_s1_analytic_nodepatch_p4_figure.py'
if (-not (Test-Path -LiteralPath $driver)) { throw "Missing APBS sweep driver: $driver" }
if (-not (Test-Path -LiteralPath $plotter)) { throw "Missing Figure S1 plotter: $plotter" }

$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$invoke = @{
    TitanExe = $TitanExe
    WorkRoot = $WorkRoot
    Fixtures = if ($Quick) { @('1vii') } else { $Fixtures }
    Dimes = $Dimes
    FineScales = @(0.55)
    FmmTimeoutSeconds = $TitanTimeoutSeconds
    ApbsTimeoutSeconds = $ApbsTimeoutSeconds
    Threads = $Threads
    TitanPrimarySurface = 'analytic-ses'
    TitanDiscretization = 'node-patch'
    TitanSesGridSpacing = 1.25
    TitanProtocolTag = 'analytic_nodepatch_p4_s1_v1'
    ApbsExe = $ApbsExe
    PlanOnly = $PlanOnly
    Force = $Force
    RerunFailed = $RerunFailed
    AllowFailedDecision = $true
}

& $driver @invoke
$driverOk = $?
$driverExit = if ($driverOk) { 0 } else { 1 }
if ($PlanOnly) { exit $driverExit }

$matrix = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_matrix.csv'
if (-not (Test-Path -LiteralPath $matrix)) {
    throw "Figure S1 matrix was not written: $matrix"
}

$python = Get-Command python.exe -ErrorAction SilentlyContinue
if ($null -eq $python) { $python = Get-Command python -ErrorAction SilentlyContinue }
if ($null -eq $python) { throw 'Python is required to generate Figure S1.' }

$png = Join-Path $WorkRoot 'fagbi_scope_s1_analytic_nodepatch_p4_apbs.png'
$svg = Join-Path $WorkRoot 'fagbi_scope_s1_analytic_nodepatch_p4_apbs.svg'
$summary = Join-Path $WorkRoot 'fagbi_scope_s1_analytic_nodepatch_p4_summary.csv'
& $python.Source $plotter --matrix $matrix --png $png --svg $svg --summary $summary
if ($LASTEXITCODE -ne 0) { throw "Figure S1 plotter failed with exit code $LASTEXITCODE" }

$rows = @(Import-Csv -LiteralPath $matrix)
$requiredFixtures = if ($Quick) { @('1vii') } else { @($Fixtures) }
$titanRows = @($rows | Where-Object {
    $_.backend -eq 'titan_guarded_p4' -and $_.fixture -in $requiredFixtures
})
$protocolPass = @($titanRows | Where-Object {
    $_.status -eq 'ok' -and
    ([string]$_.primary_surface).Replace('-','_') -eq 'titan_analytic_patch_ray_ses' -and
    ([string]$_.discretization).Replace('-','_') -eq 'node_patch'
})
$decision = if ($protocolPass.Count -eq $requiredFixtures.Count) {
    's1_analytic_nodepatch_p4_ready'
} elseif (@($titanRows | Where-Object { $_.failure_class -eq 'timeout' }).Count -gt 0) {
    's1_runtime_limited'
} else {
    's1_protocol_or_accuracy_gap'
}

@(
    "decision=$decision"
    "matrix=$matrix"
    "figure_png=$png"
    "figure_svg=$svg"
    "summary=$summary"
    'surface=analytic-ses'
    'discretization=node-patch'
    'taylor_order=4'
    'ses_grid_spacing_A=1.25'
    'threads=4'
) | Set-Content -LiteralPath (Join-Path $WorkRoot 'fagbi_scope_s1_summary.dat') -Encoding ASCII

Write-Host "Figure S1 decision: $decision"
Write-Host "Matrix: $matrix"
Write-Host "Figure: $png"
Write-Host 'Resume: rerun unchanged; add -RerunFailed for failed/timeout rows or -Force for all rows.'
if (-not $AllowFailedDecision -and $decision -ne 's1_analytic_nodepatch_p4_ready') { exit 1 }
exit $driverExit
