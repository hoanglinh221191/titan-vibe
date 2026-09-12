[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_surface_only_smoke',
    [ValidateRange(1,300)][int]$CaseTimeoutSeconds = 120
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $TitanExe))
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Titan executable not found: $TitanExe"
}

New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$pqrPath = Join-Path $WorkRoot 'surface_only.pqr'
@(
    'ATOM      1  C1  MOL A   1      -1.300   0.000   0.000   0.5000 1.5000',
    'ATOM      2  C2  MOL A   1       1.300   0.000   0.000  -0.5000 1.5000',
    'END'
) | Set-Content -LiteralPath $pqrPath -Encoding ASCII

$argPath = Write-FmmCliArgFile -RunDir $WorkRoot -RepoRoot $RepoRoot `
    -PqrPath $pqrPath -Backend tabi -Mode treecode `
    -FmmMode single-matvec-analytic -SurfacePoints 8 `
    -MaxIter 20 -Restart 10 -Tolerance 1e-4 -Threads 1 `
    -SaltM 0.15 -InternalDielectric 4.0 -ExternalDielectric 80.0 `
    -TemperatureK 298.15 -Preconditioner fmm-nearfield-ilu0 `
    -PreconditionerIluNeighborCap 32 -NearIntegration analytic_singular `
    -PrimarySurface grid-ses -SesGridSpacing 1.25
Add-Content -LiteralPath $argPath -Value '-pbe-tabi-surface-only=true' `
    -Encoding ASCII

$badArgs = @((Get-Content -LiteralPath $argPath) | Where-Object {
    -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
})
if ($badArgs.Count -ne 0) {
    throw "surface-only arg.in contains non equals-syntax options: $($badArgs -join '; ')"
}

$run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $WorkRoot `
    -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
if ($run.timed_out -or $run.exit_code -ne 0) {
    $stderr = (Get-Content (Join-Path $WorkRoot 'stderr.log') `
        -ErrorAction SilentlyContinue) -join ' '
    throw "surface-only solve failed (exit=$($run.exit_code)): $stderr"
}

$profilePath = Join-Path $WorkRoot 'pbe_tabi_profile.csv'
$memoryPath = Join-Path $WorkRoot 'pbe_memory_profile.csv'
if (-not (Test-Path -LiteralPath $profilePath)) {
    throw 'surface-only run did not write pbe_tabi_profile.csv.'
}
if (-not (Test-Path -LiteralPath $memoryPath)) {
    throw 'surface-only run did not write pbe_memory_profile.csv.'
}
$profile = @(Import-Csv -LiteralPath $profilePath) | Select-Object -Last 1
if ([string]$profile.status -ne 'ok' -or
    [string]$profile.convergence_reason -ne 'surface_only' -or
    [string]$profile.reaction_energy_source -ne 'not_computed_surface_only') {
    throw "surface-only status contract failed: status=$($profile.status) reason=$($profile.convergence_reason) energy=$($profile.reaction_energy_source)"
}
if ([int64]$profile.panels -le 0 -or [int]$profile.iterations -ne 0 -or
    [int64]$profile.matvec_calls -ne 0 -or
    [int64]$profile.fmm_local_plan_build_count -ne 0 -or
    [int64]$profile.fmm_local_cache_build_count -ne 0) {
    throw "surface-only unexpectedly entered solver setup: panels=$($profile.panels) iterations=$($profile.iterations) matvec=$($profile.matvec_calls)"
}
if ([string]$profile.matvec_mode -ne 'not_built_surface_only') {
    throw "Unexpected surface-only matvec mode: $($profile.matvec_mode)"
}

$rejectDir = Join-Path $WorkRoot 'reject_grid'
New-Item -ItemType Directory -Path $rejectDir -Force | Out-Null
$rejectArg = Join-Path $rejectDir 'arg.in'
@((Get-Content -LiteralPath $argPath) -replace '^-pbe=tabi$', '-pbe=grid') |
    Set-Content -LiteralPath $rejectArg -Encoding ASCII
$reject = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $rejectDir `
    -ArgPath $rejectArg -CaseTimeoutSeconds 30
$rejectText = ((Get-Content (Join-Path $rejectDir 'stderr.log') `
    -ErrorAction SilentlyContinue) -join ' ') + ' ' +
    ((Get-Content (Join-Path $rejectDir 'stdout.log') `
    -ErrorAction SilentlyContinue) -join ' ')
if ($reject.exit_code -eq 0 -or
    $rejectText -notlike '*surface-only=true requires -pbe=tabi*') {
    throw "Expected grid surface-only rejection; exit=$($reject.exit_code), output=$rejectText"
}

Write-Host "pbe_tabi_surface_only_smoke=pass panels=$($profile.panels)"
