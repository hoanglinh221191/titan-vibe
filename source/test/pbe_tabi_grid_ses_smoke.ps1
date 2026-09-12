[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_grid_ses_smoke',
    [ValidateRange(0.25,2.0)][double]$GridSpacing = 0.75,
    [ValidateRange(1,900)][int]$CaseTimeoutSeconds = 180
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $TitanExe))
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Titan executable not found: $TitanExe"
}

$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$pqrPath = Join-Path $WorkRoot 'narrow_gap.pqr'
@(
    'ATOM      1  C1  MOL A   1      -1.300   0.000   0.000   0.5000 1.0000',
    'ATOM      2  C2  MOL A   1       1.300   0.000   0.000  -0.5000 1.0000',
    'END'
) | Set-Content -LiteralPath $pqrPath -Encoding ASCII

$argPath = Write-FmmCliArgFile -RunDir $WorkRoot -RepoRoot $RepoRoot `
    -PqrPath $pqrPath -Backend tabi -Mode direct -FmmMode off `
    -SurfacePoints 8 -MaxIter 4 -Restart 8 -Tolerance 1e-4 -Threads 1 `
    -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 `
    -TemperatureK 298.15 -Preconditioner diagonal -NearIntegration analytic_singular `
    -PrimarySurface grid-ses -SesGridSpacing $GridSpacing

$badArgs = @((Get-Content -LiteralPath $argPath) | Where-Object {
    -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
})
if ($badArgs.Count -ne 0) {
    throw "Grid SES arg.in contains non equals-syntax options: $($badArgs -join '; ')"
}

$run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $WorkRoot `
    -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
if ($run.timed_out -or $run.exit_code -ne 0) {
    $stderr = Join-Path $WorkRoot 'stderr.log'
    throw "Grid SES solve failed (exit=$($run.exit_code)): $((Get-Content $stderr -ErrorAction SilentlyContinue) -join ' ')"
}

$qualityPath = Join-Path $WorkRoot 'pbe_tabi_surface_quality.csv'
if (-not (Test-Path -LiteralPath $qualityPath)) {
    throw 'Grid SES solve did not write pbe_tabi_surface_quality.csv.'
}
$quality = @(Import-Csv -LiteralPath $qualityPath) | Select-Object -Last 1
$profilePath = Join-Path $WorkRoot 'pbe_tabi_profile.csv'
if (-not (Test-Path -LiteralPath $profilePath)) {
    throw 'Grid SES solve did not write pbe_tabi_profile.csv.'
}
$profile = @(Import-Csv -LiteralPath $profilePath) | Select-Object -Last 1
if ([string]$profile.surface_generator -ne 'titan_grid_watertight_ses') {
    throw "Unexpected surface generator: $($profile.surface_generator)"
}
if ([string]$quality.calderon_trace_space_status -ne 'pass' -or
    [int64]$quality.calderon_trace_boundary_edges -ne 0 -or
    [int64]$quality.calderon_trace_nonmanifold_edges -ne 0) {
    throw "Grid SES topology gate failed: status=$($quality.calderon_trace_space_status) boundary=$($quality.calderon_trace_boundary_edges) nonmanifold=$($quality.calderon_trace_nonmanifold_edges)"
}
if ([int64]$quality.panel_count -le 0) {
    throw 'Grid SES produced no panels.'
}

function Assert-GridSesRejected([string]$Name, [string]$Pattern,
                                [string]$Replacement, [string]$ExpectedText) {
    $rejectDir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $rejectDir | Out-Null
    $rejectArg = Join-Path $rejectDir 'arg.in'
    @((Get-Content -LiteralPath $argPath) -replace $Pattern, $Replacement) |
        Set-Content -LiteralPath $rejectArg -Encoding ASCII
    $reject = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $rejectDir `
        -ArgPath $rejectArg -CaseTimeoutSeconds 30
    $stderr = (Get-Content (Join-Path $rejectDir 'stderr.log') -ErrorAction SilentlyContinue) -join ' '
    if ($reject.exit_code -eq 0 -or $stderr -notlike "*$ExpectedText*") {
        throw "Expected $Name rejection containing '$ExpectedText'; exit=$($reject.exit_code), stderr=$stderr"
    }
}

Assert-GridSesRejected -Name 'invalid_spacing' `
    -Pattern '^-pbe-tabi-ses-grid-spacing=.*$' `
    -Replacement '-pbe-tabi-ses-grid-spacing=0.1' `
    -ExpectedText 'must be between 0.25 and 4.0 Angstrom'
Assert-GridSesRejected -Name 'unsupported_grid_backend' `
    -Pattern '^-pbe=tabi$' -Replacement '-pbe=grid' `
    -ExpectedText 'requires triangle TABI'

Write-Host "pbe_tabi_grid_ses_smoke=pass panels=$($quality.panel_count) boundary_edges=0"
