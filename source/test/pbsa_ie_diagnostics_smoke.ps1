param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_ie_diagnostics_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$extraArgLines = @(
    '-pbsa-ie-min-samples=1',
    '-pbsa-ie-block=1',
    '-pbsa-ie-filter-sigma=0'
)
$extraCliArgs = @(
    '-pbsa-ie-min-samples', '1',
    '-pbsa-ie-block', '1',
    '-pbsa-ie-filter-sigma', '0'
)

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -ExtraArgLines $extraArgLines `
    -ExtraCliArgs $extraCliArgs `
    -ExpectedIEStatus 'ok'

foreach ($case in @('arg_equals', 'cli_spacing')) {
    $diag = Join-Path ([System.IO.Path]::GetFullPath($WorkRoot)) "$case\mm_pbsa_ie_diagnostics.csv"
    $rows = @(Import-Csv -LiteralPath $diag)
    if (-not ($rows | Where-Object { $_.scope -eq 'global' -and $_.status -eq 'ok' })) {
        throw "Missing global IE ok diagnostic row for $case"
    }
}

Write-Host "PBSA IE diagnostics smoke PASS"
