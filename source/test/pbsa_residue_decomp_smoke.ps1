param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_residue_decomp_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$extraArgLines = @(
    '-pbsa-decomp=residue',
    '-pbsa-residue-frames=all',
    '-pbsa-ie-min-samples=1',
    '-pbsa-ie-block=1',
    '-pbsa-perf=detail'
)
$extraCliArgs = @(
    '-pbsa-decomp', 'residue',
    '-pbsa-residue-frames', 'all',
    '-pbsa-ie-min-samples', '1',
    '-pbsa-ie-block', '1',
    '-pbsa-perf', 'detail'
)

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -ExtraArgLines $extraArgLines `
    -ExtraCliArgs $extraCliArgs `
    -ExpectedIEStatus 'ok'

foreach ($case in @('arg_equals', 'cli_spacing')) {
    $runDir = Join-Path ([System.IO.Path]::GetFullPath($WorkRoot)) $case
    $frameCsv = Join-Path $runDir 'mm_pbsa_residue_frame.csv'
    $perfFrames = Join-Path $runDir 'mm_pbsa_performance_frames.csv'
    foreach ($required in @($frameCsv, $perfFrames)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing residue/detail output for ${case}: $required"
        }
    }
    if (@(Import-Csv -LiteralPath $frameCsv).Count -lt 1) {
        throw "Residue frame CSV is empty for $case"
    }
}

Write-Host "PBSA residue decomposition smoke PASS"
