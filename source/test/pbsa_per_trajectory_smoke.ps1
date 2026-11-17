param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_per_trajectory_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot

foreach ($case in @('arg_equals', 'cli_spacing')) {
    $summary = Join-Path ([System.IO.Path]::GetFullPath($WorkRoot)) "$case\mm_pbsa_by_trajectory.csv"
    $rows = @(Import-Csv -LiteralPath $summary)
    if ($rows.Count -ne 1) {
        throw "Expected one trajectory row for $case, observed $($rows.Count)"
    }
    foreach ($column in @('DeltaG_bind_noIE_kJmol', 'DeltaG_bind_IE_kJmol', 'IE_status')) {
        if ($null -eq $rows[0].PSObject.Properties[$column]) {
            throw "Missing per-trajectory column $column for $case"
        }
    }
}

Write-Host "PBSA per-trajectory smoke PASS"
