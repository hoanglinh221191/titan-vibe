param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_units_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$extraArgLines = @('-pbsa-ie=false')
$extraCliArgs = @('-pbsa-ie', 'false')

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -ExtraArgLines $extraArgLines `
    -ExtraCliArgs $extraCliArgs

Write-Host "PBSA units smoke PASS"
