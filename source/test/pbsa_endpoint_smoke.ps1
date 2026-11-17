param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_endpoint_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot

Write-Host "PBSA endpoint smoke PASS"
