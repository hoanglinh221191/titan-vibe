[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\pbe_tabi_ses_quartic_gate'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$exe = Join-Path $WorkRoot 'titan_tabi_ses_quartic_smoke.exe'
$compileArgs = @(
    '-std=c++23', '-O2', '-Wall', '-Wextra',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_quartic_smoke.cpp'),
    '-o', $exe
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) { throw 'Failed to compile SES quartic smoke.' }
& $exe
if ($LASTEXITCODE -ne 0) { throw 'SES quartic smoke failed.' }
