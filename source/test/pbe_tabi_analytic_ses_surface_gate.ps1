[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\pbe_tabi_analytic_ses_surface_gate'
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
$exe = Join-Path $WorkRoot 'titan_tabi_analytic_ses_surface_smoke.exe'
$args = @(
    '-std=c++23', '-O2', '-Wall', '-Wextra',
    '-DTITAN_TABI_SES_TEST_HOOKS',
    ("-I" + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_bilevel_grid.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_analytic_ses_surface_smoke.cpp'),
    '-o', $exe
)
& $Compiler @args
if ($LASTEXITCODE -ne 0) { throw 'Failed to compile analytic SES surface smoke.' }
& $exe
if ($LASTEXITCODE -ne 0) { throw 'Analytic SES surface smoke failed.' }

