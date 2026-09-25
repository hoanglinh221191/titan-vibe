[CmdletBinding()]
param(
    [string]$Compiler = 'clang++',
    [string]$WorkRoot = 'test\tmp\bie_cartesian_full_fmm_hardening',
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$exe = Join-Path $WorkRoot 'titan_bie_cartesian_full_fmm_hardening_smoke.exe'
$sources = @(
    'test\titan_bie_cartesian_full_fmm_hardening_smoke.cpp',
    'src\bie\bie_cartesian_fmm_solver_adapter.cpp',
    'src\bie\fmm\cartesian\bie_cartesian_fmm_operator.cpp',
    'src\bie\fmm\cartesian\bie_cartesian_fmm_plan.cpp',
    'src\bie\fmm\cartesian\bie_cartesian_fmm_translations.cpp',
    'src\tabi\tabi_fmm_kifmm_hierarchy.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp',
    'src\titan_executor.cpp'
)
$sourcePaths = @($sources | ForEach-Object { Join-Path $RepoRoot $_ })
$needsBuild = $Force -or -not (Test-Path -LiteralPath $exe)
if (-not $needsBuild) {
    $exeTime = (Get-Item -LiteralPath $exe).LastWriteTimeUtc
    $needsBuild = @($sourcePaths | Where-Object {
        (Get-Item -LiteralPath $_).LastWriteTimeUtc -gt $exeTime
    }).Count -gt 0
}
if ($needsBuild) {
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "Cartesian full-FMM hardening compiler was not found: $Compiler"
    }
    $arguments = @(
        '-std=c++23', '-O2', '-mavx2', '-mfma',
        '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',
        '-DTITAN_CARTESIAN_FMM_TESTING=1',
        "-I$RepoRoot\src"
    ) + $sourcePaths + @('-o', $exe)
    Write-Host "build: Cartesian full-FMM hardening smoke with $Compiler"
    & $Compiler @arguments
    if ($LASTEXITCODE -ne 0) { throw 'Cartesian full-FMM hardening smoke build failed.' }
}

& $exe
if ($LASTEXITCODE -ne 0) { throw 'Cartesian full-FMM hardening smoke failed.' }
Write-Host 'PASS pbe_bie_cartesian_full_fmm_hardening_smoke'
