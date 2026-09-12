[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\pbe_tabi_ses_quartic_witness_replay',
    [string]$Fixture = 'example\pbe_huge\sphere_6co_full.pdb'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite artifact: $WorkRoot"
}
New-Item -ItemType Directory -Path $WorkRoot | Out-Null
$Fixture = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $Fixture))
$expected = 'FA7B83D8D9A0B6CB6A6DBF07A3F93004B76951B51E8562184139060D3752D034'
$actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $Fixture).Hash.ToUpperInvariant()
@("expected_sha256=$expected", "actual_sha256=$actual") |
    Set-Content -LiteralPath (Join-Path $WorkRoot 'input_hash.txt')
if ($actual -ne $expected) {
    throw 'Huge PDB SHA-256 lock failed.'
}

$compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$exe = Join-Path $WorkRoot 'titan_tabi_ses_quartic_witness_replay.exe'
$compileArgs = @(
    '-std=c++23', '-O3', '-march=native', '-DNDEBUG', '-Wall', '-Wextra', '-Werror',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_quartic_witness_replay.cpp'),
    '-o', $exe
)
& $compiler @compileArgs 1> (Join-Path $WorkRoot 'compile.stdout.txt') `
    2> (Join-Path $WorkRoot 'compile.stderr.txt')
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile targeted quartic witness replay.'
}
& $exe $Fixture 1> (Join-Path $WorkRoot 'stdout.txt') `
    2> (Join-Path $WorkRoot 'stderr.txt')
$exitCode = $LASTEXITCODE
Get-Content -LiteralPath (Join-Path $WorkRoot 'stdout.txt')
Get-Content -LiteralPath (Join-Path $WorkRoot 'stderr.txt')
if ($exitCode -ne 0) {
    throw "Targeted witness replay failed with exit $exitCode."
}
