param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\ptopo_kernel_smoke'),
    [string]$Python = 'python',
    [switch]$RequireOracle,
    [switch]$BuildDriver
)

# Periodic topology kernel (src/ptopo): C++ suite, external SymPy integer
# oracle, and optionally a compile/link check of tools/ptopo_xtc_audit.cpp.
# Compiled with -mavx2 -mfma like both TITAN profiles: the AVX2 contact
# kernel and the scalar path must agree bit for bit, which the suite checks.

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
$flags = @('-std=c++23', '-O2', '-mavx2', '-mfma', '-Wall', '-Wextra', '-Isrc')
$kernel = @(
    'src\ptopo\ptopo_integer_lattice.cpp',
    'src\ptopo\ptopo_periodic_graph.cpp',
    'src\ptopo\ptopo_filtration.cpp',
    'src\ptopo\ptopo_geometry.cpp',
    'src\ptopo\ptopo_gromacs_topology.cpp',
    'src\ptopo\ptopo_graph_builder.cpp',
    'src\titan_executor.cpp'
)

Push-Location $repo
try {
    $testExe = Join-Path $WorkRoot 'ptopo_kernel_test.exe'
    & $compilerPath @flags 'tests\ptopo\ptopo_kernel_test.cpp' @kernel -o $testExe
    if ($LASTEXITCODE -ne 0) { throw 'ptopo_kernel_smoke: failed to compile the C++ suite' }
    & $testExe (Join-Path $WorkRoot 'fixtures')
    if ($LASTEXITCODE -ne 0) { throw "ptopo_kernel_smoke: C++ suite failed ($LASTEXITCODE)" }

    $cliExe = Join-Path $WorkRoot 'ptopo_lattice_cli.exe'
    & $compilerPath @flags 'tests\ptopo\ptopo_lattice_cli.cpp' 'src\ptopo\ptopo_integer_lattice.cpp' -o $cliExe
    if ($LASTEXITCODE -ne 0) { throw 'ptopo_kernel_smoke: failed to compile the lattice CLI' }

    & $Python -c 'import sympy' 2>$null
    if ($LASTEXITCODE -eq 0) {
        & $Python 'tests\ptopo\sympy_oracle.py' $cliExe
        if ($LASTEXITCODE -ne 0) { throw 'ptopo_kernel_smoke: SymPy oracle failed' }
    } elseif ($RequireOracle) {
        throw "ptopo_kernel_smoke: SymPy not importable by '$Python'"
    } else {
        Write-Warning "ptopo_kernel_smoke: SymPy not importable by '$Python'; oracle skipped"
    }

    if ($BuildDriver) {
        $driverExe = Join-Path $WorkRoot 'ptopo_xtc_audit.exe'
        & $compilerPath @flags '-O3' '-Isrc\gmxtraj\include' 'tools\ptopo_xtc_audit.cpp' @kernel `
            'src\ptopo\ptopo_xtc_reader.cpp' 'src\gmxtraj\xdrfile.cpp' 'src\gmxtraj\xdrfile_xtc.cpp' -o $driverExe
        if ($LASTEXITCODE -ne 0) { throw 'ptopo_kernel_smoke: failed to build ptopo_xtc_audit' }
        Write-Host "built $driverExe"
    }
    Write-Host 'PASS ptopo_kernel_smoke'
}
finally { Pop-Location }
