param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_e1_kifmm_parser'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

function Resolve-Cxx([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        if (-not (Test-Path -LiteralPath $Requested)) {
            throw "Requested C++ compiler not found: $Requested"
        }
        return $Requested
    }
    foreach ($candidate in @(
            'C:\msys64\ucrt64\bin\clang++.exe',
            'C:\msys64\mingw64\bin\clang++.exe',
            'C:\msys64\ucrt64\bin\g++.exe',
            'C:\msys64\mingw64\bin\g++.exe')) {
        if (Test-Path -LiteralPath $candidate) { return $candidate }
    }
    $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    $cmd = Get-Command g++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    throw 'No C++ compiler found for KIFMM parser smoke.'
}

function Invoke-TitanValidation(
    [string]$Name,
    [string[]]$ArgumentList,
    [string]$ExpectedPattern
) {
    $caseDir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $caseDir | Out-Null
    $stdout = Join-Path $caseDir 'stdout.log'
    $stderr = Join-Path $caseDir 'stderr.log'
    $previousErrorActionPreference = $ErrorActionPreference
    Push-Location $caseDir
    try {
        $ErrorActionPreference = 'Continue'
        & $TitanExe @ArgumentList 1> $stdout 2> $stderr
        $exitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $previousErrorActionPreference
        Pop-Location
    }
    $combined = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
    if ($exitCode -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if ($combined -notmatch $ExpectedPattern) {
        throw "$Name error did not match '$ExpectedPattern': $combined"
    }
}

function Assert-ArgShape {
    param([string[]]$Lines)
    foreach ($line in $Lines) {
        if ($line.StartsWith('-') -and -not $line.Contains('=')) {
            throw "Generated arg.in option is not -option=value: $line"
        }
    }
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$Compiler = Resolve-Cxx $Cxx
$ExePath = Join-Path $WorkRoot 'titan_fmm_kifmm_parser_smoke.exe'
Push-Location $RepoRoot
try {
    & $Compiler -std=c++23 -O2 -mavx2 -mfma -Isrc `
        test\titan_fmm_kifmm_parser_smoke.cpp `
        src\parse_arguments.cpp `
        -o $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build titan_fmm_kifmm_parser_smoke'
    }
    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'titan_fmm_kifmm_parser_smoke failed'
    }
} finally {
    Pop-Location
}

$argLines = @(
    '-pbe=tabi',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-fmm-engine=kifmm-lite',
    '-pbe-tabi-kifmm-order=6',
    '-pbe-tabi-kifmm-equiv-scale=1.05',
    '-pbe-tabi-kifmm-check-scale=1.60'
)
Assert-ArgShape $argLines
$argDir = Join-Path $WorkRoot 'arg_shape'
New-Item -ItemType Directory -Force -Path $argDir | Out-Null
$argLines | Set-Content -LiteralPath (Join-Path $argDir 'arg.in') -Encoding ASCII

Invoke-TitanValidation -Name 'kifmm_lite_skeleton' -ArgumentList @(
    '-pbe', 'tabi',
    '-pbe-tabi-sum', 'treecode',
    '-pbe-tabi-fmm-engine', 'kifmm-lite'
) -ExpectedPattern 'KIFMM engine skeleton only.*FMM-E1'
Invoke-TitanValidation -Name 'kifmm_local_skeleton' -ArgumentList @(
    '-pbe', 'tabi',
    '-pbe-tabi-sum', 'treecode',
    '-pbe-tabi-fmm-engine', 'kifmm-local'
) -ExpectedPattern 'KIFMM engine skeleton only.*FMM-E1'
Invoke-TitanValidation -Name 'invalid_odd_order' -ArgumentList @(
    '-pbe', 'tabi',
    '-pbe-tabi-sum', 'treecode',
    '-pbe-tabi-kifmm-order', '3'
) -ExpectedPattern 'pbe-bie-kifmm-order.*even.*2.*16'
Invoke-TitanValidation -Name 'invalid_scale_order' -ArgumentList @(
    '-pbe', 'tabi',
    '-pbe-tabi-sum', 'treecode',
    '-pbe-tabi-kifmm-equiv-scale', '1.7',
    '-pbe-tabi-kifmm-check-scale', '1.6'
) -ExpectedPattern 'equiv-scale.*check-scale'

Write-Host 'PASS pbe_fagbi_fmm_e1_kifmm_parser_smoke'
