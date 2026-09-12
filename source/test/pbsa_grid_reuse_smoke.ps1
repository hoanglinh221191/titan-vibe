param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_grid_reuse_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Resolve-PbsaLibFfRoot() {
    return ((Resolve-Path -LiteralPath $ForceFieldRoot).Path.TrimEnd('\') + '\')
}

function Read-SummaryMap([string]$Path) {
    $map = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ([string]::IsNullOrWhiteSpace($line) -or $line.StartsWith('#')) {
            continue
        }
        $parts = $line -split '\s+', 2
        if ($parts.Count -eq 2) {
            $map[$parts[0]] = $parts[1]
        }
    }
    return $map
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBSA grid reuse smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-SmokePdb([string]$Path) {
    $lines = @(
        'CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Require-Key([hashtable]$Map, [string]$Key) {
    if (-not $Map.ContainsKey($Key)) {
        throw "Missing summary/performance key: $Key"
    }
    return $Map[$Key]
}

function Require-NoFile([string]$Path) {
    if (Test-Path -LiteralPath $Path) {
        throw "Unexpected file exists: $Path"
    }
}

function Require-File([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Missing expected file: $Path"
    }
}

function Invoke-PbsaRun([string]$RunDir, [string]$Mode, [string[]]$ArgLines, [string[]]$CliArgs) {
    New-CleanDir $RunDir
    $pdbPath = Join-Path $RunDir 'ala.pdb'
    Write-SmokePdb $pdbPath
    $libPath = Convert-ToTitanPath (Resolve-PbsaLibFfRoot)

    Push-Location $RunDir
    try {
        if ($Mode -eq 'arg') {
            $baseArgLines = @(
                "-f=$(Convert-ToTitanPath $pdbPath)",
                '-1pdb',
                '-ff=AMBER99SB',
                '-g=1&1',
                "-path-a=$libPath",
                "-path-n=$libPath",
                '-pbsa=true',
                '-pbsa-components=complex,receptor,ligand',
                '-pbsa-decomp=molecule',
                '-pbsa-ie=true',
                '-pbe-scale=0.35',
                '-pbe-margin=2.0',
                '-pbe-max-iter=80',
                '-pbe-rms=1e-4',
                '-pbe-dphi=1e-3',
                '-pbe-max-grid=200000',
                '-pbe-map=off',
                '-nt=1'
            )
            $argPath = Join-Path $RunDir 'arg.in'
            Set-Content -LiteralPath $argPath -Value ($baseArgLines + $ArgLines) -Encoding ASCII
            & $TitanExe -a $argPath | Out-Null
        } else {
            $argv = @(
                '-f', $pdbPath,
                '-1pdb',
                '-ff', 'AMBER99SB',
                '-g', '1&1',
                '-path-a', $libPath,
                '-path-n', $libPath,
                '-pbsa', 'true',
                '-pbsa-components', 'complex,receptor,ligand',
                '-pbsa-decomp', 'molecule',
                '-pbsa-ie', 'true',
                '-pbe-scale', '0.35',
                '-pbe-margin', '2.0',
                '-pbe-max-iter', '80',
                '-pbe-rms', '1e-4',
                '-pbe-dphi', '1e-3',
                '-pbe-max-grid', '200000',
                '-pbe-map', 'off',
                '-nt', '1'
            )
            $argv += $CliArgs
            & $TitanExe @argv | Out-Null
        }
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $RunDir"
        }
    } finally {
        Pop-Location
    }

    foreach ($required in @(
        'mm_pbsa_summary.dat',
        'mm_pbsa_energy.csv',
        'mm_pbsa_performance.dat',
        'mm_pbsa_endpoint_timing.csv',
        'pbe_summary.dat',
        'pbe_energy.csv'
    )) {
        Require-File (Join-Path $RunDir $required)
    }

    $summary = Read-SummaryMap (Join-Path $RunDir 'mm_pbsa_summary.dat')
    if ((Require-Key $summary 'pbsa_enabled') -ne 'true') {
        throw "PBSA was not enabled for $RunDir"
    }
    foreach ($key in @('DeltaG_no_entropy_kJmol_mean', 'E_MM_kJmol_mean', 'G_polar_PB_kJmol_mean')) {
        $value = Convert-ToDouble (Require-Key $summary $key)
        if ([double]::IsNaN($value) -or [double]::IsInfinity($value)) {
            throw "Non-finite PBSA summary value ${key}: $value"
        }
    }
}

function Sum-EndpointHits([string]$Path, [string]$Variant) {
    Require-File $Path
    $rows = @(Import-Csv -LiteralPath $Path | Where-Object { $_.phase -eq 'vacuum' -and $_.timing_variant -eq $Variant })
    $surface = 0
    $operator = 0
    foreach ($row in $rows) {
        $surface += [int]$row.reuse_surface_hits
        $operator += [int]$row.reuse_operator_hits
    }
    return [pscustomobject]@{
        Rows = $rows.Count
        Surface = $surface
        Operator = $operator
    }
}

function Invoke-ReuseCase(
    [string]$Name,
    [string[]]$ArgLines,
    [string[]]$CliArgs,
    [string]$ExpectedEffective,
    [bool]$ExpectValidation,
    [bool]$ExpectReuseHits
) {
    $caseRoot = Join-Path $WorkRoot $Name
    Invoke-PbsaRun -RunDir (Join-Path $caseRoot 'arg_equals') -Mode 'arg' -ArgLines $ArgLines -CliArgs $CliArgs
    Invoke-PbsaRun -RunDir (Join-Path $caseRoot 'cli_spacing') -Mode 'cli' -ArgLines $ArgLines -CliArgs $CliArgs

    foreach ($mode in @('arg_equals', 'cli_spacing')) {
        $runDir = Join-Path $caseRoot $mode
        $summary = Read-SummaryMap (Join-Path $runDir 'mm_pbsa_summary.dat')
        $perf = Read-SummaryMap (Join-Path $runDir 'mm_pbsa_performance.dat')
        $effective = Require-Key $summary 'pbsa_grid_reuse_effective'
        if ($effective -ne $ExpectedEffective) {
            throw "Expected pbsa_grid_reuse_effective=$ExpectedEffective for $Name/$mode, observed $effective"
        }
        $perfEffective = Require-Key $perf 'grid_reuse_effective'
        if ($perfEffective -ne $ExpectedEffective) {
            throw "Expected performance grid_reuse_effective=$ExpectedEffective for $Name/$mode, observed $perfEffective"
        }

        $validationCsv = Join-Path $runDir 'mm_pbsa_grid_reuse_validation.csv'
        if ($ExpectValidation) {
            Require-File $validationCsv
            $validationRows = @(Import-Csv -LiteralPath $validationCsv)
            if ($validationRows.Count -ne 3) {
                throw "Expected 3 validation rows for $Name/$mode, observed $($validationRows.Count)"
            }
            foreach ($row in $validationRows) {
                $delta = [math]::Abs((Convert-ToDouble $row.reaction_energy_delta_kJmol))
                if ($delta -gt 1e-4) {
                    throw "PBSA reuse validation delta too large for $Name/$mode/$($row.endpoint): $delta"
                }
            }
        } else {
            Require-NoFile $validationCsv
        }

        $endpointTiming = Join-Path $runDir 'mm_pbsa_endpoint_timing.csv'
        $reuseHits = Sum-EndpointHits $endpointTiming 'reuse'
        if ($ExpectReuseHits) {
            if ($reuseHits.Rows -ne 3) {
                throw "Expected 3 vacuum/reuse endpoint timing rows for $Name/$mode, observed $($reuseHits.Rows)"
            }
            if ($reuseHits.Surface -le 0 -or $reuseHits.Operator -le 0) {
                throw "Expected nonzero reuse hits for $Name/$mode, observed surface=$($reuseHits.Surface) operator=$($reuseHits.Operator)"
            }
        } else {
            if ($reuseHits.Rows -ne 0) {
                throw "Expected no vacuum/reuse endpoint timing rows for $Name/$mode, observed $($reuseHits.Rows)"
            }
            if ($reuseHits.Surface -ne 0 -or $reuseHits.Operator -ne 0) {
                throw "Expected zero reuse hits for $Name/$mode, observed surface=$($reuseHits.Surface) operator=$($reuseHits.Operator)"
            }
        }

        [void](Require-Key $perf 'grid_workspace_reuse_hits')
        [void](Require-Key $perf 'grid_workspace_reallocations')
        [void](Require-Key $perf 'grid_medium_cache_hits')
        [void](Require-Key $perf 'grid_medium_cache_misses')
        [void](Require-Key $perf 'grid_medium_cache_stores')
        if ($ExpectedEffective -ne 'ineligible' -and
            ([int64](Require-Key $perf 'grid_workspace_reuse_hits')) -le 0) {
            throw "Expected workspace reuse hits for $Name/$mode"
        }
        if ($ExpectedEffective -eq 'auto' -and ([int64](Require-Key $perf 'grid_medium_cache_hits')) -le 0) {
            throw "Expected exact medium/operator cache hits for $Name/$mode"
        }
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot

Invoke-ReuseCase `
    -Name 'auto' `
    -ArgLines @('-pbsa-grid-reuse=auto', '-pbsa-perf=detail') `
    -CliArgs @('-pbsa-grid-reuse', 'auto', '-pbsa-perf', 'detail') `
    -ExpectedEffective 'auto' `
    -ExpectValidation $false `
    -ExpectReuseHits $true

Invoke-ReuseCase `
    -Name 'validate' `
    -ArgLines @('-pbsa-grid-reuse=validate', '-pbsa-perf=detail') `
    -CliArgs @('-pbsa-grid-reuse', 'validate', '-pbsa-perf', 'detail') `
    -ExpectedEffective 'validate' `
    -ExpectValidation $true `
    -ExpectReuseHits $true

Invoke-ReuseCase `
    -Name 'ineligible_focus' `
    -ArgLines @('-pbsa-grid-reuse=auto', '-pbsa-perf=detail', '-pbe-focus=auto') `
    -CliArgs @('-pbsa-grid-reuse', 'auto', '-pbsa-perf', 'detail', '-pbe-focus', 'auto') `
    -ExpectedEffective 'ineligible' `
    -ExpectValidation $false `
    -ExpectReuseHits $false

Write-Host "PBSA grid reuse smoke PASS"
