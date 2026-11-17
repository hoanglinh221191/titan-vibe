param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_smoke'),
    [string[]]$ExtraArgLines = @(),
    [string[]]$ExtraCliArgs = @(),
    [string]$ExpectedSolverMode = '',
    [string]$ExpectedIEStatus = ''
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

function Convert-ToDouble([object]$Value) {
    # TITAN writes C++ "nan"/"inf"; .NET Framework (PowerShell 5.1) parses
    # only "NaN"/"Infinity".
    $text = ([string]$Value).Trim()
    if ($text -match '^[+-]?nan$') { return [double]::NaN }
    if ($text -match '^\+?inf(inity)?$') { return [double]::PositiveInfinity }
    if ($text -match '^-inf(inity)?$') { return [double]::NegativeInfinity }
    return [double]::Parse($text, $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBSA smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Assert-NaNValue([object]$Value, [string]$What) {
    $parsed = Convert-ToDouble $Value
    if (-not [double]::IsNaN($parsed)) {
        throw "Expected nan for $What (interaction entropy not estimated), observed $Value"
    }
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

function Assert-FiniteKey([hashtable]$Map, [string]$Key) {
    if (-not $Map.ContainsKey($Key)) {
        throw "Missing key in mm_pbsa_summary.dat: $Key"
    }
    $value = Convert-ToDouble $Map[$Key]
    if ([double]::IsNaN($value) -or [double]::IsInfinity($value)) {
        throw "Non-finite PBSA value for ${Key}: $value"
    }
    return $value
}

function Get-CsvValue([object]$Row, [string]$Name) {
    $property = $Row.PSObject.Properties[$Name]
    if ($null -eq $property) {
        throw "Missing CSV column: $Name"
    }
    return $property.Value
}

function Invoke-PbsaCase([string]$Name, [string]$Mode) {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'ala2.pdb'
    $topPath = Join-Path $runDir 'ala2.top'
    Write-PbsaTwoChainPdb $pdbPath
    Write-PbsaTwoChainTopology $topPath
    $libPath = Convert-ToTitanPath (Resolve-PbsaLibFfRoot)

    Push-Location $runDir
    try {
        if ($Mode -eq 'arg') {
            $argLines = @(
                "-f=$(Convert-ToTitanPath $pdbPath)",
                '-1pdb',
                '-ff=AMBER99SB',
                "-g=$PbsaTwoChainGroups",
                "-common_topo=$(Convert-ToTitanPath $topPath)",
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
            $argLines += $ExtraArgLines
            $argPath = Join-Path $runDir 'arg.in'
            Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
            & $TitanExe -a $argPath | Out-Null
        } else {
            $argv = @(
                '-f', $pdbPath,
                '-1pdb',
                '-ff', 'AMBER99SB',
                '-g', $PbsaTwoChainGroups,
                '-common_topo', $topPath,
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
            $argv += $ExtraCliArgs
            & $TitanExe @argv | Out-Null
        }
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $Name"
        }
    } finally {
        Pop-Location
    }

    $pbsaSummary = Join-Path $runDir 'mm_pbsa_summary.dat'
    $pbsaCsv = Join-Path $runDir 'mm_pbsa_energy.csv'
    $pbsaDecomp = Join-Path $runDir 'mm_pbsa_decomposition.dat'
    $pbsaByTraj = Join-Path $runDir 'mm_pbsa_by_trajectory.csv'
    $pbsaIeDiag = Join-Path $runDir 'mm_pbsa_ie_diagnostics.csv'
    $pbsaPerf = Join-Path $runDir 'mm_pbsa_performance.dat'
    $pbeSummary = Join-Path $runDir 'pbe_summary.dat'
    $pbeCsv = Join-Path $runDir 'pbe_energy.csv'
    foreach ($required in @($pbsaSummary, $pbsaCsv, $pbsaDecomp, $pbsaByTraj, $pbsaIeDiag, $pbsaPerf, $pbeSummary, $pbeCsv)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing expected PBSA/PBE output for ${Name}: $required"
        }
    }

    $pbsa = Read-SummaryMap $pbsaSummary
    $pbe = Read-SummaryMap $pbeSummary
    if ($pbsa['pbsa_enabled'] -ne 'true') {
        throw "PBSA was not enabled in summary for $Name"
    }
    if ($pbsa['frames_total'] -ne '1') {
        throw "Unexpected PBSA frame count for ${Name}: $($pbsa['frames_total'])"
    }
    if ($pbe['deltaG_solv_enabled'] -ne 'true') {
        throw "PBSA did not auto-enable PBE deltaG_solv for $Name"
    }
    if ($pbsa['endpoint_solves_per_frame'] -ne '6') {
        throw "PBSA summary did not record 6 endpoint solves for ${Name}: $($pbsa['endpoint_solves_per_frame'])"
    }
    if ($pbsa['energy_timeseries'] -ne 'mm_pbsa_energy.csv') {
        throw "PBSA summary did not point to mm_pbsa_energy.csv for $Name"
    }
    if ($ExpectedIEStatus.Length -gt 0 -and $pbsa['interaction_entropy_status'] -ne $ExpectedIEStatus) {
        throw "Expected interaction_entropy_status=$ExpectedIEStatus for ${Name}, observed $($pbsa['interaction_entropy_status'])"
    }

    # Interaction entropy: one frame is below the default -pbsa-ie-min-samples,
    # so the estimate is "insufficient_samples" and -TdS and DeltaG_binding_IE
    # are written as nan (since 1578e499; before that a missing entropy was
    # silently added as 0). With -pbsa-ie=false the entropy term is 0 by
    # definition, and with an "ok" estimate it is finite.
    $ieEnabled = $pbsa['interaction_entropy_enabled'] -eq 'true'
    $ieStatus = [string]$pbsa['interaction_entropy_status']
    $ieAvailable = (-not $ieEnabled) -or ($ieStatus -eq 'ok')
    $noEntropy = Assert-FiniteKey $pbsa 'DeltaG_no_entropy_kJmol_mean'
    if ($ieAvailable) {
        $minusTdS = Assert-FiniteKey $pbsa 'minus_TdS_IE_kJmol'
        $final = Assert-FiniteKey $pbsa 'DeltaG_binding_IE_kJmol'
        $finalKcal = Assert-FiniteKey $pbsa 'DeltaG_binding_IE_kcalmol'
        if (-not $ieEnabled -and $minusTdS -ne 0.0) {
            throw "PBSA with -pbsa-ie=false must report minus_TdS_IE_kJmol=0 for ${Name}, observed $minusTdS"
        }
        $expectedFinalKcal = $final / 4.184
        if ([math]::Abs($expectedFinalKcal - $finalKcal) -gt 1e-5) {
            throw "PBSA summary unit conversion mismatch for ${Name}: $expectedFinalKcal vs $finalKcal"
        }
        $delta = [math]::Abs(($noEntropy + $minusTdS) - $final)
        if ($delta -gt 1e-5) {
            throw "PBSA closure mismatch for ${Name}: no_entropy + IE != final"
        }
    } else {
        foreach ($key in @('minus_TdS_IE_kJmol', 'DeltaG_binding_IE_kJmol', 'DeltaG_binding_IE_kcalmol')) {
            if (-not $pbsa.ContainsKey($key)) {
                throw "Missing key in mm_pbsa_summary.dat: $key"
            }
            Assert-NaNValue $pbsa[$key] "$key ($Name, IE status $ieStatus)"
        }
    }

    foreach ($key in @(
        'E_elec_kJmol_mean',
        'E_vdw_kJmol_mean',
        'E_MM_kJmol_mean',
        'G_polar_PB_kJmol_mean',
        'G_nonpolar_SASA_kJmol_mean',
        'deltaSASA_A2_mean'
    )) {
        [void](Assert-FiniteKey $pbsa $key)
    }
    # MM + polar + nonpolar must add up to the no-entropy binding energy.
    $mmClosure = (Assert-FiniteKey $pbsa 'E_MM_kJmol_mean') +
        (Assert-FiniteKey $pbsa 'G_polar_PB_kJmol_mean') +
        (Assert-FiniteKey $pbsa 'G_nonpolar_SASA_kJmol_mean')
    if ([math]::Abs($mmClosure - $noEntropy) -gt 1e-5) {
        throw "PBSA no-entropy closure mismatch for ${Name}: components=$mmClosure no_entropy=$noEntropy"
    }

    $rows = @(Import-Csv -LiteralPath $pbeCsv)
    if ($rows.Count -ne 1) {
        throw "Expected one PBE energy row for ${Name}, observed $($rows.Count)"
    }
    if ($rows[0].dgsolv_status -ne 'ok') {
        throw "Expected dgsolv_status=ok for ${Name}, observed $($rows[0].dgsolv_status)"
    }

    $pbsaRows = @(Import-Csv -LiteralPath $pbsaCsv)
    if ($pbsaRows.Count -ne 1) {
        throw "Expected one MM/PBSA energy row for ${Name}, observed $($pbsaRows.Count)"
    }
    $row = $pbsaRows[0]
    if ($row.backend -ne 'grid') {
        throw "Expected backend=grid for ${Name}, observed $($row.backend)"
    }
    foreach ($statusName in @('status', 'complex_status', 'receptor_status', 'ligand_status')) {
        $statusValue = Get-CsvValue $row $statusName
        if ($statusValue -ne 'ok') {
            throw "Expected $statusName=ok for ${Name}, observed $statusValue"
        }
    }
    if ($ExpectedSolverMode.Length -gt 0) {
        $solverMode = Get-CsvValue $row 'solver_mode'
        if ($solverMode -ne $ExpectedSolverMode) {
            throw "Expected solver_mode=$ExpectedSolverMode for ${Name}, observed $solverMode"
        }
    }
    foreach ($iterName in @(
        'complex_solvent_iter',
        'complex_vacuum_iter',
        'receptor_solvent_iter',
        'receptor_vacuum_iter',
        'ligand_solvent_iter',
        'ligand_vacuum_iter'
    )) {
        $iterValue = Get-CsvValue $row $iterName
        if ([int]$iterValue -le 0) {
            throw "Expected positive endpoint iteration count $iterName for ${Name}, observed $iterValue"
        }
    }
    # A receptor and a ligand of five atoms each, not one molecule used twice.
    foreach ($countName in @('atom_count_receptor', 'atom_count_ligand')) {
        if ((Get-CsvValue $row $countName) -ne '5') {
            throw "Expected $countName=5 for ${Name}, observed $(Get-CsvValue $row $countName)"
        }
    }
    if ((Get-CsvValue $row 'atom_count_complex') -ne '10') {
        throw "Expected atom_count_complex=10 for ${Name}, observed $(Get-CsvValue $row 'atom_count_complex')"
    }
    $csvNoEntropy =
        (Convert-ToDouble (Get-CsvValue $row 'E_MM_kJmol')) +
        (Convert-ToDouble (Get-CsvValue $row 'DeltaDeltaG_solv_kJmol')) +
        (Convert-ToDouble (Get-CsvValue $row 'DeltaG_nonpolar_kJmol'))
    if ([math]::Abs($csvNoEntropy - $noEntropy) -gt 1e-5) {
        throw "PBSA CSV no-entropy closure mismatch for ${Name}: components=$csvNoEntropy summary=$noEntropy"
    }
    if ((Get-CsvValue $row 'ie_status') -ne $ieStatus) {
        throw "PBSA CSV ie_status for ${Name} is $(Get-CsvValue $row 'ie_status'), summary says $ieStatus"
    }
    if ($ieAvailable) {
        $csvBind = Convert-ToDouble (Get-CsvValue $row 'DeltaG_bind_kJmol')
        $csvClosure = $csvNoEntropy + (Convert-ToDouble (Get-CsvValue $row 'minus_TDeltaS_IE_kJmol'))
        if ([math]::Abs($csvClosure - $csvBind) -gt 1e-5) {
            throw "PBSA CSV closure mismatch for ${Name}: components=$csvClosure bind=$csvBind"
        }
        $csvBindKcal = Convert-ToDouble (Get-CsvValue $row 'DeltaG_bind_kcalmol')
        if ([math]::Abs(($csvBind / 4.184) - $csvBindKcal) -gt 1e-5) {
            throw "PBSA CSV unit conversion mismatch for ${Name}"
        }
    } else {
        foreach ($column in @('minus_TDeltaS_IE_kJmol', 'DeltaG_bind_kJmol', 'DeltaG_bind_kcalmol')) {
            Assert-NaNValue (Get-CsvValue $row $column) "CSV $column ($Name, IE status $ieStatus)"
        }
    }

    $trajRows = @(Import-Csv -LiteralPath $pbsaByTraj)
    if ($trajRows.Count -ne 1) {
        throw "Expected one MM/PBSA trajectory row for ${Name}, observed $($trajRows.Count)"
    }
    if ($trajRows[0].frames -ne '1') {
        throw "Expected one frame in MM/PBSA trajectory summary for ${Name}, observed $($trajRows[0].frames)"
    }
    $ieRows = @(Import-Csv -LiteralPath $pbsaIeDiag)
    if ($ieRows.Count -lt 1) {
        throw "Expected at least one IE diagnostic row for ${Name}"
    }

    if ($pbsa['pbsa_decomp_requested'] -eq 'residue') {
        $resSummary = Join-Path $runDir 'mm_pbsa_residue_summary.csv'
        $resByTraj = Join-Path $runDir 'mm_pbsa_residue_by_trajectory.csv'
        foreach ($required in @($resSummary, $resByTraj)) {
            if (-not (Test-Path -LiteralPath $required)) {
                throw "Missing residue PBSA output for ${Name}: $required"
            }
        }
        $resRows = @(Import-Csv -LiteralPath $resSummary)
        if ($resRows.Count -lt 1) {
            throw "Expected residue decomposition rows for ${Name}"
        }
        $decomp = Read-SummaryMap $pbsaDecomp
        foreach ($closureKey in @(
            'residue_closure_total_max_abs_kJmol',
            'residue_closure_elec_max_abs_kJmol',
            'residue_closure_vdw_max_abs_kJmol',
            'residue_closure_polar_max_abs_kJmol',
            'residue_closure_nonpolar_max_abs_kJmol'
        )) {
            $closure = Assert-FiniteKey $decomp $closureKey
            if ($closure -gt 1e-5) {
                throw "Residue closure too large for ${Name}: $closureKey=$closure"
            }
        }
    }

    Write-Host "PBSA smoke PASS: $Name"
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ForceFieldRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridForceField -ForceFieldRoot $ForceFieldRoot
Invoke-PbsaCase -Name 'arg_equals' -Mode 'arg'
Invoke-PbsaCase -Name 'cli_spacing' -Mode 'cli'
Write-Host "PBSA smoke PASS: all cases"
