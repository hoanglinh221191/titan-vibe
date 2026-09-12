param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_binding_jpcb'),
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$Pdb2PqrTarget = '',
    [switch]$RefreshFixtures
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$pdb2pqrVersion = '3.7.1'
$pdb2pqrArguments = '--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0'
$radiusFloorAngstrom = 0.6
$fixture = 'barnase_barstar_1brs'
$pdbId = '1BRS'
$assemblyId = '1'
$receptorChain = 'A'
$ligandChain = 'D'
$sourceUrl = 'https://www.rcsb.org/structure/1BRS'
$downloadUrl = 'https://files.rcsb.org/download/1BRS.pdb1'
$assemblyApiUrl = 'https://data.rcsb.org/rest/v1/core/assembly/1BRS/1'
$manifestPath = Join-Path $FixtureRoot 'binding_fixture_manifest.csv'
if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) {
    $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1'
}

function Get-JpcbSha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Get-JpcbAtomLines([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
        $_.StartsWith('ATOM') -or $_.StartsWith('HETATM')
    })
}

function Get-JpcbChain([string]$Line) {
    if ($Line.Length -gt 21) { return [string]$Line[21] }
    $fields = @($Line -split '\s+' | Where-Object { $_ -ne '' })
    if ($fields.Count -ge 5) { return [string]$fields[4] }
    return ''
}

function Get-JpcbCharge([string]$Line) {
    $culture = [Globalization.CultureInfo]::InvariantCulture
    $fields = @($Line -split '\s+' | Where-Object { $_ -ne '' })
    if ($fields.Count -lt 10) { throw "Malformed PQR atom row: $Line" }
    $charge = 0.0
    if (-not [double]::TryParse($fields[-2], [Globalization.NumberStyles]::Float,
            $culture, [ref]$charge)) {
        throw "Malformed PQR charge: $Line"
    }
    return $charge
}

function Set-JpcbContentWithRetry([string]$Path, [object]$Value) {
    for ($attempt = 1; $attempt -le 10; ++$attempt) {
        try {
            Set-Content -LiteralPath $Path -Value $Value -Encoding ASCII
            return
        } catch [IO.IOException] {
            if ($attempt -eq 10) { throw }
            Start-Sleep -Milliseconds (100 * $attempt)
        }
    }
}

function Install-JpcbPdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe)) {
        throw "Bundled Python is missing: $PythonExe"
    }
    New-Item -ItemType Directory -Path $Pdb2PqrTarget -Force | Out-Null
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        & $PythonExe -m pip install --disable-pip-version-check --target $Pdb2PqrTarget `
            ("pdb2pqr=={0}" -f $pdb2pqrVersion)
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install pinned PDB2PQR.' }
    }
}

function Invoke-JpcbPdb2Pqr([string]$InputPath, [string]$OutputPath) {
    $oldPythonPath = $env:PYTHONPATH
    try {
        $env:PYTHONPATH = if ([string]::IsNullOrWhiteSpace($oldPythonPath)) {
            $Pdb2PqrTarget
        } else {
            "$Pdb2PqrTarget;$oldPythonPath"
        }
        & $PythonExe -m pdb2pqr --ff=AMBER --ffout=AMBER --keep-chain `
            --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath
        if ($LASTEXITCODE -ne 0) { throw 'PDB2PQR failed for the 1BRS assembly.' }
    } finally {
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-JpcbPqrRadii([string]$Path, [double]$MinimumRadius) {
    $culture = [Globalization.CultureInfo]::InvariantCulture
    $output = New-Object Collections.Generic.List[string]
    $clamped = 0
    $clampedIds = New-Object Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $Path) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
            $output.Add($line) | Out-Null
            continue
        }
        $fields = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($fields.Count -lt 10) { throw "Malformed PDB2PQR atom row: $line" }
        $radius = 0.0
        if (-not [double]::TryParse($fields[-1], [Globalization.NumberStyles]::Float,
                $culture, [ref]$radius)) {
            throw "Malformed PDB2PQR radius: $line"
        }
        if ($radius -le 0.0) {
            $tokenStart = $line.LastIndexOf($fields[-1], [StringComparison]::Ordinal)
            if ($tokenStart -lt 0) { throw "Unable to locate PQR radius: $line" }
            $line = $line.Substring(0, $tokenStart) +
                $MinimumRadius.ToString('F4', $culture)
            ++$clamped
            $atomId = if ($fields.Count -ge 6) {
                "$($fields[1]):$($fields[2]):$($fields[3]):$($fields[4]):$($fields[5])"
            } else {
                "$($fields[1]):$($fields[2]):$($fields[3])"
            }
            $clampedIds.Add($atomId) | Out-Null
        }
        $output.Add($line) | Out-Null
    }
    Set-JpcbContentWithRetry -Path $Path -Value $output
    return [pscustomobject]@{
        count = $clamped
        atom_ids = ($clampedIds.ToArray() -join ';')
    }
}

function Write-JpcbPreparedAssembly([string]$RawPath, [string]$PreparedPath) {
    $output = New-Object Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $RawPath) {
        if (-not $line.StartsWith('ATOM')) { continue }
        $chain = Get-JpcbChain $line
        if ($chain -eq $receptorChain -or $chain -eq $ligandChain) {
            $output.Add($line) | Out-Null
        }
    }
    if ($output.Count -eq 0) { throw 'No A/D protein atoms found in 1BRS assembly 1.' }
    $output.Add('END') | Out-Null
    Set-JpcbContentWithRetry -Path $PreparedPath -Value $output
}

function Write-JpcbComponentPqr([string]$ComplexPath, [string]$OutputPath,
                                [string]$Chain) {
    $output = New-Object Collections.Generic.List[string]
    foreach ($line in Get-JpcbAtomLines $ComplexPath) {
        if ((Get-JpcbChain $line) -eq $Chain) { $output.Add($line) | Out-Null }
    }
    if ($output.Count -eq 0) { throw "No PQR atoms found for chain $Chain." }
    $output.Add('END') | Out-Null
    Set-JpcbContentWithRetry -Path $OutputPath -Value $output
}

if ($RefreshFixtures) {
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    Install-JpcbPdb2Pqr

    $rawPath = Join-Path $FixtureRoot '1brs_assembly1_source.pdb'
    $preparedPath = Join-Path $FixtureRoot '1brs_assembly1_chains_A_D.pdb'
    $complexPath = Join-Path $FixtureRoot '1brs_complex_amber_ph7.pqr'
    $pdb2pqrLog = [IO.Path]::ChangeExtension($complexPath, '.log')
    $receptorPath = Join-Path $FixtureRoot '1brs_barnase_chain_A_amber_ph7.pqr'
    $ligandPath = Join-Path $FixtureRoot '1brs_barstar_chain_D_amber_ph7.pqr'

    Invoke-WebRequest -UseBasicParsing -Uri $downloadUrl -OutFile $rawPath
    Write-JpcbPreparedAssembly -RawPath $rawPath -PreparedPath $preparedPath
    Invoke-JpcbPdb2Pqr -InputPath $preparedPath -OutputPath $complexPath
    if (-not (Test-Path -LiteralPath $pdb2pqrLog)) {
        throw 'Pinned PDB2PQR run did not produce its preparation log.'
    }
    $radiusNormalization = Normalize-JpcbPqrRadii -Path $complexPath `
        -MinimumRadius $radiusFloorAngstrom
    Write-JpcbComponentPqr -ComplexPath $complexPath -OutputPath $receptorPath `
        -Chain $receptorChain
    Write-JpcbComponentPqr -ComplexPath $complexPath -OutputPath $ligandPath `
        -Chain $ligandChain

    $complexAtoms = @(Get-JpcbAtomLines $complexPath)
    $receptorAtoms = @(Get-JpcbAtomLines $receptorPath)
    $ligandAtoms = @(Get-JpcbAtomLines $ligandPath)
    if ($complexAtoms.Count -ne $receptorAtoms.Count + $ligandAtoms.Count) {
        throw 'Component atom counts do not reproduce the complex PQR.'
    }
    $complexSorted = @($complexAtoms | Sort-Object)
    $componentsSorted = @($receptorAtoms + $ligandAtoms | Sort-Object)
    if (($complexSorted -join "`n") -cne ($componentsSorted -join "`n")) {
        throw 'Component PQR lines are not an exact partition of the complex.'
    }

    $rawAtoms = @(Get-JpcbAtomLines $rawPath)
    $rawHeteroCount = @($rawAtoms | Where-Object { $_.StartsWith('HETATM') }).Count
    $receptorCharge = ($receptorAtoms | ForEach-Object { Get-JpcbCharge $_ } |
        Measure-Object -Sum).Sum
    $ligandCharge = ($ligandAtoms | ForEach-Object { Get-JpcbCharge $_ } |
        Measure-Object -Sum).Sum
    [pscustomobject]@{
        fixture = $fixture
        pdb_id = $pdbId
        biological_name = 'Barnase-barstar complex; barstar C40A/C82A variant'
        assembly_id = $assemblyId
        assembly_stoichiometry = 'A1B1 heterodimer'
        receptor = 'barnase'
        receptor_author_chain = $receptorChain
        ligand = 'barstar C40A/C82A'
        ligand_author_chain = $ligandChain
        source_url = $sourceUrl
        assembly_download_url = $downloadUrl
        assembly_api_url = $assemblyApiUrl
        raw_assembly_pdb = [IO.Path]::GetFileName($rawPath)
        prepared_complex_pdb = [IO.Path]::GetFileName($preparedPath)
        complex_pqr = [IO.Path]::GetFileName($complexPath)
        pdb2pqr_log = [IO.Path]::GetFileName($pdb2pqrLog)
        receptor_pqr = [IO.Path]::GetFileName($receptorPath)
        ligand_pqr = [IO.Path]::GetFileName($ligandPath)
        raw_assembly_sha256 = Get-JpcbSha256 $rawPath
        prepared_complex_pdb_sha256 = Get-JpcbSha256 $preparedPath
        complex_pqr_sha256 = Get-JpcbSha256 $complexPath
        pdb2pqr_log_sha256 = Get-JpcbSha256 $pdb2pqrLog
        receptor_pqr_sha256 = Get-JpcbSha256 $receptorPath
        ligand_pqr_sha256 = Get-JpcbSha256 $ligandPath
        complex_atoms = $complexAtoms.Count
        receptor_atoms = $receptorAtoms.Count
        ligand_atoms = $ligandAtoms.Count
        receptor_charge_e = $receptorCharge
        ligand_charge_e = $ligandCharge
        excluded_hetero_atoms = $rawHeteroCount
        pdb2pqr_version = $pdb2pqrVersion
        pdb2pqr_arguments = $pdb2pqrArguments
        component_policy = 'PDB2PQR once on assembly; split exact output lines by author chain'
        radius_floor_angstrom = $radiusFloorAngstrom
        radius_clamped_atoms = $radiusNormalization.count
        radius_clamped_atom_ids = $radiusNormalization.atom_ids
    } | Export-Csv -LiteralPath $manifestPath -NoTypeInformation
}

if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "Missing cached binding fixture manifest. Run with -RefreshFixtures: $manifestPath"
}
$rows = @(Import-Csv -LiteralPath $manifestPath)
if ($rows.Count -ne 1) { throw "Expected one binding fixture row, found $($rows.Count)." }
$row = $rows[0]
if ($row.pdb_id -ne $pdbId -or $row.assembly_id -ne $assemblyId -or
    $row.receptor_author_chain -ne $receptorChain -or
    $row.ligand_author_chain -ne $ligandChain) {
    throw 'Binding fixture identity/chain policy drift.'
}
$checks = @{
    raw_assembly_pdb = 'raw_assembly_sha256'
    prepared_complex_pdb = 'prepared_complex_pdb_sha256'
    complex_pqr = 'complex_pqr_sha256'
    pdb2pqr_log = 'pdb2pqr_log_sha256'
    receptor_pqr = 'receptor_pqr_sha256'
    ligand_pqr = 'ligand_pqr_sha256'
}
foreach ($name in $checks.Keys) {
    $path = Join-Path $FixtureRoot $row.$name
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing cached fixture file: $path" }
    if ((Get-JpcbSha256 $path) -ne $row.($checks[$name])) {
        throw "Checksum drift for cached fixture file: $path"
    }
}
$complexAtoms = @(Get-JpcbAtomLines (Join-Path $FixtureRoot $row.complex_pqr))
$receptorAtoms = @(Get-JpcbAtomLines (Join-Path $FixtureRoot $row.receptor_pqr))
$ligandAtoms = @(Get-JpcbAtomLines (Join-Path $FixtureRoot $row.ligand_pqr))
if ($complexAtoms.Count -ne [int]$row.complex_atoms -or
    $receptorAtoms.Count -ne [int]$row.receptor_atoms -or
    $ligandAtoms.Count -ne [int]$row.ligand_atoms -or
    $complexAtoms.Count -ne $receptorAtoms.Count + $ligandAtoms.Count) {
    throw 'Binding fixture atom-count drift.'
}
if ($row.pdb2pqr_version -ne $pdb2pqrVersion -or
    $row.pdb2pqr_arguments -ne $pdb2pqrArguments) {
    throw 'Binding fixture PDB2PQR policy drift.'
}

Write-Host "PASS JPCB Barnase-barstar fixture: $manifestPath"
