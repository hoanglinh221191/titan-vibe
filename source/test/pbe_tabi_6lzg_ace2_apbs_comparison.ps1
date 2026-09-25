param(
    [string]$TitanExe = 'C:\tmp\elec_vdw-main-sync\titan.exe',
    [string]$ApbsExe = '',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_6lzg_ace2_rbd_binding_apbs',
    [string]$FixtureRoot = 'test\pbe_inputs_public_6lzg',
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$Pdb2PqrTarget = '',
    [ValidateSet('titan','apbs')]
    [string[]]$Backends = @('titan','apbs'),
    [ValidateRange(1,64)][int]$Threads = 4,
    [ValidateSet('fmm-nearfield-ilu0','fmm-nearfield-ilu0-coarse')]
    [string]$Preconditioner = 'fmm-nearfield-ilu0',
    [ValidateSet(32,64,128)][int]$IluNeighborCap = 32,
    [ValidateRange(16,8192)][int]$IluMemoryBudgetMb = 512,
    [ValidateRange(4,64)][int]$CoarseAggregateCap = 16,
    [ValidateSet('on','off')]
    [string]$NodePatchOperatorAudit = 'on',
    [ValidateRange(60,14400)][int]$TitanTimeoutSeconds = 3600,
    [ValidateRange(60,14400)][int]$ApbsTimeoutSeconds = 1800,
    [ValidateRange(0.25,1.0)][double]$TargetApbsGridSpacingA = 0.50,
    [ValidateRange(4.0,40.0)][double]$ApbsFinePaddingA = 12.0,
    [ValidateRange(4.0,60.0)][double]$ApbsCoarsePaddingA = 16.0,
    [ValidateRange(65,1025)][int]$ApbsMaxGridDimension = 385,
    [string]$ReportOutput = 'TITAN_TABI_6LZG_ACE2_RBD_BINDING_APBS_REPORT.txt',
    [switch]$RefreshFixture,
    [switch]$PrepareOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$Value) { $Value.ToString('G15', $Inv) }
function D($Value, [double]$Fallback = [double]::NaN) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
        return $Fallback
    }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float, $Inv, [ref]$parsed)) {
        return $parsed
    }
    return $Fallback
}
function I($Value, [int64]$Fallback = 0) {
    $parsed = 0L
    if ([int64]::TryParse([string]$Value, [ref]$parsed)) { return $parsed }
    return $Fallback
}
function Is-Finite([double]$Value) {
    return -not [double]::IsNaN($Value) -and
        -not [double]::IsInfinity($Value)
}
function Sha([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}
function Sha256-Text([string]$Text) {
    $sha256 = [Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [Text.Encoding]::UTF8.GetBytes($Text)
        return ([BitConverter]::ToString($sha256.ComputeHash($bytes)) -replace '-','').ToLowerInvariant()
    } finally {
        $sha256.Dispose()
    }
}
function Get-ProfileValue($Profile, [string]$Name, $Default = '') {
    if ($null -eq $Profile) { return $Default }
    $property = $Profile.PSObject.Properties[$Name]
    if ($null -eq $property) { return $Default }
    return $property.Value
}
function Normalize-List([string[]]$Values) {
    return @($Values | ForEach-Object { @([string]$_ -split ',') } |
        ForEach-Object { $_.Trim().ToLowerInvariant() } | Where-Object { $_ } |
        Sort-Object -Unique)
}
function Get-PqrAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
            $_ -match '^(ATOM|HETATM)'
        }).Count
}

$Backends = @(Normalize-List $Backends)
if ($Backends.Count -eq 0) { throw 'At least one backend is required.' }
if (@($Backends | Where-Object { $_ -notin @('titan','apbs') }).Count -gt 0) {
    throw 'Backends must contain only titan and/or apbs.'
}
if ($Threads -ne 4) {
    throw 'The manuscript comparison is locked to four threads for TITAN and APBS.'
}
if ((($ApbsMaxGridDimension - 1) % 32) -ne 0) {
    throw 'ApbsMaxGridDimension must be of the form 32*n+1.'
}

if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) {
    $FixtureRoot = Join-Path $RepoRoot $FixtureRoot
}
if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) {
    $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1'
}
New-Item -ItemType Directory -Force -Path $WorkRoot,$FixtureRoot | Out-Null

$ManifestPath = Join-Path $FixtureRoot '6lzg_binding_fixture_manifest.csv'
$RawPdbPath = Join-Path $FixtureRoot '6lzg_source.pdb'
$PreparedPdbPath = Join-Path $FixtureRoot '6lzg_model1_chains_A_B_complex.pdb'
$ComplexPqrPath = Join-Path $FixtureRoot '6lzg_complex_charmm_ph7.pqr'
$ReceptorPqrPath = Join-Path $FixtureRoot '6lzg_ace2_chain_A_charmm_ph7.pqr'
$LigandPqrPath = Join-Path $FixtureRoot '6lzg_rbd_chain_B_charmm_ph7.pqr'
$Pdb2PqrLog = Join-Path $FixtureRoot '6lzg_complex_pdb2pqr.log'
$MatrixPath = Join-Path $WorkRoot 'pbe_tabi_6lzg_ace2_rbd_binding_matrix.csv'
$ComponentPath = Join-Path $WorkRoot 'pbe_tabi_6lzg_ace2_rbd_component_energies.csv'
$ComparisonPath = Join-Path $WorkRoot 'pbe_tabi_6lzg_ace2_rbd_thermocycle.csv'
$DecisionPath = Join-Path $WorkRoot 'pbe_tabi_6lzg_ace2_rbd_binding_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'pbe_tabi_6lzg_ace2_rbd_binding_summary.dat'
$ReportPath = if ([IO.Path]::IsPathRooted($ReportOutput)) {
    $ReportOutput
} else {
    Join-Path $RepoRoot $ReportOutput
}

function Write-ModelChainProteinPdb(
    [string]$InputPath,
    [string]$OutputPath,
    [string[]]$Chains = @('A','B'),
    [int]$Model = 1
) {
    $lines = @(Get-Content -LiteralPath $InputPath)
    $hasModels = @($lines | Where-Object { $_ -match '^MODEL\s' }).Count -gt 0
    $activeModel = if ($hasModels) { 0 } else { 1 }
    $rows = [Collections.Generic.List[object]]::new()
    $ordinal = 0
    foreach ($line0 in $lines) {
        if ($line0 -match '^MODEL\s*(\d+)') {
            $activeModel = [int]$Matches[1]
            continue
        }
        if ($line0 -match '^ENDMDL') {
            if ($activeModel -eq $Model) { break }
            $activeModel = 0
            continue
        }
        if ($activeModel -ne $Model -or -not $line0.StartsWith('ATOM')) {
            continue
        }
        $line = $line0.PadRight(80)
        if ([string]$line[21] -notin $Chains) { continue }
        $atomName = $line.Substring(12,4).Trim()
        $element = $line.Substring(76,2).Trim().ToUpperInvariant()
        $nameWithoutDigits = $atomName -replace '^\d+',''
        if ($element -in @('H','D') -or $nameWithoutDigits -match '^[HD]') {
            continue
        }
        $occupancy = D $line.Substring(54,6) 0.0
        $rows.Add([pscustomobject]@{
                line = $line
                residue = $line.Substring(21,6)
                atom_name = $atomName
                alt_loc = $line.Substring(16,1)
                occupancy = $occupancy
                ordinal = $ordinal
            }) | Out-Null
        ++$ordinal
    }
    if ($rows.Count -eq 0) {
        throw "No model $Model chain(s) $($Chains -join ',') protein atoms found in $InputPath"
    }
    $selected = [Collections.Generic.List[string]]::new()
    foreach ($group in @($rows | Group-Object { "$($_.residue)|$($_.atom_name)" })) {
        $choice = @($group.Group | Sort-Object `
            @{Expression='occupancy';Descending=$true},
            @{Expression={ if ([string]::IsNullOrWhiteSpace($_.alt_loc)) { 2 }
                    elseif ($_.alt_loc -eq 'A') { 1 } else { 0 } };Descending=$true},
            ordinal | Select-Object -First 1)[0]
        $fixed = $choice.line.Remove(16,1).Insert(16,' ')
        $selected.Add($fixed.TrimEnd()) | Out-Null
    }
    $ordered = @($selected | Sort-Object {
            $line = [string]$_
            [int]($line.Substring(6,5).Trim())
        })
    $ordered + @('TER','END') | Set-Content -LiteralPath $OutputPath -Encoding ASCII
}

function Get-PqrChain([string]$Line) {
    if ($Line.Length -gt 21 -and -not [string]::IsNullOrWhiteSpace([string]$Line[21])) {
        return [string]$Line[21]
    }
    $fields = @($Line -split '\s+' | Where-Object { $_ })
    if ($fields.Count -ge 6) { return [string]$fields[4] }
    return ''
}

function Write-ComponentPqr([string]$InputPath, [string]$OutputPath,
    [string]$Chain) {
    $rows = [Collections.Generic.List[string]]::new()
    foreach ($line in Get-Content -LiteralPath $InputPath) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        if ((Get-PqrChain $line) -eq $Chain) { $rows.Add($line) | Out-Null }
    }
    if ($rows.Count -eq 0) { throw "No PQR atoms found for chain $Chain." }
    $rows.Add('END') | Out-Null
    $rows | Set-Content -LiteralPath $OutputPath -Encoding ASCII
}

function Install-Pdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe -PathType Leaf)) {
        throw "Missing bundled Python: $PythonExe"
    }
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker -PathType Leaf)) {
        New-Item -ItemType Directory -Force -Path $Pdb2PqrTarget | Out-Null
        & $PythonExe -m pip install --disable-pip-version-check `
            --target $Pdb2PqrTarget 'pdb2pqr==3.7.1'
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install PDB2PQR 3.7.1.' }
    }
}

function Invoke-Pdb2Pqr([string]$InputPath, [string]$OutputPath,
    [string]$LogPath) {
    $oldPythonPath = $env:PYTHONPATH
    $oldPreference = $ErrorActionPreference
    try {
        $env:PYTHONPATH = if ($oldPythonPath) {
            "$Pdb2PqrTarget;$oldPythonPath"
        } else { $Pdb2PqrTarget }
        $ErrorActionPreference = 'Continue'
        $output = & $PythonExe -m pdb2pqr --ff=CHARMM --ffout=CHARMM `
            --keep-chain --titration-state-method=propka --with-ph=7.0 `
            $InputPath $OutputPath 2>&1
        $exitCode = $LASTEXITCODE
        $output | Set-Content -LiteralPath $LogPath -Encoding UTF8
        if ($exitCode -ne 0 -or -not (Test-Path -LiteralPath $OutputPath)) {
            throw "PDB2PQR failed with exit code $exitCode; see $LogPath"
        }
    } finally {
        $ErrorActionPreference = $oldPreference
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-PqrRadii([string]$Path, [double]$Floor = 0.6) {
    $output = [Collections.Generic.List[string]]::new()
    $clamped = 0
    foreach ($line0 in Get-Content -LiteralPath $Path) {
        $line = $line0
        if ($line -notmatch '^(ATOM|HETATM)') {
            $output.Add($line) | Out-Null
            continue
        }
        $fields = @($line -split '\s+' | Where-Object { $_ })
        if ($fields.Count -lt 10) { throw "Malformed PQR row: $line" }
        $radius = D $fields[-1]
        if (-not (Is-Finite $radius)) { throw "Malformed PQR radius: $line" }
        if ($radius -le 0.0) {
            $start = $line.LastIndexOf($fields[-1], [StringComparison]::Ordinal)
            $line = $line.Substring(0,$start) + $Floor.ToString('F4',$Inv)
            ++$clamped
        }
        $output.Add($line) | Out-Null
    }
    $output | Set-Content -LiteralPath $Path -Encoding ASCII
    return $clamped
}

function Ensure-Fixture {
    $valid = $false
    if (-not $RefreshFixture -and (Test-Path -LiteralPath $ManifestPath)) {
        $old = @(Import-Csv -LiteralPath $ManifestPath) | Select-Object -Last 1
        if ($null -ne $old -and
            (Test-Path -LiteralPath $RawPdbPath) -and
            (Test-Path -LiteralPath $PreparedPdbPath) -and
            (Test-Path -LiteralPath $ComplexPqrPath) -and
            (Test-Path -LiteralPath $ReceptorPqrPath) -and
            (Test-Path -LiteralPath $LigandPqrPath) -and
            [string]$old.raw_pdb_sha256 -eq (Sha $RawPdbPath) -and
            [string]$old.prepared_pdb_sha256 -eq (Sha $PreparedPdbPath) -and
            [string]$old.complex_pqr_sha256 -eq (Sha $ComplexPqrPath) -and
            [string]$old.receptor_pqr_sha256 -eq (Sha $ReceptorPqrPath) -and
            [string]$old.ligand_pqr_sha256 -eq (Sha $LigandPqrPath) -and
            (I $old.complex_atom_count) -eq (Get-PqrAtomCount $ComplexPqrPath) -and
            (I $old.receptor_atom_count) -eq (Get-PqrAtomCount $ReceptorPqrPath) -and
            (I $old.ligand_atom_count) -eq (Get-PqrAtomCount $LigandPqrPath)) {
            $valid = $true
        }
    }
    if ($valid) {
        Write-Host 'fixture resume: 6LZG ACE2(A)-RBD(B) binding partition'
        return @(Import-Csv -LiteralPath $ManifestPath)[-1]
    }

    Write-Host 'fixture prepare: 6LZG ACE2(A)-RBD(B) binding partition'
    Install-Pdb2Pqr
    Invoke-WebRequest -UseBasicParsing `
        -Uri 'https://files.rcsb.org/download/6LZG.pdb' -OutFile $RawPdbPath
    Write-ModelChainProteinPdb -InputPath $RawPdbPath `
        -OutputPath $PreparedPdbPath -Chains @('A','B') -Model 1
    Invoke-Pdb2Pqr -InputPath $PreparedPdbPath -OutputPath $ComplexPqrPath `
        -LogPath $Pdb2PqrLog
    $clamped = Normalize-PqrRadii -Path $ComplexPqrPath -Floor 0.6
    Write-ComponentPqr -InputPath $ComplexPqrPath -OutputPath $ReceptorPqrPath -Chain A
    Write-ComponentPqr -InputPath $ComplexPqrPath -OutputPath $LigandPqrPath -Chain B
    $complexAtoms = @(Get-Content -LiteralPath $ComplexPqrPath | Where-Object { $_ -match '^(ATOM|HETATM)' })
    $receptorAtoms = @(Get-Content -LiteralPath $ReceptorPqrPath | Where-Object { $_ -match '^(ATOM|HETATM)' })
    $ligandAtoms = @(Get-Content -LiteralPath $LigandPqrPath | Where-Object { $_ -match '^(ATOM|HETATM)' })
    if ($complexAtoms.Count -lt 3000 -or $receptorAtoms.Count -lt 2000 -or
        $ligandAtoms.Count -lt 500) {
        throw "Unexpected 6LZG atom counts: complex=$($complexAtoms.Count) receptor=$($receptorAtoms.Count) ligand=$($ligandAtoms.Count)"
    }
    if ($complexAtoms.Count -ne $receptorAtoms.Count + $ligandAtoms.Count) {
        throw 'ACE2 and RBD component atom counts do not reproduce the complex PQR.'
    }
    if ((@($complexAtoms | Sort-Object) -join "`n") -cne
        (@($receptorAtoms + $ligandAtoms | Sort-Object) -join "`n")) {
        throw 'ACE2 and RBD PQR lines are not an exact partition of the complex.'
    }
    if (@(Get-Content -LiteralPath $ComplexPqrPath | Where-Object {
                $_ -match '^HETATM'
            }).Count -gt 0) {
        throw 'Unsupported HETATM rows remain in the prepared 6LZG complex PQR.'
    }
    $row = [pscustomobject][ordered]@{
        fixture = '6lzg_ace2_rbd_binding'
        pdb_id = '6LZG'
        biological_name = 'Human ACE2 receptor bound to SARS-CoV-2 spike receptor-binding domain'
        model = 1
        receptor_chain = 'A'
        ligand_chain = 'B'
        assembly_stoichiometry = 'A1B1 heterodimer'
        source_url = 'https://www.rcsb.org/structure/6LZG'
        download_url = 'https://files.rcsb.org/download/6LZG.pdb'
        raw_pdb = Split-Path -Leaf $RawPdbPath
        prepared_pdb = Split-Path -Leaf $PreparedPdbPath
        complex_pqr = Split-Path -Leaf $ComplexPqrPath
        receptor_pqr = Split-Path -Leaf $ReceptorPqrPath
        ligand_pqr = Split-Path -Leaf $LigandPqrPath
        raw_pdb_sha256 = Sha $RawPdbPath
        prepared_pdb_sha256 = Sha $PreparedPdbPath
        complex_pqr_sha256 = Sha $ComplexPqrPath
        receptor_pqr_sha256 = Sha $ReceptorPqrPath
        ligand_pqr_sha256 = Sha $LigandPqrPath
        complex_atom_count = $complexAtoms.Count
        receptor_atom_count = $receptorAtoms.Count
        ligand_atom_count = $ligandAtoms.Count
        pdb2pqr_version = '3.7.1'
        force_field = 'CHARMM'
        ph = '7.0'
        preparation = '--ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0; ATOM-only model 1 chains A+B protonated once as a complex; exact PQR partition; radius floor 0.6 A'
        radius_clamped_atoms = $clamped
    }
    @($row) | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation
    return $row
}

function Resolve-Executable([string]$Requested, [string[]]$Fallbacks) {
    foreach ($candidate in @($Requested) + $Fallbacks) {
        if ([string]::IsNullOrWhiteSpace($candidate)) { continue }
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
        $command = Get-Command $candidate -ErrorAction SilentlyContinue
        if ($command) { return $command.Source }
    }
    return ''
}

function Invoke-Case([string]$Exe, [string]$RunDir, [string]$Arguments,
    [int]$TimeoutSeconds, [hashtable]$Environment = @{}) {
    $psi = [Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $Exe
    $psi.WorkingDirectory = $RunDir
    $psi.Arguments = $Arguments
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    foreach ($key in $Environment.Keys) {
        $psi.EnvironmentVariables[$key] = [string]$Environment[$key]
    }
    $process = [Diagnostics.Process]::new()
    $process.StartInfo = $psi
    $watch = [Diagnostics.Stopwatch]::StartNew()
    $peak = 0L
    $started = $false
    try {
        [void]$process.Start()
        $started = $true
        $stdout = $process.StandardOutput.ReadToEndAsync()
        $stderr = $process.StandardError.ReadToEndAsync()
        $deadline = [DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
        while (-not $process.HasExited -and [DateTime]::UtcNow -lt $deadline) {
            try {
                $process.Refresh()
                $peak = [math]::Max($peak, [int64]$process.WorkingSet64)
                $peak = [math]::Max($peak, [int64]$process.PeakWorkingSet64)
            } catch {}
            Start-Sleep -Milliseconds 100
        }
        $timedOut = -not $process.HasExited
        if ($timedOut) {
            try { $process.Kill() } catch {}
            $process.WaitForExit()
        }
        try {
            $process.Refresh()
            $peak = [math]::Max($peak, [int64]$process.PeakWorkingSet64)
        } catch {}
        $stdout.Wait(); $stderr.Wait()
        Set-Content -LiteralPath (Join-Path $RunDir 'stdout.log') `
            -Value $stdout.Result -Encoding UTF8
        Set-Content -LiteralPath (Join-Path $RunDir 'stderr.log') `
            -Value $stderr.Result -Encoding UTF8
        return [pscustomobject]@{
            exit_code = if ($timedOut) { 124 } else { $process.ExitCode }
            timed_out = $timedOut
            wall_seconds = $watch.Elapsed.TotalSeconds
            peak_rss_mb = $peak / 1MB
        }
    } finally {
        if ($started -and -not $process.HasExited) {
            try { $process.Kill() } catch {}
        }
        $watch.Stop()
        $process.Dispose()
    }
}

function Get-PqrAtoms([string]$Path) {
    $atoms = [Collections.Generic.List[object]]::new()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 10) { continue }
        $atoms.Add([pscustomobject]@{
                x = D $parts[-5]
                y = D $parts[-4]
                z = D $parts[-3]
                radius = D $parts[-1] 1.5
            }) | Out-Null
    }
    return $atoms.ToArray()
}

function Get-ApbsGeometry([string]$Path) {
    $atoms = @(Get-PqrAtoms $Path)
    if ($atoms.Count -eq 0) { throw "No PQR atoms: $Path" }
    $xmin=1e300; $xmax=-1e300; $ymin=1e300; $ymax=-1e300
    $zmin=1e300; $zmax=-1e300
    foreach ($atom in $atoms) {
        $radius = [math]::Max(0.0, [double]$atom.radius)
        $xmin=[math]::Min($xmin,$atom.x-$radius)
        $xmax=[math]::Max($xmax,$atom.x+$radius)
        $ymin=[math]::Min($ymin,$atom.y-$radius)
        $ymax=[math]::Max($ymax,$atom.y+$radius)
        $zmin=[math]::Min($zmin,$atom.z-$radius)
        $zmax=[math]::Max($zmax,$atom.z+$radius)
    }
    $fgx=($xmax-$xmin)+2*$ApbsFinePaddingA
    $fgy=($ymax-$ymin)+2*$ApbsFinePaddingA
    $fgz=($zmax-$zmin)+2*$ApbsFinePaddingA
    return [pscustomobject]@{
        fgx=$fgx; fgy=$fgy; fgz=$fgz
        cgx=$fgx+2*$ApbsCoarsePaddingA
        cgy=$fgy+2*$ApbsCoarsePaddingA
        cgz=$fgz+2*$ApbsCoarsePaddingA
    }
}

function Get-ApbsDime([double]$LengthA) {
    $required = [int][math]::Ceiling($LengthA / $TargetApbsGridSpacingA) + 1
    $dimension = [int](32 * [math]::Ceiling(
            ([math]::Max(65,$required)-1)/32.0) + 1)
    if ($dimension -gt $ApbsMaxGridDimension) {
        throw "Required APBS grid dimension $dimension exceeds cap $ApbsMaxGridDimension"
    }
    return $dimension
}

function Write-ApbsInput([string]$LocalPqr, [string]$InputPath) {
    $geometry = Get-ApbsGeometry $LocalPqr
    $nx=Get-ApbsDime $geometry.fgx
    $ny=Get-ApbsDime $geometry.fgy
    $nz=Get-ApbsDime $geometry.fgz
    $pqrName = Split-Path -Leaf $LocalPqr
    $common = @(
        "    dime $nx $ny $nz",
        "    cglen $(F $geometry.cgx) $(F $geometry.cgy) $(F $geometry.cgz)",
        "    fglen $(F $geometry.fgx) $(F $geometry.fgy) $(F $geometry.fgz)",
        '    cgcent mol 1','    fgcent mol 1','    mol 1','    lpbe',
        '    bcfl sdh','    pdie 1.0','    sdie 80.0','    srfm smol',
        '    chgm spl2','    sdens 10.0','    srad 1.4','    swin 0.3',
        '    temp 298.15','    calcenergy total','    calcforce no'
    )
    $reference = @($common)
    $reference[[array]::IndexOf($reference,'    sdie 80.0')] = '    sdie 1.0'
    $lines = @(
        'read',"    mol pqr $pqrName",'end',
        'elec name solvated','    mg-auto'
    ) + $common[0..7] + @(
        '    ion charge 1 conc 0.15 radius 2.0',
        '    ion charge -1 conc 0.15 radius 2.0'
    ) + $common[8..($common.Count-1)] + @(
        'end','elec name reference','    mg-auto'
    ) + $reference + @(
        'end','print elecEnergy solvated - reference end','quit'
    )
    $lines | Set-Content -LiteralPath $InputPath -Encoding ASCII
    return [pscustomobject]@{
        nx=$nx; ny=$ny; nz=$nz
        fgx=$geometry.fgx; fgy=$geometry.fgy; fgz=$geometry.fgz
        actual_spacing_max=[math]::Max($geometry.fgx/($nx-1),
            [math]::Max($geometry.fgy/($ny-1),$geometry.fgz/($nz-1)))
    }
}

function Read-ApbsEnergy([string]$RunDir) {
    $text = ((Get-Content (Join-Path $RunDir 'stdout.log') -Raw `
            -ErrorAction SilentlyContinue) + "`n" +
        (Get-Content (Join-Path $RunDir 'stderr.log') -Raw `
            -ErrorAction SilentlyContinue))
    $matches = [regex]::Matches($text,
        'Global net ELEC energy\s*=\s*([-+0-9.eE]+)',
        [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    if ($matches.Count -eq 0) { return [double]::NaN }
    return D $matches[$matches.Count-1].Groups[1].Value
}

function New-Row([string]$CaseId, [string]$Backend, [string]$Component,
    [string]$Protocol, [string]$ExeSha, [string]$PqrSha,
    [string]$RunDir, [int64]$Atoms) {
    return [pscustomobject][ordered]@{
        case_id = $CaseId
        fixture = '6lzg_ace2_rbd_binding'
        pdb_id = '6LZG'
        biological_name = 'Human ACE2 receptor bound to SARS-CoV-2 spike RBD'
        component = $Component
        backend = $Backend
        protocol = $Protocol
        executable_sha256 = $ExeSha
        pqr_sha256 = $PqrSha
        threads = $Threads
        status = 'planned'
        failure = ''
        exit_code = ''
        timed_out = ''
        atom_count = $Atoms
        panels = ''
        energy_kj_mol = ''
        energy_kcal_mol = ''
        gpol_complex_kj_mol = ''
        gpol_receptor_kj_mol = ''
        gpol_ligand_kj_mol = ''
        delta_gpol_bind_kj_mol = ''
        e_coul_receptor_ligand_kj_mol = ''
        delta_gelec_bind_kj_mol = ''
        closure_error_kj_mol = ''
        iterations = ''
        true_residual = ''
        audit_rms = ''
        process_wall_seconds = ''
        solver_total_seconds = ''
        surface_seconds = ''
        ilu_setup_seconds = ''
        gmres_seconds = ''
        peak_rss_mb = ''
        preconditioner_effective = ''
        ilu_regularization_selected = ''
        discretization = ''
        surface_generator = ''
        apbs_dime_x = ''
        apbs_dime_y = ''
        apbs_dime_z = ''
        apbs_actual_spacing_max_A = ''
        run_dir = $RunDir
    }
}

$fixture = Ensure-Fixture
if ($PrepareOnly) {
    Write-Host "Prepared 6LZG complex: $ComplexPqrPath"
    Write-Host "Prepared ACE2 receptor: $ReceptorPqrPath"
    Write-Host "Prepared RBD ligand: $LigandPqrPath"
    exit 0
}

$resolvedTitan = Resolve-Executable $TitanExe @('.\titan.exe')
$resolvedApbs = Resolve-Executable $ApbsExe @(
    'apbs.exe','apbs',
    (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
    'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe'
)
if ('titan' -in $Backends -and [string]::IsNullOrWhiteSpace($resolvedTitan)) {
    throw "Unable to resolve TITAN executable: $TitanExe"
}
if ('apbs' -in $Backends -and [string]::IsNullOrWhiteSpace($resolvedApbs)) {
    throw 'Unable to resolve APBS executable.'
}

$pqrPaths = @{
    complex = $ComplexPqrPath
    receptor = $ReceptorPqrPath
    ligand = $LigandPqrPath
}
$pqrShas = @{
    complex = Sha $ComplexPqrPath
    receptor = Sha $ReceptorPqrPath
    ligand = Sha $LigandPqrPath
}
$atomCounts = @{
    complex = Get-PqrAtomCount $ComplexPqrPath
    receptor = Get-PqrAtomCount $ReceptorPqrPath
    ligand = Get-PqrAtomCount $LigandPqrPath
}
$partitionSha = Sha256-Text ($pqrShas.complex + '|' + $pqrShas.receptor + '|' + $pqrShas.ligand)
$titanSha = if ($resolvedTitan) { Sha $resolvedTitan } else { '' }
$apbsSha = if ($resolvedApbs) { Sha $resolvedApbs } else { '' }
$common = '6LZG-A+B|ACE2-receptor-A|RBD-ligand-B|CHARMM-ph7|salt0.15|eps1-80|T298.15|threads4'
$titanProtocol = "$common|binding-cycle-v3-large|analytic-ses-h1.25|node-patch|p1|theta0.45|asi|precond-$Preconditioner|fmm-near-ilu$IluNeighborCap|coarse-cap$CoarseAggregateCap|ilu-budget-${IluMemoryBudgetMb}MB|ilu-reg-adaptive-0.01|tol1e-4|node-patch-operator-audit-$NodePatchOperatorAudit"
$apbsProtocol = "$common|mg-auto|smol|spl2|target-grid-$(F $TargetApbsGridSpacingA)|fine-pad-$(F $ApbsFinePaddingA)|coarse-pad-$(F $ApbsCoarsePaddingA)|maxdim-$ApbsMaxGridDimension"

$checkpoint = @{}
if (Test-Path -LiteralPath $MatrixPath) {
    foreach ($row in @(Import-Csv -LiteralPath $MatrixPath)) {
        $checkpoint[[string]$row.case_id] = $row
    }
}
function Save-Matrix {
    @($checkpoint.Values | Sort-Object backend,component) |
        Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
}
function Should-Run($Old, [string]$Protocol, [string]$ExeSha, [string]$PqrSha) {
    if ($Force -or $null -eq $Old) { return $true }
    if ([string]$Old.protocol -ne $Protocol -or
        [string]$Old.executable_sha256 -ne $ExeSha -or
        [string]$Old.pqr_sha256 -ne $PqrSha) { return $true }
    if ([string]$Old.status -in @('planned','running','incomplete')) { return $true }
    if ($RerunFailed -and [string]$Old.status -ne 'ok') { return $true }
    return $false
}

$cases = [Collections.Generic.List[object]]::new()
if ('titan' -in $Backends) {
    $cases.Add([pscustomobject]@{
        case_id='6lzg_titan_binding'; backend='titan'; component='binding'
        protocol=$titanProtocol; exe=$resolvedTitan; exe_sha=$titanSha
        pqr_sha=$partitionSha; atoms=$atomCounts.complex
    }) | Out-Null
}
if ('apbs' -in $Backends) {
    foreach ($component in @('complex','receptor','ligand')) {
        $cases.Add([pscustomobject]@{
            case_id="6lzg_apbs_$component"; backend='apbs'; component=$component
            protocol=$apbsProtocol; exe=$resolvedApbs; exe_sha=$apbsSha
            pqr_sha=$pqrShas[$component]; atoms=$atomCounts[$component]
        }) | Out-Null
    }
}

foreach ($case in $cases) {
    $caseId = [string]$case.case_id
    $backend = [string]$case.backend
    $protocol = [string]$case.protocol
    $exe = [string]$case.exe
    $exeSha = [string]$case.exe_sha
    $casePqrSha = [string]$case.pqr_sha
    $old = if ($checkpoint.ContainsKey($caseId)) { $checkpoint[$caseId] } else { $null }
    if (-not (Should-Run $old $protocol $exeSha $casePqrSha)) {
        Write-Host "resume: $caseId"
        continue
    }
    $runDir = Join-Path $WorkRoot $caseId
    if (Test-Path -LiteralPath $runDir) {
        Remove-Item -LiteralPath $runDir -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $row = New-Row -CaseId $caseId -Backend $backend -Component $case.component `
        -Protocol $protocol -ExeSha $exeSha -PqrSha $casePqrSha `
        -RunDir $runDir -Atoms $case.atoms
    $row.status = 'running'
    $checkpoint[$caseId] = $row
    Save-Matrix

    if ($backend -eq 'titan') {
        $localComplex = Join-Path $runDir (Split-Path -Leaf $ComplexPqrPath)
        $localReceptor = Join-Path $runDir (Split-Path -Leaf $ReceptorPqrPath)
        $localLigand = Join-Path $runDir (Split-Path -Leaf $LigandPqrPath)
        Copy-Item -LiteralPath $ComplexPqrPath -Destination $localComplex
        Copy-Item -LiteralPath $ReceptorPqrPath -Destination $localReceptor
        Copy-Item -LiteralPath $LigandPqrPath -Destination $localLigand
        $argPath = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
            -PqrPath $localComplex -Backend tabi -Mode treecode `
            -FmmMode single-matvec-analytic -FmmLocal off `
            -SurfacePoints 8 -Threads $Threads -MaxIter 720 -Restart 128 `
            -Tolerance 1e-4 -SaltM 0.15 -InternalDielectric 1.0 `
            -ExternalDielectric 80.0 -TemperatureK 298.15 -ProbeRadius 1.4 `
            -OpeningTheta 0.45 -AuditMode sampled -AuditSampleCount 512 `
            -Preconditioner $Preconditioner `
            -PreconditionerIluNeighborCap $IluNeighborCap -PreconditionerIluDropTol 0.0 `
            -PreconditionerIluRegularization adaptive `
            -PreconditionerIluRegularizationMax 0.01 `
            -PreconditionerMemoryBudgetMb $IluMemoryBudgetMb `
            -PreconditionerCoarseAggregateCap $CoarseAggregateCap `
            -PreconditionerCoarseVerify off `
            -NearIntegration analytic_singular -TaylorOrder 4 `
            -PrimarySurface analytic-ses -SesGridSpacing 1.25 `
            -AtomPotentials csv -BindingCycle pqr-components `
            -BindingReceptorPqr $localReceptor -BindingLigandPqr $localLigand
        Add-Content -LiteralPath $argPath -Encoding ASCII `
            -Value @(
                '-pbe-bie-discretization=node-patch',
                "-pbe-bie-node-patch-operator-audit=$NodePatchOperatorAudit"
            )
        Write-Host "run: $caseId (TITAN, timeout ${TitanTimeoutSeconds}s)"
        $result = Invoke-Case -Exe $exe -RunDir $runDir `
            -Arguments ('-a "' + ($argPath -replace '"','\"') + '"') `
            -TimeoutSeconds $TitanTimeoutSeconds
        $row.exit_code = $result.exit_code
        $row.timed_out = $result.timed_out
        $row.process_wall_seconds = $result.wall_seconds
        $row.peak_rss_mb = $result.peak_rss_mb
        if ($result.timed_out) {
            $row.status = 'timeout'; $row.failure = 'titan_timeout'
        } elseif ($result.exit_code -ne 0) {
            $row.status = 'failed'; $row.failure = 'titan_exit_nonzero'
        } else {
            try {
                $components = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_component_electrostatics.csv'))
                $bindings = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_binding_electrostatics.csv'))
                $profiles = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))
                if ($components.Count -ne 3 -or $bindings.Count -ne 1) {
                    throw "Expected three TITAN component rows and one binding row; found $($components.Count)/$($bindings.Count)."
                }
                if ($profiles.Count -ne 1) {
                    throw "Expected one TITAN profile row for the integrated binding cycle; found $($profiles.Count)."
                }
                $binding = $bindings[0]
                $profile = $profiles[0]
                $profileDiscretization = ([string]$profile.discretization).Trim().ToLowerInvariant()
                $profilePanels = I $profile.panels
                $profileVertices = I $profile.mesh_vertices
                $profileUnknowns = I $profile.unknown_count
                $complexComponent = @($components | Where-Object component -eq 'complex')
                if ($profileDiscretization -ne 'node_patch') {
                    throw "Fail-closed node-patch contract: profile discretization is '$profileDiscretization'. Rebuild titan.exe and rerun."
                }
                if ($profilePanels -le 0 -or $profilePanels -ne $profileVertices -or
                    $profileUnknowns -ne (2 * $profilePanels)) {
                    throw "Fail-closed node-patch contract: panels=$profilePanels mesh_vertices=$profileVertices unknowns=$profileUnknowns."
                }
                if ($complexComponent.Count -ne 1 -or
                    (I $complexComponent[0].surface_panels) -ne $profilePanels) {
                    throw "Fail-closed node-patch contract: component/profile panel mismatch."
                }
                $row.gpol_complex_kj_mol = D $binding.Gpol_complex_kJmol
                $row.gpol_receptor_kj_mol = D $binding.Gpol_receptor_kJmol
                $row.gpol_ligand_kj_mol = D $binding.Gpol_ligand_kJmol
                $row.delta_gpol_bind_kj_mol = D $binding.Delta_Gpol_bind_kJmol
                $row.e_coul_receptor_ligand_kj_mol = D $binding.E_coul_receptor_ligand_kJmol
                $row.delta_gelec_bind_kj_mol = D $binding.Delta_Gelec_bind_kJmol
                $row.energy_kj_mol = $row.delta_gpol_bind_kj_mol
                $row.energy_kcal_mol = (D $row.delta_gpol_bind_kj_mol) / 4.184
                $row.closure_error_kj_mol = D $binding.closure_error_max_kJmol
                $row.iterations = ($components | ForEach-Object { I $_.iterations } | Measure-Object -Maximum).Maximum
                $row.true_residual = ($components | ForEach-Object { D $_.true_residual } | Measure-Object -Maximum).Maximum
                $row.audit_rms = D $binding.sampled_audit_max
                $row.panels = ($components | ForEach-Object { I $_.surface_panels } | Measure-Object -Maximum).Maximum
                $row.solver_total_seconds = $result.wall_seconds
                $row.gmres_seconds = ($components | ForEach-Object { D $_.gmres_seconds 0.0 } | Measure-Object -Sum).Sum
                $row.preconditioner_effective = [string]$profile.preconditioner_effective
                $row.ilu_regularization_selected = D $profile.preconditioner_ilu_regularization_selected 0.0
                $row.discretization = $profileDiscretization
                $row.surface_generator = (@($components.surface_generator | Sort-Object -Unique) -join ';')
                $componentOk = @($components | Where-Object {
                    $_.status -eq 'ok' -and (D $_.true_residual) -le 1e-4 -and
                    (D $_.sampled_direct_operator_rms) -le 8e-3
                }).Count -eq 3
                if (-not $componentOk -or $binding.status -ne 'ok' -or
                    -not (Is-Finite (D $row.delta_gpol_bind_kj_mol)) -or
                    (D $row.closure_error_kj_mol) -gt 1e-8 -or
                    (D $binding.sampled_audit_max) -gt 8e-3 -or
                    (D $binding.sampled_direct_true_residual_max) -gt 8e-3) {
                    $row.status='accuracy_gap'
                    $row.failure="component_gate=$componentOk binding_status=$($binding.status) closure=$($row.closure_error_kj_mol) audit=$($binding.sampled_audit_max)"
                } else { $row.status='ok'; $row.failure='' }
            } catch {
                $row.status='schema_gap'; $row.failure=$_.Exception.Message
            }
        }
    } else {
        $sourcePqr = [string]$pqrPaths[[string]$case.component]
        $localPqr = Join-Path $runDir (Split-Path -Leaf $sourcePqr)
        Copy-Item -LiteralPath $sourcePqr -Destination $localPqr
        $grid = Write-ApbsInput -LocalPqr $localPqr `
            -InputPath (Join-Path $runDir 'apbs.in')
        $row.apbs_dime_x=$grid.nx; $row.apbs_dime_y=$grid.ny
        $row.apbs_dime_z=$grid.nz
        $row.apbs_actual_spacing_max_A=$grid.actual_spacing_max
        Write-Host "run: $caseId (APBS, timeout ${ApbsTimeoutSeconds}s)"
        $environment = @{
            OMP_NUM_THREADS=[string]$Threads
            OMP_THREAD_LIMIT=[string]$Threads
            OPENBLAS_NUM_THREADS='1'
            MKL_NUM_THREADS='1'
            OMP_PROC_BIND='true'
        }
        $result = Invoke-Case -Exe $exe -RunDir $runDir `
            -Arguments 'apbs.in' -TimeoutSeconds $ApbsTimeoutSeconds `
            -Environment $environment
        $row.exit_code=$result.exit_code; $row.timed_out=$result.timed_out
        $row.process_wall_seconds=$result.wall_seconds
        $row.peak_rss_mb=$result.peak_rss_mb
        if ($result.timed_out) {
            $row.status='timeout'; $row.failure='apbs_timeout'
        } elseif ($result.exit_code -ne 0) {
            $row.status='failed'; $row.failure='apbs_exit_nonzero'
        } else {
            $energyKj = Read-ApbsEnergy $runDir
            $row.energy_kj_mol=$energyKj
            $row.energy_kcal_mol=$energyKj/4.184
            if (Is-Finite $energyKj) {
                $row.status='ok'; $row.failure=''
            } else {
                $row.status='schema_gap'; $row.failure='apbs_energy_parse_gap'
            }
        }
    }
    $checkpoint[$caseId]=$row
    Save-Matrix
}

$rows = @($checkpoint.Values | Where-Object { $_.backend -in $Backends })
$componentRows = [Collections.Generic.List[object]]::new()
$titan = @($rows | Where-Object case_id -eq '6lzg_titan_binding') | Select-Object -Last 1
if ($null -ne $titan -and (Test-Path -LiteralPath (Join-Path $titan.run_dir 'pbe_component_electrostatics.csv'))) {
    foreach ($source in @(Import-Csv -LiteralPath (Join-Path $titan.run_dir 'pbe_component_electrostatics.csv'))) {
        $componentRows.Add([pscustomobject]@{
            system='6LZG ACE2-RBD'; backend='TITAN p1'; component=$source.component
            atom_count=$source.atom_count; panels=$source.surface_panels
            energy_kj_mol=$source.reaction_energy_kJmol; iterations=$source.iterations
            true_residual=$source.true_residual; audit_rms=$source.sampled_direct_operator_rms
            wall_seconds=''; peak_rss_mb=$titan.peak_rss_mb; status=$source.status
        }) | Out-Null
    }
}
$apbsByComponent = @{}
foreach ($component in @('complex','receptor','ligand')) {
    $a = @($rows | Where-Object case_id -eq "6lzg_apbs_$component") | Select-Object -Last 1
    if ($null -ne $a) {
        $apbsByComponent[$component] = $a
        $componentRows.Add([pscustomobject]@{
            system='6LZG ACE2-RBD'; backend='APBS'; component=$component
            atom_count=$a.atom_count; panels=''; energy_kj_mol=$a.energy_kj_mol
            iterations=''; true_residual=''; audit_rms=''
            wall_seconds=$a.process_wall_seconds; peak_rss_mb=$a.peak_rss_mb; status=$a.status
        }) | Out-Null
    }
}
$componentRows.ToArray() | Export-Csv -LiteralPath $ComponentPath -NoTypeInformation

$cycles = [Collections.Generic.List[object]]::new()
if ($null -ne $titan) {
    $cycles.Add([pscustomobject][ordered]@{
        system='6LZG ACE2-RBD'; method='TITAN p1'
        gpol_complex_kj_mol=$titan.gpol_complex_kj_mol
        gpol_receptor_kj_mol=$titan.gpol_receptor_kj_mol
        gpol_ligand_kj_mol=$titan.gpol_ligand_kj_mol
        delta_gpol_bind_kj_mol=$titan.delta_gpol_bind_kj_mol
        e_coul_receptor_ligand_kj_mol=$titan.e_coul_receptor_ligand_kj_mol
        delta_gelec_bind_kj_mol=$titan.delta_gelec_bind_kj_mol
        wall_seconds=$titan.process_wall_seconds; peak_rss_mb=$titan.peak_rss_mb
        status=$titan.status; protocol=$titan.protocol
    }) | Out-Null
}
if ($apbsByComponent.Count -eq 3) {
    $gc=D $apbsByComponent.complex.energy_kj_mol
    $gr=D $apbsByComponent.receptor.energy_kj_mol
    $gl=D $apbsByComponent.ligand.energy_kj_mol
    $apbsStatus = if (@($apbsByComponent.Values | Where-Object status -ne 'ok').Count -eq 0) { 'ok' } else { 'component_failure' }
    $cycles.Add([pscustomobject][ordered]@{
        system='6LZG ACE2-RBD'; method='APBS'
        gpol_complex_kj_mol=$gc; gpol_receptor_kj_mol=$gr; gpol_ligand_kj_mol=$gl
        delta_gpol_bind_kj_mol=$gc-$gr-$gl
        e_coul_receptor_ligand_kj_mol=if($null -ne $titan){$titan.e_coul_receptor_ligand_kj_mol}else{''}
        delta_gelec_bind_kj_mol=if($null -ne $titan -and (Is-Finite (D $titan.e_coul_receptor_ligand_kj_mol))){
            $gc-$gr-$gl+(D $titan.e_coul_receptor_ligand_kj_mol)
        }else{''}
        wall_seconds=($apbsByComponent.Values | ForEach-Object { D $_.process_wall_seconds 0.0 } | Measure-Object -Sum).Sum
        peak_rss_mb=($apbsByComponent.Values | ForEach-Object { D $_.peak_rss_mb 0.0 } | Measure-Object -Maximum).Maximum
        status=$apbsStatus; protocol=$apbsProtocol
    }) | Out-Null
}
$cycles.ToArray() | Export-Csv -LiteralPath $ComparisonPath -NoTypeInformation

$decision = if (@($rows | Where-Object status -eq 'timeout').Count -gt 0) {
    '6lzg_binding_runtime_limited'
} elseif (@($rows | Where-Object status -ne 'ok').Count -gt 0 -or
    $rows.Count -ne $cases.Count) {
    '6lzg_binding_accuracy_gap'
} else {
    '6lzg_binding_titan_apbs_ready'
}
[pscustomobject]@{
    decision=$decision; fixture_manifest=$ManifestPath
    matrix=$MatrixPath; comparison=$ComparisonPath
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision",
    "fixture_manifest=$ManifestPath",
    "matrix=$MatrixPath",
    "components=$ComponentPath",
    "comparison=$ComparisonPath"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN 6LZG ACE2-RBD Polar Solvation Binding Cycle versus APBS',
    "Decision: $decision",
    'Fixture: RCSB 6LZG, model 1, ACE2 chain A as receptor and spike RBD chain B as ligand.',
    'Preparation: the A+B complex is protonated once with PDB2PQR 3.7.1, CHARMM, pH 7.0; receptor and ligand PQRs are an exact line partition used by both solvers.',
    'Reported quantities: Delta Gpol,bind = Gpol(A+B)-Gpol(A)-Gpol(B); Delta Gelec,bind additionally includes the direct inter-chain Coulomb interaction from the same PQR charges and eps_in=1.',
    "TITAN: one integrated pqr-components binding cycle; analytic SES h=1.25 A, node-patch, p1 analytical FMM/M2L, theta 0.45, analytic-singular near field, adaptive FMM-near ILU(0) cap $IluNeighborCap with ${IluMemoryBudgetMb}-MB preconditioner budget, tol 1e-4, four threads.",
    "APBS: mg-auto, smol/spl2, target spacing $(F $TargetApbsGridSpacingA) A, four threads.",
    'The comparison is contextual because the BIE surface discretization and APBS finite-difference grid are different.',
    "Matrix: $MatrixPath",
    "Components: $ComponentPath",
    "Comparison: $ComparisonPath",
    'Scope: 6LZG is used only as a public numerical test system; no biological binding interpretation is made.'
) | Set-Content -LiteralPath $ReportPath -Encoding ASCII

Write-Host "Decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Components: $ComponentPath"
Write-Host "Comparison: $ComparisonPath"
Write-Host 'Resume: rerun unchanged; add -RerunFailed for failed/timeout rows or -Force for all process rows.'
if (-not $AllowFailedDecision -and $decision -ne '6lzg_binding_titan_apbs_ready') {
    throw "6LZG ACE2-RBD binding comparison did not pass: $decision"
}
