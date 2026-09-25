param(
    [string]$TitanExe = '.\titan.exe',
    [string]$ApbsExe = '',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_fagbi27_analytic_ses_apbs',
    [string]$FixtureRoot = 'test\pbe_inputs_fagbi27',
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$Pdb2PqrTarget = '',
    [string[]]$Fixtures = @(),
    [ValidateSet('titan','apbs')]
    [string[]]$Backends = @('titan','apbs'),
    [int]$Threads = 4,
    [int]$TitanTimeoutSeconds = 3600,
    [int]$ApbsTimeoutSeconds = 1800,
    [double]$TargetApbsGridSpacingA = 0.50,
    [double]$ApbsFinePaddingA = 12.0,
    [double]$ApbsCoarsePaddingA = 16.0,
    [int]$ApbsMaxGridDimension = 385,
    [ValidateSet('off','adaptive')]
    [string]$PreconditionerIluRegularization = 'off',
    [ValidateRange(0.0,0.01)]
    [double]$PreconditionerIluRegularizationMax = 0.01,
    [string]$ReportOutput = '',
    [switch]$Quick,
    [switch]$PrepareOnly,
    [switch]$RefreshFixtures,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunFixtures = @(),
    [ValidateSet('titan','apbs')]
    [string[]]$RerunBackends = @(),
    [switch]$PlotOnly,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$Value) { $Value.ToString('G15', $Inv) }
function D($Value, [double]$Fallback = [double]::NaN) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return $Fallback }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float, $Inv, [ref]$parsed)) { return $parsed }
    return $Fallback
}
function I($Value, [int64]$Fallback = 0) {
    $parsed = 0L
    if ([int64]::TryParse([string]$Value, [ref]$parsed)) { return $parsed }
    return $Fallback
}
function Is-Finite([double]$Value) { -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value) }
function Xml([string]$Value) { [Security.SecurityElement]::Escape($Value) }
function Sha([string]$Path) { (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant() }
function Safe-Id([string]$Value) { ($Value.ToLowerInvariant() -replace '[^a-z0-9]+','_').Trim('_') }

if ($Threads -ne 4) { throw 'This manuscript gate is locked to Threads=4 for both TITAN and APBS.' }
if ($TargetApbsGridSpacingA -le 0.0) { throw 'TargetApbsGridSpacingA must be positive.' }
if ($ApbsMaxGridDimension -lt 65 -or (($ApbsMaxGridDimension - 1) % 32) -ne 0) {
    throw 'ApbsMaxGridDimension must be of the form 32*n+1 and at least 65.'
}

$PaperSpecs = @(
    [pscustomobject]@{ paper_index=1;  pdb_id='1AJJ'; paper_elements=40496;  paper_energy_kcal=-1141.15 },
    [pscustomobject]@{ paper_index=2;  pdb_id='2ERL'; paper_elements=43214;  paper_energy_kcal=-953.42 },
    [pscustomobject]@{ paper_index=3;  pdb_id='1CBN'; paper_elements=44367;  paper_energy_kcal=-305.94 },
    [pscustomobject]@{ paper_index=4;  pdb_id='1VII'; paper_elements=47070;  paper_energy_kcal=-906.11 },
    [pscustomobject]@{ paper_index=5;  pdb_id='1FCA'; paper_elements=47461;  paper_energy_kcal=-1206.48 },
    [pscustomobject]@{ paper_index=6;  pdb_id='1BBL'; paper_elements=49071;  paper_energy_kcal=-991.22 },
    [pscustomobject]@{ paper_index=7;  pdb_id='2PDE'; paper_elements=50518;  paper_energy_kcal=-829.46 },
    [pscustomobject]@{ paper_index=8;  pdb_id='1SH1'; paper_elements=51186;  paper_energy_kcal=-756.63 },
    [pscustomobject]@{ paper_index=9;  pdb_id='1VJW'; paper_elements=52536;  paper_energy_kcal=-1242.56 },
    [pscustomobject]@{ paper_index=10; pdb_id='1UXC'; paper_elements=53602;  paper_energy_kcal=-1145.38 },
    [pscustomobject]@{ paper_index=11; pdb_id='1PTQ'; paper_elements=54256;  paper_energy_kcal=-877.84 },
    [pscustomobject]@{ paper_index=12; pdb_id='1BOR'; paper_elements=54628;  paper_energy_kcal=-857.27 },
    [pscustomobject]@{ paper_index=13; pdb_id='1FXD'; paper_elements=54692;  paper_energy_kcal=-3318.14 },
    [pscustomobject]@{ paper_index=14; pdb_id='1R69'; paper_elements=57646;  paper_energy_kcal=-1094.86 },
    [pscustomobject]@{ paper_index=15; pdb_id='1MBG'; paper_elements=58473;  paper_energy_kcal=-1357.33 },
    [pscustomobject]@{ paper_index=16; pdb_id='1BPI'; paper_elements=60600;  paper_energy_kcal=-1310.02 },
    [pscustomobject]@{ paper_index=17; pdb_id='1HPT'; paper_elements=61164;  paper_energy_kcal=-817.34 },
    [pscustomobject]@{ paper_index=18; pdb_id='451C'; paper_elements=79202;  paper_energy_kcal=-1031.91 },
    [pscustomobject]@{ paper_index=19; pdb_id='1SVR'; paper_elements=88198;  paper_energy_kcal=-1718.97 },
    [pscustomobject]@{ paper_index=20; pdb_id='1FRD'; paper_elements=81792;  paper_energy_kcal=-2867.32 },
    [pscustomobject]@{ paper_index=21; pdb_id='1A2S'; paper_elements=84527;  paper_energy_kcal=-1925.24 },
    [pscustomobject]@{ paper_index=22; pdb_id='1NEQ'; paper_elements=89457;  paper_energy_kcal=-1740.49 },
    [pscustomobject]@{ paper_index=23; pdb_id='1A63'; paper_elements=132133; paper_energy_kcal=-2382.50 },
    [pscustomobject]@{ paper_index=24; pdb_id='1A7M'; paper_elements=147121; paper_energy_kcal=-2172.12 },
    [pscustomobject]@{ paper_index=25; pdb_id='2GO0'; paper_elements=111615; paper_energy_kcal=-1968.65 },
    [pscustomobject]@{ paper_index=26; pdb_id='1UV0'; paper_elements=128497; paper_energy_kcal=-2296.43 },
    [pscustomobject]@{ paper_index=27; pdb_id='4MTH'; paper_elements=123737; paper_energy_kcal=-2479.61 }
)

function Normalize-List([string[]]$Values) {
    return @($Values | ForEach-Object { @([string]$_ -split ',') } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

$Fixtures = @(Normalize-List $Fixtures)
$Backends = @(Normalize-List $Backends | ForEach-Object { $_.ToLowerInvariant() } | Sort-Object -Unique)
$RerunFixtures = @(Normalize-List $RerunFixtures | ForEach-Object { $_.ToUpperInvariant() })
$RerunBackends = @(Normalize-List $RerunBackends | ForEach-Object { $_.ToLowerInvariant() })
if ($Quick) { $Fixtures = @('1AJJ') }
if ($Fixtures.Count -eq 0) { $Fixtures = @($PaperSpecs | ForEach-Object pdb_id) }
$Fixtures = @($Fixtures | ForEach-Object { $_.ToUpperInvariant() } | Sort-Object -Unique)

$SelectedSpecs = @($PaperSpecs | Where-Object { $_.pdb_id -in $Fixtures } | Sort-Object paper_index)
if ($SelectedSpecs.Count -ne $Fixtures.Count) {
    $known = @($SelectedSpecs | ForEach-Object pdb_id)
    $missing = @($Fixtures | Where-Object { $_ -notin $known })
    throw "Unknown FAGBI Table-4 PDB IDs: $($missing -join ', ')"
}

if (-not [IO.Path]::IsPathRooted($WorkRoot)) { $WorkRoot = Join-Path $RepoRoot $WorkRoot }
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) { $FixtureRoot = Join-Path $RepoRoot $FixtureRoot }
if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) { $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1' }
New-Item -ItemType Directory -Force -Path $WorkRoot,$FixtureRoot | Out-Null

$ManifestPath = Join-Path $FixtureRoot 'fagbi27_fixture_manifest.csv'
$MatrixPath = Join-Path $WorkRoot 'fagbi27_analytic_ses_apbs_matrix.csv'
$PairPath = Join-Path $WorkRoot 'fagbi27_analytic_ses_apbs_pairs.csv'
$SupportingTablePath = Join-Path $WorkRoot 'fagbi27_supporting_table.csv'
$DecisionPath = Join-Path $WorkRoot 'fagbi27_analytic_ses_apbs_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi27_analytic_ses_apbs_summary.dat'
$FigurePath = Join-Path $WorkRoot 'fagbi27_analytic_ses_apbs_figure.svg'
$PaperSetPath = Join-Path $WorkRoot 'fagbi27_paper_table4_set.csv'
$ReportPath = if ([string]::IsNullOrWhiteSpace($ReportOutput)) {
    Join-Path $RepoRoot 'TITAN_TABI_FAGBI27_ANALYTIC_SES_APBS_REPORT.txt'
} elseif ([IO.Path]::IsPathRooted($ReportOutput)) { $ReportOutput }
else { Join-Path $RepoRoot $ReportOutput }
$PaperSpecs | Export-Csv -LiteralPath $PaperSetPath -NoTypeInformation

function Resolve-Apbs([string]$Requested) {
    if ($Requested) {
        if (Test-Path -LiteralPath $Requested -PathType Leaf) { return (Resolve-Path $Requested).Path }
        $cmd = Get-Command $Requested -ErrorAction SilentlyContinue
        if ($cmd) { return $cmd.Source }
    }
    foreach ($name in @('apbs.exe','apbs')) {
        $cmd = Get-Command $name -ErrorAction SilentlyContinue
        if ($cmd) { return $cmd.Source }
    }
    foreach ($path in @(
        (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
        'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe'
    )) {
        if (Test-Path -LiteralPath $path -PathType Leaf) { return (Resolve-Path $path).Path }
    }
    return ''
}

function Get-AtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object { $_ -match '^(ATOM|HETATM)' }).Count
}

function Write-FirstModelProteinPdb([string]$InputPath, [string]$OutputPath) {
    $lines = @(Get-Content -LiteralPath $InputPath)
    $hasModels = @($lines | Where-Object { $_ -match '^MODEL\s' }).Count -gt 0
    $inFirstModel = -not $hasModels
    $seenFirstModel = $false
    $atomRows = [Collections.Generic.List[object]]::new()
    $ordinal = 0
    foreach ($line in $lines) {
        if ($line -match '^MODEL\s') {
            if (-not $seenFirstModel) {
                $seenFirstModel = $true
                $inFirstModel = $true
            } else { $inFirstModel = $false }
            continue
        }
        if ($line -match '^ENDMDL') {
            if ($inFirstModel) { break }
            continue
        }
        if (-not $inFirstModel -or -not $line.StartsWith('ATOM')) { continue }
        $fixed = $line.PadRight(80)
        $atomName = $fixed.Substring(12,4).Trim()
        $element = $fixed.Substring(76,2).Trim().ToUpperInvariant()
        $nameWithoutLeadingDigits = $atomName -replace '^\d+',''
        if ($element -in @('H','D') -or $nameWithoutLeadingDigits -match '^[HD]') {
            continue
        }
        $altLoc = $fixed.Substring(16,1)
        $occupancy = D $fixed.Substring(54,6) 0.0
        $residueKey = @(
            $fixed.Substring(21,1), $fixed.Substring(22,4),
            $fixed.Substring(26,1)
        ) -join '|'
        $atomRows.Add([pscustomobject]@{
            line = $fixed.Remove(16,1).Insert(16,' ')
            residue_key = $residueKey
            residue_name = $fixed.Substring(17,3)
            atom_name = $atomName
            alt_loc = $altLoc
            occupancy = $occupancy
            ordinal = $ordinal
        })
        ++$ordinal
    }
    $preferredConformer = @{}
    foreach ($residueGroup in @($atomRows | Group-Object residue_key)) {
        $alternates = @($residueGroup.Group | Where-Object {
            -not [string]::IsNullOrWhiteSpace($_.alt_loc)
        })
        if ($alternates.Count -eq 0) { continue }
        $choice = @($alternates | Group-Object {
            "$($_.alt_loc)|$($_.residue_name)"
        } | ForEach-Object {
            $first = $_.Group[0]
            [pscustomobject]@{
                key = $_.Name
                occupancy = ($_.Group | Measure-Object occupancy -Maximum).Maximum
                rank = if ($first.alt_loc -eq 'A') { 1 } else { 0 }
            }
        } | Sort-Object @{Expression='occupancy';Descending=$true},
            @{Expression='rank';Descending=$true},key | Select-Object -First 1)
        $preferredConformer[$residueGroup.Name] = $choice.key
    }
    $candidates = @{}
    foreach ($row in $atomRows) {
        if (-not [string]::IsNullOrWhiteSpace($row.alt_loc)) {
            $selected = $preferredConformer[$row.residue_key]
            if ($selected -ne "$($row.alt_loc)|$($row.residue_name)") { continue }
        }
        $key = "$($row.residue_key)|$($row.atom_name)"
        $rank = if ([string]::IsNullOrWhiteSpace($row.alt_loc)) { 2 }
            elseif ($row.alt_loc -eq 'A') { 1 } else { 0 }
        $candidate = [pscustomobject]@{
            line = $row.line
            occupancy = $row.occupancy
            rank = $rank
            ordinal = $row.ordinal
        }
        $old = $candidates[$key]
        if ($null -eq $old -or $candidate.occupancy -gt $old.occupancy -or
            ($candidate.occupancy -eq $old.occupancy -and $candidate.rank -gt $old.rank)) {
            $candidates[$key] = $candidate
        }
    }
    $out = [Collections.Generic.List[string]]::new()
    foreach ($candidate in @($candidates.Values | Sort-Object ordinal)) {
        $out.Add($candidate.line.TrimEnd())
    }
    if ($out.Count -eq 0) { throw "No model-1 ATOM rows in $InputPath" }
    $out.Add('TER'); $out.Add('END')
    $out | Set-Content -LiteralPath $OutputPath -Encoding ASCII
}

function Install-Pdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe)) { throw "Missing Python: $PythonExe" }
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        New-Item -ItemType Directory -Force -Path $Pdb2PqrTarget | Out-Null
        & $PythonExe -m pip install --disable-pip-version-check --target $Pdb2PqrTarget 'pdb2pqr==3.7.1'
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install PDB2PQR 3.7.1.' }
    }
}

function Invoke-Pdb2Pqr([string]$InputPath, [string]$OutputPath, [string]$LogPath) {
    $oldPythonPath = $env:PYTHONPATH
    $oldErrorActionPreference = $ErrorActionPreference
    try {
        $env:PYTHONPATH = if ($oldPythonPath) { "$Pdb2PqrTarget;$oldPythonPath" } else { $Pdb2PqrTarget }
        # PDB2PQR writes routine INFO messages to stderr. Do not let Windows
        # PowerShell promote those lines to terminating NativeCommandError records.
        $ErrorActionPreference = 'Continue'
        $output = & $PythonExe -m pdb2pqr --ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath 2>&1
        $exitCode = $LASTEXITCODE
        $output | Set-Content -LiteralPath $LogPath -Encoding UTF8
        if ($exitCode -ne 0 -or -not (Test-Path -LiteralPath $OutputPath)) { throw "PDB2PQR failed with exit code $exitCode; see $LogPath" }
    } finally {
        $ErrorActionPreference = $oldErrorActionPreference
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-PqrRadii([string]$Path, [double]$Floor = 0.6) {
    $out = [Collections.Generic.List[string]]::new(); $clamped = 0
    foreach ($line0 in Get-Content -LiteralPath $Path) {
        $line = $line0
        if ($line -notmatch '^(ATOM|HETATM)') { $out.Add($line); continue }
        $fields = @($line -split '\s+' | Where-Object { $_ })
        if ($fields.Count -lt 10) { throw "Malformed PQR row: $line" }
        $radius = D $fields[-1]
        if (-not (Is-Finite $radius)) { throw "Malformed PQR radius: $line" }
        if ($radius -le 0.0) {
            $start = $line.LastIndexOf($fields[-1], [StringComparison]::Ordinal)
            $line = $line.Substring(0,$start) + $Floor.ToString('F4',$Inv)
            ++$clamped
        }
        $out.Add($line)
    }
    $out | Set-Content -LiteralPath $Path -Encoding ASCII
    return $clamped
}

function Save-Manifest([hashtable]$Rows) {
    @($Rows.Values | Sort-Object { I $_.paper_index }) | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation
}

function Ensure-Fixtures([object[]]$Specs) {
    $rows = @{}
    if (Test-Path $ManifestPath) {
        foreach ($row in @(Import-Csv $ManifestPath)) { $rows[[string]$row.pdb_id] = $row }
    }
    $needsPreparation = @($Specs | Where-Object {
        $id = $_.pdb_id
        if ($RefreshFixtures -or -not $rows.ContainsKey($id)) { return $true }
        $old = $rows[$id]
        $pqr = Join-Path $FixtureRoot ([string]$old.pqr_file)
        return [string]$old.status -ne 'ok' -or -not (Test-Path $pqr) -or (Sha $pqr) -ne [string]$old.pqr_sha256
    })
    if ($needsPreparation.Count -gt 0) { Install-Pdb2Pqr }
    foreach ($spec in $Specs) {
        $id = [string]$spec.pdb_id
        $key = $id.ToLowerInvariant()
        if ($spec -notin $needsPreparation) {
            Write-Host "fixture resume: $id"
            continue
        }
        $rawName = "${key}_source.pdb"
        $preparedName = "${key}_model1_all_protein.pdb"
        $pqrName = "${key}_full_charmm_ph7.pqr"
        $rawPath = Join-Path $FixtureRoot $rawName
        $preparedPath = Join-Path $FixtureRoot $preparedName
        $pqrPath = Join-Path $FixtureRoot $pqrName
        $logPath = Join-Path $FixtureRoot "${key}_pdb2pqr.log"
        $url = "https://files.rcsb.org/download/$id.pdb"
        Write-Host "fixture prepare: $id"
        try {
            Invoke-WebRequest -UseBasicParsing -Uri $url -OutFile $rawPath
            Write-FirstModelProteinPdb -InputPath $rawPath -OutputPath $preparedPath
            Invoke-Pdb2Pqr -InputPath $preparedPath -OutputPath $pqrPath -LogPath $logPath
            $clamped = Normalize-PqrRadii -Path $pqrPath
            $atoms = Get-AtomCount $pqrPath
            if ($atoms -lt 20) { throw "Unexpectedly small PQR ($atoms atoms)" }
            $rows[$id] = [pscustomobject][ordered]@{
                paper_index=$spec.paper_index; pdb_id=$id; status='ok'; source_url=$url
                model_policy='first model'; chain_policy='all protein ATOM chains'
                raw_pdb_file=$rawName; prepared_pdb_file=$preparedName; pqr_file=$pqrName
                raw_pdb_sha256=(Sha $rawPath); prepared_pdb_sha256=(Sha $preparedPath); pqr_sha256=(Sha $pqrPath)
                atom_count=$atoms; force_field='CHARMM'; ph='7.0'; pdb2pqr_version='3.7.1'
                pdb2pqr_options='--ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0'
                radius_floor_A='0.6'; radius_clamped_atoms=$clamped; failure=''
            }
        } catch {
            $rows[$id] = [pscustomobject][ordered]@{
                paper_index=$spec.paper_index; pdb_id=$id; status='failed'; source_url=$url
                model_policy='first model'; chain_policy='all protein ATOM chains'
                raw_pdb_file=$rawName; prepared_pdb_file=$preparedName; pqr_file=$pqrName
                raw_pdb_sha256=''; prepared_pdb_sha256=''; pqr_sha256=''; atom_count=0
                force_field='CHARMM'; ph='7.0'; pdb2pqr_version='3.7.1'
                pdb2pqr_options='--ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0'
                radius_floor_A='0.6'; radius_clamped_atoms=0; failure=$_.Exception.Message
            }
        }
        Save-Manifest $rows
    }
    return $rows
}

if (-not $PlotOnly) { $FixtureRows = Ensure-Fixtures $SelectedSpecs } else {
    $FixtureRows = @{}
    if (Test-Path $ManifestPath) { foreach ($row in @(Import-Csv $ManifestPath)) { $FixtureRows[[string]$row.pdb_id] = $row } }
}
if ($PrepareOnly) {
    $failed = @($SelectedSpecs | Where-Object { -not $FixtureRows.ContainsKey($_.pdb_id) -or [string]$FixtureRows[$_.pdb_id].status -ne 'ok' })
    Write-Host "Fixture preparation complete: $($SelectedSpecs.Count - $failed.Count)/$($SelectedSpecs.Count)"
    if ($failed.Count -gt 0 -and -not $AllowFailedDecision) { throw 'One or more fixture preparations failed.' }
    exit 0
}

$ResolvedApbs = Resolve-Apbs $ApbsExe
if ('apbs' -in $Backends -and -not $PlotOnly -and -not $ResolvedApbs) { throw 'APBS was requested but apbs.exe could not be resolved.' }
$ResolvedTitan = ''
if ('titan' -in $Backends -and -not $PlotOnly) {
    $candidate = if ([IO.Path]::IsPathRooted($TitanExe)) { $TitanExe } else { Join-Path $RepoRoot $TitanExe }
    $ResolvedTitan = (Resolve-Path -LiteralPath $candidate).Path
}
$TitanSha = if ($ResolvedTitan) { Sha $ResolvedTitan } else { '' }
$ApbsSha = if ($ResolvedApbs) { Sha $ResolvedApbs } else { '' }
$ProtocolCommon = 'CHARMM-ph7|salt0.15|eps1-80|T298.15|threads4'
$TitanProtocol = "$ProtocolCommon|analytic-ses-h1.25|p1|theta0.45|asi|fmm-near-ilu64|ilu-reg-$PreconditionerIluRegularization-$(F $PreconditionerIluRegularizationMax)|tol1e-4|audit512"
$ApbsProtocol = "$ProtocolCommon|mg-auto|smol|spl2|target-grid-$($TargetApbsGridSpacingA)|fine-pad-$ApbsFinePaddingA|coarse-pad-$ApbsCoarsePaddingA|maxdim-$ApbsMaxGridDimension"

function Get-PqrAtoms([string]$Path) {
    $atoms = [Collections.Generic.List[object]]::new()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 10) { continue }
        $atoms.Add([pscustomobject]@{ x=(D $parts[-5]); y=(D $parts[-4]); z=(D $parts[-3]); radius=(D $parts[-1] 1.5) })
    }
    return $atoms.ToArray()
}

function Get-ApbsGeometry([string]$PqrPath) {
    $atoms = @(Get-PqrAtoms $PqrPath)
    if ($atoms.Count -eq 0) { throw "No PQR atoms: $PqrPath" }
    $xmin=1e300; $xmax=-1e300; $ymin=1e300; $ymax=-1e300; $zmin=1e300; $zmax=-1e300
    foreach ($a in $atoms) {
        $r=[math]::Max(0.0,[double]$a.radius)
        $xmin=[math]::Min($xmin,$a.x-$r); $xmax=[math]::Max($xmax,$a.x+$r)
        $ymin=[math]::Min($ymin,$a.y-$r); $ymax=[math]::Max($ymax,$a.y+$r)
        $zmin=[math]::Min($zmin,$a.z-$r); $zmax=[math]::Max($zmax,$a.z+$r)
    }
    $fgx=($xmax-$xmin)+2*$ApbsFinePaddingA; $fgy=($ymax-$ymin)+2*$ApbsFinePaddingA; $fgz=($zmax-$zmin)+2*$ApbsFinePaddingA
    return [pscustomobject]@{ atom_count=$atoms.Count; fgx=$fgx; fgy=$fgy; fgz=$fgz; cgx=$fgx+2*$ApbsCoarsePaddingA; cgy=$fgy+2*$ApbsCoarsePaddingA; cgz=$fgz+2*$ApbsCoarsePaddingA }
}

function Get-ApbsDime([double]$LengthA) {
    $required = [int][math]::Ceiling($LengthA / $TargetApbsGridSpacingA) + 1
    $dime = [int](32 * [math]::Ceiling(([math]::Max(65,$required)-1)/32.0) + 1)
    if ($dime -gt $ApbsMaxGridDimension) { throw "Required APBS grid dimension $dime exceeds cap $ApbsMaxGridDimension" }
    return $dime
}

function Write-ApbsInput([string]$PqrPath, [string]$InputPath) {
    $g = Get-ApbsGeometry $PqrPath
    $nx=Get-ApbsDime $g.fgx; $ny=Get-ApbsDime $g.fgy; $nz=Get-ApbsDime $g.fgz
    $pqrName=Split-Path -Leaf $PqrPath
    $lines=@(
        'read', "    mol pqr $pqrName", 'end', 'elec name solvated', '    mg-auto',
        "    dime $nx $ny $nz", "    cglen $(F $g.cgx) $(F $g.cgy) $(F $g.cgz)", "    fglen $(F $g.fgx) $(F $g.fgy) $(F $g.fgz)",
        '    cgcent mol 1','    fgcent mol 1','    mol 1','    lpbe','    bcfl sdh',
        '    ion charge 1 conc 0.15 radius 2.0','    ion charge -1 conc 0.15 radius 2.0',
        '    pdie 1.0','    sdie 80.0','    srfm smol','    chgm spl2','    sdens 10.0','    srad 1.4','    swin 0.3','    temp 298.15','    calcenergy total','    calcforce no','end',
        'elec name reference','    mg-auto',"    dime $nx $ny $nz", "    cglen $(F $g.cgx) $(F $g.cgy) $(F $g.cgz)", "    fglen $(F $g.fgx) $(F $g.fgy) $(F $g.fgz)",
        '    cgcent mol 1','    fgcent mol 1','    mol 1','    lpbe','    bcfl sdh','    pdie 1.0','    sdie 1.0','    srfm smol','    chgm spl2','    sdens 10.0','    srad 1.4','    swin 0.3','    temp 298.15','    calcenergy total','    calcforce no','end',
        'print elecEnergy solvated - reference end','quit'
    )
    $lines | Set-Content -LiteralPath $InputPath -Encoding ASCII
    return [pscustomobject]@{ nx=$nx; ny=$ny; nz=$nz; fgx=$g.fgx; fgy=$g.fgy; fgz=$g.fgz; cgx=$g.cgx; cgy=$g.cgy; cgz=$g.cgz; actual_spacing_max=[math]::Max($g.fgx/($nx-1),[math]::Max($g.fgy/($ny-1),$g.fgz/($nz-1))) }
}

function Invoke-External([string]$Exe,[string]$RunDir,[string]$Arguments,[int]$Timeout,[hashtable]$Environment=@{}) {
    $psi=[Diagnostics.ProcessStartInfo]::new(); $psi.FileName=$Exe; $psi.WorkingDirectory=$RunDir; $psi.Arguments=$Arguments
    $psi.UseShellExecute=$false; $psi.CreateNoWindow=$true; $psi.RedirectStandardOutput=$true; $psi.RedirectStandardError=$true
    foreach($key in $Environment.Keys){$psi.EnvironmentVariables[$key]=[string]$Environment[$key]}
    $p=[Diagnostics.Process]::new(); $p.StartInfo=$psi; $sw=[Diagnostics.Stopwatch]::StartNew(); $peak=0L; $timedOut=$false
    try {
        [void]$p.Start(); $outTask=$p.StandardOutput.ReadToEndAsync(); $errTask=$p.StandardError.ReadToEndAsync()
        $deadline=if($Timeout-gt0){[DateTime]::UtcNow.AddSeconds($Timeout)}else{[DateTime]::MaxValue}
        while(-not $p.HasExited -and [DateTime]::UtcNow-lt$deadline){try{$p.Refresh();$peak=[math]::Max($peak,[int64]$p.WorkingSet64);$peak=[math]::Max($peak,[int64]$p.PeakWorkingSet64)}catch{};Start-Sleep -Milliseconds 100}
        $timedOut=-not $p.HasExited
        if($timedOut){try{$p.Kill()}catch{};$p.WaitForExit()}
        try{$p.Refresh();$peak=[math]::Max($peak,[int64]$p.PeakWorkingSet64)}catch{}
        $outTask.Wait();$errTask.Wait()
        Set-Content (Join-Path $RunDir 'stdout.log') $outTask.Result -Encoding UTF8
        Set-Content (Join-Path $RunDir 'stderr.log') $errTask.Result -Encoding UTF8
        return [pscustomobject]@{exit_code=if($timedOut){124}else{$p.ExitCode};timed_out=$timedOut;wall_seconds=$sw.Elapsed.TotalSeconds;peak_rss_mb=$peak/1MB}
    } finally {
        if(-not $p.HasExited){try{$p.Kill()}catch{}}
        $sw.Stop();$p.Dispose()
    }
}

function Invoke-Fagbi27TitanCase([string]$Exe,[string]$RunDir,[string]$ArgPath,[int]$Timeout) {
    $psi=[Diagnostics.ProcessStartInfo]::new(); $psi.FileName=$Exe; $psi.WorkingDirectory=$RunDir
    $psi.Arguments='-a "' + ($ArgPath -replace '"','\"') + '"'
    $psi.UseShellExecute=$false; $psi.CreateNoWindow=$true; $psi.RedirectStandardOutput=$true; $psi.RedirectStandardError=$true
    $p=[Diagnostics.Process]::new(); $p.StartInfo=$psi
    $sw=[Diagnostics.Stopwatch]::StartNew(); $peak=0L; $timedOut=$false; $started=$false
    try {
        [void]$p.Start(); $started=$true
        $outTask=$p.StandardOutput.ReadToEndAsync(); $errTask=$p.StandardError.ReadToEndAsync()
        $deadline=if($Timeout-gt0){[DateTime]::UtcNow.AddSeconds($Timeout)}else{[DateTime]::MaxValue}
        while(-not $p.HasExited -and [DateTime]::UtcNow-lt$deadline){
            try{$p.Refresh();$peak=[math]::Max($peak,[int64]$p.WorkingSet64);$peak=[math]::Max($peak,[int64]$p.PeakWorkingSet64)}catch{}
            Start-Sleep -Milliseconds 100
        }
        $timedOut=-not $p.HasExited
        if($timedOut){try{$p.Kill()}catch{};$p.WaitForExit()}
        try{$p.Refresh();$peak=[math]::Max($peak,[int64]$p.PeakWorkingSet64)}catch{}
        $outTask.Wait();$errTask.Wait()
        Set-Content (Join-Path $RunDir 'stdout.log') $outTask.Result -Encoding UTF8
        Set-Content (Join-Path $RunDir 'stderr.log') $errTask.Result -Encoding UTF8
        return [pscustomobject]@{exit_code=if($timedOut){124}else{$p.ExitCode};timed_out=$timedOut;wall_seconds=$sw.Elapsed.TotalSeconds;peak_working_set_mb=$peak/1MB}
    } finally {
        if($started){
            try{if(-not $p.HasExited){$p.Kill();$p.WaitForExit(5000)|Out-Null}}catch{}
        }
        $sw.Stop();$p.Dispose()
    }
}

function Read-ApbsEnergy([string]$RunDir) {
    $text=((Get-Content (Join-Path $RunDir 'stdout.log') -Raw -ErrorAction SilentlyContinue)+"`n"+(Get-Content (Join-Path $RunDir 'stderr.log') -Raw -ErrorAction SilentlyContinue))
    $matches=[regex]::Matches($text,'Global net ELEC energy\s*=\s*([-+0-9.eE]+)',[Text.RegularExpressions.RegexOptions]::IgnoreCase)
    if($matches.Count -gt 0){return D $matches[$matches.Count-1].Groups[1].Value}
    return [double]::NaN
}

$Configs=[Collections.Generic.List[object]]::new()
foreach($spec in $SelectedSpecs){foreach($backend in $Backends){$Configs.Add([pscustomobject]@{case_id="$(Safe-Id $spec.pdb_id)_$backend";paper_index=$spec.paper_index;pdb_id=$spec.pdb_id;backend=$backend})}}
$Checkpoint=@{}
if(Test-Path $MatrixPath){foreach($row in @(Import-Csv $MatrixPath)){$Checkpoint[[string]$row.case_id]=$row}}
function Save-Matrix {
    @($Checkpoint.Values | Sort-Object {I $_.paper_index},backend) | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
}
function Get-ProfileValue($Profile,[string]$Name,$Default=''){if($null-eq$Profile){return $Default};$p=$Profile.PSObject.Properties[$Name];if($null-eq$p){return $Default};return $p.Value}
function New-Row($Config,[string]$Protocol,[string]$ExeSha,[string]$PqrSha,[string]$RunDir){
    [pscustomobject][ordered]@{case_id=$Config.case_id;paper_index=$Config.paper_index;pdb_id=$Config.pdb_id;backend=$Config.backend;protocol=$Protocol;executable_sha256=$ExeSha;pqr_sha256=$PqrSha;threads_requested=4;status='planned';failure='';exit_code='';timed_out='';atom_count=0;panels='';energy_kj_mol='';energy_kcal_mol='';iterations='';true_residual='';audit_rms='';process_wall_seconds='';solver_total_seconds='';surface_seconds='';ilu_setup_seconds='';gmres_seconds='';peak_rss_mb='';preconditioner_effective='';ilu_regularization_mode='';ilu_regularization_selected='';ilu_regularization_attempts='';ilu_regularization_rejected_attempts='';ilu_regularization_shifted_rows='';ilu_quality_status='';fmm_mode='';apbs_dime_x='';apbs_dime_y='';apbs_dime_z='';apbs_fglen_x='';apbs_fglen_y='';apbs_fglen_z='';apbs_actual_spacing_max_A='';run_dir=$RunDir}
}

function Should-Run($Config,[string]$Protocol,[string]$ExeSha,[string]$PqrSha){
    if($Force){return $true}
    if(-not $Checkpoint.ContainsKey($Config.case_id)){return $true}
    $old=$Checkpoint[$Config.case_id]
    $matches=[string]$old.protocol-eq$Protocol -and [string]$old.executable_sha256-eq$ExeSha -and [string]$old.pqr_sha256-eq$PqrSha
    if(-not $matches){return $true}
    if($Config.pdb_id -in $RerunFixtures -or $Config.backend -in $RerunBackends){return $true}
    if($RerunFailed -and [string]$old.status-ne'ok'){return $true}
    return $false
}

if(-not $PlotOnly){
foreach($config in $Configs){
    if(-not $FixtureRows.ContainsKey($config.pdb_id) -or [string]$FixtureRows[$config.pdb_id].status-ne'ok'){
        $runDir=Join-Path $WorkRoot $config.case_id;$row=New-Row $config '' '' '' $runDir;$row.status='fixture_gap';$row.failure=if($FixtureRows.ContainsKey($config.pdb_id)){[string]$FixtureRows[$config.pdb_id].failure}else{'missing manifest row'};$Checkpoint[$config.case_id]=$row;Save-Matrix;continue
    }
    $fixture=$FixtureRows[$config.pdb_id];$pqr=Join-Path $FixtureRoot ([string]$fixture.pqr_file);$pqrSha=Sha $pqr
    $protocol=if($config.backend-eq'titan'){$TitanProtocol}else{$ApbsProtocol};$exeSha=if($config.backend-eq'titan'){$TitanSha}else{$ApbsSha}
    if(-not (Should-Run $config $protocol $exeSha $pqrSha)){Write-Host "resume: $($config.case_id)";continue}
    $runDir=Join-Path $WorkRoot $config.case_id
    if(Test-Path $runDir){Remove-Item -LiteralPath $runDir -Recurse -Force}
    New-Item -ItemType Directory -Force -Path $runDir|Out-Null
    $localPqr=Join-Path $runDir (Split-Path -Leaf $pqr);Copy-Item $pqr $localPqr
    $row=New-Row $config $protocol $exeSha $pqrSha $runDir;$row.atom_count=I $fixture.atom_count
    if($config.backend-eq'titan'){
        $arg=Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir -PqrPath $localPqr -Backend tabi -Mode treecode -FmmMode single-matvec-analytic -FmmLocal off -SurfacePoints 8 -Threads 4 -MaxIter 720 -Restart 128 -Tolerance 1e-4 -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 -ProbeRadius 1.4 -OpeningTheta 0.45 -AuditMode sampled -AuditSampleCount 512 -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0 -PreconditionerIluRegularization $PreconditionerIluRegularization -PreconditionerIluRegularizationMax $PreconditionerIluRegularizationMax -NearIntegration analytic_singular -TaylorOrder 4 -PrimarySurface analytic-ses -SesGridSpacing 1.25 -DgSolv $true
        Write-Host "run: $($config.case_id) (TITAN, 4 cores, timeout ${TitanTimeoutSeconds}s)"
        $result=Invoke-Fagbi27TitanCase -Exe $ResolvedTitan -RunDir $runDir -ArgPath $arg -Timeout $TitanTimeoutSeconds
        $row.exit_code=$result.exit_code;$row.timed_out=$result.timed_out;$row.process_wall_seconds=$result.wall_seconds;$row.peak_rss_mb=$result.peak_working_set_mb
        if($result.timed_out){$row.status='timeout';$row.failure='titan_timeout'}elseif($result.exit_code-ne0){$row.status='failed';$row.failure='titan_exit_nonzero'}else{
            try{$profile=@(Import-Csv (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1];$energy=@(Import-Csv (Join-Path $runDir 'pbe_energy.csv'))[-1];$e=D $energy.E_rxn_kJmol;$row.energy_kj_mol=$e;$row.energy_kcal_mol=$e/4.184;$row.panels=I $profile.panels;$row.iterations=I $profile.iterations;$row.true_residual=D $profile.true_residual_rel;$row.audit_rms=D $profile.tree_vs_direct_matvec_rel;$row.solver_total_seconds=D $profile.total_seconds;$row.surface_seconds=D $profile.surface_build_seconds;$row.ilu_setup_seconds=D $profile.preconditioner_setup_seconds;$row.gmres_seconds=D $profile.gmres_seconds;$row.preconditioner_effective=[string]$profile.preconditioner_effective;$row.ilu_regularization_mode=Get-ProfileValue $profile 'preconditioner_ilu_regularization';$row.ilu_regularization_selected=Get-ProfileValue $profile 'preconditioner_ilu_regularization_selected';$row.ilu_regularization_attempts=Get-ProfileValue $profile 'preconditioner_ilu_regularization_attempts';$row.ilu_regularization_rejected_attempts=Get-ProfileValue $profile 'preconditioner_ilu_regularization_rejected_attempts';$row.ilu_regularization_shifted_rows=Get-ProfileValue $profile 'preconditioner_ilu_regularization_shifted_rows';$row.ilu_quality_status=Get-ProfileValue $profile 'preconditioner_quality_status';$row.fmm_mode=[string]$profile.fmm_m2l_mode;$row.status=if([string]$profile.status-eq'ok' -and (Is-Finite $e)){'ok'}else{'accuracy_gap'};$row.failure=if($row.status-eq'ok'){''}else{'titan_profile_or_energy_gap'}}catch{$row.status='schema_gap';$row.failure=$_.Exception.Message}
        }
    }else{
        Write-Host "run: $($config.case_id) (APBS, OMP threads=4, timeout ${ApbsTimeoutSeconds}s)"
        try{$grid=Write-ApbsInput -PqrPath $localPqr -InputPath (Join-Path $runDir 'apbs.in');$row.apbs_dime_x=$grid.nx;$row.apbs_dime_y=$grid.ny;$row.apbs_dime_z=$grid.nz;$row.apbs_fglen_x=$grid.fgx;$row.apbs_fglen_y=$grid.fgy;$row.apbs_fglen_z=$grid.fgz;$row.apbs_actual_spacing_max_A=$grid.actual_spacing_max
            $env4=@{OMP_NUM_THREADS='4';OMP_THREAD_LIMIT='4';OPENBLAS_NUM_THREADS='4';MKL_NUM_THREADS='4'}
            $result=Invoke-External -Exe $ResolvedApbs -RunDir $runDir -Arguments 'apbs.in' -Timeout $ApbsTimeoutSeconds -Environment $env4
            $row.exit_code=$result.exit_code;$row.timed_out=$result.timed_out;$row.process_wall_seconds=$result.wall_seconds;$row.peak_rss_mb=$result.peak_rss_mb
            if($result.timed_out){$row.status='timeout';$row.failure='apbs_timeout'}elseif($result.exit_code-ne0){$row.status='failed';$row.failure='apbs_exit_nonzero'}else{$e=Read-ApbsEnergy $runDir;$row.energy_kj_mol=$e;$row.energy_kcal_mol=$e/4.184;$row.status=if(Is-Finite $e){'ok'}else{'schema_gap'};$row.failure=if($row.status-eq'ok'){''}else{'apbs_energy_parse_gap'}}
        }catch{$row.status='failed';$row.failure=$_.Exception.Message}
    }
    $Checkpoint[$config.case_id]=$row;Save-Matrix
}
}

function Build-Pairs {
    $pairs=[Collections.Generic.List[object]]::new()
    foreach($spec in $PaperSpecs){
        $t=$Checkpoint["$(Safe-Id $spec.pdb_id)_titan"];$a=$Checkpoint["$(Safe-Id $spec.pdb_id)_apbs"]
        if($null-eq$t -and $null-eq$a){continue}
        $te=if($t){D $t.energy_kcal_mol}else{[double]::NaN};$ae=if($a){D $a.energy_kcal_mol}else{[double]::NaN}
        $diff=if((Is-Finite $te)-and(Is-Finite $ae)){$te-$ae}else{[double]::NaN};$pct=if((Is-Finite $diff)-and[math]::Abs($ae)-gt1e-12){100*$diff/[math]::Abs($ae)}else{[double]::NaN}
        $titanPass=$t -and [string]$t.status-eq'ok' -and (D $t.true_residual 1)-le1e-4 -and (D $t.audit_rms 1)-le8e-3
        $apbsPass=$a -and [string]$a.status-eq'ok' -and (Is-Finite $ae)
        $pairs.Add([pscustomobject][ordered]@{paper_index=$spec.paper_index;pdb_id=$spec.pdb_id;paper_elements=$spec.paper_elements;paper_energy_kcal_mol=$spec.paper_energy_kcal;atom_count=if($t){$t.atom_count}elseif($a){$a.atom_count}else{''};panels=if($t){$t.panels}else{''};titan_deltaG_solv_kcal_mol=$te;apbs_deltaG_solv_kcal_mol=$ae;titan_minus_apbs_kcal_mol=$diff;titan_minus_apbs_percent=$pct;titan_iterations=if($t){$t.iterations}else{''};titan_true_residual=if($t){$t.true_residual}else{''};titan_audit_rms=if($t){$t.audit_rms}else{''};titan_wall_seconds=if($t){$t.process_wall_seconds}else{''};apbs_wall_seconds=if($a){$a.process_wall_seconds}else{''};titan_peak_rss_mb=if($t){$t.peak_rss_mb}else{''};apbs_peak_rss_mb=if($a){$a.peak_rss_mb}else{''};apbs_grid=if($a){"$($a.apbs_dime_x)x$($a.apbs_dime_y)x$($a.apbs_dime_z)"}else{''};status=if($titanPass-and$apbsPass){'ok'}elseif(-not$titanPass){'titan_gap'}else{'apbs_gap'}})
    }
    return $pairs.ToArray()
}

$Pairs=@(Build-Pairs)
$Pairs|Export-Csv -LiteralPath $PairPath -NoTypeInformation
$Pairs|Select-Object paper_index,pdb_id,atom_count,panels,titan_deltaG_solv_kcal_mol,apbs_deltaG_solv_kcal_mol,titan_minus_apbs_percent,titan_iterations,titan_true_residual,titan_audit_rms,titan_wall_seconds,apbs_wall_seconds,titan_peak_rss_mb,apbs_peak_rss_mb,apbs_grid,status|Export-Csv -LiteralPath $SupportingTablePath -NoTypeInformation

function LogY([double]$v,[double]$min,[double]$max,[double]$top,[double]$bottom){$bottom-([math]::Log10($v)-[math]::Log10($min))/([math]::Log10($max)-[math]::Log10($min))*($bottom-$top)}
function Write-Figure([object[]]$Rows){
    $ok = @($Rows | Where-Object { [string]$_.status -eq 'ok' })
    if ($ok.Count -eq 0) { return }
    $w=1800;$h=650;$top=135;$bottom=555;$panels=@(@{l=110;r=575},@{l=690;r=1155},@{l=1270;r=1735});$blue='#0072B2';$orange='#D55E00';$gray='#444444'
    $svg=[Collections.Generic.List[string]]::new();$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="65mm" viewBox="0 0 1800 650"><rect width="100%" height="100%" fill="white"/><style>text{font-family:Arial;fill:#222}.axis{stroke:#222;stroke-width:2}.grid{stroke:#dedede;stroke-width:1}</style>')
    $svg.Add('<text x="100" y="108" font-size="34" font-weight="bold">A</text><text x="680" y="108" font-size="34" font-weight="bold">B</text><text x="1260" y="108" font-size="34" font-weight="bold">C</text>')
    $svg.Add("<line x1=`"690`" y1=`"42`" x2=`"745`" y2=`"42`" stroke=`"$blue`" stroke-width=`"5`"/><circle cx=`"718`" cy=`"42`" r=`"7`" fill=`"$blue`"/><text x=`"760`" y=`"51`" font-size=`"25`">TITAN analytic SES</text><line x1=`"1110`" y1=`"42`" x2=`"1165`" y2=`"42`" stroke=`"$orange`" stroke-width=`"5`"/><rect x=`"1138`" y=`"35`" width=`"14`" height=`"14`" fill=`"white`" stroke=`"$orange`" stroke-width=`"3`"/><text x=`"1180`" y=`"51`" font-size=`"25`">APBS grid</text>")
    foreach($p in $panels){$svg.Add("<line class=`"axis`" x1=`"$($p.l)`" y1=`"$top`" x2=`"$($p.l)`" y2=`"$bottom`"/><line class=`"axis`" x1=`"$($p.l)`" y1=`"$bottom`" x2=`"$($p.r)`" y2=`"$bottom`"/>")}
    $energies=@($ok|ForEach-Object{D $_.titan_deltaG_solv_kcal_mol;D $_.apbs_deltaG_solv_kcal_mol});$emin=($energies|Measure-Object -Minimum).Minimum;$emax=($energies|Measure-Object -Maximum).Maximum;$pad=0.06*($emax-$emin);$emin-=$pad;$emax+=$pad
    function EX([double]$v){$panels[0].l+($v-$emin)/($emax-$emin)*($panels[0].r-$panels[0].l)};function EY([double]$v){$bottom-($v-$emin)/($emax-$emin)*($bottom-$top)}
    for($i=0;$i-le4;$i++){$v=$emin+$i*($emax-$emin)/4;$x=EX $v;$y=EY $v;$svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/><line class=`"grid`" x1=`"$($panels[0].l)`" y1=`"$y`" x2=`"$($panels[0].r)`" y2=`"$y`"/><text x=`"$x`" y=`"590`" text-anchor=`"middle`" font-size=`"20`">$([math]::Round($v))</text><text x=`"$($panels[0].l-12)`" y=`"$($y+7)`" text-anchor=`"end`" font-size=`"20`">$([math]::Round($v))</text>")}
    $svg.Add("<line x1=`"$(EX $emin)`" y1=`"$(EY $emin)`" x2=`"$(EX $emax)`" y2=`"$(EY $emax)`" stroke=`"$gray`" stroke-width=`"2`" stroke-dasharray=`"7 6`"/>")
    foreach($r in $ok){$x=EX (D $r.apbs_deltaG_solv_kcal_mol);$y=EY (D $r.titan_deltaG_solv_kcal_mol);$svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"6`" fill=`"$blue`" fill-opacity=`"0.82`"/><text x=`"$($x+7)`" y=`"$($y-7)`" font-size=`"13`">$($r.paper_index)</text>")}
    $svg.Add('<text x="342" y="625" text-anchor="middle" font-size="24">APBS &#916;Gsolv (kcal mol&#8722;1)</text><text transform="translate(32,338) rotate(-90)" text-anchor="middle" font-size="24">TITAN &#916;Gsolv (kcal mol&#8722;1)</text>')
    foreach($metric in @(@{panel=1;t='titan_peak_rss_mb';a='apbs_peak_rss_mb';label='Peak RSS (MB, log scale)'},@{panel=2;t='titan_wall_seconds';a='apbs_wall_seconds';label='Total wall time (s, log scale)'})){
        $p=$panels[$metric.panel];$vals=@($ok|ForEach-Object{D $_.($metric.t);D $_.($metric.a)}|Where-Object{$_-gt0});$vmin=[math]::Max(1,($vals|Measure-Object -Minimum).Minimum*0.75);$vmax=($vals|Measure-Object -Maximum).Maximum*1.35
        function IX([double]$idx){$p.l+($idx-1)/26.0*($p.r-$p.l)}
        $ticks=@(1,5,10,15,20,25,27);foreach($idx in $ticks){$x=IX $idx;$svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/><text x=`"$x`" y=`"590`" text-anchor=`"middle`" font-size=`"20`">$idx</text>")}
        foreach($frac in 0,0.25,0.5,0.75,1){$v=[math]::Pow(10,[math]::Log10($vmin)+$frac*([math]::Log10($vmax)-[math]::Log10($vmin)));$y=LogY $v $vmin $vmax $top $bottom;$svg.Add("<line class=`"grid`" x1=`"$($p.l)`" y1=`"$y`" x2=`"$($p.r)`" y2=`"$y`"/><text x=`"$($p.l-10)`" y=`"$($y+7)`" text-anchor=`"end`" font-size=`"19`">$([math]::Round($v))</text>")}
        foreach($series in @(@{field=$metric.t;color=$blue;marker='circle'},@{field=$metric.a;color=$orange;marker='square'})){$pts=@();foreach($r in $ok|Sort-Object paper_index){$x=IX (D $r.paper_index);$y=LogY (D $r.($series.field)) $vmin $vmax $top $bottom;$pts+="$x,$y"};$svg.Add("<polyline points=`"$($pts-join' ')`" fill=`"none`" stroke=`"$($series.color)`" stroke-width=`"3`"/>");foreach($r in $ok){$x=IX (D $r.paper_index);$y=LogY (D $r.($series.field)) $vmin $vmax $top $bottom;if($series.marker-eq'circle'){$svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"5`" fill=`"$($series.color)`"/>")}else{$svg.Add("<rect x=`"$($x-5)`" y=`"$($y-5)`" width=`"10`" height=`"10`" fill=`"white`" stroke=`"$($series.color)`" stroke-width=`"2`"/>")}}}
        $cx=($p.l+$p.r)/2;$svg.Add("<text x=`"$cx`" y=`"625`" text-anchor=`"middle`" font-size=`"24`">Protein index</text><text transform=`"translate($($p.l-62),338) rotate(-90)`" text-anchor=`"middle`" font-size=`"24`">$(Xml $metric.label)</text>")
    }
    $svg.Add('</svg>');$svg|Set-Content -LiteralPath $FigurePath -Encoding UTF8
    $ink=Get-Command inkscape.exe -ErrorAction SilentlyContinue
    if($ink){
        $inkCli=Join-Path (Split-Path -Parent $ink.Source) 'inkscape.com'
        if(-not(Test-Path -LiteralPath $inkCli)){$inkCli=$ink.Source}
        & $inkCli $FigurePath "--export-filename=$([IO.Path]::ChangeExtension($FigurePath,'.png'))" '--export-dpi=600'|Out-Null
        if($LASTEXITCODE-ne0){throw 'Inkscape PNG export failed.'}
        & $inkCli $FigurePath "--export-filename=$([IO.Path]::ChangeExtension($FigurePath,'.pdf'))"|Out-Null
        if($LASTEXITCODE-ne0){throw 'Inkscape PDF export failed.'}
    }
}
Write-Figure $Pairs

$required=$SelectedSpecs.Count;$paired=@($Pairs|Where-Object{$_.pdb_id-in$Fixtures});$okPairs=@($paired|Where-Object { [string]$_.status -eq 'ok' })
$fixtureGaps=@($SelectedSpecs|Where-Object{-not$FixtureRows.ContainsKey($_.pdb_id)-or[string]$FixtureRows[$_.pdb_id].status-ne'ok'})
$matrixSelected=@($Checkpoint.Values|Where-Object{$_.pdb_id-in$Fixtures-and$_.backend-in$Backends})
$decision=if($fixtureGaps.Count-gt0){'fagbi27_fixture_preparation_gap'}elseif($matrixSelected.Count-lt($required*$Backends.Count)){'fagbi27_incomplete'}elseif(@($matrixSelected|Where-Object { [string]$_.status -eq 'timeout' }).Count-gt0){'fagbi27_runtime_limited'}elseif(@($matrixSelected|Where-Object{$_.backend-eq'titan'-and$_.status-ne'ok'}).Count-gt0){'fagbi27_titan_accuracy_gap'}elseif(@($matrixSelected|Where-Object{$_.backend-eq'apbs'-and$_.status-ne'ok'}).Count-gt0){'fagbi27_apbs_gap'}elseif($okPairs.Count-eq$required){if($Quick){'fagbi27_quick_ready'}else{'fagbi27_comparison_ready'}}else{'fagbi27_pairing_gap'}
[pscustomobject]@{decision=$decision;required_proteins=$required;paired_ok=$okPairs.Count;threads=4;matrix=$MatrixPath;supporting_table=$SupportingTablePath;figure=$FigurePath}|Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision","required_proteins=$required","paired_ok=$($okPairs.Count)","threads=4","matrix=$MatrixPath","supporting_table=$SupportingTablePath","figure=$FigurePath")|Set-Content $SummaryPath -Encoding ASCII
@('TITAN FAGBI-27 Analytic-SES versus APBS Report',"Decision: $decision","Proteins: $required; paired rows: $($okPairs.Count)",'Manuscript target: replace the current panel-count preconditioner plot (Figure 5); retain the dual-SES h-sensitivity Figure 4.',"TITAN: analytic SES h=1.25 A, degree-one analytical FMM/M2L, theta 0.45, FMM-near ILU(0) cap 64, ILU regularization $PreconditionerIluRegularization (max $(F $PreconditionerIluRegularizationMax)), tolerance 1e-4, four cores.",'APBS: same CHARMM PQR, mg-auto, smol/spl2, target grid spacing 0.50 A, four OpenMP threads requested.','The APBS comparison is contextual because the surface and volume discretizations differ; APBS is not a dense BIE oracle.',"Matrix: $MatrixPath","Supporting table: $SupportingTablePath","Figure: $FigurePath")|Set-Content $ReportPath -Encoding ASCII
Write-Host "Decision: $decision";Write-Host "Matrix: $MatrixPath";Write-Host "Supporting table: $SupportingTablePath";Write-Host "Figure: $FigurePath";Write-Host 'Resume: rerun the same command; add -RerunFailed for failed/timeout rows or use -RerunFixtures/-RerunBackends.'
if(-not$AllowFailedDecision-and$decision-notin@('fagbi27_comparison_ready','fagbi27_quick_ready')){throw "FAGBI-27 gate did not pass: $decision"}
