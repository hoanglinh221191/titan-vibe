param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_grid_1a63_same_pqr_gate'),
    [string]$PqrPath = (Join-Path $PSScriptRoot 'pbe_inputs_public\1a63_full_amber_ph7.pqr'),
    [double[]]$GridSpacings = @(1.0, 0.75, 0.5),
    [double]$TabiSurfaceSpacing = 1.25,
    [int]$Threads = 4,
    [string]$ApbsExe = '',
    [int]$TimeoutSeconds = 1800,
    [double]$MaxFinestRelativeDifference = 0.15,
    [double]$MaxAdjacentGridDifference = 0.20
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Culture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_grid_fixture_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$Value) { return $Value.ToString('G15', $Culture) }
function D([object]$Value) { return [double]::Parse([string]$Value, $Culture) }
function Is-Finite([double]$Value) {
    return -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}
function To-TitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function Safe-H([double]$Value) { return (F $Value).Replace('.', 'p') }

function Resolve-Apbs([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        if (Test-Path -LiteralPath $Requested -PathType Leaf) {
            return (Resolve-Path -LiteralPath $Requested).Path
        }
        $command = Get-Command $Requested -ErrorAction SilentlyContinue
        if ($null -ne $command) { return $command.Source }
        throw "Unable to resolve APBS executable: $Requested"
    }
    foreach ($name in @('apbs.exe', 'apbs')) {
        $command = Get-Command $name -ErrorAction SilentlyContinue
        if ($null -ne $command) { return $command.Source }
    }
    foreach ($candidate in @(
        (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
        'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe'
    )) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return ''
}

function Get-PqrBounds([string]$Path) {
    $atoms = @()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 10) { continue }
        $atoms += [pscustomobject]@{ x=D $parts[5]; y=D $parts[6]; z=D $parts[7] }
    }
    if ($atoms.Count -eq 0) { throw "No atoms found in PQR: $Path" }
    return [pscustomobject]@{
        x = (($atoms | Measure-Object x -Maximum).Maximum - ($atoms | Measure-Object x -Minimum).Minimum)
        y = (($atoms | Measure-Object y -Maximum).Maximum - ($atoms | Measure-Object y -Minimum).Minimum)
        z = (($atoms | Measure-Object z -Maximum).Maximum - ($atoms | Measure-Object z -Minimum).Minimum)
    }
}

function Invoke-ApbsReference([string]$Executable, [string]$RunDir,
                              [string]$SourcePqr, [int]$CaseTimeoutSeconds) {
    New-Item -ItemType Directory -Path $RunDir -Force | Out-Null
    $localPqr = Join-Path $RunDir (Split-Path -Leaf $SourcePqr)
    Copy-Item -LiteralPath $SourcePqr -Destination $localPqr -Force
    $bounds = Get-PqrBounds $SourcePqr
    $cg = @([math]::Max(30.0, $bounds.x + 32.0),
            [math]::Max(30.0, $bounds.y + 32.0),
            [math]::Max(30.0, $bounds.z + 32.0))
    $fg = @([math]::Max(20.0, $bounds.x + 16.0),
            [math]::Max(20.0, $bounds.y + 16.0),
            [math]::Max(20.0, $bounds.z + 16.0))
    $common = @(
        '    mg-auto', '    dime 161 161 161',
        "    cglen $(F $cg[0]) $(F $cg[1]) $(F $cg[2])",
        "    fglen $(F $fg[0]) $(F $fg[1]) $(F $fg[2])",
        '    cgcent mol 1', '    fgcent mol 1', '    mol 1', '    lpbe',
        '    bcfl sdh', '    pdie 1.0', '    srfm smol', '    chgm spl2',
        '    sdens 10.0', '    srad 1.4', '    swin 0.3', '    temp 298.15',
        '    calcenergy total', '    calcforce no'
    )
    $input = Join-Path $RunDir 'apbs.in'
    Set-Content -LiteralPath $input -Encoding ASCII -Value (@(
        'read', "    mol pqr $(Split-Path -Leaf $localPqr)", 'end',
        'elec name solvated'
    ) + $common + @(
        '    ion charge 1 conc 0.15 radius 2.0',
        '    ion charge -1 conc 0.15 radius 2.0', '    sdie 80.0', 'end',
        'elec name reference'
    ) + $common + @('    sdie 1.0', 'end',
        'print elecEnergy solvated - reference end', 'quit'))

    $stdout = Join-Path $RunDir 'stdout.log'
    $stderr = Join-Path $RunDir 'stderr.log'
    $process = Start-Process -FilePath $Executable -ArgumentList 'apbs.in' `
        -WorkingDirectory $RunDir -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
        Stop-Process -Id $process.Id -Force -ErrorAction SilentlyContinue
        return [pscustomobject]@{ status='timeout'; energy=[double]::NaN }
    }
    $process.Refresh()
    $text = (Get-Content -LiteralPath $stdout -Raw) + "`n" +
            (Get-Content -LiteralPath $stderr -Raw)
    $matches = [regex]::Matches(
        $text, 'Global net ELEC energy\s*=\s*([-+0-9.eE]+)\s*kJ/mol',
        [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    if ($matches.Count -eq 0) {
        $status = if ($process.ExitCode -eq 0) { 'parse_gap' } else { 'failed' }
        return [pscustomobject]@{ status=$status; energy=[double]::NaN }
    }
    return [pscustomobject]@{
        status='ok'; energy=D $matches[$matches.Count - 1].Groups[1].Value
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $PqrPath)) { throw "Missing 1A63 PQR: $PqrPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$PqrPath = (Resolve-Path -LiteralPath $PqrPath).Path
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$workDriveRoot = [IO.Path]::GetPathRoot($WorkRoot).TrimEnd('\')
if ($WorkRoot.TrimEnd('\') -eq $workDriveRoot -or
    $WorkRoot.TrimEnd('\') -eq [IO.Path]::GetFullPath($RepoRoot).TrimEnd('\')) {
    throw "Refusing to use an unsafe 1A63 work root: $WorkRoot"
}
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$ffRoot = Join-Path $WorkRoot 'lib_ff'
Initialize-PbeGridForceField -ForceFieldRoot $ffRoot
$pqrHash = (Get-FileHash -LiteralPath $PqrPath -Algorithm SHA256).Hash.ToLowerInvariant()
$rows = [Collections.Generic.List[object]]::new()

$tabiDir = Join-Path $WorkRoot 'tabi_p1'
New-Item -ItemType Directory -Path $tabiDir -Force | Out-Null
$tabiArg = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $tabiDir `
    -PqrPath $PqrPath -ForceFieldRoot $ffRoot -Backend tabi -Mode treecode `
    -FmmMode single-matvec-analytic -FmmLocal off -SurfacePoints 8 `
    -Threads $Threads -MaxIter 720 -Restart 128 -Tolerance 1e-4 -SaltM 0.15 `
    -InternalDielectric 1.0 -ExternalDielectric 80.0 -TemperatureK 298.15 `
    -ProbeRadius 1.4 -OpeningTheta 0.45 -AuditMode sampled `
    -AuditSampleCount 512 -Preconditioner fmm-nearfield-ilu0 `
    -PreconditionerIluNeighborCap 64 -PreconditionerIluRegularization adaptive `
    -PreconditionerIluRegularizationMax 0.01 -NearIntegration analytic_singular `
    -TaylorOrder 4 -PrimarySurface analytic-ses -SesGridSpacing $TabiSurfaceSpacing `
    -DgSolv $true
Add-Content -LiteralPath $tabiArg -Encoding ASCII -Value '-pbe-tabi-discretization=node-patch'
$tabiRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $tabiDir `
    -ArgPath $tabiArg -CaseTimeoutSeconds $TimeoutSeconds
if ($tabiRun.exit_code -ne 0 -or $tabiRun.timed_out) { throw '1A63 TABI p1 reference failed.' }
$tabiEnergyRow = @(Import-Csv (Join-Path $tabiDir 'pbe_energy.csv'))[-1]
$tabiProfile = @(Import-Csv (Join-Path $tabiDir 'pbe_tabi_profile.csv'))[-1]
$tabiEnergy = D $tabiEnergyRow.deltaG_solv_kJmol
if ($tabiEnergyRow.dgsolv_status -ne 'ok' -or -not (Is-Finite $tabiEnergy)) {
    throw '1A63 TABI p1 reference did not produce a converged finite dgsolv.'
}
$rows.Add([pscustomobject][ordered]@{
    backend='tabi_p1'; spacing_A=F $TabiSurfaceSpacing; energy_kJmol=F $tabiEnergy
    relative_vs_tabi=F 0.0; iterations=$tabiProfile.iterations
    residual=$tabiProfile.true_residual_rel; solute_points=$tabiProfile.panels
    status='ok'; pqr_sha256=$pqrHash
}) | Out-Null

$gridResults = @()
foreach ($spacing in @($GridSpacings | Sort-Object -Descending)) {
    $runDir = Join-Path $WorkRoot "grid_h$(Safe-H $spacing)"
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $driver = Join-Path $runDir 'driver.pdb'
    Write-FmmCliDummyPdb $driver
    $lib = To-TitanPath ((Resolve-Path -LiteralPath $ffRoot).Path.TrimEnd('\') + '\')
    $arg = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $arg -Encoding ASCII -Value @(
        "-f=$(To-TitanPath $driver)", '-1pdb=true', '-ff=AMBER99SB', '-g=1&1',
        "-path-a=$lib", "-path-n=$lib", '-pbe=grid',
        "-pbe-test-pqr=$(To-TitanPath $PqrPath)", '-pbe-dgsolv=true',
        '-pbe-model=linear', '-pbe-biomodel=pbe', '-pbe-dielectric=two-dielectric',
        '-pbe-surface=delphi-ses', '-pbe-charge-assignment=delphi',
        '-pbe-radius-convention=auto', '-pbe-boundary=coulombic',
        '-pbe-ion-access=dielectric', '-pbe-probe=1.4',
        "-pbe-scale=$(F (1.0 / $spacing))", '-pbe-margin=15',
        '-pbe-max-grid=6000000', '-pbe-max-iter=1200', '-pbe-rms=1e-6',
        '-pbe-dphi=1e-5', '-pbe-int-diel=1', '-pbe-ext-diel=80',
        '-pbe-salt=0.15', '-pbe-temp=298.15', '-pbe-focus=off',
        '-pbe-multigrid=off', '-pbe-map=off', "-nt=$Threads"
    )
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $TimeoutSeconds
    if ($run.exit_code -ne 0 -or $run.timed_out) { throw "1A63 Grid h=$spacing failed." }
    $energyRow = @(Import-Csv (Join-Path $runDir 'pbe_energy.csv'))[-1]
    $energy = D $energyRow.deltaG_solv_kJmol
    if ($energyRow.dgsolv_status -ne 'ok' -or -not (Is-Finite $energy)) {
        throw "1A63 Grid h=$spacing did not converge to finite dgsolv."
    }
    $relative = [math]::Abs($energy - $tabiEnergy) / [math]::Max(1.0, [math]::Abs($tabiEnergy))
    $record = [pscustomobject][ordered]@{
        backend='grid'; spacing_A=F $spacing; energy_kJmol=F $energy
        relative_vs_tabi=F $relative; iterations="$($energyRow.iterations_solvent)/$($energyRow.iterations_vacuum)"
        residual=$energyRow.linear_residual; solute_points=$energyRow.solute_points
        status='ok'; pqr_sha256=$pqrHash
    }
    $gridResults += $record
    $rows.Add($record) | Out-Null
}

$resolvedApbs = Resolve-Apbs $ApbsExe
if ([string]::IsNullOrWhiteSpace($resolvedApbs)) {
    $rows.Add([pscustomobject][ordered]@{
        backend='apbs'; spacing_A=''; energy_kJmol=''; relative_vs_tabi=''
        iterations=''; residual=''; solute_points=''; status='skipped_missing'
        pqr_sha256=$pqrHash
    }) | Out-Null
} else {
    $apbs = Invoke-ApbsReference -Executable $resolvedApbs `
        -RunDir (Join-Path $WorkRoot 'apbs') -SourcePqr $PqrPath `
        -CaseTimeoutSeconds $TimeoutSeconds
    $relative = if (Is-Finite $apbs.energy) {
        [math]::Abs($apbs.energy - $tabiEnergy) / [math]::Max(1.0, [math]::Abs($tabiEnergy))
    } else { [double]::NaN }
    $rows.Add([pscustomobject][ordered]@{
        backend='apbs'; spacing_A='auto'; energy_kJmol=F $apbs.energy
        relative_vs_tabi=F $relative; iterations=''; residual=''; solute_points=''
        status=$apbs.status; pqr_sha256=$pqrHash
    }) | Out-Null
}

$rows | Export-Csv -LiteralPath (Join-Path $WorkRoot 'comparison.csv') -NoTypeInformation
$orderedGrid = @($gridResults | Sort-Object { D $_.spacing_A })
$finest = $orderedGrid[0]
$adjacentMax = 0.0
for ($i = 1; $i -lt $orderedGrid.Count; ++$i) {
    $coarse = D $orderedGrid[$i].energy_kJmol
    $fine = D $orderedGrid[$i - 1].energy_kJmol
    $adjacentMax = [math]::Max($adjacentMax,
        [math]::Abs($coarse - $fine) / [math]::Max(1.0, [math]::Abs($fine)))
}
$finestRelative = D $finest.relative_vs_tabi
$decision = if ($finestRelative -le $MaxFinestRelativeDifference -and
                $adjacentMax -le $MaxAdjacentGridDifference) {
    'grid_1a63_same_pqr_ready'
} else { 'grid_1a63_same_pqr_accuracy_gap' }
Set-Content -LiteralPath (Join-Path $WorkRoot 'decision.dat') -Encoding ASCII -Value @(
    "decision $decision", "pqr_sha256 $pqrHash",
    "tabi_energy_kJmol $(F $tabiEnergy)",
    "grid_finest_spacing_A $($finest.spacing_A)",
    "grid_finest_energy_kJmol $($finest.energy_kJmol)",
    "grid_finest_relative_vs_tabi $(F $finestRelative)",
    "grid_adjacent_max_relative_difference $(F $adjacentMax)",
    "apbs_status $($rows[-1].status)"
)
Write-Host "PBE 1A63 same-PQR decision: $decision"
if ($decision -ne 'grid_1a63_same_pqr_ready') { throw $decision }
