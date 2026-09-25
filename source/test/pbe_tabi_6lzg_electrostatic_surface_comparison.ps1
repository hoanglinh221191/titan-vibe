[CmdletBinding()]
param(
    [string]$TitanExe = 'C:\tmp\elec_vdw-main-sync\titan.exe',
    [string]$ApbsExe = '',
    [string]$PythonExe = '',
    [string]$FixtureRoot = 'test\pbe_inputs_public_6lzg',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_6lzg_electrostatic_surface',
    [ValidateSet('titan','apbs')][string[]]$Backends = @('titan','apbs'),
    [ValidateRange(1,64)][int]$Threads = 4,
    [ValidateRange(60,14400)][int]$TitanTimeoutSeconds = 3600,
    [ValidateRange(60,14400)][int]$ApbsTimeoutSeconds = 1800,
    [ValidateRange(0.25,1.0)][double]$TargetApbsGridSpacingA = 0.50,
    [ValidateRange(65,1025)][int]$ApbsMaxGridDimension = 385,
    [ValidateRange(0.5,20.0)][double]$ColorLimitKtE = 5.0,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$Value) { $Value.ToString('G15', $Inv) }
function Resolve-Exe([string]$Requested, [string[]]$Fallbacks) {
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
function Normalize-List([string[]]$Values) {
    return @($Values | ForEach-Object { @([string]$_ -split ',') } |
        ForEach-Object { $_.Trim().ToLowerInvariant() } | Where-Object { $_ } |
        Sort-Object -Unique)
}
function Get-PqrAtoms([string]$Path) {
    $atoms = [Collections.Generic.List[object]]::new()
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 10) { throw "Malformed PQR row: $line" }
        $atoms.Add([pscustomobject]@{
                x=[double]::Parse($parts[-5],$Inv); y=[double]::Parse($parts[-4],$Inv)
                z=[double]::Parse($parts[-3],$Inv); radius=[double]::Parse($parts[-1],$Inv)
            }) | Out-Null
    }
    return $atoms.ToArray()
}
function Get-ApbsGeometry([string]$Path) {
    $atoms = @(Get-PqrAtoms $Path)
    if ($atoms.Count -eq 0) { throw "No PQR atoms: $Path" }
    $xmin=1e300; $xmax=-1e300; $ymin=1e300; $ymax=-1e300; $zmin=1e300; $zmax=-1e300
    foreach ($atom in $atoms) {
        $r=[math]::Max(0.0,[double]$atom.radius)
        $xmin=[math]::Min($xmin,$atom.x-$r); $xmax=[math]::Max($xmax,$atom.x+$r)
        $ymin=[math]::Min($ymin,$atom.y-$r); $ymax=[math]::Max($ymax,$atom.y+$r)
        $zmin=[math]::Min($zmin,$atom.z-$r); $zmax=[math]::Max($zmax,$atom.z+$r)
    }
    $finePadding=12.0; $coarsePadding=16.0
    $fgx=($xmax-$xmin)+2*$finePadding; $fgy=($ymax-$ymin)+2*$finePadding
    $fgz=($zmax-$zmin)+2*$finePadding
    return [pscustomobject]@{ fgx=$fgx; fgy=$fgy; fgz=$fgz
        cgx=$fgx+2*$coarsePadding; cgy=$fgy+2*$coarsePadding; cgz=$fgz+2*$coarsePadding }
}
function Get-ApbsDime([double]$LengthA) {
    $required=[int][math]::Ceiling($LengthA/$TargetApbsGridSpacingA)+1
    $dimension=[int](32*[math]::Ceiling(([math]::Max(65,$required)-1)/32.0)+1)
    if ($dimension -gt $ApbsMaxGridDimension) {
        throw "Required APBS dimension $dimension exceeds cap $ApbsMaxGridDimension."
    }
    return $dimension
}
function Write-ApbsPotentialInput([string]$LocalPqr,[string]$InputPath) {
    $g=Get-ApbsGeometry $LocalPqr
    $nx=Get-ApbsDime $g.fgx; $ny=Get-ApbsDime $g.fgy; $nz=Get-ApbsDime $g.fgz
    $gridLines=@(
        "    dime $nx $ny $nz",
        "    cglen $(F $g.cgx) $(F $g.cgy) $(F $g.cgz)",
        "    fglen $(F $g.fgx) $(F $g.fgy) $(F $g.fgz)",
        '    cgcent mol 1','    fgcent mol 1','    mol 1','    lpbe','    bcfl sdh',
        '    pdie 1.0','    srfm smol','    chgm spl2','    sdens 10.0',
        '    srad 1.4','    swin 0.3','    temp 298.15',
        '    calcenergy total','    calcforce no'
    )
    @(
        'read',"    mol pqr $(Split-Path -Leaf $LocalPqr)",'end',
        'elec name solvated','    mg-auto'
    ) + $gridLines[0..8] + @(
        '    sdie 80.0','    ion charge 1 conc 0.15 radius 2.0',
        '    ion charge -1 conc 0.15 radius 2.0'
    ) + $gridLines[9..($gridLines.Count-1)] + @(
        '    write pot dx apbs_6lzg_potential','end',
        'elec name homogeneous_reference','    mg-auto'
    ) + $gridLines[0..8] + @(
        '    sdie 1.0'
    ) + $gridLines[9..($gridLines.Count-1)] + @(
        '    write pot dx apbs_6lzg_reference_potential','end','quit'
    ) | Set-Content -LiteralPath $InputPath -Encoding ASCII
    return [pscustomobject]@{ nx=$nx; ny=$ny; nz=$nz
        spacing_max_A=[math]::Max($g.fgx/($nx-1),[math]::Max($g.fgy/($ny-1),$g.fgz/($nz-1))) }
}
function Invoke-Case([string]$Exe,[string]$RunDir,[string]$Arguments,[int]$TimeoutSeconds) {
    $psi=[Diagnostics.ProcessStartInfo]::new(); $psi.FileName=$Exe
    $psi.WorkingDirectory=$RunDir; $psi.Arguments=$Arguments
    $psi.UseShellExecute=$false; $psi.CreateNoWindow=$true
    $psi.RedirectStandardOutput=$true; $psi.RedirectStandardError=$true
    $psi.EnvironmentVariables['OMP_NUM_THREADS']=[string]$Threads
    $psi.EnvironmentVariables['OPENBLAS_NUM_THREADS']=[string]$Threads
    $process=[Diagnostics.Process]::new(); $process.StartInfo=$psi
    $watch=[Diagnostics.Stopwatch]::StartNew(); $peak=0L
    try {
        [void]$process.Start(); $stdout=$process.StandardOutput.ReadToEndAsync()
        $stderr=$process.StandardError.ReadToEndAsync()
        $deadline=[DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
        while (-not $process.HasExited -and [DateTime]::UtcNow -lt $deadline) {
            try { $process.Refresh(); $peak=[math]::Max($peak,[int64]$process.PeakWorkingSet64) } catch {}
            Start-Sleep -Milliseconds 100
        }
        $timedOut=-not $process.HasExited
        if ($timedOut) { try { $process.Kill() } catch {}; $process.WaitForExit() }
        $stdout.Wait(); $stderr.Wait()
        Set-Content (Join-Path $RunDir 'stdout.log') $stdout.Result -Encoding UTF8
        Set-Content (Join-Path $RunDir 'stderr.log') $stderr.Result -Encoding UTF8
        return [pscustomobject]@{ exit_code=if($timedOut){124}else{$process.ExitCode}
            timed_out=$timedOut; wall_seconds=$watch.Elapsed.TotalSeconds; peak_rss_mb=$peak/1MB }
    } finally { $watch.Stop(); $process.Dispose() }
}

$Backends=@(Normalize-List $Backends)
if (@($Backends | Where-Object { $_ -notin @('titan','apbs') }).Count -gt 0) {
    throw 'Backends must contain only titan and/or apbs.'
}
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) { $FixtureRoot=Join-Path $RepoRoot $FixtureRoot }
if (-not [IO.Path]::IsPathRooted($WorkRoot)) { $WorkRoot=Join-Path $RepoRoot $WorkRoot }
$FixtureRoot=[IO.Path]::GetFullPath($FixtureRoot); $WorkRoot=[IO.Path]::GetFullPath($WorkRoot)
$pqr=Join-Path $FixtureRoot '6lzg_complex_charmm_ph7.pqr'
if (-not (Test-Path -LiteralPath $pqr)) {
    throw "Missing cached 6LZG complex PQR: $pqr. Run pbe_tabi_6lzg_ace2_apbs_comparison.ps1 -PrepareOnly first."
}
if ($Force -and (Test-Path $WorkRoot)) { Remove-Item $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$statusPath=Join-Path $WorkRoot 'pbe_tabi_6lzg_surface_status.csv'
$status=@(if(Test-Path $statusPath){Import-Csv $statusPath}else{@()})
function Save-Status([object]$Row) {
    $script:status=@($script:status | Where-Object { $_.backend -ne $Row.backend }) + @($Row)
    $script:status | Export-Csv $statusPath -NoTypeInformation
}
function Need-Run([string]$Backend,[string]$Artifact) {
    $old=@($script:status | Where-Object backend -eq $Backend | Select-Object -Last 1)
    if ($Force -or $old.Count -eq 0 -or -not (Test-Path $Artifact)) { return $true }
    if ($RerunFailed -and $old[0].status -ne 'ok') { return $true }
    return $false
}

$resolvedTitan=Resolve-Exe $TitanExe @('.\titan.exe')
$resolvedApbs=Resolve-Exe $ApbsExe @('apbs.exe','apbs',
    (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
    'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe')

if ('titan' -in $Backends) {
    if (-not $resolvedTitan) { throw "Unable to resolve TITAN executable: $TitanExe" }
    $runDir=Join-Path $WorkRoot 'titan'; New-Item -ItemType Directory $runDir -Force | Out-Null
    $vtkPath=Join-Path $runDir 'pbe_tabi_surface_solution.vtk'
    if (Need-Run 'titan' $vtkPath) {
        Write-Host "run: TITAN 6LZG analytic-SES node-patch surface (timeout ${TitanTimeoutSeconds}s)"
        Get-ChildItem $runDir -File -ErrorAction SilentlyContinue | Remove-Item -Force
        $argPath=Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
            -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
            -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1e-4 -Threads $Threads `
            -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 -TemperatureK 298.15 `
            -OpeningTheta 0.45 -AuditMode sampled -AuditSampleCount 512 `
            -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
            -NearIntegration analytic_singular -PrimarySurface analytic-ses -SesGridSpacing 1.25
        Add-Content $argPath @('-pbe-tabi-discretization=node-patch',
            '-pbe-tabi-surface-output=all','-pbe-tabi-precond-memory-budget-mb=1024') -Encoding ASCII
        $run=Invoke-FmmCliTitanCase -TitanExe $resolvedTitan -RunDir $runDir `
            -ArgPath $argPath -CaseTimeoutSeconds $TitanTimeoutSeconds
        $row=[ordered]@{backend='titan';status=if($run.exit_code -eq 0){'ok'}else{'failed'}
            exit_code=$run.exit_code;wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb;detail=''}
        if ($run.exit_code -eq 0) {
            $profile=@(Import-Csv (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
            $csv=@(Import-Csv (Join-Path $runDir 'pbe_tabi_surface_solution.csv'))
            $vtk=Get-Content $vtkPath -Raw
            if ($profile.surface_generator -ne 'titan_analytic_patch_ray_ses' -or
                $profile.discretization -ne 'node_patch' -or
                [int64]$profile.panels -ne [int64]$profile.mesh_vertices -or
                $csv[0].data_association -ne 'point' -or
                -not ($csv[0].PSObject.Properties.Name -contains 'raw_boundary_unknown') -or
                -not ($csv[0].PSObject.Properties.Name -contains 'total_potential_kT_e') -or
                -not ($csv[0].PSObject.Properties.Name -contains 'coulomb_potential_kT_e') -or
                -not ($csv[0].PSObject.Properties.Name -contains 'reaction_potential_kT_e') -or
                $vtk -notmatch "POINT_DATA $($profile.panels)" -or
                $vtk -notmatch 'SCALARS total_potential_kT_e double 1' -or
                $vtk -notmatch 'SCALARS reaction_potential_kT_e double 1') {
                $row.status='failed'; $row.detail='binary_did_not_emit_connected_analytic_ses_node_patch_output'
            } else { $row.detail="points=$($profile.panels);triangles=$($profile.mesh_triangles)" }
        }
        Save-Status ([pscustomobject]$row)
    } else { Write-Host 'resume: TITAN surface output' }
}

if ('apbs' -in $Backends) {
    if (-not $resolvedApbs) { throw 'Unable to resolve APBS executable.' }
    $runDir=Join-Path $WorkRoot 'apbs'; New-Item -ItemType Directory $runDir -Force | Out-Null
    $dxPath=Join-Path $runDir 'apbs_6lzg_potential.dx'
    $referenceDxPath=Join-Path $runDir 'apbs_6lzg_reference_potential.dx'
    if ((Need-Run 'apbs' $dxPath) -or -not (Test-Path $referenceDxPath)) {
        Write-Host "run: APBS 6LZG OpenDX potential (timeout ${ApbsTimeoutSeconds}s)"
        Get-ChildItem $runDir -File -ErrorAction SilentlyContinue | Remove-Item -Force
        $localPqr=Join-Path $runDir '6lzg_complex_charmm_ph7.pqr'; Copy-Item $pqr $localPqr
        $input=Join-Path $runDir 'apbs_6lzg_potential.in'; $grid=Write-ApbsPotentialInput $localPqr $input
        $run=Invoke-Case $resolvedApbs $runDir ('"'+$input+'"') $ApbsTimeoutSeconds
        $ok=$run.exit_code -eq 0 -and (Test-Path $dxPath) -and
            (Test-Path $referenceDxPath)
        Save-Status ([pscustomobject][ordered]@{backend='apbs';status=if($ok){'ok'}else{'failed'}
            exit_code=$run.exit_code;wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_rss_mb
            detail="dime=$($grid.nx)x$($grid.ny)x$($grid.nz);spacing_max_A=$(F $grid.spacing_max_A)"})
    } else { Write-Host 'resume: APBS OpenDX output' }
}

$titanVtk=Join-Path $WorkRoot 'titan\pbe_tabi_surface_solution.vtk'
$apbsDx=Join-Path $WorkRoot 'apbs\apbs_6lzg_potential.dx'
$apbsReferenceDx=Join-Path $WorkRoot 'apbs\apbs_6lzg_reference_potential.dx'
$ready=@($status | Where-Object status -eq 'ok' | Select-Object -ExpandProperty backend)
if ('titan' -in $ready -and 'apbs' -in $ready -and
    (Test-Path $titanVtk) -and (Test-Path $apbsDx) -and
    (Test-Path $apbsReferenceDx)) {
    if (-not $PythonExe) {
        $PythonExe=Resolve-Exe '' @('python.exe','python',
            'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe')
    } else { $PythonExe=Resolve-Exe $PythonExe @() }
    if (-not $PythonExe) { throw 'Unable to resolve Python for matched-surface asset generation.' }
    $assetDir=Join-Path $WorkRoot 'comparison'
    & $PythonExe (Join-Path $PSScriptRoot 'pbe_tabi_6lzg_electrostatic_surface_assets.py') `
        --titan-vtk $titanVtk --apbs-dx $apbsDx --apbs-reference-dx $apbsReferenceDx `
        --pqr $pqr --output-dir $assetDir `
        --color-limit-kt-e (F $ColorLimitKtE)
    if ($LASTEXITCODE -ne 0) { throw 'Matched TITAN/APBS surface asset generation failed.' }
}

$failed=@($status | Where-Object status -ne 'ok')
$haveBoth=@($status | Where-Object status -eq 'ok' |
    Select-Object -ExpandProperty backend -Unique)
$decision=if($failed.Count -eq 0 -and 'titan' -in $haveBoth -and 'apbs' -in $haveBoth -and
    (Test-Path $titanVtk) -and (Test-Path $apbsDx) -and
    (Test-Path $apbsReferenceDx)) {
    'surface_comparison_ready'
} else { 'surface_comparison_incomplete' }
Set-Content (Join-Path $WorkRoot 'pbe_tabi_6lzg_surface_decision.dat') "decision $decision" -Encoding ASCII
Write-Host "Decision: $decision"
$pymolPath=Join-Path $WorkRoot 'comparison\view_6lzg_titan_vs_apbs.py'
if (Test-Path $pymolPath) { Write-Host "PyMOL: $pymolPath" }
if (-not $AllowFailedDecision -and $decision -ne 'surface_comparison_ready') { throw $decision }
