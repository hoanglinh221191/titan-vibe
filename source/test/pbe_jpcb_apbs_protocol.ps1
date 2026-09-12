param(
    [string]$TitanExe = '.\titan.exe',
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_binding_jpcb'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\jpcb_apbs_protocol'),
    [string]$TitanWorkRoot = (Join-Path $PSScriptRoot 'tmp\jpcb_binding_salt'),
    [string]$Salts = '0.05,0.15,0.30,0.50',
    [string]$Dimes = '129,161,193,225',
    [ValidateRange(1,64)][int]$TitanThreads = 4,
    [ValidateRange(1,64)][int]$ApbsThreads = 4,
    [string]$WslDistribution = 'Ubuntu',
    [string]$ApbsWslExe = '/home/linh/.local/apbs-3.4.1-openmp/bin/apbs',
    [ValidateRange(60,7200)][int]$TitanTimeoutSeconds = 3600,
    [ValidateRange(60,7200)][int]$ApbsTimeoutSeconds = 1800,
    [ValidateRange(0.001,0.10)][double]$GridStabilityRelativeThreshold = 0.02,
    [ValidateRange(0.10,0.90)][double]$FineScale = 0.55,
    [ValidateRange(8.0,100.0)][double]$CoarsePaddingAngstrom = 32.0,
    [switch]$RunTitan,
    [switch]$Quick,
    [switch]$ContractOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$culture = [Globalization.CultureInfo]::InvariantCulture
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TitanWorkRoot = [IO.Path]::GetFullPath($TitanWorkRoot)
$TitanExe = [IO.Path]::GetFullPath($TitanExe)

function Convert-JpcbApbsDouble([object]$Value, [double]$Default = [double]::NaN) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $Default }
    return [double]::Parse($text, [Globalization.CultureInfo]::InvariantCulture)
}

function Test-JpcbApbsFinite([double]$Value) {
    return -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}

function Convert-JpcbApbsList([string]$Text, [string]$Kind) {
    $values = New-Object System.Collections.Generic.List[object]
    foreach ($token in ($Text -split ',')) {
        $trim = $token.Trim()
        if ([string]::IsNullOrWhiteSpace($trim)) { continue }
        if ($Kind -eq 'salt') {
            $value = 0.0
            if (-not [double]::TryParse($trim, [Globalization.NumberStyles]::Float,
                    $culture, [ref]$value) -or $value -lt 0.0 -or $value -gt 1.0) {
                throw "Invalid ionic strength: $trim"
            }
            $values.Add($value) | Out-Null
        } else {
            $value = 0
            if (-not [int]::TryParse($trim, [ref]$value) -or $value -lt 33) {
                throw "Invalid APBS dime: $trim"
            }
            $values.Add($value) | Out-Null
        }
    }
    return @($values.ToArray() | Sort-Object -Unique)
}

function Get-JpcbApbsAtoms([string]$PqrPath) {
    $atoms = New-Object System.Collections.Generic.List[object]
    foreach ($line in Get-Content -LiteralPath $PqrPath) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($parts.Count -lt 10) { continue }
        $atoms.Add([pscustomobject]@{
            x = Convert-JpcbApbsDouble $parts[5]
            y = Convert-JpcbApbsDouble $parts[6]
            z = Convert-JpcbApbsDouble $parts[7]
        }) | Out-Null
    }
    if ($atoms.Count -eq 0) { throw "No PQR atoms found: $PqrPath" }
    return $atoms.ToArray()
}

function Get-JpcbApbsBox([string]$PqrPath, [double]$Padding) {
    $atoms = @(Get-JpcbApbsAtoms $PqrPath)
    $xmin = ($atoms | Measure-Object x -Minimum).Minimum
    $xmax = ($atoms | Measure-Object x -Maximum).Maximum
    $ymin = ($atoms | Measure-Object y -Minimum).Minimum
    $ymax = ($atoms | Measure-Object y -Maximum).Maximum
    $zmin = ($atoms | Measure-Object z -Minimum).Minimum
    $zmax = ($atoms | Measure-Object z -Maximum).Maximum
    return [pscustomobject]@{
        x = [Math]::Max(20.0, ($xmax - $xmin) + $Padding)
        y = [Math]::Max(20.0, ($ymax - $ymin) + $Padding)
        z = [Math]::Max(20.0, ($zmax - $zmin) + $Padding)
        atoms = $atoms.Count
    }
}

function Write-JpcbApbsInput(
    [string]$PqrPath,
    [string]$InputPath,
    [int]$Dime,
    [double]$SaltM
) {
    $box = Get-JpcbApbsBox -PqrPath $PqrPath -Padding $CoarsePaddingAngstrom
    $fmt = { param([double]$v) $v.ToString('G12', $culture) }
    $fineX = [Math]::Max(12.0, $box.x * $FineScale)
    $fineY = [Math]::Max(12.0, $box.y * $FineScale)
    $fineZ = [Math]::Max(12.0, $box.z * $FineScale)
    $pqrName = Split-Path -Leaf $PqrPath
    $dimeLine = "    dime $Dime $Dime $Dime"
    $lengths = @(
        "    cglen $(& $fmt $box.x) $(& $fmt $box.y) $(& $fmt $box.z)",
        "    fglen $(& $fmt $fineX) $(& $fmt $fineY) $(& $fmt $fineZ)",
        '    cgcent mol 1',
        '    fgcent mol 1',
        '    mol 1',
        '    lpbe',
        '    bcfl sdh',
        "    ion charge 1 conc $(& $fmt $SaltM) radius 2.0",
        "    ion charge -1 conc $(& $fmt $SaltM) radius 2.0",
        '    pdie 1.0',
        '    sdie 80.0',
        '    srfm smol',
        '    chgm spl2',
        '    sdens 10.0',
        '    srad 1.4',
        '    swin 0.3',
        '    temp 298.15',
        '    calcenergy total',
        '    calcforce no'
    )
    $reference = @(
        '    mg-auto',
        $dimeLine,
        "    cglen $(& $fmt $box.x) $(& $fmt $box.y) $(& $fmt $box.z)",
        "    fglen $(& $fmt $fineX) $(& $fmt $fineY) $(& $fmt $fineZ)",
        '    cgcent mol 1',
        '    fgcent mol 1',
        '    mol 1',
        '    lpbe',
        '    bcfl sdh',
        '    pdie 1.0',
        '    sdie 1.0',
        '    srfm smol',
        '    chgm spl2',
        '    sdens 10.0',
        '    srad 1.4',
        '    swin 0.3',
        '    temp 298.15',
        '    calcenergy total',
        '    calcforce no'
    )
    $lines = @(
        'read', "    mol pqr $pqrName", 'end',
        'elec name solvated', '    mg-auto', $dimeLine
    ) + $lengths + @(
        'end', 'elec name reference'
    ) + $reference + @(
        'end', 'print elecEnergy solvated - reference end', 'quit'
    )
    Set-Content -LiteralPath $InputPath -Value $lines -Encoding ASCII
    return $box
}

function Convert-JpcbApbsWslPath([string]$Path) {
    return (& wsl.exe -d $WslDistribution -- wslpath -a ($Path -replace '\\','/')).Trim()
}

function Invoke-JpcbApbsWslCase([string]$RunDir, [int]$TimeoutSeconds) {
    $wslRunDir = Convert-JpcbApbsWslPath $RunDir
    if ([string]::IsNullOrWhiteSpace($wslRunDir)) { throw "Cannot map run directory into WSL: $RunDir" }
    $runnerPath = Join-Path $RunDir 'run_apbs_wsl.sh'
    $runner = @(
        '#!/usr/bin/env bash',
        'set -euo pipefail',
        "cd '$wslRunDir'",
        "export OMP_NUM_THREADS=$ApbsThreads",
        'export OPENBLAS_NUM_THREADS=1',
        'export OMP_PROC_BIND=true',
        'export OMP_PLACES=cores',
        "/usr/bin/time -v -o apbs_resource.log '$ApbsWslExe' apbs.in > stdout.log 2> stderr.log"
    )
    # Bash under WSL rejects CRLF after shell options, so keep this transport
    # wrapper explicitly LF-only even when the harness runs from Windows.
    [IO.File]::WriteAllText($runnerPath, ($runner -join "`n") + "`n",
        [Text.UTF8Encoding]::new($false))
    $wslRunnerPath = Convert-JpcbApbsWslPath $runnerPath
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = 'wsl.exe'
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = "-d $WslDistribution -- bash $wslRunnerPath"
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $psi
    $watch = [Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdout = $process.StandardOutput.ReadToEndAsync()
    $stderr = $process.StandardError.ReadToEndAsync()
    $deadline = [DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
    while (-not $process.HasExited -and [DateTime]::UtcNow -lt $deadline) {
        Start-Sleep -Milliseconds 100
    }
    $timedOut = -not $process.HasExited
    if ($timedOut) {
        try { $process.Kill() } catch {}
        $process.WaitForExit()
    }
    $stdout.Wait(); $stderr.Wait(); $watch.Stop()
    if (-not (Test-Path -LiteralPath (Join-Path $RunDir 'stdout.log'))) {
        Set-Content -LiteralPath (Join-Path $RunDir 'stdout.log') -Value $stdout.Result -Encoding UTF8
    }
    if (-not (Test-Path -LiteralPath (Join-Path $RunDir 'stderr.log'))) {
        Set-Content -LiteralPath (Join-Path $RunDir 'stderr.log') -Value $stderr.Result -Encoding UTF8
    }
    $rssMb = [double]::NaN
    $resourcePath = Join-Path $RunDir 'apbs_resource.log'
    if (Test-Path -LiteralPath $resourcePath) {
        $match = [regex]::Match((Get-Content -LiteralPath $resourcePath -Raw),
            'Maximum resident set size \(kbytes\):\s*([0-9]+)')
        if ($match.Success) { $rssMb = [double]$match.Groups[1].Value / 1024.0 }
    }
    return [pscustomobject]@{
        exit_code = if ($timedOut) { 124 } else { $process.ExitCode }
        timed_out = $timedOut
        wall_seconds = $watch.Elapsed.TotalSeconds
        peak_rss_mb = $rssMb
    }
}

function Read-JpcbApbsEnergy([string]$RunDir) {
    $text = @('stdout.log','stderr.log' | ForEach-Object {
        $path = Join-Path $RunDir $_
        if (Test-Path -LiteralPath $path) { Get-Content -LiteralPath $path -Raw }
    }) -join "`n"
    foreach ($pattern in @(
        'Global net ELEC energy\s*=\s*([-+0-9.eE]+)',
        'Total electrostatic energy\s*=\s*([-+0-9.eE]+)',
        'elecEnergy\s*=\s*([-+0-9.eE]+)')) {
        $match = [regex]::Match($text, $pattern,
            [Text.RegularExpressions.RegexOptions]::IgnoreCase)
        if ($match.Success) { return Convert-JpcbApbsDouble $match.Groups[1].Value }
    }
    return [double]::NaN
}

function Get-JpcbApbsCaseId([string]$Phase, [string]$Component, [double]$SaltM, [int]$Dime) {
    $saltId = $SaltM.ToString('0.00', $culture).Replace('.','p')
    return "$Phase`_$Component`_s$saltId`_d$dime"
}

function Invoke-JpcbApbsCase([string]$Phase, [string]$Component, [string]$PqrPath,
                              [double]$SaltM, [int]$Dime) {
    $caseId = Get-JpcbApbsCaseId $Phase $Component $SaltM $Dime
    $runDir = Join-Path $WorkRoot $caseId
    $resultPath = Join-Path $runDir 'case_result.csv'
    if (-not $Force -and (Test-Path -LiteralPath $resultPath)) {
        $prior = @(Import-Csv -LiteralPath $resultPath)
        if ($prior.Count -eq 1 -and -not ($RerunFailed -and $prior[0].status -ne 'ok')) {
            Write-Host "resume: $caseId"
            return
        }
    }
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    Copy-Item -LiteralPath $PqrPath -Destination (Join-Path $runDir (Split-Path -Leaf $PqrPath))
    $box = Write-JpcbApbsInput -PqrPath (Join-Path $runDir (Split-Path -Leaf $PqrPath)) `
        -InputPath (Join-Path $runDir 'apbs.in') -Dime $Dime -SaltM $SaltM
    Write-Host "run: $caseId (APBS OpenMP=$ApbsThreads timeout ${ApbsTimeoutSeconds}s)"
    $run = Invoke-JpcbApbsWslCase -RunDir $runDir -TimeoutSeconds $ApbsTimeoutSeconds
    $energy = if ($run.exit_code -eq 0) { Read-JpcbApbsEnergy $runDir } else { [double]::NaN }
    $status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) {
        'failed'
    } elseif (-not (Test-JpcbApbsFinite $energy)) { 'missing_energy' } else { 'ok' }
    [pscustomobject]@{
        protocol_id = 'jpcb-apbs-v1-openmp'
        case_id = $caseId; phase = $Phase; component = $Component; salt_M = $SaltM
        dime = $Dime; fine_scale = $FineScale; coarse_padding_angstrom = $CoarsePaddingAngstrom
        atom_count = $box.atoms; apbs_threads = $ApbsThreads; openblas_threads = 1
        threading_policy = "OMP_NUM_THREADS=$ApbsThreads;OPENBLAS_NUM_THREADS=1;OMP_PROC_BIND=true"
        apbs_wsl_exe = $ApbsWslExe; status = $status; exit_code = $run.exit_code
        timed_out = $run.timed_out; reaction_energy_kJmol = $energy
        wall_seconds = $run.wall_seconds; peak_rss_mb = $run.peak_rss_mb
    } | Export-Csv -LiteralPath $resultPath -NoTypeInformation
}

function Write-JpcbApbsPlot([object[]]$Rows, [string]$Path) {
    $valid = @($Rows | Where-Object { $_.status -eq 'ok' } |
        Sort-Object { Convert-JpcbApbsDouble $_.salt_M })
    if ($valid.Count -eq 0) { return }
    $width=850; $height=520; $left=100; $right=50; $top=45; $bottom=85
    $plotW=$width-$left-$right; $plotH=$height-$top-$bottom
    $values=@($valid | ForEach-Object { Convert-JpcbApbsDouble $_.delta_g_pol_bind_kJmol })
    $min=($values | Measure-Object -Minimum).Minimum; $max=($values | Measure-Object -Maximum).Maximum
    if ($max -le $min) { $max=$min+1.0 }
    $xMax=[Math]::Max(0.5,($valid | ForEach-Object { Convert-JpcbApbsDouble $_.salt_M } | Measure-Object -Maximum).Maximum)
    function X([double]$v) { return $left+$plotW*$v/$xMax }
    function Y([double]$v) { return $top+$plotH*(1.0-($v-$min)/($max-$min)) }
    $points=($valid | ForEach-Object { "$(X (Convert-JpcbApbsDouble $_.salt_M)),$(Y (Convert-JpcbApbsDouble $_.delta_g_pol_bind_kJmol))" }) -join ' '
    $svg=@(
        "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<line x1='$left' y1='$top' x2='$left' y2='$($top+$plotH)' stroke='black' stroke-width='2'/>",
        "<line x1='$left' y1='$($top+$plotH)' x2='$($left+$plotW)' y2='$($top+$plotH)' stroke='black' stroke-width='2'/>",
        "<polyline points='$points' fill='none' stroke='#0072B2' stroke-width='3'/>",
        "<text x='$($width/2)' y='$($height-25)' text-anchor='middle' font-family='Arial' font-size='22'>Ionic strength (M)</text>",
        "<text x='28' y='$($height/2)' text-anchor='middle' font-family='Arial' font-size='22' transform='rotate(-90 28 $($height/2))'>APBS polar binding cycle (kJ mol-1)</text>",
        '</svg>'
    )
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$manifestPath = Join-Path $FixtureRoot 'binding_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw 'jpcb_apbs_fixture_gap: run test\pbe_jpcb_prepare_barnase_barstar.ps1 -RefreshFixtures'
}
$fixtureRows = @(Import-Csv -LiteralPath $manifestPath)
if ($fixtureRows.Count -ne 1) { throw 'Expected exactly one binding fixture manifest row.' }
$fixture = $fixtureRows[0]
$components = @(
    [pscustomobject]@{ name='complex'; pqr=(Join-Path $FixtureRoot $fixture.complex_pqr); sha=$fixture.complex_pqr_sha256 },
    [pscustomobject]@{ name='barnase'; pqr=(Join-Path $FixtureRoot $fixture.receptor_pqr); sha=$fixture.receptor_pqr_sha256 },
    [pscustomobject]@{ name='barstar'; pqr=(Join-Path $FixtureRoot $fixture.ligand_pqr); sha=$fixture.ligand_pqr_sha256 }
)
foreach ($component in $components) {
    if (-not (Test-Path -LiteralPath $component.pqr)) { throw "Missing PQR: $($component.pqr)" }
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $component.pqr).Hash.ToLowerInvariant() -ne $component.sha) {
        throw "jpcb_apbs_fixture_gap: checksum drift for $($component.name)"
    }
}
& wsl.exe -d $WslDistribution -- test -x $ApbsWslExe
if ($LASTEXITCODE -ne 0) { throw "APBS OpenMP executable is unavailable in WSL: $ApbsWslExe" }

$saltValues = Convert-JpcbApbsList $Salts 'salt'
$dimeValues = Convert-JpcbApbsList $Dimes 'dime'
if ($RunTitan -and @($saltValues | Where-Object { [double]$_ -le 0.0 }).Count -gt 0) {
    throw 'Paired TITAN/APBS protocol requires positive ionic strength because the guarded FMM/ILU0 TITAN path operates on the coupled JCP system.'
}
if ($Quick) { $saltValues=@(0.15); $dimeValues=@($dimeValues | Select-Object -First 2) }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
if ($ContractOnly) {
    $contractDir=Join-Path $WorkRoot 'contract_only'
    New-FmmCliCleanDir -Path $contractDir -WorkRoot $WorkRoot
    Copy-Item -LiteralPath $components[0].pqr -Destination (Join-Path $contractDir (Split-Path -Leaf $components[0].pqr))
    Write-JpcbApbsInput -PqrPath (Join-Path $contractDir (Split-Path -Leaf $components[0].pqr)) `
        -InputPath (Join-Path $contractDir 'apbs.in') -Dime $dimeValues[0] -SaltM 0.15 | Out-Null
    $text=Get-Content -LiteralPath (Join-Path $contractDir 'apbs.in') -Raw
    if ($text -notmatch 'mg-auto' -or $text -notmatch 'pdie 1.0' -or $text -notmatch 'sdie 80.0') {
        throw 'APBS input contract missing required PB settings.'
    }
    Write-Host "PASS JPCB APBS contract: $contractDir"
    exit 0
}

# Calibrate only the complex at 0.15 M, then use the least expensive stable grid for every component/salt row.
foreach ($dime in $dimeValues) { Invoke-JpcbApbsCase -Phase calibration -Component complex -PqrPath $components[0].pqr -SaltM 0.15 -Dime $dime }
$caseFiles=@(Get-ChildItem -LiteralPath $WorkRoot -Directory | ForEach-Object {
    $path=Join-Path $_.FullName 'case_result.csv'; if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$calibration=@($caseFiles | Where-Object { $_.phase -eq 'calibration' -and $_.component -eq 'complex' } |
    Sort-Object { [int]$_.dime })
$calibration | ForEach-Object {
    $_ | Add-Member -NotePropertyName adjacent_relative_energy_delta -NotePropertyValue '' -Force
}
$selectedDime=0; $selection='unresolved'
for($i=1;$i -lt $calibration.Count;$i++) {
    $previous=$calibration[$i-1]; $current=$calibration[$i]
    if($previous.status -eq 'ok' -and $current.status -eq 'ok') {
        $e0=Convert-JpcbApbsDouble $previous.reaction_energy_kJmol; $e1=Convert-JpcbApbsDouble $current.reaction_energy_kJmol
        $delta=[Math]::Abs($e1-$e0)/[Math]::Max(1.0,[Math]::Abs($e1))
        $current | Add-Member -NotePropertyName adjacent_relative_energy_delta -NotePropertyValue $delta -Force
        if($selectedDime -eq 0 -and $delta -le $GridStabilityRelativeThreshold) { $selectedDime=[int]$current.dime; $selection='adjacent_energy_plateau' }
    }
}
if($selectedDime -eq 0) {
    $available=@($calibration | Where-Object status -eq 'ok' | Sort-Object { [int]$_.dime })
    if($available.Count -gt 0) { $selectedDime=[int]$available[-1].dime; $selection='finest_available_no_plateau' }
}
if($selectedDime -eq 0) { throw 'jpcb_apbs_grid_calibration_failed: no finite complex APBS calibration row.' }

foreach($salt in $saltValues) {
    foreach($component in $components) { Invoke-JpcbApbsCase -Phase production -Component $component.name -PqrPath $component.pqr -SaltM $salt -Dime $selectedDime }
}
$caseFiles=@(Get-ChildItem -LiteralPath $WorkRoot -Directory | ForEach-Object {
    $path=Join-Path $_.FullName 'case_result.csv'; if(Test-Path -LiteralPath $path){ Import-Csv -LiteralPath $path }
})
$calibration | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_grid_calibration.csv') -NoTypeInformation
$production=@($caseFiles | Where-Object { $_.phase -eq 'production' })
$production | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_component_matrix.csv') -NoTypeInformation
$bindingRows=New-Object System.Collections.Generic.List[object]
foreach($salt in $saltValues) {
    $selected=@($production | Where-Object { (Convert-JpcbApbsDouble $_.salt_M) -eq $salt })
    $complex=@($selected | Where-Object component -eq 'complex'); $barnase=@($selected | Where-Object component -eq 'barnase'); $barstar=@($selected | Where-Object component -eq 'barstar')
    $allOk=$complex.Count -eq 1 -and $barnase.Count -eq 1 -and $barstar.Count -eq 1 -and @($complex+$barnase+$barstar | Where-Object status -ne 'ok').Count -eq 0
    $delta=[double]::NaN
    if($allOk){ $delta=(Convert-JpcbApbsDouble $complex[0].reaction_energy_kJmol)-(Convert-JpcbApbsDouble $barnase[0].reaction_energy_kJmol)-(Convert-JpcbApbsDouble $barstar[0].reaction_energy_kJmol) }
    $bindingRows.Add([pscustomobject]@{ protocol_id='jpcb-apbs-v1-openmp'; salt_M=$salt; dime=$selectedDime; grid_selection=$selection; apbs_threads=$ApbsThreads; delta_g_pol_bind_kJmol=$delta; status=if($allOk){'ok'}else{'component_failure'}; complex_energy_kJmol=if($complex.Count){$complex[0].reaction_energy_kJmol}else{''}; barnase_energy_kJmol=if($barnase.Count){$barnase[0].reaction_energy_kJmol}else{''}; barstar_energy_kJmol=if($barstar.Count){$barstar[0].reaction_energy_kJmol}else{''} }) | Out-Null
}
$bindingArray=@($bindingRows.ToArray())
$bindingArray | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_thermocycle.csv') -NoTypeInformation
Write-JpcbApbsPlot -Rows $bindingArray -Path (Join-Path $WorkRoot 'jpcb_apbs_binding_plot.svg')

if($RunTitan) {
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'pbe_jpcb_binding_salt_sweep.ps1') `
        -TitanExe $TitanExe -WorkRoot $TitanWorkRoot -Salts $Salts -Threads $TitanThreads `
        -CaseTimeoutSeconds $TitanTimeoutSeconds -AllowFailedDecision
    if($LASTEXITCODE -ne 0) { throw 'TITAN binding salt wrapper failed.' }
}
$titanThermo=Join-Path $TitanWorkRoot 'jpcb_binding_thermocycle.csv'
if(Test-Path -LiteralPath $titanThermo) {
    $titanRows=@(Import-Csv -LiteralPath $titanThermo)
    $context=foreach($apbsRow in $bindingArray) {
        $titan=@($titanRows | Where-Object { [Math]::Abs((Convert-JpcbApbsDouble $_.salt_M)-(Convert-JpcbApbsDouble $apbsRow.salt_M)) -lt 1e-12 })
        [pscustomobject]@{salt_M=$apbsRow.salt_M; apbs_delta_g_pol_bind_kJmol=$apbsRow.delta_g_pol_bind_kJmol; apbs_status=$apbsRow.status; titan_delta_g_pol_bind_kJmol=if($titan.Count){$titan[0].delta_g_pol_bind_kJmol}else{''}; titan_status=if($titan.Count){$titan[0].status}else{'not_run'}; titan_threads=$TitanThreads; apbs_threads=$ApbsThreads; timing_comparison_policy='cross_os_context_only'}
    }
    $context | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_titan_binding_context.csv') -NoTypeInformation
}
$productionOk=@($production | Where-Object status -eq 'ok').Count
$expected=$saltValues.Count*3
$decision=if($productionOk -ne $expected){'jpcb_apbs_runtime_or_component_gap'}elseif($selection -eq 'adjacent_energy_plateau'){'jpcb_apbs_protocol_pass'}else{'jpcb_apbs_protocol_pass_grid_conservative'}
[pscustomobject]@{decision=$decision; fixture=$fixture.fixture; pdb_id=$fixture.pdb_id; selected_dime=$selectedDime; grid_selection=$selection; apbs_threads=$ApbsThreads; titan_threads=$TitanThreads; production_rows=$expected; production_ok=$productionOk; generated_utc=[DateTime]::UtcNow.ToString('o')} | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_decision.csv') -NoTypeInformation
Set-Content -LiteralPath (Join-Path $WorkRoot 'jpcb_apbs_summary.dat') -Encoding ASCII -Value @("decision=$decision","selected_dime=$selectedDime","grid_selection=$selection","apbs_threads=$ApbsThreads","titan_threads=$TitanThreads","calibration=$(Join-Path $WorkRoot 'jpcb_apbs_grid_calibration.csv')","thermocycle=$(Join-Path $WorkRoot 'jpcb_apbs_thermocycle.csv')")
Set-Content -LiteralPath (Join-Path $RepoRoot 'TITAN_JPCB_APBS_PROTOCOL_REPORT.txt') -Encoding ASCII -Value @('TITAN JPCB APBS Protocol Gate',"Decision: $decision","Fixture: $($fixture.biological_name) ($($fixture.pdb_id))", "APBS executable: $ApbsWslExe", "APBS thread policy: OMP_NUM_THREADS=$ApbsThreads; OPENBLAS_NUM_THREADS=1; OMP_PROC_BIND=true", "TITAN thread policy: -nt=$TitanThreads", 'Wall-time comparisons between native Windows TITAN and WSL APBS are context-only because they use distinct OS environments.', "Canonical APBS dime: $selectedDime ($selection)", "Thermocycle: $(Join-Path $WorkRoot 'jpcb_apbs_thermocycle.csv')")
Write-Host "JPCB APBS decision: $decision"
Write-Host "Canonical APBS dime: $selectedDime ($selection)"
if(-not $AllowFailedDecision -and $decision -notlike 'jpcb_apbs_protocol_pass*'){throw "JPCB APBS protocol gate did not pass: $decision"}
