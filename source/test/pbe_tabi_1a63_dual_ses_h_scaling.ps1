param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_1a63_dual_ses_h_scaling_node_patch',
    [string]$PqrPath = 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr',
    [ValidateSet('grid-ses','analytic-ses')]
    [string[]]$Surfaces = @('grid-ses','analytic-ses'),
    [double[]]$SurfaceSpacings = @(1.00,1.25,1.50,1.75,2.00),
    [int[]]$Threads = @(1,4),
    [ValidateSet('centroid','node-patch')]
    [string]$Discretization = 'node-patch',
    [string]$ReportOutput = 'TITAN_TABI_1A63_DUAL_SES_H_SCALING_NODE_PATCH_REPORT.txt',
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [double[]]$RerunSpacings = @(),
    [int[]]$RerunThreads = @(),
    [ValidateSet('grid-ses','analytic-ses')]
    [string[]]$RerunSurfaces = @(),
    [switch]$PlotOnly,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$v) { $v.ToString('G15', $Inv) }
function D($v, [double]$fallback = [double]::NaN) {
    if ($null -eq $v -or [string]::IsNullOrWhiteSpace([string]$v)) { return $fallback }
    $parsed = 0.0
    if ([double]::TryParse([string]$v, [Globalization.NumberStyles]::Float, $Inv, [ref]$parsed)) {
        return $parsed
    }
    return $fallback
}
function I($v, [int64]$fallback = 0) {
    $parsed = 0L
    if ([int64]::TryParse([string]$v, [ref]$parsed)) { return $parsed }
    return $fallback
}
function Safe([double]$v) { (F $v).Replace('-', 'm').Replace('.', 'p') }
function Xml([string]$v) { [Security.SecurityElement]::Escape($v) }
function Is-Finite([double]$v) { -not [double]::IsNaN($v) -and -not [double]::IsInfinity($v) }

if ($Quick) {
    $SurfaceSpacings = @(1.25)
}
$Surfaces = @($Surfaces | Sort-Object -Unique)
$SurfaceSpacings = @($SurfaceSpacings | Sort-Object -Unique)
$Threads = @($Threads | Sort-Object -Unique)
foreach ($threadCount in $Threads) {
    if ($threadCount -notin @(1,4)) { throw 'This gate supports Threads=1,4 only.' }
}
if (4 -notin $Threads) { throw 'Threads must include 4 because panel A uses the matched four-core energy rows.' }

$resolvedTitan = Resolve-Path -LiteralPath (Join-Path $RepoRoot $TitanExe) -ErrorAction SilentlyContinue
$TitanExe = if ($null -ne $resolvedTitan) { $resolvedTitan.Path } else { (Resolve-Path $TitanExe).Path }
if (-not [IO.Path]::IsPathRooted($PqrPath)) { $PqrPath = Join-Path $RepoRoot $PqrPath }
$PqrPath = (Resolve-Path -LiteralPath $PqrPath).Path
$WorkRoot = if ([IO.Path]::IsPathRooted($WorkRoot)) { $WorkRoot } else { Join-Path $RepoRoot $WorkRoot }
New-Item -ItemType Directory -Force $WorkRoot | Out-Null
$TitanSha = (Get-FileHash -LiteralPath $TitanExe -Algorithm SHA256).Hash.ToLowerInvariant()
$Protocol = "1a63|amber-ph7|p1|theta0.45|asi|ilu64|ilu-regularization=adaptive|tol1e-4|salt0.15|probe1.4|audit512|discretization=$Discretization"

$MatrixPath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling_matrix.csv'
$DecisionPath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling_summary.dat'
$FigurePath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling.svg'
$ReportPath = if ([IO.Path]::IsPathRooted($ReportOutput)) {
    $ReportOutput
} else {
    Join-Path $RepoRoot $ReportOutput
}

$configs = [Collections.Generic.List[object]]::new()
foreach ($surface in $Surfaces) {
    foreach ($spacing in $SurfaceSpacings) {
        foreach ($threadCount in $Threads) {
            $configs.Add([pscustomobject]@{
                surface = $surface
                spacing = [double]$spacing
                threads = [int]$threadCount
                row_id = ('1a63_{0}_h{1}_nt{2}' -f
                    ($surface.Replace('-','_')),(Safe $spacing),$threadCount)
            })
        }
    }
}

$checkpoint = @{}
if (Test-Path $MatrixPath) {
    foreach ($row in @(Import-Csv $MatrixPath)) { $checkpoint[[string]$row.row_id] = $row }
}
function Write-Checkpoint {
    $ordered = foreach ($config in $configs) {
        if ($checkpoint.ContainsKey([string]$config.row_id)) { $checkpoint[[string]$config.row_id] }
    }
    @($ordered) | Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
}

$rerunSpacingSet = @($RerunSpacings)
$rerunThreadSet = @($RerunThreads)
$rerunSurfaceSet = @($RerunSurfaces)
if ($PlotOnly -and -not (Test-Path $MatrixPath)) {
    throw "PlotOnly requires an existing matrix: $MatrixPath"
}
if (-not $PlotOnly) {
foreach ($config in $configs) {
    $old = if ($checkpoint.ContainsKey($config.row_id)) { $checkpoint[$config.row_id] } else { $null }
    $hashMatches = $null -ne $old -and
        $old.PSObject.Properties.Name -contains 'titan_sha256' -and
        [string]$old.titan_sha256 -eq $TitanSha -and
        $old.PSObject.Properties.Name -contains 'protocol' -and
        [string]$old.protocol -eq $Protocol
    $targeted =
        (@($rerunSpacingSet).Count -gt 0 -and $config.spacing -in $rerunSpacingSet) -or
        (@($rerunThreadSet).Count -gt 0 -and $config.threads -in $rerunThreadSet) -or
        (@($rerunSurfaceSet).Count -gt 0 -and $config.surface -in $rerunSurfaceSet)
    $rerun = $Force -or $targeted -or ($RerunFailed -and $null -ne $old -and [string]$old.status -ne 'ok')
    if ($hashMatches -and -not $rerun) {
        Write-Host "resume: $($config.row_id)"
        continue
    }

    $runDir = Join-Path $WorkRoot $config.row_id
    if (Test-Path $runDir) { Remove-Item -LiteralPath $runDir -Recurse -Force }
    New-Item -ItemType Directory -Force $runDir | Out-Null
    $arg = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
        -PqrPath $PqrPath -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -FmmLocal off `
        -SurfacePoints 8 -Threads $config.threads -MaxIter 720 -Restart 128 `
        -Tolerance 1.0e-4 -SaltM 0.15 -InternalDielectric 1.0 `
        -ExternalDielectric 80.0 -TemperatureK 298.15 -ProbeRadius 1.4 `
        -OpeningTheta 0.45 -AuditMode sampled -AuditSampleCount 512 `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -PreconditionerIluRegularization adaptive `
        -PreconditionerIluRegularizationMax 0.01 -NearIntegration analytic_singular `
        -TaylorOrder 4 -PrimarySurface $config.surface `
        -SesGridSpacing $config.spacing -DgSolv $true
    Add-Content -LiteralPath $arg -Encoding ASCII `
        -Value "-pbe-tabi-discretization=$Discretization"

    Write-Host "run: $($config.row_id) (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) { 'failed' } else { 'ok' }
    $profile = $null; $energy = $null; $quality = $null
    if ($status -eq 'ok') {
        try {
            $profile = @(Import-Csv (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
            $energy = @(Import-Csv (Join-Path $runDir 'pbe_energy.csv'))[-1]
            $quality = @(Import-Csv (Join-Path $runDir 'pbe_tabi_surface_quality.csv'))[-1]
            if ([string]$profile.status -ne 'ok') { $status = [string]$profile.status }
        } catch {
            $status = 'profile_schema_gap'
        }
    }
    $energyKj = if ($null -ne $energy) { D $energy.E_rxn_kJmol } else { [double]::NaN }
    $checkpoint[$config.row_id] = [pscustomobject][ordered]@{
        row_id = $config.row_id
        protocol = $Protocol
        titan_sha256 = $TitanSha
        fixture = '1A63'
        pqr_path = $PqrPath
        surface = $config.surface
        discretization = $Discretization
        grid_spacing_A = F $config.spacing
        threads = $config.threads
        status = $status
        exit_code = $run.exit_code
        timed_out = $run.timed_out
        atoms = if ($profile) { I $profile.atom_count } else { 0 }
        panels = if ($profile) { I $profile.panels } else { 0 }
        mesh_vertices = if ($profile) { I $profile.mesh_vertices } else { 0 }
        topology_status = if ($quality) { [string]$quality.calderon_trace_space_status } else { '' }
        boundary_edges = if ($quality) { I $quality.calderon_trace_boundary_edges -1 } else { -1 }
        nonmanifold_edges = if ($quality) { I $quality.calderon_trace_nonmanifold_edges -1 } else { -1 }
        deltaG_solv_kJ_mol = $energyKj
        deltaG_solv_kcal_mol = $energyKj / 4.184
        iterations = if ($profile) { I $profile.iterations } else { 0 }
        true_residual = if ($profile) { D $profile.true_residual_rel } else { [double]::NaN }
        audit_rms = if ($profile) { D $profile.tree_vs_direct_matvec_rel } else { [double]::NaN }
        surface_seconds = if ($profile) { D $profile.surface_build_seconds } else { [double]::NaN }
        ilu_setup_seconds = if ($profile) { D $profile.preconditioner_setup_seconds } else { [double]::NaN }
        gmres_seconds = if ($profile) { D $profile.gmres_seconds } else { [double]::NaN }
        total_seconds = if ($profile) { D $profile.total_seconds } else { D $run.wall_seconds }
        peak_rss_mb = D $run.peak_working_set_mb
        preconditioner_effective = if ($profile) { [string]$profile.preconditioner_effective } else { '' }
        ilu_regularization_mode = if ($profile) { [string]$profile.preconditioner_ilu_regularization } else { '' }
        ilu_regularization_selected = if ($profile) { D $profile.preconditioner_ilu_regularization_selected 0.0 } else { 0.0 }
        fmm_mode = if ($profile) { [string]$profile.fmm_m2l_mode } else { '' }
    }
    Write-Checkpoint
}
}

$rows = @($configs | ForEach-Object { if ($checkpoint.ContainsKey($_.row_id)) { $checkpoint[$_.row_id] } })
$requiredCount = $configs.Count
$valid = @($rows | Where-Object {
    [string]$_.status -eq 'ok' -and
    (Is-Finite (D $_.deltaG_solv_kcal_mol)) -and
    (D $_.true_residual 1) -le 1.0e-4 -and
    (D $_.audit_rms 1) -le 8.0e-3 -and
    [string]$_.preconditioner_effective -like 'fmm_nearfield_ilu0*' -and
    [string]$_.discretization -eq $Discretization
})
$decision = if ($rows.Count -lt $requiredCount) {
    'dual_ses_h_scaling_incomplete'
} elseif ($valid.Count -eq $requiredCount) {
    'dual_ses_h_scaling_ready'
} else {
    'dual_ses_h_scaling_accuracy_gap'
}

function Write-Figure([object[]]$data) {
    $ok = @($data | Where-Object { [string]$_.status -eq 'ok' })
    if ($ok.Count -eq 0) { return }
    $energyRows = @($ok | Where-Object { [int]$_.threads -eq 4 -and (Is-Finite (D $_.deltaG_solv_kcal_mol)) })
    $timeRows = @($ok | Where-Object { (Is-Finite (D $_.total_seconds)) -and (D $_.total_seconds) -gt 0 })
    if ($energyRows.Count -eq 0 -or $timeRows.Count -eq 0) { return }

    $width=1800; $height=680; $top=130; $bottom=555
    $aLeft=150; $aRight=800; $bLeft=1040; $bRight=1690
    $hMin=($SurfaceSpacings | Measure-Object -Minimum).Minimum
    $hMax=($SurfaceSpacings | Measure-Object -Maximum).Maximum
    $energies=@($energyRows | ForEach-Object { D $_.deltaG_solv_kcal_mol })
    $eMin=($energies|Measure-Object -Minimum).Minimum; $eMax=($energies|Measure-Object -Maximum).Maximum
    $ePad=[math]::Max(25.0,0.08*($eMax-$eMin)); $eMin-=$ePad; $eMax+=$ePad
    $times=@($timeRows|ForEach-Object { D $_.total_seconds }); $tMin=[math]::Max(1.0,($times|Measure-Object -Minimum).Minimum*0.75); $tMax=($times|Measure-Object -Maximum).Maximum*1.35
    function X([double]$h,[double]$left,[double]$right){$left+($h-$hMin)/[math]::Max(1e-12,$hMax-$hMin)*($right-$left)}
    function YE([double]$e){$bottom-($e-$eMin)/[math]::Max(1e-12,$eMax-$eMin)*($bottom-$top)}
    function YT([double]$t){$bottom-([math]::Log10($t)-[math]::Log10($tMin))/[math]::Max(1e-12,[math]::Log10($tMax)-[math]::Log10($tMin))*($bottom-$top)}
    $colors=@{'grid-ses'='#0072B2';'analytic-ses'='#D55E00'}
    $labels=@{'grid-ses'='Watertight Grid-SES';'analytic-ses'='Analytic patch-ray SES'}
    $svg=[Collections.Generic.List[string]]::new()
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="68mm" viewBox="0 0 1800 680"><rect width="100%" height="100%" fill="white"/>')
    $svg.Add('<style>text{font-family:Arial;fill:#222}.axis{stroke:#222;stroke-width:2}.grid{stroke:#dedede;stroke-width:1}.guide{stroke:#777;stroke-width:2;stroke-dasharray:7 7}</style>')
    $svg.Add('<text x="150" y="78" font-size="34" font-weight="bold">A</text><text x="1040" y="78" font-size="34" font-weight="bold">B</text>')
    $legendX=420
    foreach($surface in @('grid-ses','analytic-ses')){$c=$colors[$surface];$svg.Add("<line x1=`"$legendX`" y1=`"42`" x2=`"$($legendX+55)`" y2=`"42`" stroke=`"$c`" stroke-width=`"5`"/><circle cx=`"$($legendX+28)`" cy=`"42`" r=`"7`" fill=`"$c`"/><text x=`"$($legendX+68)`" y=`"51`" font-size=`"27`">$(Xml $labels[$surface])</text>");$legendX+=390}
    $svg.Add('<line x1="1220" y1="42" x2="1275" y2="42" stroke="#333" stroke-width="4"/><text x="1288" y="51" font-size="26">4 cores</text><line x1="1430" y1="42" x2="1485" y2="42" stroke="#333" stroke-width="4" stroke-dasharray="10 7"/><rect x="1450" y="34" width="16" height="16" fill="white" stroke="#333" stroke-width="3"/><text x="1498" y="51" font-size="26">1 core</text>')
    foreach($left in @($aLeft,$bLeft)){$right=if($left-eq$aLeft){$aRight}else{$bRight};$svg.Add("<line class=`"axis`" x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$bottom`"/><line class=`"axis`" x1=`"$left`" y1=`"$bottom`" x2=`"$right`" y2=`"$bottom`"/>");foreach($h in $SurfaceSpacings){$x=X $h $left $right;$svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/><text x=`"$x`" y=`"590`" text-anchor=`"middle`" font-size=`"25`">$(F $h)</text>")};$hx=X 1.25 $left $right;$svg.Add("<line class=`"guide`" x1=`"$hx`" y1=`"$top`" x2=`"$hx`" y2=`"$bottom`"/>")}
    for($i=0;$i-le4;$i++){ $e=$eMin+$i*($eMax-$eMin)/4;$y=YE $e;$svg.Add("<line class=`"grid`" x1=`"$aLeft`" y1=`"$y`" x2=`"$aRight`" y2=`"$y`"/><text x=`"$($aLeft-15)`" y=`"$($y+8)`" text-anchor=`"end`" font-size=`"24`">$($e.ToString('F0',$Inv))</text>") }
    $timeTicks=@(1,2,5,10,20,50,100,200,500,1000,2000,5000)|Where-Object{$_-ge$tMin-and$_-le$tMax}
    foreach($tick in $timeTicks){$y=YT $tick;$svg.Add("<line class=`"grid`" x1=`"$bLeft`" y1=`"$y`" x2=`"$bRight`" y2=`"$y`"/><text x=`"$($bLeft-15)`" y=`"$($y+8)`" text-anchor=`"end`" font-size=`"24`">$tick</text>")}
    foreach($surface in @('grid-ses','analytic-ses')){
        $c=$colors[$surface]
        $series=@($energyRows|Where-Object{$_.surface-eq$surface}|Sort-Object{D $_.grid_spacing_A});if($series.Count-gt0){$points=@($series|ForEach-Object{"$(X (D $_.grid_spacing_A) $aLeft $aRight),$(YE (D $_.deltaG_solv_kcal_mol))"})-join' ';$svg.Add("<polyline points=`"$points`" fill=`"none`" stroke=`"$c`" stroke-width=`"5`"/>");foreach($r in $series){$x=X (D $r.grid_spacing_A) $aLeft $aRight;$y=YE (D $r.deltaG_solv_kcal_mol);$svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"8`" fill=`"$c`" stroke=`"white`" stroke-width=`"2`"/>")}}
        foreach($nt in @(1,4)){$series=@($timeRows|Where-Object{$_.surface -eq $surface -and [int]$_.threads -eq $nt}|Sort-Object{D $_.grid_spacing_A});if($series.Count-eq0){continue};$dash=if($nt-eq1){' stroke-dasharray="11 8"'}else{''};$points=@($series|ForEach-Object{"$(X (D $_.grid_spacing_A) $bLeft $bRight),$(YT (D $_.total_seconds))"})-join' ';$svg.Add("<polyline points=`"$points`" fill=`"none`" stroke=`"$c`" stroke-width=`"5`"$dash/>");foreach($r in $series){$x=X (D $r.grid_spacing_A) $bLeft $bRight;$y=YT (D $r.total_seconds);if($nt-eq1){$svg.Add("<rect x=`"$($x-7)`" y=`"$($y-7)`" width=`"14`" height=`"14`" fill=`"white`" stroke=`"$c`" stroke-width=`"3`"/>")}else{$svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"8`" fill=`"$c`" stroke=`"white`" stroke-width=`"2`"/>")}}}
    }
    $svg.Add('<text x="475" y="640" text-anchor="middle" font-size="29">Surface spacing, h (&#197;)</text><text x="1365" y="640" text-anchor="middle" font-size="29">Surface spacing, h (&#197;)</text>')
    $svg.Add('<text transform="translate(38,360) rotate(-90)" text-anchor="middle" font-size="29">&#916;Gsolv (kcal mol&#8722;1)</text><text transform="translate(925,360) rotate(-90)" text-anchor="middle" font-size="29">Total wall time (s, log scale)</text>')
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $FigurePath -Encoding UTF8
    $inkscape = Get-Command inkscape.exe -ErrorAction SilentlyContinue
    if ($null -ne $inkscape) {
        $pngPath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling.png'
        $pdfPath = Join-Path $WorkRoot 'tabi_1a63_dual_ses_h_scaling.pdf'
        & $inkscape.Source $FigurePath '--export-type=png' '--export-dpi=600' "--export-filename=$pngPath" | Out-Null
        & $inkscape.Source $FigurePath '--export-type=pdf' "--export-filename=$pdfPath" | Out-Null
    }
}

Write-Figure $rows
[pscustomobject]@{decision=$decision;required_rows=$requiredCount;valid_rows=$valid.Count;matrix=$MatrixPath;figure=$FigurePath} |
    Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@("decision=$decision","required_rows=$requiredCount","valid_rows=$($valid.Count)","matrix=$MatrixPath","figure=$FigurePath") |
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII
@(
    'TITAN TABI 1A63 Dual-SES Surface-Spacing Report',
    "Decision: $decision",
    "Protocol: full 1A63 AMBER PQR; $Discretization discretization; 0.15 M; p1 analytical FMM/M2L; theta 0.45; analytic-singular near field; FMM-near ILU(0) cap 64 with adaptive regularization; tolerance 1e-4.",
    'Panel A uses matched four-core rows. Panel B compares one- and four-core wall time.',
    'The spacing sweep is a sensitivity study; no convergence order is claimed.',
    "Matrix: $MatrixPath",
    "Figure: $FigurePath"
) | Set-Content -LiteralPath $ReportPath -Encoding UTF8
Write-Host "Decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Figure: $FigurePath"
Write-Host 'Resume: rerun unchanged; use -RerunFailed, -RerunSpacings, -RerunThreads, or -Force as needed.'
if (-not $AllowFailedDecision -and $decision -ne 'dual_ses_h_scaling_ready') { exit 1 }
