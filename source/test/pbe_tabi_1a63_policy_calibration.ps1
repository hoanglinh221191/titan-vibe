param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_1a63_policy_calibration',
    [string]$PqrPath = 'test\pbe_inputs_public\1a63_full_amber_ph7.pqr',
    [string]$ReportPath = 'TITAN_TABI_1A63_POLICY_CALIBRATION_REPORT.txt',
    [string[]]$Stages = @('surface','probe','theta','ilu'),
    [double[]]$GridSpacings = @(1.50,1.25,1.00,0.80),
    [double[]]$ProbeRadii = @(1.20,1.40,1.60),
    [double[]]$OpeningThetas = @(0.35,0.45,0.55,0.65),
    [int[]]$IluNeighborCaps = @(32,64,128),
    [double]$PolicyGridSpacing = 1.00,
    [double]$PolicyProbeRadius = 1.40,
    [double]$PolicyOpeningTheta = 0.45,
    [int]$PolicyIluNeighborCap = 64,
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$IncludeDirect,
    [switch]$IncludePreconditionerComparators,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function F([double]$v) { $v.ToString('G15', $Inv) }
function Tick([double]$v,[int]$digits=2) { $v.ToString("F$digits", $Inv) }
function Safe([double]$v) { (F $v).Replace('-', 'm').Replace('.', 'p') }
function D($v, [double]$fallback = [double]::NaN) {
    if ($null -eq $v -or [string]::IsNullOrWhiteSpace([string]$v)) { return $fallback }
    $parsed = 0.0
    if ([double]::TryParse([string]$v, [Globalization.NumberStyles]::Float,
            $Inv, [ref]$parsed)) { return $parsed }
    return $fallback
}
function I($v, [int64]$fallback = 0) {
    $parsed = 0L
    if ([int64]::TryParse([string]$v, [ref]$parsed)) { return $parsed }
    return $fallback
}
function B($v) { ([string]$v).ToLowerInvariant() -eq 'true' }
function Split-Values([string[]]$values) {
    @($values | ForEach-Object { ([string]$_).Split(',') } |
        ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
}
function Xml([string]$v) { [Security.SecurityElement]::Escape($v) }

$stageSet = @(Split-Values $Stages)
foreach ($stage in $stageSet) {
    if ($stage -notin @('surface','probe','theta','ilu','comparators')) {
        throw "Unknown stage '$stage'."
    }
}
if ($Quick) {
    $stageSet = @('surface')
    $GridSpacings = @($GridSpacings | Sort-Object -Descending | Select-Object -First 2)
}
$repoTitan = Resolve-Path -LiteralPath (Join-Path $RepoRoot $TitanExe) -ErrorAction SilentlyContinue
$TitanExe = if ($null -ne $repoTitan) { $repoTitan.Path } else { (Resolve-Path $TitanExe).Path }
if (-not [IO.Path]::IsPathRooted($PqrPath)) { $PqrPath = Join-Path $RepoRoot $PqrPath }
$PqrPath = (Resolve-Path -LiteralPath $PqrPath).Path
$WorkRoot = if ([IO.Path]::IsPathRooted($WorkRoot)) { $WorkRoot } else { Join-Path $RepoRoot $WorkRoot }
New-Item -ItemType Directory -Force $WorkRoot | Out-Null
$TitanSha = (Get-FileHash $TitanExe -Algorithm SHA256).Hash.ToLowerInvariant()

$MatrixPath = Join-Path $WorkRoot 'tabi_1a63_policy_matrix.csv'
$SurfaceTablePath = Join-Path $WorkRoot 'tabi_1a63_surface_table.csv'
$PerformanceTablePath = Join-Path $WorkRoot 'tabi_1a63_performance_table.csv'
$ConvergencePath = Join-Path $WorkRoot 'tabi_1a63_surface_convergence.csv'
$DecisionPath = Join-Path $WorkRoot 'tabi_1a63_policy_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'tabi_1a63_policy_summary.dat'
$EnergyFigurePath = Join-Path $WorkRoot 'tabi_1a63_energy_convergence.svg'
$RuntimeFigurePath = Join-Path $WorkRoot 'tabi_1a63_runtime_accuracy.svg'
$BreakdownFigurePath = Join-Path $WorkRoot 'tabi_1a63_runtime_breakdown.svg'
$ReportPath = if ([IO.Path]::IsPathRooted($ReportPath)) {
    $ReportPath
} else {
    Join-Path $RepoRoot $ReportPath
}

$configs = @{}
function Add-Config([string]$purpose, [string]$backend, [double]$spacing,
                    [double]$probe, [double]$theta, [int]$iluCap,
                    [string]$preconditioner = 'fmm-nearfield-ilu0') {
    $key = '{0}|{1}|{2}|{3}|{4}|{5}' -f $backend,(F $spacing),(F $probe),
        (F $theta),$iluCap,$preconditioner
    if ($configs.ContainsKey($key)) {
        $configs[$key].purpose = @($configs[$key].purpose.Split(';') + $purpose |
            Sort-Object -Unique) -join ';'
        return
    }
    $configs[$key] = [pscustomobject]@{
        purpose=$purpose; backend=$backend; spacing=$spacing; probe=$probe
        theta=$theta; ilu_cap=$iluCap; preconditioner=$preconditioner
    }
}
if ('surface' -in $stageSet) {
    foreach ($h in $GridSpacings) { Add-Config surface fmm $h $PolicyProbeRadius $PolicyOpeningTheta $PolicyIluNeighborCap }
}
if ('probe' -in $stageSet) {
    foreach ($probe in $ProbeRadii) { Add-Config probe fmm $PolicyGridSpacing $probe $PolicyOpeningTheta $PolicyIluNeighborCap }
}
if ('theta' -in $stageSet) {
    foreach ($theta in $OpeningThetas) { Add-Config theta fmm $PolicyGridSpacing $PolicyProbeRadius $theta $PolicyIluNeighborCap }
}
if ('ilu' -in $stageSet) {
    foreach ($cap in $IluNeighborCaps) { Add-Config ilu fmm $PolicyGridSpacing $PolicyProbeRadius $PolicyOpeningTheta $cap }
}
if ($IncludeDirect) {
    $coarsest = $GridSpacings | Sort-Object -Descending | Select-Object -First 1
    Add-Config comparator direct $coarsest $PolicyProbeRadius $PolicyOpeningTheta 0 diagonal
}
if ($IncludePreconditionerComparators -or 'comparators' -in $stageSet) {
    Add-Config comparator fmm $PolicyGridSpacing $PolicyProbeRadius $PolicyOpeningTheta 0 diagonal
    Add-Config comparator fmm $PolicyGridSpacing $PolicyProbeRadius $PolicyOpeningTheta 32 fmm-block-jacobi
}

$existing = @{}
if (Test-Path $MatrixPath) {
    foreach ($row in @(Import-Csv $MatrixPath)) { $existing[[string]$row.row_id] = $row }
}
$rows = [Collections.Generic.List[object]]::new()
foreach ($config in @($configs.Values | Sort-Object purpose,spacing,probe,theta,ilu_cap)) {
    $rowId = '1a63_{0}_h{1}_probe{2}_theta{3}_cap{4}_{5}' -f
        $config.backend,(Safe $config.spacing),(Safe $config.probe),
        (Safe $config.theta),$config.ilu_cap,$config.preconditioner.Replace('-','_')
    $old = $existing[$rowId]
    $hashMatches = $null -ne $old -and
        ($old.PSObject.Properties.Name -contains 'titan_sha256') -and
        [string]$old.titan_sha256 -eq $TitanSha
    $reuse = -not $Force -and $hashMatches -and
        (-not $RerunFailed -or [string]$old.status -eq 'ok')
    if ($reuse) {
        Write-Host "resume: $rowId"
        $rows.Add($old)
        continue
    }
    $runDir = Join-Path $WorkRoot $rowId
    if (Test-Path $runDir) { Remove-Item $runDir -Recurse -Force }
    New-Item -ItemType Directory -Force $runDir | Out-Null
    $isDirect = $config.backend -eq 'direct'
    $arg = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
        -PqrPath $PqrPath -Backend tabi -Mode $(if($isDirect){'direct'}else{'treecode'}) `
        -FmmMode $(if($isDirect){'off'}else{'single-matvec-analytic'}) `
        -SurfacePoints 8 -Threads $Threads -MaxIter 720 -Restart 128 `
        -Tolerance 1.0e-4 -SaltM 0.15 -InternalDielectric 1.0 `
        -ExternalDielectric 80.0 -TemperatureK 298.15 `
        -ProbeRadius $config.probe -OpeningTheta $config.theta `
        -AuditMode $(if($isDirect){'off'}else{'sampled'}) `
        -AuditSampleCount $(if($isDirect){0}else{512}) `
        -Preconditioner $config.preconditioner `
        -PreconditionerBlockPanelCap 32 `
        -PreconditionerIluNeighborCap $(if($config.ilu_cap -gt 0){$config.ilu_cap}else{64}) `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -PrimarySurface grid-ses `
        -SesGridSpacing $config.spacing -DgSolv $true
    Write-Host "run: $rowId (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $status = if ($run.timed_out) {'timeout'} elseif ($run.exit_code -ne 0) {'failed'} else {'ok'}
    $profile=$null; $quality=$null; $energy=$null; $memory=$null
    if ($status -eq 'ok') {
        try {
            $profile = @(Import-Csv (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
            $quality = @(Import-Csv (Join-Path $runDir 'pbe_tabi_surface_quality.csv'))[-1]
            $energy = @(Import-Csv (Join-Path $runDir 'pbe_energy.csv'))[-1]
            $memory = @(Import-Csv (Join-Path $runDir 'pbe_memory_profile.csv'))[-1]
            if ([string]$profile.status -ne 'ok') { $status = [string]$profile.status }
        } catch { $status = 'profile_schema_gap' }
    }
    $energyKj = if ($energy) { D $energy.E_rxn_kJmol } else {[double]::NaN}
    $total = if($profile){D $profile.total_seconds}else{D $run.wall_seconds}
    $known = if($profile){
        (D $profile.surface_build_seconds 0)+(D $profile.rhs_seconds 0)+
        (D $profile.preconditioner_setup_seconds 0)+(D $profile.gmres_seconds 0)+
        (D $profile.matvec_audit_seconds 0)
    }else{0}
    $rows.Add([pscustomobject][ordered]@{
        row_id=$rowId; purpose=$config.purpose; backend=$config.backend
        titan_sha256=$TitanSha; fixture='1A63'; pqr_path=$PqrPath; threads=$Threads
        grid_spacing_A=F $config.spacing; probe_radius_A=F $config.probe
        opening_theta=F $config.theta; ilu_neighbor_cap=$config.ilu_cap
        preconditioner_requested=$config.preconditioner
        status=$status; exit_code=$run.exit_code; timed_out=$run.timed_out
        atoms=if($profile){I $profile.atom_count}else{0}
        panels=if($profile){I $profile.panels}else{0}
        mesh_vertices=if($profile){I $profile.mesh_vertices}else{0}
        topology_status=if($quality){[string]$quality.calderon_trace_space_status}else{''}
        boundary_edges=if($quality){I $quality.calderon_trace_boundary_edges}else{-1}
        nonmanifold_edges=if($quality){I $quality.calderon_trace_nonmanifold_edges}else{-1}
        orientation_mismatches=if($quality){I $quality.calderon_trace_orientation_mismatches}else{-1}
        triangle_area_cv=if($quality){D $quality.triangle_area_cv}else{[double]::NaN}
        deltaG_solv_kJ_mol=$energyKj; deltaG_solv_kcal_mol=$energyKj/4.184
        iterations=if($profile){I $profile.iterations}else{0}
        convergence_reason=if($profile){[string]$profile.convergence_reason}else{''}
        true_residual=if($profile){D $profile.true_residual_rel}else{[double]::NaN}
        sampled_direct_residual=if($profile){D $profile.sampled_true_residual_rel}else{[double]::NaN}
        audit_rms=if($profile){D $profile.tree_vs_direct_matvec_rel}else{[double]::NaN}
        audit_rows=if($profile){I $profile.matvec_audit_rows}else{0}
        surface_seconds=if($profile){D $profile.surface_build_seconds}else{[double]::NaN}
        solid_angle_seconds=if($profile){D $profile.solid_angle_seconds}else{[double]::NaN}
        rhs_seconds=if($profile){D $profile.rhs_seconds}else{[double]::NaN}
        fmm_setup_seconds=if($profile){D $profile.fmm_m2l_setup_seconds}else{[double]::NaN}
        ilu_setup_seconds=if($profile){D $profile.preconditioner_setup_seconds}else{[double]::NaN}
        ilu_apply_seconds=if($profile){D $profile.preconditioner_apply_seconds}else{[double]::NaN}
        gmres_seconds=if($profile){D $profile.gmres_seconds}else{[double]::NaN}
        matvec_seconds_total=if($profile){D $profile.matvec_seconds_total}else{[double]::NaN}
        matvec_seconds_avg=if($profile){D $profile.matvec_seconds_avg}else{[double]::NaN}
        matvec_calls=if($profile){I $profile.matvec_calls}else{0}
        orthogonalization_seconds=if($profile){D $profile.orthogonalization_seconds}else{[double]::NaN}
        true_residual_probe_seconds=if($profile){D $profile.true_residual_probe_seconds}else{[double]::NaN}
        audit_seconds=if($profile){D $profile.matvec_audit_seconds}else{[double]::NaN}
        total_seconds=$total; other_seconds=[math]::Max(0,$total-$known)
        peak_rss_mb=D $run.peak_working_set_mb
        estimated_memory_mb=if($memory){(D $memory.estimated_bytes_total)/1MB}else{[double]::NaN}
        fmm_near_pairs=if($profile){I $profile.fmm_m2l_near_pairs}else{0}
        fmm_far_pairs=if($profile){I $profile.fmm_m2l_far_pairs}else{0}
        ilu_nnz_blocks=if($profile){I $profile.preconditioner_sparse_nnz_blocks}else{0}
        ilu_avg_degree=if($profile){D $profile.preconditioner_sparse_avg_degree}else{[double]::NaN}
        preconditioner_effective=if($profile){[string]$profile.preconditioner_effective}else{''}
        near_kernel=if($profile){[string]$profile.near_kernel_policy}else{''}
    })
    $checkpoint = [Collections.Generic.List[object]]::new()
    $checkpointIds = @{}
    foreach ($item in $rows) {
        $checkpoint.Add($item)
        $checkpointIds[[string]$item.row_id] = $true
    }
    foreach ($prior in $existing.Values) {
        if (-not $checkpointIds.ContainsKey([string]$prior.row_id) -and
            ($prior.PSObject.Properties.Name -contains 'titan_sha256') -and
            [string]$prior.titan_sha256 -eq $TitanSha) {
            $checkpoint.Add($prior)
        }
    }
    $checkpoint | Sort-Object purpose,row_id |
        Export-Csv $MatrixPath -NoTypeInformation
}

# Preserve completed rows from stages not requested in this invocation. This
# lets a long calibration be extended with direct/comparator rows without
# discarding the surface, theta, probe, or ILU evidence already collected.
$currentIds = @{}
foreach ($row in $rows) { $currentIds[[string]$row.row_id] = $true }
foreach ($old in $existing.Values) {
    if (-not $currentIds.ContainsKey([string]$old.row_id) -and
        ($old.PSObject.Properties.Name -contains 'titan_sha256') -and
        [string]$old.titan_sha256 -eq $TitanSha) {
        $rows.Add($old)
    }
}
$rows | Sort-Object purpose,row_id | Export-Csv $MatrixPath -NoTypeInformation

$surfaceRows = @($rows | Where-Object {
    $_.purpose -match 'surface' -and $_.backend -eq 'fmm' -and $_.status -eq 'ok'
} | Sort-Object {[double]$_.grid_spacing_A} -Descending)
$energyInfinity = [double]::NaN
$energyFitR2 = [double]::NaN
$finestExtrapolationGapPercent = [double]::NaN
if ($surfaceRows.Count -ge 2) {
    $sx=0.0;$sy=0.0;$sxx=0.0;$sxy=0.0
    foreach($r in $surfaceRows){$x=1.0/[double]$r.panels;$y=[double]$r.deltaG_solv_kcal_mol;$sx+=$x;$sy+=$y;$sxx+=$x*$x;$sxy+=$x*$y}
    $den=$surfaceRows.Count*$sxx-$sx*$sx
    if([math]::Abs($den)-gt 1e-300){
        $slope=($surfaceRows.Count*$sxy-$sx*$sy)/$den
        $energyInfinity=($sy-$slope*$sx)/$surfaceRows.Count
        $meanEnergy=$sy/$surfaceRows.Count;$ssTotal=0.0;$ssResidual=0.0
        foreach($row in $surfaceRows){$x=1.0/[double]$row.panels;$y=[double]$row.deltaG_solv_kcal_mol;$fit=$energyInfinity+$slope*$x;$ssTotal+=($y-$meanEnergy)*($y-$meanEnergy);$ssResidual+=($y-$fit)*($y-$fit)}
        if($ssTotal-gt1e-300){$energyFitR2=1.0-$ssResidual/$ssTotal}
        $finest=$surfaceRows|Sort-Object {[double]$_.grid_spacing_A}|Select-Object -First 1
        $finestExtrapolationGapPercent=100*[math]::Abs([double]$finest.deltaG_solv_kcal_mol-$energyInfinity)/[math]::Max(1e-300,[math]::Abs($energyInfinity))
    }
}
$surfaceTable=[Collections.Generic.List[object]]::new()
$performance=[Collections.Generic.List[object]]::new()
$convergence=[Collections.Generic.List[object]]::new()
for($i=0;$i-lt$surfaceRows.Count;++$i){
    $r=$surfaceRows[$i]
    $finiteReference = -not [double]::IsNaN($energyInfinity) -and
        -not [double]::IsInfinity($energyInfinity)
    $errorPct=if($finiteReference){100*[math]::Abs([double]$r.deltaG_solv_kcal_mol-$energyInfinity)/[math]::Max(1e-300,[math]::Abs($energyInfinity))}else{[double]::NaN}
    $surfaceTable.Add([pscustomobject]@{grid_spacing_A=$r.grid_spacing_A;probe_radius_A=$r.probe_radius_A;panels=$r.panels;vertices=$r.mesh_vertices;area_cv=$r.triangle_area_cv;topology=$r.topology_status})
    $performance.Add([pscustomobject]@{grid_spacing_A=$r.grid_spacing_A;panels=$r.panels;deltaG_solv_kcal_mol=$r.deltaG_solv_kcal_mol;extrapolated_deltaG_kcal_mol=$energyInfinity;relative_error_percent=$errorPct;iterations=$r.iterations;true_residual=$r.true_residual;audit_rms=$r.audit_rms;surface_seconds=$r.surface_seconds;ilu_setup_seconds=$r.ilu_setup_seconds;gmres_seconds=$r.gmres_seconds;total_seconds=$r.total_seconds;peak_rss_mb=$r.peak_rss_mb})
    if($i+1-lt$surfaceRows.Count){$fine=$surfaceRows[$i+1];$rel=[math]::Abs([double]$r.deltaG_solv_kcal_mol-[double]$fine.deltaG_solv_kcal_mol)/[math]::Max(1e-300,[math]::Abs([double]$fine.deltaG_solv_kcal_mol));$convergence.Add([pscustomobject]@{coarse_spacing_A=$r.grid_spacing_A;fine_spacing_A=$fine.grid_spacing_A;coarse_panels=$r.panels;fine_panels=$fine.panels;energy_delta_kcal_mol=[math]::Abs([double]$r.deltaG_solv_kcal_mol-[double]$fine.deltaG_solv_kcal_mol);relative_energy_delta=$rel;practical_stability_pass=($rel-le0.03-and[double]$r.peak_rss_mb-le512)})}
}
$surfaceTable|Export-Csv $SurfaceTablePath -NoTypeInformation
$performance|Export-Csv $PerformanceTablePath -NoTypeInformation
$convergence|Export-Csv $ConvergencePath -NoTypeInformation

function Map([double]$v,[double]$lo,[double]$hi,[double]$a,[double]$b){if($hi-le$lo){return($a+$b)/2};$a+($v-$lo)*($b-$a)/($hi-$lo)}
function Svg-Plot([string]$path,[object[]]$data,[string]$kind){
    if($data.Count-lt2){return}
    $w=1080;$h=460;$leftA=92;$rightA=492;$leftB=640;$rightB=1040;$top=38;$bottom=368
    $svg=[Collections.Generic.List[string]]::new()
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$w`" height=`"$h`" viewBox=`"0 0 $w $h`">")
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    if($kind-eq'energy'){
        $x=[double[]]@($data|%{1e5/[double]$_.panels})
        $yA=[double[]]@($data|%{[double]$_.deltaG_solv_kcal_mol})
        $yB=[double[]]@($data|%{[math]::Max(1e-6,[double]$_.relative_error_percent)})
        $xB=$x;$xlab='10^5 / N panels';$xBlab=$xlab
        $ylA='Delta G_solv (kcal/mol)';$ylB='Relative error (%)'
    }else{
        $x=[double[]]@($data|%{[double]$_.panels/1000})
        $yA=[double[]]@($data|%{[double]$_.total_seconds})
        $yB=$yA
        $xB=[double[]]@($data|%{[math]::Max(1e-6,[double]$_.relative_error_percent)})
        $xlab='Panels (10^3)';$xBlab='Relative error (%)'
        $ylA='Total wall time (s)';$ylB='Total wall time (s)'
    }
    $xmin=($x|Measure-Object -Minimum).Minimum;$xmax=($x|Measure-Object -Maximum).Maximum
    $xbmin=($xB|Measure-Object -Minimum).Minimum;$xbmax=($xB|Measure-Object -Maximum).Maximum
    $aMin=($yA|Measure-Object -Minimum).Minimum;$aMax=($yA|Measure-Object -Maximum).Maximum
    $bMin=($yB|Measure-Object -Minimum).Minimum;$bMax=($yB|Measure-Object -Maximum).Maximum
    if($kind-eq'energy'){$yBMap=[double[]]@($yB|%{[math]::Log10($_)});$bMin=($yBMap|Measure-Object -Minimum).Minimum;$bMax=($yBMap|Measure-Object -Maximum).Maximum}else{$yBMap=$yB}
    foreach($panel in @(@($leftA,$rightA,'A'),@($leftB,$rightB,'B'))){$l=$panel[0];$r=$panel[1];$svg.Add("<line x1=`"$l`" y1=`"$bottom`" x2=`"$r`" y2=`"$bottom`" stroke=`"#111`" stroke-width=`"1.4`"/>");$svg.Add("<line x1=`"$l`" y1=`"$top`" x2=`"$l`" y2=`"$bottom`" stroke=`"#111`" stroke-width=`"1.4`"/>");$svg.Add("<text x=`"$($l+12)`" y=`"$($top+24)`" font-family=`"Arial`" font-size=`"18`" font-weight=`"bold`">$($panel[2])</text>")}
    for($tick=0;$tick-le4;++$tick){
        $f=$tick/4.0;$y=$bottom-$f*($bottom-$top)
        $va=$aMin+$f*($aMax-$aMin);$vb=$bMin+$f*($bMax-$bMin)
        $labelA=if($kind-eq'energy'){Tick $va 0}else{Tick $va 1}
        $labelB=if($kind-eq'energy'){Tick ([math]::Pow(10,$vb)) 2}else{Tick $vb 1}
        $svg.Add("<line x1=`"$leftA`" y1=`"$y`" x2=`"$rightA`" y2=`"$y`" stroke=`"#dddddd`"/><text x=`"$($leftA-8)`" y=`"$($y+5)`" text-anchor=`"end`" font-family=`"Arial`" font-size=`"12`">$labelA</text>")
        $svg.Add("<line x1=`"$leftB`" y1=`"$y`" x2=`"$rightB`" y2=`"$y`" stroke=`"#dddddd`"/><text x=`"$($leftB-8)`" y=`"$($y+5)`" text-anchor=`"end`" font-family=`"Arial`" font-size=`"12`">$labelB</text>")
    }
    $pointsA=[Collections.Generic.List[string]]::new();$pointsB=[Collections.Generic.List[string]]::new()
    for($i=0;$i-lt$data.Count;++$i){
        $xPosA=Map $x[$i] $xmin $xmax $leftA $rightA;$yPosA=Map $yA[$i] $aMin $aMax $bottom $top
        $xPosB=Map $xB[$i] $xbmin $xbmax $leftB $rightB;$yPosB=Map $yBMap[$i] $bMin $bMax $bottom $top
        $pointsA.Add("$xPosA,$yPosA");$pointsB.Add("$xPosB,$yPosB")
        $svg.Add("<line x1=`"$xPosA`" y1=`"$bottom`" x2=`"$xPosA`" y2=`"$($bottom+5)`" stroke=`"#111`"/><text x=`"$xPosA`" y=`"$($bottom+20)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"11`">$(Tick $x[$i] 2)</text>")
        $svg.Add("<line x1=`"$xPosB`" y1=`"$bottom`" x2=`"$xPosB`" y2=`"$($bottom+5)`" stroke=`"#111`"/><text x=`"$xPosB`" y=`"$($bottom+20)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"11`">$(Tick $xB[$i] 2)</text>")
    }
    $svg.Add("<polyline points=`"$($pointsA -join ' ')`" fill=`"none`" stroke=`"#1769aa`" stroke-width=`"2.5`"/>");$svg.Add("<polyline points=`"$($pointsB -join ' ')`" fill=`"none`" stroke=`"#b33c2e`" stroke-width=`"2.5`"/>")
    foreach($p in $pointsA){$q=$p.Split(',');$svg.Add("<circle cx=`"$($q[0])`" cy=`"$($q[1])`" r=`"4.5`" fill=`"#1769aa`"/>")};foreach($p in $pointsB){$q=$p.Split(',');$svg.Add("<circle cx=`"$($q[0])`" cy=`"$($q[1])`" r=`"4.5`" fill=`"#b33c2e`"/>")}
    $svg.Add("<text x=`"$((($leftA+$rightA)/2)-45)`" y=`"440`" font-family=`"Arial`" font-size=`"15`">$(Xml $xlab)</text><text x=`"$((($leftB+$rightB)/2)-45)`" y=`"440`" font-family=`"Arial`" font-size=`"15`">$(Xml $xBlab)</text>")
    $svg.Add("<text transform=`"translate(20,285) rotate(-90)`" font-family=`"Arial`" font-size=`"15`">$(Xml $ylA)</text><text transform=`"translate(568,285) rotate(-90)`" font-family=`"Arial`" font-size=`"15`">$(Xml $ylB)</text></svg>")
    $svg|Set-Content $path -Encoding UTF8
}
Svg-Plot $EnergyFigurePath @($performance) energy
Svg-Plot $RuntimeFigurePath @($performance) runtime

if($surfaceRows.Count-gt0){$barW=76;$plotL=92;$plotB=350;$plotH=270;$svg=[Collections.Generic.List[string]]::new();$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="430" viewBox="0 0 980 430"><rect width="100%" height="100%" fill="white"/>');$svg.Add("<line x1=`"$plotL`" y1=`"$plotB`" x2=`"690`" y2=`"$plotB`" stroke=`"#111`"/><line x1=`"$plotL`" y1=`"80`" x2=`"$plotL`" y2=`"$plotB`" stroke=`"#111`"/>");for($tick=0;$tick-le4;++$tick){$f=$tick/4;$y=$plotB-$f*$plotH;$svg.Add("<line x1=`"$plotL`" y1=`"$y`" x2=`"690`" y2=`"$y`" stroke=`"#ddd`"/><text x=`"82`" y=`"$($y+4)`" text-anchor=`"end`" font-family=`"Arial`" font-size=`"12`">$(Tick $f 2)</text>")};$spacing=500/[math]::Max(1,$surfaceRows.Count);$parts=@(@('surface_seconds','#4c78a8','Surface'),@('ilu_setup_seconds','#f58518','ILU setup'),@('gmres_seconds','#54a24b','GMRES'),@('audit_seconds','#e45756','Direct audit'),@('other_seconds','#9d9da1','Other'));for($i=0;$i-lt$surfaceRows.Count;++$i){$r=$surfaceRows[$i];$x=$plotL+35+$i*$spacing;$y=$plotB;$accounted=0.0;foreach($part in $parts){$accounted+=D $r.($part[0]) 0};foreach($part in $parts){$fraction=(D $r.($part[0]) 0)/[math]::Max(1e-300,$accounted);$bh=$plotH*$fraction;$y-=$bh;$svg.Add("<rect x=`"$x`" y=`"$y`" width=`"$barW`" height=`"$bh`" fill=`"$($part[1])`"/>")};$svg.Add("<text x=`"$($x+$barW/2)`" y=`"375`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"14`">h=$($r.grid_spacing_A)</text>")};for($i=0;$i-lt$parts.Count;++$i){$y=105+$i*34;$svg.Add("<rect x=`"740`" y=`"$($y-14)`" width=`"22`" height=`"14`" fill=`"$($parts[$i][1])`"/><text x=`"772`" y=`"$y`" font-family=`"Arial`" font-size=`"14`">$($parts[$i][2])</text>")};$svg.Add('<text transform="translate(24,285) rotate(-90)" font-family="Arial" font-size="15">Fraction of accounted wall time</text></svg>');$svg|Set-Content $BreakdownFigurePath -Encoding UTF8}

$valid=@($rows|?{$_.status-eq'ok'-and($_.backend-eq'direct'-or([double]$_.true_residual-le1e-4-and[double]$_.audit_rms-le8e-3))})
$stable=@($convergence|?{$_.practical_stability_pass-eq$true})
$selectedSpacing=if($stable.Count){($stable|Sort-Object {[double]$_.coarse_spacing_A} -Descending|Select-Object -First 1).coarse_spacing_A}else{''}
$thetaRows=@($rows|?{$_.purpose-match'theta'-and$_.status-eq'ok'-and[double]$_.audit_rms-le8e-3-and[double]$_.true_residual-le1e-4});$selectedTheta=if($thetaRows.Count){($thetaRows|Sort-Object {[double]$_.total_seconds}|Select-Object -First 1).opening_theta}else{''}
$iluRows=@($rows|?{$_.purpose-match'ilu'-and$_.status-eq'ok'-and[double]$_.audit_rms-le8e-3-and[double]$_.true_residual-le1e-4-and[double]$_.peak_rss_mb-le512});$selectedCap=if($iluRows.Count){($iluRows|Sort-Object {[double]$_.total_seconds}|Select-Object -First 1).ilu_neighbor_cap}else{''}
$probePass=@($rows|?{$_.purpose-match'probe'-and$_.status-eq'ok'-and[double]$_.audit_rms-le8e-3-and[double]$_.true_residual-le1e-4}).Count-gt0
$axesComplete=$surfaceRows.Count-ge2-and$probePass-and$thetaRows.Count-gt0-and$iluRows.Count-gt0
$decision=if(-not$axesComplete){'calibration_incomplete'}elseif($stable.Count-eq0){'surface_not_stable'}elseif($valid.Count-lt$rows.Count){'calibration_policy_ready_with_rejected_alternatives'}else{'calibration_policy_ready'}
[pscustomobject]@{decision=$decision;selected_grid_spacing_A=$selectedSpacing;probe_radius_policy_A=$PolicyProbeRadius;selected_opening_theta=$selectedTheta;selected_ilu_neighbor_cap=$selectedCap;extrapolated_deltaG_kcal_mol=$energyInfinity;energy_fit_r_squared=$energyFitR2;finest_extrapolation_gap_percent=$finestExtrapolationGapPercent;valid_rows=$valid.Count;rejected_alternatives=$rows.Count-$valid.Count;total_rows=$rows.Count}|Export-Csv $DecisionPath -NoTypeInformation
@("decision=$decision","selected_grid_spacing_A=$selectedSpacing","probe_radius_policy_A=$PolicyProbeRadius","selected_opening_theta=$selectedTheta","selected_ilu_neighbor_cap=$selectedCap","extrapolated_deltaG_kcal_mol=$energyInfinity","energy_fit_r_squared=$energyFitR2","finest_extrapolation_gap_percent=$finestExtrapolationGapPercent","matrix=$MatrixPath")|Set-Content $SummaryPath -Encoding ASCII
@('TITAN TABI 1A63 Policy Calibration Report',"Decision: $decision",'Protocol: watertight grid SES; degree-one analytical FMM/M2L; analytic-singular self treatment; sparse near-field ILU(0); tol 1e-4; salt 0.15 M.',"Selected practical spacing: $selectedSpacing A","Probe-radius policy: $PolicyProbeRadius A (physical sensitivity, not speed tuning)","Selected theta: $selectedTheta","Selected ILU neighbor cap: $selectedCap","Empirical N->infinity energy: $energyInfinity kcal/mol","Linear 1/N fit R^2: $energyFitR2","Finest-mesh gap to extrapolation: $finestExtrapolationGapPercent percent",'The extrapolation is descriptive only; the finest mesh remains outside a strict asymptotic-convergence claim.','Outputs:',"  $MatrixPath","  $SurfaceTablePath","  $PerformanceTablePath","  $EnergyFigurePath","  $RuntimeFigurePath","  $BreakdownFigurePath")|Set-Content $ReportPath -Encoding UTF8
Write-Host "1A63 calibration decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Figures: $EnergyFigurePath ; $RuntimeFigurePath ; $BreakdownFigurePath"
if($decision-notin@('calibration_policy_ready','calibration_policy_ready_with_rejected_alternatives')-and-not$AllowFailedDecision){exit 2}
