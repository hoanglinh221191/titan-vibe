[CmdletBinding()]
param(
    [string]$PairMatrix = 'test/tmp/pbe_tabi_fagbi27_publication_figure5/fagbi27_publication_pairs.csv',
    [string]$PaperMatrix = 'test/tmp/pbe_tabi_fagbi27_analytic_ses_apbs/fagbi27_paper_table4_set.csv',
    [string]$WorkRoot = 'test/tmp/pbe_tabi_fagbi27_threeway_comparison',
    [string]$ReportOutput = 'TITAN_TABI_FAGBI27_THREEWAY_COMPARISON_REPORT.txt'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Invariant = [Globalization.CultureInfo]::InvariantCulture

function Resolve-PathValue([string]$Value) {
    if ([IO.Path]::IsPathRooted($Value)) { return [IO.Path]::GetFullPath($Value) }
    return [IO.Path]::GetFullPath((Join-Path (Get-Location) $Value))
}

function D([object]$Value) {
    $parsed = 0.0
    if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float, $Invariant, [ref]$parsed)) { return $parsed }
    throw "Expected a finite numeric value, got '$Value'."
}

function F([double]$Value, [string]$Pattern = '0.###') {
    return $Value.ToString($Pattern, $Invariant)
}

function Xml([object]$Value) { return [Security.SecurityElement]::Escape([string]$Value) }

function Is-Finite([double]$Value) {
    return -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}

function Mean([double[]]$Values) { return ($Values | Measure-Object -Average).Average }

function Median([double[]]$Values) {
    $sorted = @($Values | Sort-Object)
    if (($sorted.Count % 2) -eq 1) { return $sorted[[int](($sorted.Count - 1) / 2)] }
    return 0.5 * ($sorted[$sorted.Count / 2 - 1] + $sorted[$sorted.Count / 2])
}

function Pearson([double[]]$X, [double[]]$Y) {
    $mx = Mean $X
    $my = Mean $Y
    $numerator = 0.0
    $dx2 = 0.0
    $dy2 = 0.0
    for ($i = 0; $i -lt $X.Count; $i++) {
        $dx = $X[$i] - $mx
        $dy = $Y[$i] - $my
        $numerator += $dx * $dy
        $dx2 += $dx * $dx
        $dy2 += $dy * $dy
    }
    return $numerator / [math]::Sqrt($dx2 * $dy2)
}

function LX([double]$Value, [double]$Min, [double]$Max, [double]$Left, [double]$Right) {
    return $Left + ($Value - $Min) / ($Max - $Min) * ($Right - $Left)
}

function LY([double]$Value, [double]$Min, [double]$Max, [double]$Top, [double]$Bottom) {
    return $Bottom - ($Value - $Min) / ($Max - $Min) * ($Bottom - $Top)
}

function LogX([double]$Value, [double]$Min, [double]$Max, [double]$Left, [double]$Right) {
    return $Left + ([math]::Log10($Value) - [math]::Log10($Min)) / ([math]::Log10($Max) - [math]::Log10($Min)) * ($Right - $Left)
}

function LogY([double]$Value, [double]$Min, [double]$Max, [double]$Top, [double]$Bottom) {
    return $Bottom - ([math]::Log10($Value) - [math]::Log10($Min)) / ([math]::Log10($Max) - [math]::Log10($Min)) * ($Bottom - $Top)
}

function LogTicks([double]$Min, [double]$Max) {
    $ticks = [Collections.Generic.List[double]]::new()
    foreach ($power in ([int][math]::Floor([math]::Log10($Min)))..([int][math]::Ceiling([math]::Log10($Max)))) {
        foreach ($factor in @(1.0, 2.0, 5.0)) {
            $value = $factor * [math]::Pow(10.0, $power)
            if ($value -ge $Min -and $value -le $Max) { $ticks.Add($value) }
        }
    }
    return $ticks.ToArray()
}

function Add-Axes([Collections.Generic.List[string]]$Svg, [hashtable]$Panel, [double]$Top, [double]$Bottom, [string]$Letter, [double]$LetterY = 66.0) {
    $Svg.Add("<text class=`"panel`" x=`"$($Panel.left - 22)`" y=`"$LetterY`">$Letter</text>")
    $Svg.Add("<line class=`"axis`" x1=`"$($Panel.left)`" y1=`"$Top`" x2=`"$($Panel.left)`" y2=`"$Bottom`"/>")
    $Svg.Add("<line class=`"axis`" x1=`"$($Panel.left)`" y1=`"$Bottom`" x2=`"$($Panel.right)`" y2=`"$Bottom`"/>")
}

function Write-SvgExports([Collections.Generic.List[string]]$Svg, [string]$SvgPath) {
    $Svg.Add('</svg>')
    $Svg | Set-Content -LiteralPath $SvgPath -Encoding UTF8
    $inkscape = Get-Command inkscape.exe -ErrorAction SilentlyContinue
    if ($null -eq $inkscape) { throw 'inkscape.exe is required for publication exports.' }
    $inkscapeCli = Join-Path (Split-Path -Parent $inkscape.Source) 'inkscape.com'
    if (-not (Test-Path -LiteralPath $inkscapeCli)) { $inkscapeCli = $inkscape.Source }
    $png = [IO.Path]::ChangeExtension($SvgPath, '.png')
    $pdf = [IO.Path]::ChangeExtension($SvgPath, '.pdf')
    & $inkscapeCli $SvgPath "--export-filename=$png" '--export-dpi=600' | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "PNG export failed for $SvgPath" }
    & $inkscapeCli $SvgPath "--export-filename=$pdf" | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "PDF export failed for $SvgPath" }
}

$pairPath = Resolve-PathValue $PairMatrix
$paperPath = Resolve-PathValue $PaperMatrix
$workPath = Resolve-PathValue $WorkRoot
$reportPath = Resolve-PathValue $ReportOutput
foreach ($path in @($pairPath, $paperPath)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing required comparison input: $path" }
}
New-Item -ItemType Directory -Force -Path $workPath | Out-Null

$paperById = @{}
foreach ($row in @(Import-Csv -LiteralPath $paperPath)) { $paperById[[string]$row.pdb_id] = $row }
$comparison = [Collections.Generic.List[object]]::new()
foreach ($row in @(Import-Csv -LiteralPath $pairPath | Sort-Object { [int]$_.paper_index })) {
    $id = [string]$row.pdb_id
    if (-not $paperById.ContainsKey($id)) { throw "Missing FAGBI value for $id." }
    $paper = D $paperById[$id].paper_energy_kcal
    $titan = D $row.titan_deltaG_solv_kcal_mol
    $apbs = D $row.apbs_deltaG_solv_kcal_mol
    $comparison.Add([pscustomobject][ordered]@{
        paper_index = [int]$row.paper_index
        pdb_id = $id
        paper_elements = [int](D $paperById[$id].paper_elements)
        titan_panels = [int](D $row.panels)
        fagbi_2110_kcal_mol = $paper
        titan_kcal_mol = $titan
        apbs_kcal_mol = $apbs
        titan_peak_rss_mb = D $row.titan_peak_rss_mb
        apbs_peak_rss_mb = D $row.apbs_peak_rss_mb
        titan_wall_seconds = D $row.titan_wall_seconds
        apbs_wall_seconds = D $row.apbs_wall_seconds
        titan_minus_fagbi_kcal_mol = $titan - $paper
        apbs_minus_fagbi_kcal_mol = $apbs - $paper
        titan_minus_fagbi_percent = 100.0 * ($titan - $paper) / [math]::Abs($paper)
        apbs_minus_fagbi_percent = 100.0 * ($apbs - $paper) / [math]::Abs($paper)
    })
}
$rows = @($comparison.ToArray())
if ($rows.Count -ne 27) { throw "Expected 27 comparison rows; found $($rows.Count)." }
$comparisonPath = Join-Path $workPath 'fagbi27_threeway_comparison.csv'
$statisticsPath = Join-Path $workPath 'fagbi27_threeway_statistics.csv'
$outlierPath = Join-Path $workPath 'fagbi27_threeway_outliers.csv'
$rows | Export-Csv -LiteralPath $comparisonPath -NoTypeInformation

$paperValues = [double[]]@($rows | ForEach-Object { $_.fagbi_2110_kcal_mol })
$stats = [Collections.Generic.List[object]]::new()
foreach ($method in @(
    @{ name = 'TITAN'; value = 'titan_kcal_mol'; delta = 'titan_minus_fagbi_kcal_mol'; percent = 'titan_minus_fagbi_percent' },
    @{ name = 'APBS'; value = 'apbs_kcal_mol'; delta = 'apbs_minus_fagbi_kcal_mol'; percent = 'apbs_minus_fagbi_percent' }
)) {
    $values = [double[]]@($rows | ForEach-Object { $_.$($method.value) })
    $deltas = [double[]]@($rows | ForEach-Object { $_.$($method.delta) })
    $absDeltas = [double[]]@($deltas | ForEach-Object { [math]::Abs($_) })
    $absPercents = [double[]]@($rows | ForEach-Object { [math]::Abs($_.$($method.percent)) })
    $stats.Add([pscustomobject][ordered]@{
        method = $method.name
        n = $rows.Count
        pearson_r = Pearson $paperValues $values
        mean_signed_difference_kcal_mol = Mean $deltas
        mean_absolute_difference_kcal_mol = Mean $absDeltas
        root_mean_square_difference_kcal_mol = [math]::Sqrt((Mean ([double[]]@($deltas | ForEach-Object { $_ * $_ }))))
        median_absolute_percent_difference = Median $absPercents
        mean_absolute_percent_difference = Mean $absPercents
        maximum_absolute_percent_difference = ($absPercents | Measure-Object -Maximum).Maximum
    })
}
$stats.ToArray() | Export-Csv -LiteralPath $statisticsPath -NoTypeInformation
$rows | ForEach-Object {
    [pscustomobject]@{ pdb_id = $_.pdb_id; method = 'TITAN'; difference_kcal_mol = $_.titan_minus_fagbi_kcal_mol; absolute_percent_difference = [math]::Abs($_.titan_minus_fagbi_percent) }
    [pscustomobject]@{ pdb_id = $_.pdb_id; method = 'APBS'; difference_kcal_mol = $_.apbs_minus_fagbi_kcal_mol; absolute_percent_difference = [math]::Abs($_.apbs_minus_fagbi_percent) }
} | Sort-Object absolute_percent_difference -Descending | Export-Csv -LiteralPath $outlierPath -NoTypeInformation

$blue = '#1769AA'
$orange = '#D97706'
$black = '#30343B'
$gray = '#6B7280'
$grid = '#D9DDE3'
$baseStyle = '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.axis{stroke:#202124;stroke-width:2.3}.grid{stroke:#D9DDE3;stroke-width:1.2}.tick{font-size:23px}.label{font-size:28px}.panel{font-size:38px;font-weight:700}.legend{font-size:27px}</style>'

# Figure 1: direct parity against the values reported in 2110.13778.
$paritySvg = Join-Path $workPath 'fagbi27_threeway_parity.svg'
$svg = [Collections.Generic.List[string]]::new()
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="78mm" viewBox="0 0 1800 780"><rect width="100%" height="100%" fill="white"/>')
$svg.Add($baseStyle)
$top = 105.0; $bottom = 650.0
$parityPanels = @(@{left=115.0;right=825.0}, @{left=990.0;right=1700.0})
$allEnergy = [double[]]@($rows | ForEach-Object { $_.fagbi_2110_kcal_mol; $_.titan_kcal_mol; $_.apbs_kcal_mol })
$minEnergy = [double](($allEnergy | Measure-Object -Minimum).Minimum) - 120.0
$maxEnergy = [double](($allEnergy | Measure-Object -Maximum).Maximum) + 120.0
$step = 600.0
foreach ($i in 0..1) {
    $p = $parityPanels[$i]
    Add-Axes $svg $p $top $bottom ([string][char](65 + $i))
    for ($tick = [math]::Ceiling($minEnergy / $step) * $step; $tick -le $maxEnergy; $tick += $step) {
        $x = LX $tick $minEnergy $maxEnergy $p.left $p.right
        $y = LY $tick $minEnergy $maxEnergy $top $bottom
        $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/><line class=`"grid`" x1=`"$($p.left)`" y1=`"$y`" x2=`"$($p.right)`" y2=`"$y`"/>")
        $svg.Add("<text class=`"tick`" x=`"$x`" y=`"689`" text-anchor=`"middle`">$(F $tick '0')</text><text class=`"tick`" x=`"$($p.left - 14)`" y=`"$($y + 8)`" text-anchor=`"end`">$(F $tick '0')</text>")
    }
    $svg.Add("<line x1=`"$(LX $minEnergy $minEnergy $maxEnergy $p.left $p.right)`" y1=`"$(LY $minEnergy $minEnergy $maxEnergy $top $bottom)`" x2=`"$(LX $maxEnergy $minEnergy $maxEnergy $p.left $p.right)`" y2=`"$(LY $maxEnergy $minEnergy $maxEnergy $top $bottom)`" stroke=`"$gray`" stroke-width=`"2.2`" stroke-dasharray=`"10 8`"/>")
    $methodField = if ($i -eq 0) { 'titan_kcal_mol' } else { 'apbs_kcal_mol' }
    $percentField = if ($i -eq 0) { 'titan_minus_fagbi_percent' } else { 'apbs_minus_fagbi_percent' }
    $color = if ($i -eq 0) { $blue } else { $orange }
    $outlierIds = @($rows | Sort-Object { [math]::Abs($_.$percentField) } -Descending | Select-Object -First 3 -ExpandProperty pdb_id)
    foreach ($row in $rows) {
        $x = LX $row.fagbi_2110_kcal_mol $minEnergy $maxEnergy $p.left $p.right
        $y = LY $row.$methodField $minEnergy $maxEnergy $top $bottom
        if ($i -eq 0) { $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"8`" fill=`"$color`" fill-opacity=`"0.88`"/>") }
        else { $svg.Add("<rect x=`"$($x - 8)`" y=`"$($y - 8)`" width=`"16`" height=`"16`" fill=`"white`" stroke=`"$color`" stroke-width=`"3.2`"/>") }
        if ($row.pdb_id -in $outlierIds) { $svg.Add("<text x=`"$($x + 11)`" y=`"$($y - 10)`" font-size=`"21`" font-weight=`"600`">$(Xml $row.pdb_id)</text>") }
    }
    $methodLabel = if ($i -eq 0) { 'TITAN' } else { 'APBS' }
    $svg.Add("<text class=`"label`" x=`"$(($p.left + $p.right) / 2)`" y=`"748`" text-anchor=`"middle`">FAGBI &#916;G<tspan baseline-shift=`"sub`" font-size=`"19`">solv</tspan> (kcal mol&#8722;1)</text>")
    $svg.Add("<text class=`"label`" transform=`"translate($($p.left - 78),$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">$methodLabel &#916;G<tspan baseline-shift=`"sub`" font-size=`"19`">solv</tspan> (kcal mol&#8722;1)</text>")
}
Write-SvgExports $svg $paritySvg

# Combined three-way parity figure with identical scales in all panels.
$combinedSvg = Join-Path $workPath 'fagbi27_threeway_combined_parity.svg'
$svg = [Collections.Generic.List[string]]::new()
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="59mm" viewBox="0 0 2400 790"><rect width="100%" height="100%" fill="white"/>')
$svg.Add($baseStyle)
$top = 105.0; $bottom = 650.0
$combinedPanels = @(
    @{left=105.0;right=735.0},
    @{left=905.0;right=1535.0},
    @{left=1705.0;right=2335.0}
)
$combinedSpecs = @(
    @{ x='apbs_kcal_mol'; y='titan_kcal_mol'; xlabel='APBS'; ylabel='TITAN'; color=$blue; marker='circle'; outlier_metric={ param($r) [math]::Abs($r.titan_kcal_mol-$r.apbs_kcal_mol) } },
    @{ x='fagbi_2110_kcal_mol'; y='titan_kcal_mol'; xlabel='FAGBI'; ylabel='TITAN'; color=$blue; marker='circle'; outlier_metric={ param($r) [math]::Abs($r.titan_minus_fagbi_percent) } },
    @{ x='fagbi_2110_kcal_mol'; y='apbs_kcal_mol'; xlabel='FAGBI'; ylabel='APBS'; color=$orange; marker='square'; outlier_metric={ param($r) [math]::Abs($r.apbs_minus_fagbi_percent) } }
)
foreach ($i in 0..2) {
    $p = $combinedPanels[$i]
    $spec = $combinedSpecs[$i]
    Add-Axes $svg $p $top $bottom ([string][char](65 + $i))
    for ($tick = [math]::Ceiling($minEnergy / $step) * $step; $tick -le $maxEnergy; $tick += $step) {
        $x = LX $tick $minEnergy $maxEnergy $p.left $p.right
        $y = LY $tick $minEnergy $maxEnergy $top $bottom
        $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/><line class=`"grid`" x1=`"$($p.left)`" y1=`"$y`" x2=`"$($p.right)`" y2=`"$y`"/>")
        $svg.Add("<text class=`"tick`" x=`"$x`" y=`"689`" text-anchor=`"middle`">$(F $tick '0')</text><text class=`"tick`" x=`"$($p.left - 13)`" y=`"$($y + 8)`" text-anchor=`"end`">$(F $tick '0')</text>")
    }
    $svg.Add("<line x1=`"$(LX $minEnergy $minEnergy $maxEnergy $p.left $p.right)`" y1=`"$(LY $minEnergy $minEnergy $maxEnergy $top $bottom)`" x2=`"$(LX $maxEnergy $minEnergy $maxEnergy $p.left $p.right)`" y2=`"$(LY $maxEnergy $minEnergy $maxEnergy $top $bottom)`" stroke=`"$gray`" stroke-width=`"2.2`" stroke-dasharray=`"10 8`"/>")
    $outlier = @($rows | Sort-Object { & $spec.outlier_metric $_ } -Descending | Select-Object -First 1)[0]
    foreach ($row in $rows) {
        $x = LX $row.$($spec.x) $minEnergy $maxEnergy $p.left $p.right
        $y = LY $row.$($spec.y) $minEnergy $maxEnergy $top $bottom
        if ($spec.marker -eq 'circle') {
            $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"7.5`" fill=`"$($spec.color)`" fill-opacity=`"0.87`"/>")
        } else {
            $svg.Add("<rect x=`"$($x - 7.5)`" y=`"$($y - 7.5)`" width=`"15`" height=`"15`" fill=`"white`" stroke=`"$($spec.color)`" stroke-width=`"3`"/>")
        }
        if ($row.pdb_id -eq $outlier.pdb_id) {
            $svg.Add("<text x=`"$($x + 10)`" y=`"$($y - 10)`" font-size=`"21`" font-weight=`"600`">$(Xml $row.pdb_id)</text>")
        }
    }
    $center = ($p.left + $p.right) / 2
    $svg.Add("<text class=`"label`" x=`"$center`" y=`"750`" text-anchor=`"middle`">$($spec.xlabel) &#916;G<tspan baseline-shift=`"sub`" font-size=`"19`">solv</tspan> (kcal mol&#8722;1)</text>")
    $svg.Add("<text class=`"label`" transform=`"translate($($p.left - 73),$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">$($spec.ylabel) &#916;G<tspan baseline-shift=`"sub`" font-size=`"19`">solv</tspan> (kcal mol&#8722;1)</text>")
}
Write-SvgExports $svg $combinedSvg

# Five-panel manuscript figure: three energy comparisons plus memory and wall time.
$fivePanelSvg = Join-Path $workPath 'fagbi27_threeway_five_panel_comparison.svg'
$svg = [Collections.Generic.List[string]]::new()
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="112mm" viewBox="0 0 2400 1500"><rect width="100%" height="100%" fill="white"/>')
$svg.Add($baseStyle)
$svg.Add('<style>.tick{font-size:33px}.label{font-size:39px}.panel{font-size:46px}.legend{font-size:36px}</style>')
$energyTop = 92.0; $energyBottom = 650.0
$topPanels = @(
    @{left=155.0;right=745.0},
    @{left=955.0;right=1545.0},
    @{left=1755.0;right=2345.0}
)
foreach ($i in 0..2) {
    $p = $topPanels[$i]
    $spec = $combinedSpecs[$i]
    Add-Axes $svg $p $energyTop $energyBottom ([string][char](65 + $i))
    for ($tick = [math]::Ceiling($minEnergy / $step) * $step; $tick -le $maxEnergy; $tick += $step) {
        $x = LX $tick $minEnergy $maxEnergy $p.left $p.right
        $y = LY $tick $minEnergy $maxEnergy $energyTop $energyBottom
        $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$energyTop`" x2=`"$x`" y2=`"$energyBottom`"/><line class=`"grid`" x1=`"$($p.left)`" y1=`"$y`" x2=`"$($p.right)`" y2=`"$y`"/>")
        $svg.Add("<text class=`"tick`" x=`"$x`" y=`"689`" text-anchor=`"middle`">$(F $tick '0')</text><text class=`"tick`" x=`"$($p.left - 13)`" y=`"$($y + 8)`" text-anchor=`"end`">$(F $tick '0')</text>")
    }
    $svg.Add("<line x1=`"$(LX $minEnergy $minEnergy $maxEnergy $p.left $p.right)`" y1=`"$(LY $minEnergy $minEnergy $maxEnergy $energyTop $energyBottom)`" x2=`"$(LX $maxEnergy $minEnergy $maxEnergy $p.left $p.right)`" y2=`"$(LY $maxEnergy $minEnergy $maxEnergy $energyTop $energyBottom)`" stroke=`"$gray`" stroke-width=`"2.2`" stroke-dasharray=`"10 8`"/>")
    $outlier = @($rows | Sort-Object { & $spec.outlier_metric $_ } -Descending | Select-Object -First 1)[0]
    foreach ($row in $rows) {
        $x = LX $row.$($spec.x) $minEnergy $maxEnergy $p.left $p.right
        $y = LY $row.$($spec.y) $minEnergy $maxEnergy $energyTop $energyBottom
        if ($spec.marker -eq 'circle') {
            $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"7.5`" fill=`"$($spec.color)`" fill-opacity=`"0.87`"/>")
        } else {
            $svg.Add("<rect x=`"$($x - 7.5)`" y=`"$($y - 7.5)`" width=`"15`" height=`"15`" fill=`"white`" stroke=`"$($spec.color)`" stroke-width=`"3`"/>")
        }
        if ($row.pdb_id -eq $outlier.pdb_id) {
            $svg.Add("<text x=`"$($x + 12)`" y=`"$($y - 12)`" font-size=`"30`" font-weight=`"600`">$(Xml $row.pdb_id)</text>")
        }
    }
    $center = ($p.left + $p.right) / 2
    $svg.Add("<text class=`"label`" x=`"$center`" y=`"750`" text-anchor=`"middle`">$($spec.xlabel) &#916;G<tspan baseline-shift=`"sub`" font-size=`"27`">solv</tspan> (kcal mol&#8722;1)</text>")
    $svg.Add("<text class=`"label`" transform=`"translate($($p.left - 112),$(($energyTop + $energyBottom) / 2)) rotate(-90)`" text-anchor=`"middle`">$($spec.ylabel) &#916;G<tspan baseline-shift=`"sub`" font-size=`"27`">solv</tspan> (kcal mol&#8722;1)</text>")
}

$performanceTop = 860.0; $performanceBottom = 1370.0
$performancePanels = @(
    @{left=235.0;right=1085.0;letter='D';x='apbs_peak_rss_mb';y='titan_peak_rss_mb';xlabel='APBS peak RSS (MB)';ylabel='TITAN peak RSS (MB)'},
    @{left=1415.0;right=2265.0;letter='E';x='apbs_wall_seconds';y='titan_wall_seconds';xlabel='APBS total wall time (s)';ylabel='TITAN total wall time (s)'}
)
foreach ($spec in $performancePanels) {
    Add-Axes $svg $spec $performanceTop $performanceBottom $spec.letter 824.0
    $values = [double[]]@($rows | ForEach-Object { $_.$($spec.x); $_.$($spec.y) })
    $minimum = [double](($values | Measure-Object -Minimum).Minimum) * 0.78
    $maximum = [double](($values | Measure-Object -Maximum).Maximum) * 1.28
    foreach ($tick in @(LogTicks $minimum $maximum)) {
        $x = LogX $tick $minimum $maximum $spec.left $spec.right
        $y = LogY $tick $minimum $maximum $performanceTop $performanceBottom
        $label = if ($tick -lt 10) { F $tick '0.#' } else { F $tick '0' }
        $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$performanceTop`" x2=`"$x`" y2=`"$performanceBottom`"/><line class=`"grid`" x1=`"$($spec.left)`" y1=`"$y`" x2=`"$($spec.right)`" y2=`"$y`"/>")
        $svg.Add("<text class=`"tick`" x=`"$x`" y=`"1409`" text-anchor=`"middle`">$label</text><text class=`"tick`" x=`"$($spec.left - 14)`" y=`"$($y + 8)`" text-anchor=`"end`">$label</text>")
    }
    $svg.Add("<line x1=`"$(LogX $minimum $minimum $maximum $spec.left $spec.right)`" y1=`"$(LogY $minimum $minimum $maximum $performanceTop $performanceBottom)`" x2=`"$(LogX $maximum $minimum $maximum $spec.left $spec.right)`" y2=`"$(LogY $maximum $minimum $maximum $performanceTop $performanceBottom)`" stroke=`"$gray`" stroke-width=`"2.2`" stroke-dasharray=`"10 8`"/>")
    foreach ($row in $rows) {
        $x = LogX $row.$($spec.x) $minimum $maximum $spec.left $spec.right
        $y = LogY $row.$($spec.y) $minimum $maximum $performanceTop $performanceBottom
        $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"7.5`" fill=`"$blue`" fill-opacity=`"0.86`"/>")
    }
    $center = ($spec.left + $spec.right) / 2
    $svg.Add("<text class=`"label`" x=`"$center`" y=`"1470`" text-anchor=`"middle`">$(Xml $spec.xlabel)</text>")
    $svg.Add("<text class=`"label`" transform=`"translate($($spec.left - 108),$(($performanceTop + $performanceBottom) / 2)) rotate(-90)`" text-anchor=`"middle`">$(Xml $spec.ylabel)</text>")
}
Write-SvgExports $svg $fivePanelSvg

# Figure 2: signed deviation from FAGBI for every protein.
$deviationSvg = Join-Path $workPath 'fagbi27_threeway_signed_deviation.svg'
$svg = [Collections.Generic.List[string]]::new()
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="78mm" viewBox="0 0 1900 820"><rect width="100%" height="100%" fill="white"/>')
$svg.Add($baseStyle)
$left = 105.0; $right = 1860.0; $top = 105.0; $bottom = 635.0
$svg.Add("<line class=`"axis`" x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$bottom`"/><line class=`"axis`" x1=`"$left`" y1=`"$bottom`" x2=`"$right`" y2=`"$bottom`"/>")
$maxAbsPercent = [math]::Ceiling(([double](($rows | ForEach-Object { [math]::Abs($_.titan_minus_fagbi_percent); [math]::Abs($_.apbs_minus_fagbi_percent) } | Measure-Object -Maximum).Maximum) + 2.0) / 5.0) * 5.0
$minPercent = -$maxAbsPercent; $maxPercent = $maxAbsPercent
for ($tick = $minPercent; $tick -le $maxPercent; $tick += 10.0) {
    $y = LY $tick $minPercent $maxPercent $top $bottom
    $stroke = if ([math]::Abs($tick) -lt 1e-9) { $gray } else { $grid }
    $dash = if ([math]::Abs($tick) -lt 1e-9) { ' stroke-dasharray="9 7"' } else { '' }
    $svg.Add("<line x1=`"$left`" y1=`"$y`" x2=`"$right`" y2=`"$y`" stroke=`"$stroke`" stroke-width=`"1.7`"$dash/>")
    $svg.Add("<text class=`"tick`" x=`"$($left - 15)`" y=`"$($y + 8)`" text-anchor=`"end`">$(F $tick '0')</text>")
}
foreach ($row in $rows) {
    $x = LX $row.paper_index 0.5 27.5 $left $right
    $ty = LY $row.titan_minus_fagbi_percent $minPercent $maxPercent $top $bottom
    $ay = LY $row.apbs_minus_fagbi_percent $minPercent $maxPercent $top $bottom
    $svg.Add("<line x1=`"$x`" y1=`"$ty`" x2=`"$x`" y2=`"$ay`" stroke=`"#C3C7CE`" stroke-width=`"1.4`"/>")
    $svg.Add("<circle cx=`"$x`" cy=`"$ty`" r=`"7`" fill=`"$blue`"/>")
    $svg.Add("<rect x=`"$($x - 7)`" y=`"$($ay - 7)`" width=`"14`" height=`"14`" fill=`"white`" stroke=`"$orange`" stroke-width=`"3`"/>")
    $svg.Add("<text class=`"tick`" transform=`"translate($($x + 7),$($bottom + 28)) rotate(58)`" text-anchor=`"start`">$(Xml $row.pdb_id)</text>")
}
$svg.Add("<circle cx=`"715`" cy=`"57`" r=`"8`" fill=`"$blue`"/><text class=`"legend`" x=`"735`" y=`"67`">TITAN</text><rect x=`"960`" y=`"49`" width=`"16`" height=`"16`" fill=`"white`" stroke=`"$orange`" stroke-width=`"3`"/><text class=`"legend`" x=`"991`" y=`"67`">APBS</text>")
$svg.Add("<text class=`"label`" transform=`"translate(36,$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">Difference from FAGBI (%)</text>")
Write-SvgExports $svg $deviationSvg

# Figure 3: all three reported energy series across the 27 systems.
$profileSvg = Join-Path $workPath 'fagbi27_threeway_energy_profile.svg'
$svg = [Collections.Generic.List[string]]::new()
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="78mm" viewBox="0 0 1900 820"><rect width="100%" height="100%" fill="white"/>')
$svg.Add($baseStyle)
$left = 105.0; $right = 1860.0; $top = 105.0; $bottom = 635.0
$svg.Add("<line class=`"axis`" x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$bottom`"/><line class=`"axis`" x1=`"$left`" y1=`"$bottom`" x2=`"$right`" y2=`"$bottom`"/>")
$minEnergy = [math]::Floor(([double](($allEnergy | Measure-Object -Minimum).Minimum) - 100.0) / 500.0) * 500.0
$maxEnergy = 0.0
for ($tick = $minEnergy; $tick -le $maxEnergy; $tick += 500.0) {
    $y = LY $tick $minEnergy $maxEnergy $top $bottom
    $svg.Add("<line class=`"grid`" x1=`"$left`" y1=`"$y`" x2=`"$right`" y2=`"$y`"/><text class=`"tick`" x=`"$($left - 15)`" y=`"$($y + 8)`" text-anchor=`"end`">$(F $tick '0')</text>")
}
$series = @(
    @{ field='fagbi_2110_kcal_mol'; color=$black; marker='diamond'; label='FAGBI' },
    @{ field='titan_kcal_mol'; color=$blue; marker='circle'; label='TITAN' },
    @{ field='apbs_kcal_mol'; color=$orange; marker='square'; label='APBS' }
)
foreach ($s in $series) {
    $points = [Collections.Generic.List[string]]::new()
    foreach ($row in $rows) {
        $x = LX $row.paper_index 0.5 27.5 $left $right
        $y = LY $row.$($s.field) $minEnergy $maxEnergy $top $bottom
        $points.Add("$x,$y")
    }
    $svg.Add("<polyline points=`"$($points -join ' ')`" fill=`"none`" stroke=`"$($s.color)`" stroke-width=`"2.5`" stroke-opacity=`"0.76`"/>")
    foreach ($row in $rows) {
        $x = LX $row.paper_index 0.5 27.5 $left $right
        $y = LY $row.$($s.field) $minEnergy $maxEnergy $top $bottom
        if ($s.marker -eq 'circle') { $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"6.5`" fill=`"$($s.color)`"/>") }
        elseif ($s.marker -eq 'square') { $svg.Add("<rect x=`"$($x - 6)`" y=`"$($y - 6)`" width=`"12`" height=`"12`" fill=`"white`" stroke=`"$($s.color)`" stroke-width=`"2.8`"/>") }
        else { $svg.Add("<path d=`"M $x $($y - 7) L $($x + 7) $y L $x $($y + 7) L $($x - 7) $y Z`" fill=`"$($s.color)`"/>") }
    }
}
foreach ($row in $rows) {
    $x = LX $row.paper_index 0.5 27.5 $left $right
    $svg.Add("<text class=`"tick`" transform=`"translate($($x + 7),$($bottom + 28)) rotate(58)`" text-anchor=`"start`">$(Xml $row.pdb_id)</text>")
}
$svg.Add("<path d=`"M 700 57 L 708 49 L 716 57 L 708 65 Z`" fill=`"$black`"/><text class=`"legend`" x=`"729`" y=`"67`">FAGBI</text><circle cx=`"959`" cy=`"57`" r=`"7`" fill=`"$blue`"/><text class=`"legend`" x=`"979`" y=`"67`">TITAN</text><rect x=`"1204`" y=`"50`" width=`"14`" height=`"14`" fill=`"white`" stroke=`"$orange`" stroke-width=`"3`"/><text class=`"legend`" x=`"1233`" y=`"67`">APBS</text>")
$svg.Add("<text class=`"label`" transform=`"translate(36,$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">&#916;G<tspan baseline-shift=`"sub`" font-size=`"19`">solv</tspan> (kcal mol&#8722;1)</text>")
Write-SvgExports $svg $profileSvg

$captionPath = Join-Path $workPath 'fagbi27_threeway_figure_captions.txt'
@(
    'Figure. Three-way energy comparison and computational footprint over the 27-protein set. (A) TITAN versus APBS, (B) TITAN versus FAGBI, and (C) APBS versus FAGBI electrostatic solvation energies. (D) Peak resident-set size and (E) total wall time for TITAN versus APBS; both programs used four CPU threads. All pairwise panels use matched axis limits and dashed equality lines; panels D and E use logarithmic axes. Labels in A-C identify the largest deviation in each comparison. FAGBI values are from arXiv:2110.13778. Cross-implementation energy comparisons are contextual because molecular-surface construction and discretization differ.'
    'Figure. Three-way pairwise comparison of electrostatic solvation energies for the 27-protein set. (A) TITAN versus APBS, (B) TITAN versus FAGBI, and (C) APBS versus FAGBI. All panels use identical limits; dashed lines denote equality. Labels identify the largest deviation in each comparison. The FAGBI values are those reported in arXiv:2110.13778. The comparison is contextual because molecular-surface construction and discretization differ among implementations.'
    'Figure A. Pairwise comparison against the electrostatic solvation energies reported for the 27-protein FAGBI set in arXiv:2110.13778. The dashed line denotes equality. Labels identify the three largest relative deviations in each panel.'
    'Figure B. Signed percentage differences of TITAN and APBS energies from the FAGBI values for each protein. Grey segments connect the two calculations for the same system; zero denotes agreement with FAGBI.'
    'Figure C. Electrostatic solvation-energy profiles for FAGBI, TITAN, and APBS across the 27 proteins. These are contextual cross-implementation comparisons because molecular-surface construction and discretization differ among the calculations.'
) | Set-Content -LiteralPath $captionPath -Encoding UTF8

$decision = if (@($rows | Where-Object { -not (Is-Finite ([double]$_.titan_kcal_mol)) -or -not (Is-Finite ([double]$_.apbs_kcal_mol)) }).Count -gt 0) { 'fagbi27_threeway_data_gap' } else { 'fagbi27_threeway_figures_ready' }
$decisionPath = Join-Path $workPath 'fagbi27_threeway_decision.csv'
[pscustomobject]@{decision=$decision;systems=$rows.Count;comparison_csv=$comparisonPath;statistics_csv=$statisticsPath;five_panel_figure=$fivePanelSvg;combined_parity_figure=$combinedSvg;parity_figure=$paritySvg;signed_deviation_figure=$deviationSvg;energy_profile_figure=$profileSvg} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation

$titanStats = @($stats | Where-Object method -eq 'TITAN')[0]
$apbsStats = @($stats | Where-Object method -eq 'APBS')[0]
@(
    'TITAN FAGBI-27 Three-Way Comparison Report'
    "Decision: $decision"
    'Reference values: the 27 electrostatic solvation energies reported in arXiv:2110.13778.'
    "TITAN mean absolute percent difference: $(F $titanStats.mean_absolute_percent_difference '0.###')%"
    "APBS mean absolute percent difference: $(F $apbsStats.mean_absolute_percent_difference '0.###')%"
    "TITAN Pearson r: $(F $titanStats.pearson_r '0.#####')"
    "APBS Pearson r: $(F $apbsStats.pearson_r '0.#####')"
    'The comparison is contextual rather than an oracle test because the published FAGBI, TITAN, and APBS calculations do not share one identical surface discretization.'
    "Comparison table: $comparisonPath"
    "Statistics: $statisticsPath"
    "Five-panel manuscript figure: $fivePanelSvg"
    "Combined three-way parity figure: $combinedSvg"
    "Parity figure: $paritySvg"
    "Signed-deviation figure: $deviationSvg"
    "Energy-profile figure: $profileSvg"
) | Set-Content -LiteralPath $reportPath -Encoding UTF8

Write-Host "Decision: $decision"
Write-Host "Comparison: $comparisonPath"
Write-Host "Five-panel figure: $fivePanelSvg"
Write-Host "Combined figure: $combinedSvg"
Write-Host "Parity figure: $paritySvg"
Write-Host "Signed-deviation figure: $deviationSvg"
Write-Host "Energy-profile figure: $profileSvg"

if ($decision -ne 'fagbi27_threeway_figures_ready') { throw "Three-way comparison did not pass: $decision" }
