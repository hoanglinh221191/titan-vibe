[CmdletBinding()]
param(
    [string]$OriginalMatrix = 'test/tmp/pbe_tabi_fagbi27_analytic_ses_apbs/fagbi27_analytic_ses_apbs_matrix.csv',
    [string]$RecoveryRoot = 'test/tmp/pbe_tabi_fagbi27_recovery_v2',
    [string]$WorkRoot = 'test/tmp/pbe_tabi_fagbi27_publication_figure5',
    [string]$ReportOutput = 'TITAN_TABI_FAGBI27_PUBLICATION_FIGURE5_REPORT.txt'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Resolve-RepoPath([string]$PathValue) {
    if ([IO.Path]::IsPathRooted($PathValue)) { return [IO.Path]::GetFullPath($PathValue) }
    return [IO.Path]::GetFullPath((Join-Path (Get-Location) $PathValue))
}

function Parse-Double([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return $Default }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) { return $parsed }
    return $Default
}

function Is-Finite([double]$Value) {
    return -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}

function Format-Number([double]$Value, [string]$Pattern = '0.###') {
    return $Value.ToString($Pattern, [Globalization.CultureInfo]::InvariantCulture)
}

function Xml-Escape([object]$Value) {
    return [Security.SecurityElement]::Escape([string]$Value)
}

function Copy-RowWithSource([object]$Row, [string]$Source) {
    $ordered = [ordered]@{}
    foreach ($property in $Row.PSObject.Properties) { $ordered[$property.Name] = $property.Value }
    $ordered['publication_data_source'] = $Source
    return [pscustomobject]$ordered
}

function Get-SingleRow([object[]]$Rows, [string]$PdbId, [string]$Backend, [string]$Description) {
    $matched = @($Rows | Where-Object { [string]$_.pdb_id -eq $PdbId -and [string]$_.backend -eq $Backend })
    if ($matched.Count -ne 1) { throw "Expected exactly one $Description row for $PdbId/$Backend; found $($matched.Count)." }
    return $matched[0]
}

$originalPath = Resolve-RepoPath $OriginalMatrix
$recoveryRootPath = Resolve-RepoPath $RecoveryRoot
$workRootPath = Resolve-RepoPath $WorkRoot
$reportPath = Resolve-RepoPath $ReportOutput
$cleanMatrixPath = Join-Path $recoveryRootPath '1cbn_clean_fixture/fagbi27_analytic_ses_apbs_matrix.csv'
$regularizedMatrixPath = Join-Path $recoveryRootPath 'regularized_original_pqr/fagbi27_analytic_ses_apbs_matrix.csv'

foreach ($requiredPath in @($originalPath, $cleanMatrixPath, $regularizedMatrixPath)) {
    if (-not (Test-Path -LiteralPath $requiredPath)) { throw "Missing required FAGBI-27 artifact: $requiredPath" }
}

New-Item -ItemType Directory -Force -Path $workRootPath | Out-Null
$mergedPath = Join-Path $workRootPath 'fagbi27_publication_merged_matrix.csv'
$pairsPath = Join-Path $workRootPath 'fagbi27_publication_pairs.csv'
$supportingPath = Join-Path $workRootPath 'fagbi27_publication_supporting_table.csv'
$decisionPath = Join-Path $workRootPath 'fagbi27_publication_decision.csv'
$summaryPath = Join-Path $workRootPath 'fagbi27_publication_summary.dat'
$captionPath = Join-Path $workRootPath 'fagbi27_publication_figure5_caption.txt'
$figureSvgPath = Join-Path $workRootPath 'fagbi27_publication_figure5.svg'
$figurePngPath = Join-Path $workRootPath 'fagbi27_publication_figure5.png'
$figurePdfPath = Join-Path $workRootPath 'fagbi27_publication_figure5.pdf'

$originalRows = @(Import-Csv -LiteralPath $originalPath)
$cleanRows = @(Import-Csv -LiteralPath $cleanMatrixPath)
$regularizedRows = @(Import-Csv -LiteralPath $regularizedMatrixPath)
if ($originalRows.Count -ne 54) { throw "Expected 54 rows in the original matrix; found $($originalRows.Count)." }

$replacementRows = @{
    '1CBN|titan' = Get-SingleRow $cleanRows '1CBN' 'titan' 'clean-fixture'
    '1CBN|apbs' = Get-SingleRow $cleanRows '1CBN' 'apbs' 'clean-fixture'
    '2ERL|titan' = Get-SingleRow $regularizedRows '2ERL' 'titan' 'recovered'
    '1HPT|titan' = Get-SingleRow $regularizedRows '1HPT' 'titan' 'recovered'
}

$merged = [Collections.Generic.List[object]]::new()
foreach ($row in ($originalRows | Sort-Object { [int]$_.paper_index }, backend)) {
    $key = "$([string]$row.pdb_id)|$([string]$row.backend)"
    if ($replacementRows.ContainsKey($key)) {
        $merged.Add((Copy-RowWithSource $replacementRows[$key] 'recovery_verification'))
    } else {
        $merged.Add((Copy-RowWithSource $row 'original_full_matrix'))
    }
}
$merged.ToArray() | Export-Csv -LiteralPath $mergedPath -NoTypeInformation

$pairs = [Collections.Generic.List[object]]::new()
foreach ($index in 1..27) {
    $titan = @($merged | Where-Object { [int]$_.paper_index -eq $index -and [string]$_.backend -eq 'titan' })
    $apbs = @($merged | Where-Object { [int]$_.paper_index -eq $index -and [string]$_.backend -eq 'apbs' })
    if ($titan.Count -ne 1 -or $apbs.Count -ne 1) { throw "Incomplete merged pair at paper index $index." }
    $t = $titan[0]
    $a = $apbs[0]
    $tEnergy = Parse-Double $t.energy_kcal_mol
    $aEnergy = Parse-Double $a.energy_kcal_mol
    $delta = $tEnergy - $aEnergy
    $relative = if ([math]::Abs($aEnergy) -gt 1e-14) { 100.0 * $delta / [math]::Abs($aEnergy) } else { [double]::NaN }
    $tPass = [string]$t.status -eq 'ok' -and (Is-Finite $tEnergy) -and (Parse-Double $t.true_residual 1.0) -le 1e-4
    $aPass = [string]$a.status -eq 'ok' -and (Is-Finite $aEnergy)
    $pairs.Add([pscustomobject][ordered]@{
        paper_index = $index
        pdb_id = [string]$t.pdb_id
        atom_count = [int](Parse-Double $t.atom_count 0)
        panels = [int](Parse-Double $t.panels 0)
        titan_deltaG_solv_kcal_mol = $tEnergy
        apbs_deltaG_solv_kcal_mol = $aEnergy
        titan_minus_apbs_kcal_mol = $delta
        titan_minus_apbs_percent = $relative
        titan_iterations = [int](Parse-Double $t.iterations 0)
        titan_true_residual = Parse-Double $t.true_residual
        titan_audit_rms = Parse-Double $t.audit_rms
        titan_wall_seconds = Parse-Double $t.process_wall_seconds
        apbs_wall_seconds = Parse-Double $a.process_wall_seconds
        titan_peak_rss_mb = Parse-Double $t.peak_rss_mb
        apbs_peak_rss_mb = Parse-Double $a.peak_rss_mb
        status = if ($tPass -and $aPass) { 'ok' } elseif (-not $tPass) { 'titan_gap' } else { 'apbs_gap' }
        titan_data_source = [string]$t.publication_data_source
        apbs_data_source = [string]$a.publication_data_source
    })
}
$pairRows = @($pairs.ToArray())
$pairRows | Export-Csv -LiteralPath $pairsPath -NoTypeInformation
$pairRows | Select-Object paper_index,pdb_id,atom_count,panels,titan_deltaG_solv_kcal_mol,apbs_deltaG_solv_kcal_mol,titan_minus_apbs_kcal_mol,titan_minus_apbs_percent,titan_iterations,titan_true_residual,titan_audit_rms,titan_wall_seconds,apbs_wall_seconds,titan_peak_rss_mb,apbs_peak_rss_mb,status | Export-Csv -LiteralPath $supportingPath -NoTypeInformation

function Linear-X([double]$Value, [double]$Minimum, [double]$Maximum, [double]$Left, [double]$Right) {
    return $Left + ($Value - $Minimum) / ($Maximum - $Minimum) * ($Right - $Left)
}

function Linear-Y([double]$Value, [double]$Minimum, [double]$Maximum, [double]$Top, [double]$Bottom) {
    return $Bottom - ($Value - $Minimum) / ($Maximum - $Minimum) * ($Bottom - $Top)
}

function Log-X([double]$Value, [double]$Minimum, [double]$Maximum, [double]$Left, [double]$Right) {
    return $Left + ([math]::Log10($Value) - [math]::Log10($Minimum)) / ([math]::Log10($Maximum) - [math]::Log10($Minimum)) * ($Right - $Left)
}

function Log-Y([double]$Value, [double]$Minimum, [double]$Maximum, [double]$Top, [double]$Bottom) {
    return $Bottom - ([math]::Log10($Value) - [math]::Log10($Minimum)) / ([math]::Log10($Maximum) - [math]::Log10($Minimum)) * ($Bottom - $Top)
}

function Nice-Log-Ticks([double]$Minimum, [double]$Maximum) {
    $ticks = [Collections.Generic.List[double]]::new()
    $start = [math]::Floor([math]::Log10($Minimum))
    $end = [math]::Ceiling([math]::Log10($Maximum))
    foreach ($power in $start..$end) {
        foreach ($factor in @(1.0, 2.0, 5.0)) {
            $value = $factor * [math]::Pow(10.0, $power)
            if ($value -ge $Minimum -and $value -le $Maximum) { $ticks.Add($value) }
        }
    }
    return $ticks.ToArray()
}

function Write-Publication-Figure([object[]]$Rows) {
    $blue = '#1769AA'
    $orange = '#D97706'
    $ink = '#202124'
    $grid = '#D9DDE3'
    $width = 2280
    $height = 790
    $top = 128
    $bottom = 650
    $panels = @(
        @{ left = 112; right = 720 },
        @{ left = 860; right = 1468 },
        @{ left = 1610; right = 2218 }
    )
    $svg = [Collections.Generic.List[string]]::new()
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"180mm`" height=`"62.4mm`" viewBox=`"0 0 $width $height`">")
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    $svg.Add('<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.axis{stroke:#202124;stroke-width:2.4}.grid{stroke:#D9DDE3;stroke-width:1.25}.tick{font-size:24px}.label{font-size:29px}.panel{font-size:38px;font-weight:700}</style>')
    foreach ($panelIndex in 0..2) {
        $panel = $panels[$panelIndex]
        $letter = [char](65 + $panelIndex)
        $svg.Add("<text class=`"panel`" x=`"$($panel.left - 22)`" y=`"74`">$letter</text>")
        $svg.Add("<line class=`"axis`" x1=`"$($panel.left)`" y1=`"$top`" x2=`"$($panel.left)`" y2=`"$bottom`"/>")
        $svg.Add("<line class=`"axis`" x1=`"$($panel.left)`" y1=`"$bottom`" x2=`"$($panel.right)`" y2=`"$bottom`"/>")
    }
    $svg.Add("<circle cx=`"965`" cy=`"62`" r=`"9`" fill=`"$blue`"/><text x=`"986`" y=`"72`" font-size=`"29`">TITAN</text>")
    $svg.Add("<rect x=`"1207`" y=`"53`" width=`"18`" height=`"18`" fill=`"white`" stroke=`"$orange`" stroke-width=`"4`"/><text x=`"1240`" y=`"72`" font-size=`"29`">APBS</text>")

    $energyValues = @($Rows | ForEach-Object { [double]$_.titan_deltaG_solv_kcal_mol; [double]$_.apbs_deltaG_solv_kcal_mol })
    $energyMin = [double](($energyValues | Measure-Object -Minimum).Minimum)
    $energyMax = [double](($energyValues | Measure-Object -Maximum).Maximum)
    $energyPad = [math]::Max(50.0, 0.06 * ($energyMax - $energyMin))
    $energyMin -= $energyPad
    $energyMax += $energyPad
    $energyStep = [math]::Ceiling(($energyMax - $energyMin) / 4.0 / 100.0) * 100.0
    $energyTickMin = [math]::Floor($energyMin / $energyStep) * $energyStep
    $energyTickMax = [math]::Ceiling($energyMax / $energyStep) * $energyStep
    $p = $panels[0]
    for ($value = $energyTickMin; $value -le $energyTickMax + 0.5 * $energyStep; $value += $energyStep) {
        if ($value -lt $energyMin -or $value -gt $energyMax) { continue }
        $x = Linear-X $value $energyMin $energyMax $p.left $p.right
        $y = Linear-Y $value $energyMin $energyMax $top $bottom
        $label = Format-Number $value '0'
        $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/>")
        $svg.Add("<line class=`"grid`" x1=`"$($p.left)`" y1=`"$y`" x2=`"$($p.right)`" y2=`"$y`"/>")
        $svg.Add("<text class=`"tick`" x=`"$x`" y=`"690`" text-anchor=`"middle`">$label</text>")
        $svg.Add("<text class=`"tick`" x=`"$($p.left - 15)`" y=`"$($y + 8)`" text-anchor=`"end`">$label</text>")
    }
    $x1 = Linear-X $energyMin $energyMin $energyMax $p.left $p.right
    $y1 = Linear-Y $energyMin $energyMin $energyMax $top $bottom
    $x2 = Linear-X $energyMax $energyMin $energyMax $p.left $p.right
    $y2 = Linear-Y $energyMax $energyMin $energyMax $top $bottom
    $svg.Add("<line x1=`"$x1`" y1=`"$y1`" x2=`"$x2`" y2=`"$y2`" stroke=`"#6B7280`" stroke-width=`"2.3`" stroke-dasharray=`"10 8`"/>")
    foreach ($row in $Rows) {
        $x = Linear-X ([double]$row.apbs_deltaG_solv_kcal_mol) $energyMin $energyMax $p.left $p.right
        $y = Linear-Y ([double]$row.titan_deltaG_solv_kcal_mol) $energyMin $energyMax $top $bottom
        $svg.Add("<circle cx=`"$x`" cy=`"$y`" r=`"8`" fill=`"$blue`" fill-opacity=`"0.86`"/>")
    }
    $svg.Add("<text class=`"label`" x=`"$(($p.left + $p.right) / 2)`" y=`"752`" text-anchor=`"middle`">APBS &#916;G<tspan baseline-shift=`"sub`" font-size=`"20`">solv</tspan> (kcal mol&#8722;1)</text>")
    $svg.Add("<text class=`"label`" transform=`"translate(35,$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">TITAN &#916;G<tspan baseline-shift=`"sub`" font-size=`"20`">solv</tspan> (kcal mol&#8722;1)</text>")

    $panelSpecs = @(
        @{ index = 1; titan = 'titan_peak_rss_mb'; apbs = 'apbs_peak_rss_mb'; label = 'Peak RSS (MB)' },
        @{ index = 2; titan = 'titan_wall_seconds'; apbs = 'apbs_wall_seconds'; label = 'Total wall time (s)' }
    )
    $panelMin = [double](($Rows.panels | Measure-Object -Minimum).Minimum) * 0.85
    $panelMax = [double](($Rows.panels | Measure-Object -Maximum).Maximum) * 1.18
    foreach ($spec in $panelSpecs) {
        $p = $panels[$spec.index]
        $values = @($Rows | ForEach-Object { [double]$_.$($spec.titan); [double]$_.$($spec.apbs) })
        $yMin = [math]::Max(0.5, [double](($values | Measure-Object -Minimum).Minimum) * 0.72)
        $yMax = [double](($values | Measure-Object -Maximum).Maximum) * 1.42
        foreach ($tick in @(Nice-Log-Ticks $panelMin $panelMax)) {
            $x = Log-X $tick $panelMin $panelMax $p.left $p.right
            $svg.Add("<line class=`"grid`" x1=`"$x`" y1=`"$top`" x2=`"$x`" y2=`"$bottom`"/>")
            $svg.Add("<text class=`"tick`" x=`"$x`" y=`"690`" text-anchor=`"middle`">$(Format-Number $tick '0')</text>")
        }
        foreach ($tick in @(Nice-Log-Ticks $yMin $yMax)) {
            $y = Log-Y $tick $yMin $yMax $top $bottom
            $svg.Add("<line class=`"grid`" x1=`"$($p.left)`" y1=`"$y`" x2=`"$($p.right)`" y2=`"$y`"/>")
            $svg.Add("<text class=`"tick`" x=`"$($p.left - 15)`" y=`"$($y + 8)`" text-anchor=`"end`">$(Format-Number $tick '0.#')</text>")
        }
        foreach ($row in $Rows) {
            $x = Log-X ([double]$row.panels) $panelMin $panelMax $p.left $p.right
            $titanY = Log-Y ([double]$row.$($spec.titan)) $yMin $yMax $top $bottom
            $apbsY = Log-Y ([double]$row.$($spec.apbs)) $yMin $yMax $top $bottom
            $svg.Add("<line x1=`"$x`" y1=`"$titanY`" x2=`"$x`" y2=`"$apbsY`" stroke=`"#B9BEC7`" stroke-width=`"1.4`"/>")
            $svg.Add("<circle cx=`"$x`" cy=`"$titanY`" r=`"7`" fill=`"$blue`"/>")
            $svg.Add("<rect x=`"$($x - 7)`" y=`"$($apbsY - 7)`" width=`"14`" height=`"14`" fill=`"white`" stroke=`"$orange`" stroke-width=`"3`"/>")
        }
        $center = ($p.left + $p.right) / 2
        $svg.Add("<text class=`"label`" x=`"$center`" y=`"752`" text-anchor=`"middle`">TITAN boundary panels</text>")
        $svg.Add("<text class=`"label`" transform=`"translate($($p.left - 67),$(($top + $bottom) / 2)) rotate(-90)`" text-anchor=`"middle`">$(Xml-Escape $spec.label)</text>")
    }
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $figureSvgPath -Encoding UTF8

    $inkscape = Get-Command inkscape.exe -ErrorAction SilentlyContinue
    if ($null -eq $inkscape) { throw 'inkscape.exe is required to export the publication figure.' }
    $inkscapeCli = Join-Path (Split-Path -Parent $inkscape.Source) 'inkscape.com'
    if (-not (Test-Path -LiteralPath $inkscapeCli)) { $inkscapeCli = $inkscape.Source }
    & $inkscapeCli $figureSvgPath "--export-filename=$figurePngPath" '--export-dpi=600' | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Inkscape PNG export failed.' }
    & $inkscapeCli $figureSvgPath "--export-filename=$figurePdfPath" | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Inkscape PDF export failed.' }
}

$okPairs = @($pairRows | Where-Object { [string]$_.status -eq 'ok' })
if ($okPairs.Count -eq 27) { Write-Publication-Figure $okPairs }

$titanRows = @($merged | Where-Object { [string]$_.backend -eq 'titan' })
$diagonalFallbackRows = @($titanRows | Where-Object { [string]$_.preconditioner_effective -match 'diagonal|fallback' })
$auditWarnings = @($pairRows | Where-Object { (Parse-Double $_.titan_audit_rms 1.0) -gt 8e-3 })
$missingOutputs = @($figureSvgPath, $figurePngPath, $figurePdfPath | Where-Object { -not (Test-Path -LiteralPath $_) })
$decision = if ($merged.Count -ne 54 -or $pairRows.Count -ne 27) {
    'fagbi27_publication_merge_gap'
} elseif ($okPairs.Count -ne 27) {
    'fagbi27_publication_pair_gap'
} elseif ($diagonalFallbackRows.Count -gt 0) {
    'fagbi27_publication_solver_gap'
} elseif ($missingOutputs.Count -gt 0) {
    'fagbi27_publication_figure_gap'
} else {
    'fagbi27_publication_figure_ready'
}

$caption = 'Figure 5. Solvation-energy, memory, and timing comparison over 27 proteins. (A) TITAN and APBS electrostatic solvation energies; the dashed line denotes equality. (B) Peak resident-set size and (C) total wall time plotted against the TITAN boundary-panel count, used here as a common measure of molecular-system size. TITAN and APBS calculations used the same CHARMM PQR inputs and four CPU threads. Vertical grey segments connect results for the same protein. APBS is a contextual finite-difference comparison rather than a dense boundary-integral oracle.'
$caption | Set-Content -LiteralPath $captionPath -Encoding UTF8

[pscustomobject][ordered]@{
    decision = $decision
    merged_rows = $merged.Count
    paired_rows = $pairRows.Count
    paired_ok = $okPairs.Count
    audit_warning_rows = $auditWarnings.Count
    merged_matrix = $mergedPath
    supporting_table = $supportingPath
    figure_svg = $figureSvgPath
    figure_png = $figurePngPath
    figure_pdf = $figurePdfPath
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation

@(
    "decision=$decision"
    "merged_rows=$($merged.Count)"
    "paired_rows=$($pairRows.Count)"
    "paired_ok=$($okPairs.Count)"
    "audit_warning_rows=$($auditWarnings.Count)"
    "figure_png=$figurePngPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

@(
    'TITAN FAGBI-27 Publication Figure Report'
    "Decision: $decision"
    "Merged matrix: $mergedPath"
    "Paired systems: $($okPairs.Count)/27"
    "Rows above the 8e-3 sampled-audit reporting gate: $($auditWarnings.Count)"
    if ($auditWarnings.Count -gt 0) { "Audit warning systems: $((@($auditWarnings.pdb_id) -join ', '))" }
    'The original full matrix and recovery matrices were read only and were not modified.'
    'The visual uses only the method labels TITAN and APBS; solver recovery categories are intentionally not encoded in the figure.'
    "SVG: $figureSvgPath"
    "PNG (600 dpi): $figurePngPath"
    "PDF: $figurePdfPath"
    "Caption: $captionPath"
) | Set-Content -LiteralPath $reportPath -Encoding UTF8

Write-Host "Decision: $decision"
Write-Host "Merged matrix: $mergedPath"
Write-Host "Supporting table: $supportingPath"
Write-Host "Figure: $figurePngPath"

if ($decision -ne 'fagbi27_publication_figure_ready') { throw "Publication figure gate did not pass: $decision" }
