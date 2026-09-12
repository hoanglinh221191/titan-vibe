param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$BaselineRoot = 'D:\data_MD\vibe\pbctopo_long_5case_101frames_20260702_113315',
    [string]$OutputRoot = 'D:\data_MD\vibe\pbctopo_review_12point_validation_20260716'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$BaselineRoot = (Resolve-Path -LiteralPath $BaselineRoot).Path
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null

function Get-CoordinateHash {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][int]$MaxModels
    )
    $sha = [System.Security.Cryptography.SHA256]::Create()
    $reader = [System.IO.StreamReader]::new($Path)
    $encoding = [System.Text.Encoding]::ASCII
    $model = 0
    try {
        while (($line = $reader.ReadLine()) -ne $null) {
            if ($line.StartsWith('MODEL')) {
                $model++
                if ($model -gt $MaxModels) { break }
                continue
            }
            if (($model -eq 0 -or $model -le $MaxModels) -and
                ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
                $bytes = $encoding.GetBytes($line + "`n")
                [void]$sha.TransformBlock($bytes, 0, $bytes.Length, $bytes, 0)
            }
        }
        [void]$sha.TransformFinalBlock([byte[]]::new(0), 0, 0)
        return (($sha.Hash | ForEach-Object { $_.ToString('X2') }) -join '')
    }
    finally {
        $reader.Dispose()
        $sha.Dispose()
    }
}

function Get-LogInteger {
    param([string]$Text, [string]$Name)
    $matches = [regex]::Matches($Text, [regex]::Escape($Name) + '=(\d+)')
    if ($matches.Count -eq 0) { return 0L }
    return [long]$matches[$matches.Count - 1].Groups[1].Value
}

$cases = @(
    [pscustomobject]@{
        Name='7wbl_protein_glycan'; Base='7wbl_ba275_protein'; End=8000
        Group='28-28'; Extra=@()
    },
    [pscustomobject]@{ Name='amylin_llps_local'; Base='amylin_128_martini'; End=40000; Extra=@() },
    [pscustomobject]@{ Name='dimer_ab42'; Base='dimerab42'; End=8000; Extra=@() },
    [pscustomobject]@{ Name='npro_6mer'; Base='npro_6mer_wt'; End=4000; Extra=@() },
    [pscustomobject]@{ Name='npro_tetramer'; Base='npro_tetramer_wt'; End=8000; Extra=@() },
    [pscustomobject]@{
        Name='npro_tetramer_viterbi'; Base='npro_tetramer_wt'; End=8000
        Extra=@('-pbctopo-temporal=viterbi','-pbctopo-temporal-k=4','-pbctopo-temporal-penalty=2.0')
    },
    [pscustomobject]@{ Name='npro_6mer_anchor_full'; Base='npro_6mer_wt'; End=0; Extra=@('-pbctopo-anchor-metrics=full') },
    [pscustomobject]@{ Name='npro_6mer_anchor_mst_only'; Base='npro_6mer_wt'; End=0; Extra=@('-pbctopo-anchor-metrics=mst-only') },
    [pscustomobject]@{ Name='npro_6mer_anchor_no_order'; Base='npro_6mer_wt'; End=0; Extra=@('-pbctopo-anchor-metrics=no-order') },
    [pscustomobject]@{ Name='npro_6mer_anchor_no_contact'; Base='npro_6mer_wt'; End=0; Extra=@('-pbctopo-anchor-metrics=no-contact') }
)

$summary = [System.Collections.Generic.List[object]]::new()
foreach ($case in $cases) {
    $baselineDir = Join-Path $BaselineRoot $case.Base
    $baselineArg = Join-Path $baselineDir 'arg.in'
    $baseText = Get-Content -Raw -LiteralPath $baselineArg
    $baselineOutputMatch = [regex]::Match($baseText, '-opdb=([^\s]+)')
    if (-not $baselineOutputMatch.Success) {
        throw "Missing baseline opdb in $baselineArg"
    }
    $baselineOutput = Join-Path $baselineDir $baselineOutputMatch.Groups[1].Value

    $runDir = Join-Path $OutputRoot $case.Name
    if (Test-Path -LiteralPath $runDir) {
        Remove-Item -LiteralPath $runDir -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $outputName = $case.Name + '.pdb'
    $argText = [regex]::Replace($baseText, '(?<!\S)-e\s+\d+', "-e $($case.End)")
    $argText = [regex]::Replace($argText, '-opdb=([^\s]+)', "-opdb=$outputName")
    $groupProperty = $case.PSObject.Properties['Group']
    if ($null -ne $groupProperty -and -not [string]::IsNullOrWhiteSpace($groupProperty.Value)) {
        $argText = [regex]::Replace($argText, '-g=\d+-\d+', "-g=$($groupProperty.Value)")
    }
    if ($case.Extra.Count -gt 0) {
        $argText += "`r`n" + ($case.Extra -join "`r`n") + "`r`n"
    }
    $argPath = Join-Path $runDir 'arg.in'
    [System.IO.File]::WriteAllText($argPath, $argText, [System.Text.Encoding]::ASCII)

    Write-Host "RUN $($case.Name)"
    $console = Join-Path $runDir 'titan_console.txt'
    $consoleError = Join-Path $runDir 'titan_stderr.txt'
    $process = Start-Process -FilePath $TitanExe `
        -ArgumentList @('-a', ('"' + $argPath + '"')) `
        -WorkingDirectory $runDir -WindowStyle Hidden -Wait -PassThru `
        -RedirectStandardOutput $console -RedirectStandardError $consoleError
    if ($process.ExitCode -ne 0) {
        throw "TITAN failed for $($case.Name): $($process.ExitCode)"
    }

    $output = Join-Path $runDir $outputName
    $report = Join-Path $runDir 'pbctopo_frame_report.csv'
    $debug = Join-Path $runDir 'pbctopo_debug.log'
    foreach ($required in @($output, $report, $debug)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing output for $($case.Name): $required"
        }
    }

    $rows = @(Import-Csv -LiteralPath $report)
    if ($rows.Count -eq 0) { throw "Empty report for $($case.Name)" }
    $fallback = @($rows | Where-Object { $_.status -ne 'certified' }).Count
    if ($fallback -ne 0) {
        throw "$($case.Name) has $fallback non-certified frames"
    }
    $hardResiduals = ($rows | Measure-Object -Property hard_edge_residuals -Sum).Sum
    $cycleMismatches = ($rows | Measure-Object -Property hard_cycle_mismatches -Sum).Sum
    if ($case.Base -ne 'amylin_128_martini' -and
        ([double]$hardResiduals -ne 0 -or [double]$cycleMismatches -ne 0)) {
        throw "$($case.Name) has hard residuals or cycle mismatches"
    }

    $modelCount = $rows.Count
    $currentHash = Get-CoordinateHash -Path $output -MaxModels $modelCount
    $baselineHash = Get-CoordinateHash -Path $baselineOutput -MaxModels $modelCount
    $baselineComparable = $null -eq $groupProperty
    $baselineCoordinateMatch = if ($baselineComparable) {
        $currentHash -eq $baselineHash
    } else {
        $null
    }
    $debugText = Get-Content -Raw -LiteralPath $debug
    $modeCounts = ($rows.mode | Group-Object | ForEach-Object {
        "$($_.Name):$($_.Count)"
    }) -join ';'
    $componentCounts = (($rows.component_count | Sort-Object -Unique) -join ';')
    $signatures = @($rows.assignment_signature | Where-Object {
        $_ -and $_ -ne '0'
    } | Sort-Object -Unique).Count

    $summary.Add([pscustomobject]@{
        case = $case.Name
        frames = $rows.Count
        certified = $rows.Count - $fallback
        fallback = $fallback
        modes = $modeCounts
        component_counts = $componentCounts
        hard_edge_residuals = $hardResiduals
        hard_cycle_mismatches = $cycleMismatches
        assignment_signatures = $signatures
        baseline_comparable = $baselineComparable
        baseline_coordinate_match = $baselineCoordinateMatch
        coordinate_sha256 = $currentHash
        baseline_coordinate_sha256 = $baselineHash
        topology_source = if ($debugText -match 'pbctopo_topology_source=([^,\r\n]+)') { $Matches[1] } else { 'na' }
        pbci_fallback_frames = Get-LogInteger $debugText 'pbci_fallback_frames'
        component_graph_failed = Get-LogInteger $debugText 'component_graph_failed'
        temporal_replay_applied = Get-LogInteger $debugText 'temporal_replay_applied'
        sweep_workspace_reuse_frames = Get-LogInteger $debugText 'pbctopo_virtual_sweep_workspace_reuse_frames'
        sweep_workspace_peak_bytes = Get-LogInteger $debugText 'pbctopo_virtual_sweep_workspace_peak_bytes'
        point_buffer_peak_bytes = Get-LogInteger $debugText 'pbctopo_frame_point_buffer_peak_bytes_estimate'
    })
}

$summaryPath = Join-Path $OutputRoot 'validation_summary.csv'
$summary | Export-Csv -LiteralPath $summaryPath -NoTypeInformation -Encoding UTF8
$summary | Format-Table case,frames,certified,fallback,modes,component_counts,baseline_coordinate_match -AutoSize
Write-Host "PASS pbctopo_review_12point_validation: $summaryPath"
