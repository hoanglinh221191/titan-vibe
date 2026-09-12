param(
    [string]$TitanExe =
        (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot =
        (Join-Path $PSScriptRoot 'tmp\pbctopo_temporal_report_integration_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$titanPath = (Resolve-Path -LiteralPath $TitanExe).Path
$confidenceSmoke =
    Join-Path $PSScriptRoot 'pbctopo_confidence_report_smoke.ps1'

& $confidenceSmoke -TitanExe $titanPath -WorkRoot $WorkRoot | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "PBCTOPO confidence fixture setup failed: $LASTEXITCODE"
}

$argPath = Join-Path $WorkRoot 'arg.in'
$requiredContactsPath = Join-Path $WorkRoot 'temporal_required.contacts'
@'
1 5 30.0 1.0
3 5 45.0 1.0
'@ | Set-Content -LiteralPath $requiredContactsPath -Encoding ASCII
$argText = [System.IO.File]::ReadAllText($argPath)
if (-not $argText.EndsWith("`n")) {
    $argText += "`r`n"
}
$argText += @"
-pbctopo-required-contacts=$requiredContactsPath
-pbctopo-temporal=viterbi
-pbctopo-temporal-k=8
-pbctopo-temporal-penalty=1.0
"@
[System.IO.File]::WriteAllText(
    $argPath, $argText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a arg.in | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN Viterbi integration run failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

$frameRows = @(
    Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbctopo_frame_report.csv')
)
$assignmentRows = @(
    Import-Csv -LiteralPath (
        Join-Path $WorkRoot 'pbctopo_assignment_candidates.csv')
)
$temporalPath = Join-Path $WorkRoot 'pbctopo_temporal_path.csv'
if (-not (Test-Path -LiteralPath $temporalPath)) {
    throw 'Missing pbctopo_temporal_path.csv'
}
$temporalRows = @(Import-Csv -LiteralPath $temporalPath)
if ($temporalRows.Count -ne $frameRows.Count -or $temporalRows.Count -eq 0) {
    throw 'Temporal path does not cover every certified fixture frame'
}

foreach ($temporalRow in $temporalRows) {
    if ($temporalRow.observation_scope -ne 'experimental_temporal_path' -or
        $temporalRow.temporal_policy_status -ne
            'experimental_viterbi_diagnostic' -or
        [int]$temporalRow.selected_rank -lt 1) {
        throw 'Temporal row does not expose the experimental 1-based contract'
    }
    $frameRow = $frameRows |
        Where-Object frame -eq $temporalRow.frame |
        Select-Object -First 1
    if ($null -eq $frameRow -or $frameRow.temporal_selected -ne '1' -or
        $frameRow.temporal_assignment_signature -ne
            $temporalRow.assignment_signature -or
        $frameRow.time_ps -ne $temporalRow.time_ps) {
        throw 'Temporal row does not match its frame report'
    }
    $frameCandidates = @(
        $assignmentRows |
            Where-Object frame -eq $temporalRow.frame
    )
    $candidateSignatures = @(
        $frameCandidates | Select-Object -ExpandProperty assignment_signature)
    if ($candidateSignatures -notcontains $temporalRow.assignment_signature) {
        throw 'Temporal state is absent from the exported candidate set'
    }
    $selectedCandidate = $frameCandidates |
        Where-Object assignment_signature -eq $temporalRow.assignment_signature |
        Select-Object -First 1
    if ($null -eq $selectedCandidate -or
        $selectedCandidate.time_ps -ne $temporalRow.time_ps) {
        throw 'Temporal time does not exactly match the selected candidate'
    }
}

$nonRoundFrame = $frameRows | Where-Object frame -eq '1' | Select-Object -First 1
$nonRoundTemporal = $temporalRows |
    Where-Object frame -eq '1' | Select-Object -First 1
if ($null -eq $nonRoundFrame -or $null -eq $nonRoundTemporal -or
    $nonRoundFrame.time_ps -ne '1.1234567890123457' -or
    $nonRoundTemporal.time_ps -ne '1.1234567890123457') {
    throw 'Non-round temporal time_ps did not round-trip across artifacts'
}

$canonicalNegativeZero = '(?m)(?:^|,)-0(?:\.0*)?(?=,|$)'
foreach ($artifact in @(
        'pbctopo_frame_report.csv',
        'pbctopo_assignment_candidates.csv',
        'pbctopo_temporal_path.csv')) {
    $text = [System.IO.File]::ReadAllText((Join-Path $WorkRoot $artifact))
    if ($text -match $canonicalNegativeZero) {
        throw "Artifact violates positive-zero serialization: $artifact"
    }
}

$pythonPath = (Get-Command python -ErrorAction Stop).Source
$validator = Join-Path $repo 'tools\validate_vibe_observation.py'
& $pythonPath $validator $WorkRoot | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "VIBE temporal package validation failed: $LASTEXITCODE"
}

$tamperRoot = "$WorkRoot-validator-tamper"
if (Test-Path -LiteralPath $tamperRoot) {
    Remove-Item -LiteralPath $tamperRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $tamperRoot | Out-Null
foreach ($artifact in @(
        'pbctopo_component_layout.csv',
        'pbctopo_hard_edges.csv',
        'pbctopo_frame_report.csv',
        'pbctopo_assignment_candidates.csv',
        'pbctopo_temporal_path.csv')) {
    Copy-Item -LiteralPath (Join-Path $WorkRoot $artifact) `
        -Destination (Join-Path $tamperRoot $artifact)
}
$tamperedTemporalPath = Join-Path $tamperRoot 'pbctopo_temporal_path.csv'
$tamperedTemporal = @(Import-Csv -LiteralPath $tamperedTemporalPath)
$tamperedTemporal[0].assignment_signature = '1'
$tamperedTemporal | Export-Csv -LiteralPath $tamperedTemporalPath `
    -NoTypeInformation -Encoding ASCII
$savedErrorActionPreference = $ErrorActionPreference
try {
    $ErrorActionPreference = 'Continue'
    & $pythonPath $validator $tamperRoot 2>&1 | Out-Null
    $validatorExitCode = $LASTEXITCODE
}
finally {
    $ErrorActionPreference = $savedErrorActionPreference
}
if ($validatorExitCode -eq 0) {
    throw 'VIBE validator accepted a temporal signature outside candidate set'
}
$global:LASTEXITCODE = 0

Write-Host 'PASS pbctopo_temporal_report_integration_smoke'
