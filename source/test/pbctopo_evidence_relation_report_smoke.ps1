param(
    [string]$Compiler = 'g++',
    [string]$Python = 'python',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_evidence_relation_report_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$source = Join-Path $repo 'test\titan_pbctopo_evidence_relation_report_smoke.cpp'
$validator = Join-Path $repo 'tools\validate_pbctopo_evidence_relations.py'
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$exe = Join-Path $WorkRoot 'smoke.exe'

Push-Location $repo
try {
    & $Compiler -std=c++23 -O2 -Wall -Wextra -pedantic -Isrc $source `
        src/pbctopo/evidence_relation_report.cpp `
        src/pbctopo/evidence_connectivity_margin.cpp `
        -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile evidence relation-report smoke'
    }

    $runA = Join-Path $WorkRoot 'run_a'
    $runB = Join-Path $WorkRoot 'run_b'
    & $exe $runA
    if ($LASTEXITCODE -ne 0) {
        throw 'Evidence relation-report fixture run A failed'
    }
    & $exe $runB
    if ($LASTEXITCODE -ne 0) {
        throw 'Evidence relation-report fixture run B failed'
    }

    $artifactNames = @(
        'pbctopo_candidate_evidence_summary.csv',
        'pbctopo_candidate_evidence_relations.csv',
        'pbctopo_frame_report.csv'
    )
    foreach ($artifact in $artifactNames) {
        $hashA = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $runA $artifact)).Hash
        $hashB = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $runB $artifact)).Hash
        if ($hashA -ne $hashB) {
            throw "Nondeterministic evidence relation artifact: $artifact"
        }
    }

    $summary = Join-Path $runA 'pbctopo_candidate_evidence_summary.csv'
    $relations = Join-Path $runA 'pbctopo_candidate_evidence_relations.csv'
    $frames = Join-Path $runA 'pbctopo_frame_report.csv'
    $validation = Join-Path $runA 'replay_validation.json'
    & $Python $validator --summary $summary --relations $relations `
        --frame-report $frames --output $validation
    if ($LASTEXITCODE -ne 0) {
        throw 'Independent evidence relation replay failed'
    }

    $result = Get-Content -LiteralPath $validation -Raw | ConvertFrom-Json
    if ($result.status -ne 'PASS' -or $result.candidate_groups -ne 3 -or `
            $result.relation_rows -ne 4 -or $result.frame_report_matches -ne 2) {
        throw 'Independent replay returned unexpected fixture counts'
    }

    $tampered = Join-Path $WorkRoot 'tampered_relations.csv'
    $relationText = [System.IO.File]::ReadAllText($relations)
    $tamperedText = [System.Text.RegularExpressions.Regex]::Replace(
        $relationText, ',1,decision_score,', ',0,decision_score,', 1)
    if ($tamperedText -eq $relationText) {
        throw 'Could not construct the tampered relation fixture'
    }
    [System.IO.File]::WriteAllText(
        $tampered, $tamperedText, [System.Text.UTF8Encoding]::new($false))
    $tamperOutput = & $Python $validator --summary $summary `
        --relations $tampered --frame-report $frames 2>&1
    if ($LASTEXITCODE -eq 0) {
        throw 'Independent replay accepted a tampered relation row'
    }
    if (($tamperOutput -join "`n") -notmatch 'row payload hash mismatch') {
        throw "Tamper rejection occurred for the wrong reason: $tamperOutput"
    }

    # The rejected validator call above is intentional. Clear its native exit
    # status so callers can safely compose this script as a package gate.
    $global:LASTEXITCODE = 0
    Write-Output 'PASS pbctopo_evidence_relation_report_smoke determinism=3/3 replay=PASS tamper=REJECT'
}
finally {
    Pop-Location
}
