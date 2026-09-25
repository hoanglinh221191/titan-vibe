[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$ApbsExe = '',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_fagbi27_recovery_v2',
    [string]$RecoveryFixtureRoot = 'test\pbe_inputs_fagbi27_recovery',
    [int]$Threads = 4,
    [int]$TitanTimeoutSeconds = 1800,
    [int]$ApbsTimeoutSeconds = 900,
    [switch]$RefreshFixtures,
    [switch]$RerunFailed,
    [switch]$Force,
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$mainScript = Join-Path $PSScriptRoot 'pbe_tabi_fagbi27_analytic_ses_apbs.ps1'
$originalWorkRoot = Join-Path $RepoRoot 'test\tmp\pbe_tabi_fagbi27_analytic_ses_apbs'
$originalFixtureRoot = Join-Path $RepoRoot 'test\pbe_inputs_fagbi27'
$originalMatrix = Join-Path $originalWorkRoot 'fagbi27_analytic_ses_apbs_matrix.csv'
$originalManifest = Join-Path $originalFixtureRoot 'fagbi27_fixture_manifest.csv'

function Full-Path([string]$Path) {
    if ([IO.Path]::IsPathRooted($Path)) { return [IO.Path]::GetFullPath($Path) }
    return [IO.Path]::GetFullPath((Join-Path $RepoRoot $Path))
}
function Hash-Or-Empty([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return '' }
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash
}
function As-Double($Value, [double]$Fallback = [double]::NaN) {
    $parsed = 0.0
    if ([double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,
            [ref]$parsed)) { return $parsed }
    return $Fallback
}

$resolvedWorkRoot = Full-Path $WorkRoot
$resolvedRecoveryFixtureRoot = Full-Path $RecoveryFixtureRoot
if ($resolvedWorkRoot -eq [IO.Path]::GetFullPath($originalWorkRoot) -or
    $resolvedRecoveryFixtureRoot -eq [IO.Path]::GetFullPath($originalFixtureRoot)) {
    throw 'Recovery roots must differ from the original FAGBI-27 roots.'
}
New-Item -ItemType Directory -Force -Path $resolvedWorkRoot | Out-Null

$matrixHashBefore = Hash-Or-Empty $originalMatrix
$manifestHashBefore = Hash-Or-Empty $originalManifest
$common = @{
    TitanExe = $TitanExe
    ApbsExe = $ApbsExe
    Threads = $Threads
    TitanTimeoutSeconds = $TitanTimeoutSeconds
    ApbsTimeoutSeconds = $ApbsTimeoutSeconds
    PreconditionerIluRegularization = 'adaptive'
    PreconditionerIluRegularizationMax = 0.01
    AllowFailedDecision = $true
}
if ($RerunFailed) { $common.RerunFailed = $true }
if ($Force) { $common.Force = $true }

# Recover the malformed alternate-residue fixture in an isolated fixture cache.
$cleanInvoke = @{} + $common
$cleanInvoke.WorkRoot = Join-Path $resolvedWorkRoot '1cbn_clean_fixture'
$cleanInvoke.FixtureRoot = $resolvedRecoveryFixtureRoot
$cleanInvoke.Fixtures = @('1CBN')
$cleanInvoke.Backends = @('titan','apbs')
$cleanInvoke.ReportOutput = 'TITAN_TABI_FAGBI27_RECOVERY_1CBN_REPORT.txt'
if ($RefreshFixtures) { $cleanInvoke.RefreshFixtures = $true }
& $mainScript @cleanInvoke

# Isolate regularization using exactly the PQR files from the original matrix.
$regularizedInvoke = @{} + $common
$regularizedInvoke.WorkRoot = Join-Path $resolvedWorkRoot 'regularized_original_pqr'
$regularizedInvoke.FixtureRoot = $originalFixtureRoot
$regularizedInvoke.Fixtures = @('2ERL','1HPT')
$regularizedInvoke.Backends = @('titan')
$regularizedInvoke.ReportOutput = 'TITAN_TABI_FAGBI27_RECOVERY_REGULARIZED_REPORT.txt'
& $mainScript @regularizedInvoke

$matrixHashAfter = Hash-Or-Empty $originalMatrix
$manifestHashAfter = Hash-Or-Empty $originalManifest
if ($matrixHashAfter -ne $matrixHashBefore -or
    $manifestHashAfter -ne $manifestHashBefore) {
    throw 'Recovery isolation failure: an original FAGBI-27 artifact changed.'
}

$cleanMatrix = Join-Path $cleanInvoke.WorkRoot 'fagbi27_analytic_ses_apbs_matrix.csv'
$regularizedMatrix = Join-Path $regularizedInvoke.WorkRoot 'fagbi27_analytic_ses_apbs_matrix.csv'
foreach ($path in @($originalMatrix,$cleanMatrix,$regularizedMatrix)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing recovery input: $path" }
}
$originalRows = @(Import-Csv -LiteralPath $originalMatrix)
$cleanRows = @(Import-Csv -LiteralPath $cleanMatrix)
$regularizedRows = @(Import-Csv -LiteralPath $regularizedMatrix)
$combined = [Collections.Generic.List[object]]::new()
foreach ($row in @($originalRows | Where-Object {
        $_.backend -eq 'titan' -and $_.pdb_id -in @('2ERL','1HPT')
    })) {
    $combined.Add(($row | Select-Object @{Name='recovery_stage';Expression={'original_unregularized'}},*))
}
foreach ($row in $regularizedRows) {
    $combined.Add(($row | Select-Object @{Name='recovery_stage';Expression={'adaptive_regularized_original_pqr'}},*))
}
foreach ($row in $cleanRows) {
    $combined.Add(($row | Select-Object @{Name='recovery_stage';Expression={'clean_altloc_fixture'}},*))
}
$combinedPath = Join-Path $resolvedWorkRoot 'fagbi27_recovery_matrix.csv'
$combined | Sort-Object {[int]$_.paper_index},recovery_stage,backend |
    Export-Csv -LiteralPath $combinedPath -NoTypeInformation

$cleanTitan = @($cleanRows | Where-Object {
    $_.pdb_id -eq '1CBN' -and $_.backend -eq 'titan'
})[-1]
$cleanApbs = @($cleanRows | Where-Object {
    $_.pdb_id -eq '1CBN' -and $_.backend -eq 'apbs'
})[-1]
$regularizedTitan = @($regularizedRows | Where-Object {
    $_.backend -eq 'titan' -and $_.pdb_id -in @('2ERL','1HPT')
})
$regularizedPass = @($regularizedTitan | Where-Object {
    $_.status -eq 'ok' -and
    (As-Double $_.true_residual 1.0) -le 1.0e-4 -and
    (As-Double $_.audit_rms 1.0) -le 8.0e-3 -and
    $_.preconditioner_effective -notmatch 'diagonal|fallback'
}).Count -eq 2
$cleanPass = $null -ne $cleanTitan -and $null -ne $cleanApbs -and
    $cleanTitan.status -eq 'ok' -and $cleanApbs.status -eq 'ok' -and
    (As-Double $cleanTitan.true_residual 1.0) -le 1.0e-4 -and
    (As-Double $cleanTitan.audit_rms 1.0) -le 8.0e-3
$decision = if (-not $cleanPass) { 'fagbi27_recovery_1cbn_gap' }
    elseif (-not $regularizedPass) { 'fagbi27_recovery_regularization_gap' }
    else { 'fagbi27_recovery_pass' }

$decisionPath = Join-Path $resolvedWorkRoot 'fagbi27_recovery_decision.csv'
[pscustomobject]@{
    decision = $decision
    clean_1cbn_pass = $cleanPass
    regularized_original_pqr_pass = $regularizedPass
    original_artifacts_unchanged = $true
    matrix = $combinedPath
} | Export-Csv -LiteralPath $decisionPath -NoTypeInformation
$summaryPath = Join-Path $resolvedWorkRoot 'fagbi27_recovery_summary.dat'
@(
    "decision=$decision",
    "clean_1cbn_pass=$cleanPass",
    "regularized_original_pqr_pass=$regularizedPass",
    'original_artifacts_unchanged=true',
    "matrix=$combinedPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
$reportPath = Join-Path $RepoRoot 'TITAN_TABI_FAGBI27_RECOVERY_REPORT.txt'
@(
    'TITAN FAGBI-27 Recovery Report',
    "Decision: $decision",
    '1CBN: first model, input hydrogens removed, one occupancy-ranked residue conformer retained, CHARMM PDB2PQR regenerated.',
    '2ERL and 1HPT: adaptive regularized ILU(0) evaluated with the exact original PQR fixtures.',
    'All recovery work used isolated work directories; original matrix and fixture manifest hashes were unchanged.',
    "Matrix: $combinedPath"
) | Set-Content -LiteralPath $reportPath -Encoding ASCII

Write-Host "Recovery decision: $decision"
Write-Host "Recovery matrix: $combinedPath"
Write-Host 'Original FAGBI-27 matrix and fixture manifest were hash-verified unchanged.'
Write-Host 'Resume: rerun the same command; add -RerunFailed for failed/timeout rows.'
if (-not $AllowFailedDecision -and $decision -ne 'fagbi27_recovery_pass') {
    throw "FAGBI-27 recovery did not pass: $decision"
}
