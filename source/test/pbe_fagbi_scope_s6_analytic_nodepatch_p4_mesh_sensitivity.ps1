param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_s6_analytic_nodepatch_p4_mesh'),
    [string]$PqrPath = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q\1vii_full_amber_ph7.pqr'),
    [double[]]$SesGridSpacings = @(1.50,1.25,1.00),
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$Force,
    [switch]$RerunFailed,
    [double[]]$RerunSpacings = @(),
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$PqrPath = (Resolve-Path -LiteralPath $PqrPath).Path
$TitanSha = (Get-FileHash -LiteralPath $TitanExe -Algorithm SHA256).Hash.ToLowerInvariant()
$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $gitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($gitDir)) {
        $gitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $gitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $gitDir
}
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    throw "Missing lib_ff under runtime repo root: $RuntimeRepoRoot"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$MatrixPath = Join-Path $WorkRoot 'fagbi_scope_s6_analytic_nodepatch_p4_mesh_matrix.csv'
$TablePath = Join-Path $WorkRoot 'fagbi_scope_s6_analytic_nodepatch_p4_mesh_table.csv'
$DecisionPath = Join-Path $WorkRoot 'fagbi_scope_s6_analytic_nodepatch_p4_mesh_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_s6_analytic_nodepatch_p4_mesh_summary.dat'
$Protocol = '1vii-full|analytic-ses|node-patch|guarded-p4|theta0.6|static-tail0.002|strict-geometry0.5|ilu64|tol1e-6|salt0.15|nt4'

function Convert-Double([object]$Value, [double]$Default = [double]::NaN) {
    $parsed = 0.0
    if ($null -ne $Value -and [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float, $Inv, [ref]$parsed)) { return $parsed }
    return $Default
}
function Convert-Int([object]$Value, [int64]$Default = 0) {
    $parsed = 0L
    if ($null -ne $Value -and [int64]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Integer, $Inv, [ref]$parsed)) { return $parsed }
    return $Default
}
function Get-Property([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    $property = $Object.PSObject.Properties[$Name]
    if ($null -eq $property) { return $Default }
    return $property.Value
}
function Safe-Spacing([double]$Value) {
    $Value.ToString('0.00', $Inv).Replace('.', 'p')
}
function Is-Finite([double]$Value) {
    -not [double]::IsNaN($Value) -and -not [double]::IsInfinity($Value)
}

$SesGridSpacings = @($SesGridSpacings | Sort-Object -Descending -Unique)
if ($SesGridSpacings.Count -lt 2) { throw 'At least two SES grid spacings are required.' }
foreach ($spacing in $SesGridSpacings) {
    if ($spacing -lt 0.25 -or $spacing -gt 4.0) { throw "Invalid SES grid spacing: $spacing" }
}

$rows = [Collections.Generic.List[object]]::new()
foreach ($spacing in $SesGridSpacings) {
    $runId = '1vii_analytic_nodepatch_p4_h{0}' -f (Safe-Spacing $spacing)
    $runDir = Join-Path $WorkRoot $runId
    $rowPath = Join-Path $runDir 'fagbi_scope_s6_row.csv'
    $old = $null
    if (Test-Path -LiteralPath $rowPath) {
        try { $old = @(Import-Csv -LiteralPath $rowPath)[0] } catch { $old = $null }
    }
    $hashMatches = $null -ne $old -and
        [string](Get-Property $old 'titan_sha256') -eq $TitanSha -and
        [string](Get-Property $old 'protocol') -eq $Protocol
    $targeted = @($RerunSpacings).Count -gt 0 -and $spacing -in @($RerunSpacings)
    $failed = $null -ne $old -and [string](Get-Property $old 'failure_class') -ne 'pass'
    if ($hashMatches -and -not $Force -and -not $targeted -and -not ($RerunFailed -and $failed)) {
        Write-Host "resume: $runId"
        $rows.Add($old) | Out-Null
        continue
    }

    if (Test-Path -LiteralPath $runDir) { Remove-Item -LiteralPath $runDir -Recurse -Force }
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $PqrPath -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -FmmLocal gmres-p4 `
        -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1.0e-6 -Threads $Threads `
        -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalGeometryGuard strict-geometry-mac -FmmLocalGeometryTheta 0.5 `
        -PrimarySurface analytic-ses -SesGridSpacing $spacing -DgSolv $true
    Add-Content -LiteralPath $arg -Value '-pbe-tabi-discretization=node-patch' -Encoding ASCII
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count -gt 0) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }

    Write-Host "run: $runId (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg `
        -CaseTimeoutSeconds $CaseTimeoutSeconds
    $status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) { 'failed' } else { 'ok' }
    $failure = if ($run.timed_out) { 'runtime_limited' } elseif ($run.exit_code -ne 0) { "exit_code_$($run.exit_code)" } else { 'pass' }
    $profile = $null; $energy = $null
    if ($status -eq 'ok') {
        try {
            $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
            $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[-1]
        } catch {
            $status = 'failed'; $failure = 'profile_schema_gap'
        }
    }
    $energyKj = if ($energy) { Convert-Double (Get-Property $energy 'E_rxn_kJmol') } else { [double]::NaN }
    $audit = if ($profile) { Convert-Double (Get-Property $profile 'tree_vs_direct_matvec_rel') } else { [double]::NaN }
    $sampledResidual = if ($profile) { Convert-Double (Get-Property $profile 'fmm_local_sampled_direct_true_residual') } else { [double]::NaN }
    $trueResidual = if ($profile) { Convert-Double (Get-Property $profile 'true_residual_rel') } else { [double]::NaN }
    $surface = if ($profile) { [string](Get-Property $profile 'surface_generator') } else { '' }
    $discretization = if ($profile) { [string](Get-Property $profile 'discretization') } else { '' }
    $protocolOk = $surface.Replace('-','_') -eq 'titan_analytic_patch_ray_ses' -and
        $discretization.Replace('-','_') -eq 'node_patch'
    $accuracyOk = (Is-Finite $energyKj) -and (Is-Finite $audit) -and $audit -le 8.0e-3 -and
        (Is-Finite $sampledResidual) -and $sampledResidual -le 8.0e-3 -and
        (Is-Finite $trueResidual) -and $trueResidual -le 1.0e-6
    if ($failure -eq 'pass' -and -not $protocolOk) { $status = 'failed'; $failure = 'surface_or_discretization_mismatch' }
    elseif ($failure -eq 'pass' -and -not $accuracyOk) { $status = 'failed'; $failure = 'accuracy_gap' }

    $row = [pscustomobject][ordered]@{
        run_id = $runId
        protocol = $Protocol
        titan_sha256 = $TitanSha
        fixture = '1vii'
        pdb_id = '1VII'
        protein_name = 'Villin headpiece subdomain'
        ses_grid_spacing_A = $spacing.ToString('0.00', $Inv)
        surface_generator = $surface
        discretization = $discretization
        status = $status
        failure_class = $failure
        exit_code = $run.exit_code
        timed_out = [string]$run.timed_out
        atom_count = if ($profile) { Convert-Int (Get-Property $profile 'atom_count') } else { 0 }
        panels = if ($profile) { Convert-Int (Get-Property $profile 'panels') } else { 0 }
        energy_kj_mol = $energyKj
        energy_kcal_mol = $energyKj / 4.184
        sampled_audit_rms = $audit
        sampled_direct_true_residual = $sampledResidual
        true_residual = $trueResidual
        iterations = if ($profile) { Convert-Int (Get-Property $profile 'iterations') } else { 0 }
        wall_seconds = $run.wall_seconds
        peak_rss_mb = $run.peak_working_set_mb
        plan_build_count = if ($profile) { Convert-Int (Get-Property $profile 'fmm_local_plan_build_count') } else { 0 }
        cache_build_count = if ($profile) { Convert-Int (Get-Property $profile 'fmm_local_cache_build_count') } else { 0 }
        run_dir = $runDir
    }
    $row | Export-Csv -LiteralPath $rowPath -NoTypeInformation
    $rows.Add($row) | Out-Null
}

@($rows | Sort-Object { Convert-Double $_.ses_grid_spacing_A } -Descending) |
    Export-Csv -LiteralPath $MatrixPath -NoTypeInformation
@($rows | Sort-Object { Convert-Double $_.ses_grid_spacing_A } -Descending) |
    Select-Object pdb_id,ses_grid_spacing_A,panels,energy_kcal_mol,sampled_audit_rms,
        sampled_direct_true_residual,true_residual,iterations,wall_seconds,peak_rss_mb,status,failure_class |
    Export-Csv -LiteralPath $TablePath -NoTypeInformation

$passes = @($rows | Where-Object { $_.failure_class -eq 'pass' })
$decision = if ($rows.Count -ne $SesGridSpacings.Count) { 's6_mesh_incomplete' }
elseif (@($rows | Where-Object { $_.failure_class -eq 'runtime_limited' }).Count -gt 0) { 's6_mesh_runtime_limited' }
elseif ($passes.Count -eq $SesGridSpacings.Count) { 's6_analytic_nodepatch_p4_mesh_ready' }
else { 's6_mesh_protocol_or_accuracy_gap' }
[pscustomobject]@{ decision=$decision; rows=$rows.Count; pass_rows=$passes.Count; titan_sha256=$TitanSha } |
    Export-Csv -LiteralPath $DecisionPath -NoTypeInformation
@(
    "decision=$decision"
    "matrix=$MatrixPath"
    "table=$TablePath"
    "surface=analytic-ses"
    "discretization=node-patch"
    "taylor_order=4"
    "spacings_A=$($SesGridSpacings -join ',')"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

Write-Host "Table S6 decision: $decision"
Write-Host "Matrix: $MatrixPath"
Write-Host "Table: $TablePath"
Write-Host 'Resume: rerun unchanged; use -RerunFailed or -RerunSpacings 1.00 for targeted reruns.'
if (-not $AllowFailedDecision -and $decision -ne 's6_analytic_nodepatch_p4_mesh_ready') { exit 1 }
