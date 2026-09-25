[CmdletBinding()]
param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$Compiler = 'clang++',
    [string]$WorkRoot = 'test\tmp\fagbi_calderon_c5d1r_smoke',
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$smokeExe = Join-Path $WorkRoot 'c5d1r_regularized_ilu_smoke.exe'
$sources = @(
    'test\titan_fmm_m2l_scope_u_preconditioner_smoke.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_fmm_m2l_preconditioner.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp',
    'src\tabi\tabi_preconditioner.cpp',
    'src\titan_executor.cpp'
)
$sourcePaths = @($sources | ForEach-Object { Join-Path $RepoRoot $_ })
$needsBuild = $Force -or -not (Test-Path -LiteralPath $smokeExe)
if (-not $needsBuild) {
    $exeTime = (Get-Item -LiteralPath $smokeExe).LastWriteTimeUtc
    $needsBuild = @($sourcePaths | Where-Object {
        (Get-Item -LiteralPath $_).LastWriteTimeUtc -gt $exeTime
    }).Count -gt 0
}
if ($needsBuild) {
    if ($null -eq (Get-Command $Compiler -ErrorAction SilentlyContinue)) {
        throw "C5D.1R compiler was not found: $Compiler"
    }
    $arguments = @('-std=c++23','-O2','-mavx2','-mfma',
        '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',"-I$RepoRoot\src") +
        $sourcePaths + @('-o',$smokeExe)
    Write-Host "build: C5D.1R regularized ILU smoke with $Compiler"
    & $Compiler @arguments
    if ($LASTEXITCODE -ne 0) { throw 'C5D.1R C++ smoke build failed.' }
}
& $smokeExe
if ($LASTEXITCODE -ne 0) { throw 'C5D.1R C++ smoke failed.' }

$pqr = Join-Path $WorkRoot 'parser_fixture.pqr'
@(
    'ATOM      1  C   GLY A   1       0.000   0.000   0.000  1.0000  1.7000',
    'ATOM      2  O   GLY A   1       4.000   0.000   0.000 -1.0000  1.5000'
) | Set-Content -LiteralPath $pqr -Encoding ASCII

function New-C5D1rArgs([string]$Name) {
    $dir = Join-Path $WorkRoot $Name
    New-FmmCliCleanDir -Path $dir -WorkRoot $WorkRoot
    return Write-FmmCliArgFile -RunDir $dir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter 20 -Restart 16 -Tolerance 1e-4 -Threads 1 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditSampleCount 0 `
        -AuditMode off -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -PreconditionerIluRegularization adaptive `
        -PreconditionerIluRegularizationMax 0.01 `
        -NearIntegration analytic_singular -TaylorOrder 4 -DgSolv $true
}

function Invoke-C5D1rParserCase([string]$Name, [string[]]$Lines,
        [bool]$ExpectSuccess, [string]$ExpectedMessage = '') {
    $dir = Join-Path $WorkRoot $Name
    New-FmmCliCleanDir -Path $dir -WorkRoot $WorkRoot
    $arg = Join-Path $dir 'arg.in'
    $Lines | Set-Content -LiteralPath $arg -Encoding ASCII
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $dir `
        -ArgPath $arg -CaseTimeoutSeconds 180
    $log = (Get-Content -Raw (Join-Path $dir 'stdout.log')) +
        (Get-Content -Raw (Join-Path $dir 'stderr.log'))
    if ($ExpectSuccess -and $run.exit_code -ne 0) {
        throw "$Name unexpectedly failed: $log"
    }
    if (-not $ExpectSuccess -and $run.exit_code -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if (-not $ExpectSuccess -and $log -notmatch $ExpectedMessage) {
        throw "$Name error was unclear: $log"
    }
}

$seedArg = New-C5D1rArgs 'seed'
$accepted = @(Get-Content -LiteralPath $seedArg | Where-Object {
    $_ -notmatch '^-pbe-tabi-precond-ilu-regularization(-max)?='
})
Invoke-C5D1rParserCase 'accepted' $accepted $true
$profile = @(Import-Csv -LiteralPath (Join-Path $WorkRoot `
    'accepted\pbe_tabi_profile.csv'))[0]
$energy = @(Import-Csv -LiteralPath (Join-Path $WorkRoot `
    'accepted\pbe_energy.csv'))[0]
foreach ($column in @(
        'preconditioner_ilu_regularization',
        'preconditioner_ilu_regularization_max',
        'preconditioner_ilu_regularization_selected',
        'preconditioner_ilu_regularization_absolute_max',
        'preconditioner_ilu_regularization_attempts',
        'preconditioner_ilu_regularization_rejected_attempts',
        'preconditioner_ilu_regularization_shifted_rows',
        'preconditioner_ilu_unshifted_failure_stage')) {
    if (-not ($profile.PSObject.Properties.Name -contains $column)) {
        throw "accepted profile is missing C5D.1R column $column"
    }
}
if ([string]$profile.preconditioner_ilu_regularization -ne 'adaptive') {
    throw 'automatic ILU fallback pipeline was not enabled.'
}
if ([string]$profile.convergence_reason -ne 'converged_early_true_residual') {
    throw "accepted solve did not use early true-residual confirmation: $($profile.convergence_reason)"
}
if ([int]$profile.iterations -ge [int]$profile.gmres_restart_effective) {
    throw 'accepted solve did not stop before the end of its GMRES restart window.'
}
if ([int64]$profile.full_true_residual_confirm_calls -le 0) {
    throw 'accepted solve did not record a full true-residual confirmation.'
}
if ([double]::Parse([string]$profile.true_residual_rel,
        [Globalization.CultureInfo]::InvariantCulture) -gt 1.0e-4) {
    throw "accepted solve stopped above the requested true-residual tolerance: $($profile.true_residual_rel)"
}
if ([string]$energy.dgsolv_status -ne 'ok') {
    throw "accepted TABI solvation-energy endpoint failed: $($energy.dgsolv_status)"
}
$rxnKcal = [double]::Parse([string]$energy.E_rxn_kcalmol,
    [Globalization.CultureInfo]::InvariantCulture)
$dgKcal = [double]::Parse([string]$energy.deltaG_solv_kcalmol,
    [Globalization.CultureInfo]::InvariantCulture)
if ([Math]::Abs($rxnKcal - $dgKcal) -gt 1.0e-10 *
        [Math]::Max(1.0, [Math]::Abs($rxnKcal))) {
    throw "single-system TABI deltaG_solv does not match reaction-field energy: $rxnKcal vs $dgKcal"
}

$legacyOff = @($accepted) + @(
    '-pbe-tabi-precond-ilu-regularization=off',
    '-pbe-tabi-precond-ilu-regularization-max=0'
)
Invoke-C5D1rParserCase 'explicit_off' $legacyOff $true
$legacyProfile = @(Import-Csv -LiteralPath (Join-Path $WorkRoot `
    'explicit_off\pbe_tabi_profile.csv'))[0]
if ([string]$legacyProfile.preconditioner_ilu_regularization -ne 'off' -or
    [double]$legacyProfile.preconditioner_ilu_regularization_max -ne 0.0) {
    throw 'explicit off/max options did not disable regularized ILU.'
}

$badMode = @($accepted) + '-pbe-tabi-precond-ilu-regularization=aggressive'
Invoke-C5D1rParserCase 'bad_mode' $badMode $false `
    'pbe-(tabi|bie)-precond-ilu-regularization.*off.*adaptive'

$badMax = @($accepted) + '-pbe-tabi-precond-ilu-regularization-max=0.02'
Invoke-C5D1rParserCase 'bad_max' $badMax $false `
    'pbe-(tabi|bie)-precond-ilu-regularization-max.*0.*1e-2'

$bareOptions = @($accepted | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
if ($bareOptions.Count -ne 0) {
    throw "generated arg.in contains bare options: $($bareOptions -join '; ')"
}
Write-Host 'PASS pbe_fagbi_calderon_c5d1r_regularized_ilu_smoke'
