[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_analytic_ses_smoke',
    [ValidateRange(30,600)][int]$CaseTimeoutSeconds = 180,
    [ValidateSet('centroid','node-patch')][string]$Discretization = 'centroid',
    [ValidateSet('current','patch-ray-experimental')][string]$Engine = 'current'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $TitanExe))
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Titan executable not found: $TitanExe"
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$pqr = Join-Path $WorkRoot 'analytic_triple.pqr'
$root3 = [Math]::Sqrt(3.0)
@(
    ('ATOM      1  C1  MOL A   1      {0,8:F3} {1,8:F3}   0.000   0.2500 1.0000' -f -1.5,(-0.5*$root3)),
    ('ATOM      2  C2  MOL A   1      {0,8:F3} {1,8:F3}   0.000  -0.2500 1.0000' -f 1.5,(-0.5*$root3)),
    ('ATOM      3  C3  MOL A   1       0.000 {0,8:F3}   0.000   0.0000 1.0000' -f $root3),
    'END'
) | Set-Content -LiteralPath $pqr -Encoding ASCII

$argPath = Write-FmmCliArgFile -RunDir $WorkRoot -RepoRoot $RepoRoot `
    -PqrPath $pqr -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
    -SurfacePoints 8 -MaxIter 4 -Restart 8 -Tolerance 1e-4 -Threads 2 `
    -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 `
    -TemperatureK 298.15 -Preconditioner diagonal `
    -NearIntegration analytic_singular -PrimarySurface analytic-ses `
    -SesGridSpacing 0.25
Add-Content -LiteralPath $argPath -Value '-pbe-tabi-surface-only=true' -Encoding ASCII
Add-Content -LiteralPath $argPath -Value "-pbe-tabi-discretization=$Discretization" -Encoding ASCII
Add-Content -LiteralPath $argPath -Value "-pbe-tabi-analytic-ses-engine=$Engine" -Encoding ASCII

$run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $WorkRoot `
    -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
if ($run.timed_out -or $run.exit_code -ne 0) {
    $stderr = (Get-Content (Join-Path $WorkRoot 'stderr.log') -ErrorAction SilentlyContinue) -join ' '
    throw "Analytic SES smoke failed (exit=$($run.exit_code)): $stderr"
}
$profile = @(Import-Csv (Join-Path $WorkRoot 'pbe_tabi_profile.csv')) | Select-Object -Last 1
$quality = @(Import-Csv (Join-Path $WorkRoot 'pbe_tabi_surface_quality.csv')) | Select-Object -Last 1
if ([string]$profile.surface_generator -ne 'titan_analytic_patch_ray_ses' -or
    [string]$profile.convergence_reason -ne 'surface_only') {
    throw "Analytic SES profile contract failed: generator=$($profile.surface_generator) reason=$($profile.convergence_reason)"
}
if ([string]$quality.calderon_trace_space_status -ne 'pass' -or
    [int64]$quality.calderon_trace_boundary_edges -ne 0 -or
    [int64]$quality.calderon_trace_nonmanifold_edges -ne 0) {
    throw "Analytic SES topology failed: status=$($quality.calderon_trace_space_status) boundary=$($quality.calderon_trace_boundary_edges) nonmanifold=$($quality.calderon_trace_nonmanifold_edges)"
}
if ([double]$quality.triangle_aspect_max -gt 500.0) {
    throw "Analytic SES retained a solver-hostile aspect ratio: $($quality.triangle_aspect_max)"
}
if ([string]$quality.analytic_ses_engine -ne $Engine -or
    [string]$quality.surface_status_detail -ne 'pass' -or
    [string]$quality.cavity_policy -ne 'external-only') {
    throw "Analytic SES telemetry contract failed: engine=$($quality.analytic_ses_engine) status=$($quality.surface_status_detail) cavity=$($quality.cavity_policy)"
}
foreach ($timing in @('initial_collapse_ms','stable_sliver_ms','full_projection_ms',
                      'orientation_repair_ms','wide_repair_ms','self_intersection_audit_ms')) {
    if ($quality.PSObject.Properties.Name -notcontains $timing) {
        throw "Analytic SES telemetry column is missing: $timing"
    }
}
if ($Engine -eq 'patch-ray-experimental' -and
    ([int64]$quality.odd_ray_rows -ne 0 -or
     [int64]$quality.patch_ray_unresolved_uncertified_rows -ne 0 -or
     [int64]$quality.morphology_fallback_vertices -ne 0 -or
     [string]$quality.analytic_orientation_status -ne 'pass' -or
     [int64]$quality.analytic_inverted_faces -ne 0 -or
     [int64]$quality.analytic_orientation_uncertified_faces -ne 0 -or
     [int64]$quality.analytic_self_intersections -ne 0)) {
    throw "Patch-ray SES certificate telemetry failed: odd=$($quality.odd_ray_rows) unresolved=$($quality.patch_ray_unresolved_uncertified_rows) morphology=$($quality.morphology_fallback_vertices) orientation=$($quality.analytic_orientation_status) inverted=$($quality.analytic_inverted_faces) uncertified=$($quality.analytic_orientation_uncertified_faces) intersections=$($quality.analytic_self_intersections)"
}
if ($Discretization -eq 'node-patch' -and
    ([int64]$profile.panels -ne [int64]$profile.mesh_vertices -or
     [int64]$profile.panels -ge [int64]$profile.mesh_triangles)) {
    throw "Analytic SES node-patch size contract failed: panels=$($profile.panels) vertices=$($profile.mesh_vertices) triangles=$($profile.mesh_triangles)"
}
Write-Host "pbe_tabi_analytic_ses_smoke=pass engine=$Engine panels=$($quality.panel_count) aspect_max=$($quality.triangle_aspect_max)"
