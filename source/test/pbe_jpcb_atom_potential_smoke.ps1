param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\jpcb_atom_potential_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$pqrPath = Join-Path $WorkRoot 'single_ion.pqr'
Set-Content -LiteralPath $pqrPath -Encoding ASCII -Value @(
    'ATOM      1 C1   MOL A   1       0.000   0.000   0.000  1.0000 2.0000 C',
    'END'
)
$pqrSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $pqrPath).Hash.ToLowerInvariant()

function Invoke-JpcbAtomCase([string]$Backend) {
    $runDir = Join-Path $WorkRoot $Backend
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $argPath = Write-FmmCliArgFile `
        -RunDir $runDir `
        -RepoRoot $RepoRoot `
        -PqrPath $pqrPath `
        -Backend $Backend `
        -Mode direct `
        -FmmMode off `
        -SurfacePoints 8 `
        -MaxIter 2000 `
        -Restart 32 `
        -Tolerance 1.0e-6 `
        -Threads 1 `
        -SaltM 0.15 `
        -InternalDielectric 1.0 `
        -ExternalDielectric 80.0 `
        -TemperatureK 298.15 `
        -NearIntegration analytic_singular `
        -Preconditioner diagonal `
        -AtomPotentials csv

    $badOptions = @(Get-Content -LiteralPath $argPath | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "$Backend generated non-equals options: $($badOptions -join '; ')"
    }

    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds 180
    if ($run.exit_code -ne 0) {
        $stderr = Get-Content -Raw -LiteralPath (Join-Path $runDir 'stderr.log')
        throw "$Backend atom-potential smoke failed with exit $($run.exit_code): $stderr"
    }

    $atomPath = Join-Path $runDir 'pbe_atom_reaction_potential.csv'
    $residuePath = Join-Path $runDir 'pbe_residue_electrostatics.csv'
    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    foreach ($path in @($atomPath, $residuePath, $energyPath)) {
        if (-not (Test-Path -LiteralPath $path)) {
            throw "$Backend missing expected output: $path"
        }
    }

    $atomRows = @(Import-Csv -LiteralPath $atomPath)
    $residueRows = @(Import-Csv -LiteralPath $residuePath)
    $energyRows = @(Import-Csv -LiteralPath $energyPath)
    if ($atomRows.Count -ne 1 -or $residueRows.Count -ne 1 -or
        $energyRows.Count -ne 1) {
        throw "$Backend unexpected output row counts: atom=$($atomRows.Count), residue=$($residueRows.Count), energy=$($energyRows.Count)"
    }

    $atom = $atomRows[0]
    $energy = $energyRows[0]
    if ($atom.schema_version -ne '2' -or
        $atom.pqr_sha256 -ne $pqrSha256 -or
        $atom.atom_serial -ne '1') {
        throw "$Backend atom identity/SHA256 provenance mismatch."
    }
    if ($atom.surface_model -ne 'vdw' -or
        [string]::IsNullOrWhiteSpace($atom.surface_generator) -or
        [string]::IsNullOrWhiteSpace($atom.probe_radius_A)) {
        throw "$Backend surface provenance is incomplete."
    }
    if ($Backend -eq 'tabi' -and
        (([int]$atom.surface_points_requested) -ne 8 -or
         ([int]$atom.surface_panels) -le 0)) {
        throw 'TABI atom output did not report the requested/generated surface.'
    }
    if ($Backend -eq 'grid' -and
        (([int]$atom.surface_points_requested) -ne 0 -or
         ([int]$atom.surface_panels) -ne 0)) {
        throw 'Grid atom output reported unexpected panel-surface metadata.'
    }
    foreach ($name in @('reaction_potential_kJmol_e',
                         'reaction_energy_contribution_kJmol',
                         'closure_energy_kJmol',
                         'closure_error_kJmol')) {
        $value = Convert-FmmCliDouble $atom.$name
        if (-not (Test-FmmCliFinite $value)) {
            throw "$Backend non-finite atom field ${name}: $($atom.$name)"
        }
    }

    $closure = Convert-FmmCliDouble $atom.closure_energy_kJmol
    $contribution = Convert-FmmCliDouble $atom.reaction_energy_contribution_kJmol
    $closureError = [Math]::Abs((Convert-FmmCliDouble $atom.closure_error_kJmol))
    $selectedEnergy = Convert-FmmCliDouble $energy.E_rxn_kJmol
    if ([Math]::Abs($closure - $contribution) -gt 1.0e-10) {
        throw "$Backend single-atom closure does not equal its contribution."
    }
    if ($closureError -gt 1.0e-9) {
        throw "$Backend closure error $closureError exceeds 1e-9 kJ/mol."
    }
    if ([Math]::Abs($closure - $selectedEnergy) -gt
        [Math]::Max(1.0e-8, 1.0e-10 * [Math]::Abs($selectedEnergy))) {
        throw "$Backend closure $closure does not reproduce selected energy $selectedEnergy."
    }
    if ($atom.backend -ne $Backend -or $atom.component -ne 'complex') {
        throw "$Backend output identity mismatch."
    }
    return [pscustomobject]@{
        backend = $Backend
        energy_kj_mol = $selectedEnergy
        closure_kj_mol = $closure
        closure_error_kj_mol = $closureError
        atom_rows = $atomRows.Count
        residue_rows = $residueRows.Count
        wall_seconds = $run.wall_seconds
    }
}

$rows = @(
    Invoke-JpcbAtomCase -Backend 'tabi'
    Invoke-JpcbAtomCase -Backend 'grid'
)
$rows | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_atom_potential_smoke.csv') `
    -NoTypeInformation

$offRunDir = Join-Path $WorkRoot 'tabi_output_off'
New-FmmCliCleanDir -Path $offRunDir -WorkRoot $WorkRoot
$offArgPath = Write-FmmCliArgFile `
    -RunDir $offRunDir -RepoRoot $RepoRoot -PqrPath $pqrPath `
    -Backend tabi -Mode direct -FmmMode off -SurfacePoints 8 `
    -MaxIter 2000 -Restart 32 -Tolerance 1.0e-6 -Threads 1 `
    -SaltM 0.15 -InternalDielectric 1.0 -ExternalDielectric 80.0 `
    -TemperatureK 298.15 -NearIntegration analytic_singular `
    -Preconditioner diagonal -AtomPotentials off
$offRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $offRunDir `
    -ArgPath $offArgPath -CaseTimeoutSeconds 180
if ($offRun.exit_code -ne 0) {
    throw "TABI atom-output-off control failed with exit $($offRun.exit_code)."
}
if (Test-Path -LiteralPath (Join-Path $offRunDir 'pbe_atom_reaction_potential.csv')) {
    throw 'Atom-potential CSV was written despite -pbe-atom-potentials=off.'
}
$offEnergy = Convert-FmmCliDouble (@(Import-Csv -LiteralPath `
    (Join-Path $offRunDir 'pbe_energy.csv'))[0].E_rxn_kJmol)
$onEnergy = [double](@($rows | Where-Object backend -eq 'tabi')[0].energy_kj_mol)
if ([Math]::Abs($offEnergy - $onEnergy) -gt 1.0e-10) {
    throw "Atom-potential opt-in changed TABI energy: off=$offEnergy on=$onEnergy"
}

Write-Host 'JPCB atom-potential smoke passed.'
Write-Host "Summary: $(Join-Path $WorkRoot 'jpcb_atom_potential_smoke.csv')"
