param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\jpcb_binding_energy_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$receptorLine = 'ATOM      1 C1   REC A   1      -3.000   0.000   0.000  1.0000 1.5000 C'
$ligandLine = 'ATOM      2 C2   LIG B   2       3.000   0.000   0.000 -1.0000 1.5000 C'
$complexPqr = Join-Path $WorkRoot 'complex.pqr'
$receptorPqr = Join-Path $WorkRoot 'receptor.pqr'
$ligandPqr = Join-Path $WorkRoot 'ligand.pqr'
$badLigandPqr = Join-Path $WorkRoot 'ligand_bad_partition.pqr'
Set-Content -LiteralPath $complexPqr -Encoding ASCII -Value @($receptorLine, $ligandLine, 'END')
Set-Content -LiteralPath $receptorPqr -Encoding ASCII -Value @($receptorLine, 'END')
Set-Content -LiteralPath $ligandPqr -Encoding ASCII -Value @($ligandLine, 'END')
Set-Content -LiteralPath $badLigandPqr -Encoding ASCII -Value @(
    'ATOM      2 C2   LIG B   2       3.000   0.000   0.000 -0.5000 1.5000 C',
    'END'
)
$complexSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $complexPqr).Hash.ToLowerInvariant()
$receptorSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $receptorPqr).Hash.ToLowerInvariant()
$ligandSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $ligandPqr).Hash.ToLowerInvariant()

function New-JpcbBindingArgFile(
    [string]$RunDir,
    [string]$Backend,
    [string]$ComplexPath = $complexPqr,
    [string]$ReceptorPath = $receptorPqr,
    [string]$LigandPath = $ligandPqr
) {
    return Write-FmmCliArgFile `
        -RunDir $RunDir `
        -RepoRoot $RepoRoot `
        -PqrPath $ComplexPath `
        -Backend $Backend `
        -Mode direct `
        -FmmMode off `
        -SurfacePoints 8 `
        -MaxIter 4000 `
        -Restart 64 `
        -Tolerance 1.0e-6 `
        -Threads 1 `
        -SaltM 0.15 `
        -InternalDielectric 1.0 `
        -ExternalDielectric 80.0 `
        -TemperatureK 298.15 `
        -NearIntegration analytic_singular `
        -Preconditioner diagonal `
        -AtomPotentials csv `
        -BindingCycle pqr-components `
        -BindingReceptorPqr $ReceptorPath `
        -BindingLigandPqr $LigandPath
}

function Assert-Near([double]$Actual, [double]$Expected, [double]$Tolerance,
                     [string]$Label) {
    if (-not (Test-FmmCliFinite $Actual) -or
        [Math]::Abs($Actual - $Expected) -gt $Tolerance) {
        throw "$Label mismatch: actual=$Actual expected=$Expected tolerance=$Tolerance"
    }
}

function Invoke-JpcbBindingCase([string]$Backend) {
    $runDir = Join-Path $WorkRoot $Backend
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $argPath = New-JpcbBindingArgFile -RunDir $runDir -Backend $Backend
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds 300
    if ($run.exit_code -ne 0) {
        $stderr = Get-Content -Raw -LiteralPath (Join-Path $runDir 'stderr.log')
        $stdout = Get-Content -Raw -LiteralPath (Join-Path $runDir 'stdout.log')
        throw "$Backend binding smoke failed with exit $($run.exit_code): $stderr $stdout"
    }

    $atomRows = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_atom_reaction_potential.csv'))
    $componentRows = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_component_electrostatics.csv'))
    $bindingRows = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_binding_electrostatics.csv'))
    if ($atomRows.Count -ne 4 -or $componentRows.Count -ne 3 -or
        $bindingRows.Count -ne 1) {
        throw "$Backend row-count mismatch: atoms=$($atomRows.Count), components=$($componentRows.Count), binding=$($bindingRows.Count)"
    }
    foreach ($component in @('complex','receptor','ligand')) {
        if (@($componentRows | Where-Object component -eq $component).Count -ne 1) {
            throw "$Backend missing unique $component component row."
        }
    }
    if (@($atomRows | Where-Object component -eq 'complex').Count -ne 2 -or
        @($atomRows | Where-Object component -eq 'receptor').Count -ne 1 -or
        @($atomRows | Where-Object component -eq 'ligand').Count -ne 1) {
        throw "$Backend atomwise component coverage mismatch."
    }

    $binding = $bindingRows[0]
    if ($binding.schema_version -ne '3' -or
        $binding.complex_pqr_sha256 -ne $complexSha256 -or
        $binding.receptor_pqr_sha256 -ne $receptorSha256 -or
        $binding.ligand_pqr_sha256 -ne $ligandSha256) {
        throw "$Backend binding SHA256 provenance mismatch."
    }
    $expectedHashes = @{
        complex = $complexSha256
        receptor = $receptorSha256
        ligand = $ligandSha256
    }
    foreach ($componentRow in $componentRows) {
        if ($componentRow.schema_version -ne '6' -or
            $componentRow.pqr_sha256 -ne $expectedHashes[$componentRow.component]) {
            throw "$Backend component SHA256 provenance mismatch for $($componentRow.component)."
        }
        foreach ($timerColumn in @('surface_seconds','rhs_seconds',
                'operator_audit_seconds','ilu_setup_seconds')) {
            if ($null -eq $componentRow.$timerColumn -or
                [double]$componentRow.$timerColumn -lt 0.0) {
                throw "$Backend invalid component timer $timerColumn for $($componentRow.component)."
            }
        }
    }
    foreach ($atomRow in $atomRows) {
        if ($atomRow.schema_version -ne '2' -or
            $atomRow.pqr_sha256 -ne $expectedHashes[$atomRow.component]) {
            throw "$Backend atom SHA256 provenance mismatch for $($atomRow.component)."
        }
    }
    $complexSerials = @($atomRows | Where-Object component -eq 'complex' |
        ForEach-Object { [int]$_.atom_serial } | Sort-Object)
    if (($complexSerials -join ',') -ne '1,2' -or
        [int]($atomRows | Where-Object component -eq 'receptor')[0].atom_serial -ne 1 -or
        [int]($atomRows | Where-Object component -eq 'ligand')[0].atom_serial -ne 2) {
        throw "$Backend did not preserve original PQR atom serials."
    }
    if ($binding.surface_model -ne 'vdw' -or
        [string]::IsNullOrWhiteSpace($binding.surface_generator)) {
        throw "$Backend binding surface provenance is incomplete."
    }
    if ($Backend -eq 'tabi' -and
        (([int]$binding.surface_points_requested) -ne 8 -or
         ([int]$binding.complex_surface_panels) -le 0 -or
         ([int]$binding.receptor_surface_panels) -le 0 -or
         ([int]$binding.ligand_surface_panels) -le 0)) {
        throw 'TABI binding output did not report all generated surfaces.'
    }
    if ($Backend -eq 'grid' -and
        (([int]$binding.complex_surface_panels) -ne 0 -or
         ([int]$binding.receptor_surface_panels) -ne 0 -or
         ([int]$binding.ligand_surface_panels) -ne 0)) {
        throw 'Grid binding output reported unexpected panel-surface metadata.'
    }
    $complexEnergy = Convert-FmmCliDouble $binding.Gpol_complex_kJmol
    $receptorEnergy = Convert-FmmCliDouble $binding.Gpol_receptor_kJmol
    $ligandEnergy = Convert-FmmCliDouble $binding.Gpol_ligand_kJmol
    $deltaPol = Convert-FmmCliDouble $binding.Delta_Gpol_bind_kJmol
    $crossCoulomb = Convert-FmmCliDouble $binding.E_coul_receptor_ligand_kJmol
    $deltaElec = Convert-FmmCliDouble $binding.Delta_Gelec_bind_kJmol
    $expectedDeltaPol = $complexEnergy - $receptorEnergy - $ligandEnergy
    $expectedCross = 332.063713299 * 4.184 * -1.0 / 6.0
    Assert-Near $deltaPol $expectedDeltaPol 1.0e-8 "$Backend polar binding cycle"
    Assert-Near $crossCoulomb $expectedCross 1.0e-8 "$Backend cross Coulomb"
    Assert-Near $deltaElec ($deltaPol + $crossCoulomb) 1.0e-8 "$Backend total electrostatics"
    if ((Convert-FmmCliDouble $binding.closure_error_max_kJmol) -gt 1.0e-8) {
        throw "$Backend component atom-energy closure exceeds 1e-8 kJ/mol."
    }
    if ($binding.status -ne 'ok' -or $binding.backend -ne $Backend) {
        throw "$Backend binding status or backend label mismatch."
    }

    return [pscustomobject]@{
        backend = $Backend
        delta_g_pol_kj_mol = $deltaPol
        cross_coulomb_kj_mol = $crossCoulomb
        delta_g_elec_kj_mol = $deltaElec
        atom_rows = $atomRows.Count
        component_rows = $componentRows.Count
        closure_error_max_kj_mol = Convert-FmmCliDouble $binding.closure_error_max_kJmol
        wall_seconds = $run.wall_seconds
    }
}

$rows = @(
    Invoke-JpcbBindingCase -Backend 'tabi'
    Invoke-JpcbBindingCase -Backend 'grid'
)
$rows | Export-Csv -LiteralPath (Join-Path $WorkRoot 'jpcb_binding_energy_smoke.csv') `
    -NoTypeInformation

$swapRunDir = Join-Path $WorkRoot 'swap_invariance'
New-FmmCliCleanDir -Path $swapRunDir -WorkRoot $WorkRoot
$swapArgPath = New-JpcbBindingArgFile -RunDir $swapRunDir -Backend 'tabi' `
    -ReceptorPath $ligandPqr -LigandPath $receptorPqr
$swapRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $swapRunDir `
    -ArgPath $swapArgPath -CaseTimeoutSeconds 300
if ($swapRun.exit_code -ne 0) {
    throw "Receptor/ligand swap control failed with exit $($swapRun.exit_code)."
}
$swapBinding = @(Import-Csv -LiteralPath `
    (Join-Path $swapRunDir 'pbe_binding_electrostatics.csv'))[0]
$tabiBinding = @(Import-Csv -LiteralPath `
    (Join-Path $WorkRoot 'tabi\pbe_binding_electrostatics.csv'))[0]
Assert-Near (Convert-FmmCliDouble $swapBinding.E_coul_receptor_ligand_kJmol) `
    (Convert-FmmCliDouble $tabiBinding.E_coul_receptor_ligand_kJmol) 1.0e-10 `
    'Receptor/ligand swap Coulomb invariance'
Assert-Near (Convert-FmmCliDouble $swapBinding.Delta_Gpol_bind_kJmol) `
    (Convert-FmmCliDouble $tabiBinding.Delta_Gpol_bind_kJmol) 1.0e-10 `
    'Receptor/ligand swap polar-cycle invariance'

$zeroLigandLine = 'ATOM      2 C2   LIG B   2       3.000   0.000   0.000  0.0000 1.5000 C'
$zeroComplexPqr = Join-Path $WorkRoot 'complex_zero_ligand_charge.pqr'
$zeroLigandPqr = Join-Path $WorkRoot 'ligand_zero_charge.pqr'
Set-Content -LiteralPath $zeroComplexPqr -Encoding ASCII `
    -Value @($receptorLine, $zeroLigandLine, 'END')
Set-Content -LiteralPath $zeroLigandPqr -Encoding ASCII `
    -Value @($zeroLigandLine, 'END')
$zeroRunDir = Join-Path $WorkRoot 'zero_charge_component'
New-FmmCliCleanDir -Path $zeroRunDir -WorkRoot $WorkRoot
$zeroArgPath = New-JpcbBindingArgFile -RunDir $zeroRunDir -Backend 'tabi' `
    -ComplexPath $zeroComplexPqr -LigandPath $zeroLigandPqr
$zeroRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $zeroRunDir `
    -ArgPath $zeroArgPath -CaseTimeoutSeconds 300
if ($zeroRun.exit_code -ne 0) {
    throw "Zero-charge component control failed with exit $($zeroRun.exit_code)."
}
$zeroBinding = @(Import-Csv -LiteralPath `
    (Join-Path $zeroRunDir 'pbe_binding_electrostatics.csv'))[0]
Assert-Near (Convert-FmmCliDouble $zeroBinding.E_coul_receptor_ligand_kJmol) `
    0.0 1.0e-12 'Zero-charge component cross Coulomb'

$badRunDir = Join-Path $WorkRoot 'partition_rejection'
New-FmmCliCleanDir -Path $badRunDir -WorkRoot $WorkRoot
$badArgPath = New-JpcbBindingArgFile -RunDir $badRunDir -Backend 'tabi' `
    -LigandPath $badLigandPqr
$badRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $badRunDir `
    -ArgPath $badArgPath -CaseTimeoutSeconds 60
if ($badRun.exit_code -eq 0) {
    throw 'Malformed receptor/ligand partition was accepted unexpectedly.'
}
$badText = ((Get-Content -Raw -LiteralPath (Join-Path $badRunDir 'stdout.log')) +
            (Get-Content -Raw -LiteralPath (Join-Path $badRunDir 'stderr.log')))
if ($badText -notmatch 'atomwise partition') {
    throw "Malformed partition rejection did not report the expected reason: $badText"
}

Write-Host 'JPCB binding-energy smoke passed.'
Write-Host "Summary: $(Join-Path $WorkRoot 'jpcb_binding_energy_smoke.csv')"
