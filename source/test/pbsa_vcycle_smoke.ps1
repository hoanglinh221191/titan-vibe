param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbsa_vcycle_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

# -pbe-max-iter is a budget in smoother sweeps: one V-cycle here costs
# (2 + 2) + 120 = 124, so the base smoke's -pbe-max-iter=80 allows exactly one
# cycle. Convergence is judged on the change between two cycles, so one cycle
# can never converge on a charged solute (the old one-ALA fixture had zero
# charges, phi stayed 0 and that hid it). 1240 allows 10 cycles; 3 are used.
$extraArgLines = @(
    '-pbe-multigrid=vcycle',
    '-pbe-mg-levels=2',
    '-pbe-mg-pre=2',
    '-pbe-mg-post=2',
    '-pbe-mg-coarse-iter=120',
    '-pbe-max-iter=1240'
)
$extraCliArgs = @(
    '-pbe-multigrid', 'vcycle',
    '-pbe-mg-levels', '2',
    '-pbe-mg-pre', '2',
    '-pbe-mg-post', '2',
    '-pbe-mg-coarse-iter', '120',
    '-pbe-max-iter', '1240'
)

& (Join-Path $PSScriptRoot 'pbsa_smoke.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -ExtraArgLines $extraArgLines `
    -ExtraCliArgs $extraCliArgs `
    -ExpectedSolverMode 'vcycle'

foreach ($case in @('arg_equals', 'cli_spacing')) {
    $csv = Join-Path ([System.IO.Path]::GetFullPath($WorkRoot)) "$case\mm_pbsa_energy.csv"
    $row = @(Import-Csv -LiteralPath $csv)[0]
    foreach ($column in @('mg_cycles_complex', 'mg_cycles_receptor', 'mg_cycles_ligand')) {
        if ([int]$row.$column -lt 2) {
            throw "Expected at least two V-cycles in $column for $case, observed $($row.$column)"
        }
    }
}

Write-Host "PBSA vcycle smoke PASS"
