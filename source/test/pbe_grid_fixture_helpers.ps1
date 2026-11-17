$PbeGridFixtureCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Format-PbeGridPqrNumber([double]$Value) {
    return $Value.ToString('F4', $PbeGridFixtureCulture)
}

function New-PbeGridPqrLine(
    [int]$Serial,
    [string]$AtomName,
    [string]$ResName,
    [string]$Chain,
    [int]$ResId,
    [double]$X,
    [double]$Y,
    [double]$Z,
    [double]$Charge,
    [double]$Radius
) {
    $xText = Format-PbeGridPqrNumber $X
    $yText = Format-PbeGridPqrNumber $Y
    $zText = Format-PbeGridPqrNumber $Z
    $qText = Format-PbeGridPqrNumber $Charge
    $rText = Format-PbeGridPqrNumber $Radius
    return ('ATOM  {0,5} {1,-4} {2,3} {3}{4,4}    {5,8}{6,8}{7,8}{8,8}{9,8}' -f
        $Serial, $AtomName, $ResName, $Chain, $ResId,
        $xText, $yText, $zText, $qText, $rText)
}

function Get-PbeGridFixtureRows([string]$FixtureName) {
    switch ($FixtureName) {
        'mini-protein' {
            return @(
                @(1, 'N',  'ALA', 'A', 1, -4.0000,  0.0000,  0.0000,  0.2500, 1.5500),
                @(2, 'CA', 'ALA', 'A', 1, -2.4000,  0.2000,  0.0000,  0.1000, 1.7000),
                @(3, 'C',  'ALA', 'A', 1, -0.6000,  0.5000, -1.2000,  0.4500, 1.7000),
                @(4, 'O',  'ALA', 'A', 1, -1.1000,  1.4000, -2.0000, -0.5500, 1.5000),
                @(5, 'CB', 'ALA', 'A', 1, -2.8000, -1.2000,  1.4000, -0.1500, 1.7000),
                @(6, 'NZ', 'LYS', 'A', 2,  2.8000,  0.1000,  0.4000,  0.3500, 1.6000)
            )
        }
        'charged-peptide' {
            return @(
                @(1, 'N',  'LYS', 'A', 1, -3.3000, -0.4000,  0.0000,  0.3500, 1.5500),
                @(2, 'CA', 'LYS', 'A', 1, -1.8000,  0.1000,  0.1000,  0.1000, 1.7000),
                @(3, 'C',  'LYS', 'A', 1, -0.4000,  0.4000, -0.8000,  0.4500, 1.7000),
                @(4, 'O',  'LYS', 'A', 1, -0.8000,  1.2000, -1.6000, -0.5000, 1.5000),
                @(5, 'NZ', 'LYS', 'A', 1,  1.7000, -0.9000,  1.1000,  0.6000, 1.6000),
                @(6, 'OD1','ASP', 'A', 2,  3.4000,  0.7000, -0.5000, -0.3000, 1.5000)
            )
        }
        'compact-mini-protein' {
            return @(
                @(1, 'N',  'GLY', 'A', 1, -1.2000,  0.0000,  0.0000,  0.2500, 1.5500),
                @(2, 'CA', 'GLY', 'A', 1,  0.0000,  0.1000,  0.2000,  0.0500, 1.7000),
                @(3, 'C',  'GLY', 'A', 1,  1.2000,  0.0000, -0.1000,  0.4000, 1.7000),
                @(4, 'O',  'GLY', 'A', 1,  1.8000,  0.9000, -0.2000, -0.5000, 1.5000),
                @(5, 'CB', 'ALA', 'A', 2, -0.2000, -1.1000,  0.9000, -0.2000, 1.7000)
            )
        }
        'sphere' {
            return @(, @(
                @(1, 'CA', 'GLY', 'A', 1, 0.0000, 0.0000, 0.0000, 1.0000, 2.0000)
            ))
        }
        'two-atoms' {
            return @(
                @(1, 'CA', 'GLY', 'A', 1, -1.7000, 0.0000, 0.0000,  0.8000, 1.7000),
                @(2, 'O',  'GLY', 'A', 1,  1.7000, 0.0000, 0.0000, -0.8000, 1.5000)
            )
        }
        default {
            throw "Unknown PBE grid fixture: $FixtureName"
        }
    }
}

function Write-PbeGridFixture([string]$FixtureRoot, [string]$FixtureName) {
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    $path = Join-Path $FixtureRoot "${FixtureName}.pqr"
    $lines = @()
    foreach ($row in Get-PbeGridFixtureRows $FixtureName) {
        $lines += New-PbeGridPqrLine `
            -Serial $row[0] -AtomName $row[1] -ResName $row[2] -Chain $row[3] `
            -ResId $row[4] -X $row[5] -Y $row[6] -Z $row[7] `
            -Charge $row[8] -Radius $row[9]
    }
    $lines += 'END'
    Set-Content -LiteralPath $path -Value $lines -Encoding ASCII
    return $path
}

function Initialize-PbeGridFixtures([string]$FixtureRoot, [string[]]$FixtureNames) {
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    foreach ($fixture in $FixtureNames) {
        [void](Write-PbeGridFixture -FixtureRoot $FixtureRoot -FixtureName $fixture)
    }
}

function Initialize-PbeGridForceField([string]$ForceFieldRoot) {
    New-Item -ItemType Directory -Path $ForceFieldRoot -Force | Out-Null
    Set-Content -LiteralPath (Join-Path $ForceFieldRoot 'ffAMBER99SBnonbonded.dat') -Value @(
        '; atom type  symbol  mass  charge  ptype  c6  c12',
        'C C 12.011 0.0 A 0.001 0.000001',
        'N N 14.007 0.0 A 0.001 0.000001',
        'O O 15.999 0.0 A 0.001 0.000001'
    ) -Encoding ASCII
    Set-Content -LiteralPath (Join-Path $ForceFieldRoot 'AMBER99SB_aminoacids.rtp') -Value @(
        '[ ALA ]',
        ' [ atoms ]',
        ' N  N  0.0  1',
        ' CA C  0.0  1',
        ' C  C  0.0  1',
        ' O  O  0.0  1',
        ' CB C  0.0  1',
        ' [ bonds ]'
    ) -Encoding ASCII
}

function Get-PbeGridFixturePath([string]$FixtureRoot, [string]$FixtureName) {
    $path = Join-Path $FixtureRoot "${FixtureName}.pqr"
    if (-not (Test-Path -LiteralPath $path)) {
        [void](Write-PbeGridFixture -FixtureRoot $FixtureRoot -FixtureName $FixtureName)
    }
    return $path
}

# MM/PBSA smoke fixture: a receptor and a ligand, each a five-atom ALA
# (N CA C O CB), chain A and chain B, 4 A apart along z.
#
# The PBSA smokes used to run one ALA with -g=1&1. That makes the receptor and
# the ligand the same atoms: the "complex" PB solve then holds two superposed
# copies of the molecule and E_MM is the molecule's intramolecular energy,
# so the binding numbers mean nothing (MM/GBSA refuses that input with
# invalid_gbsa_identical_selections). Since ENE-02 the run also stops with
# unsupported_topology_self_energy, because -g=1&1 needs topology exclusions.
#
# With -1pdb the index groups are 0 System, 1 protein, 2 chain_A, 3 chain_B,
# so -g=2-3 is receptor = chain A, ligand = chain B. The topology below is
# passed as -common_topo: TITAN then takes charges and LJ from it (log line
# "Charge Source Mol1 : topology"), and its bonds and [ pairs ] are complete,
# so no exclusion or 1-4 is inferred from distances. The charges are synthetic
# (net 0 per ALA); the smokes check the PBSA contract, not ALA chemistry.
$PbsaTwoChainGroups = '2-3'

function Write-PbsaTwoChainPdb([string]$Path) {
    $lines = @(
        'CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'TER',
        'ATOM      6  N   ALA B   2       0.000   0.000   4.000  1.00  0.00           N',
        'ATOM      7  CA  ALA B   2       1.450   0.000   4.000  1.00  0.00           C',
        'ATOM      8  C   ALA B   2       2.050   1.350   4.000  1.00  0.00           C',
        'ATOM      9  O   ALA B   2       1.450   2.350   4.000  1.00  0.00           O',
        'ATOM     10  CB  ALA B   2       1.950  -0.750   5.250  1.00  0.00           C',
        'TER',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Write-PbsaTwoChainTopology([string]$Path) {
    $lines = @(
        '[ defaults ]',
        '; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ',
        '1 2 yes 0.5 0.8333333333',
        '',
        '[ atomtypes ]',
        '; name mass charge ptype sigma(nm) epsilon(kJ/mol)',
        'N 14.007 0.0 A 0.325000 0.711280',
        'C 12.011 0.0 A 0.339967 0.457730',
        'O 15.999 0.0 A 0.295992 0.878640',
        '',
        '[ moleculetype ]',
        'ALA 3',
        '',
        '[ atoms ]',
        '; nr type resnr res atom cgnr charge mass',
        '1 N 1 ALA N  1 -0.30 14.007',
        '2 C 1 ALA CA 1  0.10 12.011',
        '3 C 1 ALA C  2  0.50 12.011',
        '4 O 1 ALA O  2 -0.50 15.999',
        '5 C 1 ALA CB 3  0.20 12.011',
        '',
        '[ bonds ]',
        '1 2 1',
        '2 3 1',
        '3 4 1',
        '2 5 1',
        '',
        '; 1-4 pairs N-CA-C-O and CB-CA-C-O; parameters from gen-pairs',
        '[ pairs ]',
        '1 4 1',
        '5 4 1',
        '',
        '[ system ]',
        'PBSA smoke: ALA receptor (chain A) and ALA ligand (chain B)',
        '',
        '[ molecules ]',
        'ALA 2'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}
