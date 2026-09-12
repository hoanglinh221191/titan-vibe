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
