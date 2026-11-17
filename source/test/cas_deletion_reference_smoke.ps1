param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\cas_deletion_reference'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'lib_ff'),
    [string]$ComplexPqr = (Join-Path $PSScriptRoot 'pbe_inputs_binding_jpcb\1brs_complex_amber_ph7.pqr'),
    [double]$Tolerance = 5.0e-4,
    [switch]$Alpb
)

# CAS against explicit side-chain deletion on barnase-barstar (1BRS).
#
# CAS neutralises a side chain (charge, epsilon, sigma = 0; radius = the
# 0.001 A sentinel that removes it from descreening and SASA).
# That is the same as deleting those atoms, so every ddG term from CAS must
# equal GBSA(complex with the side chain removed) - GBSA(wild type), with
# every endpoint recomputed from scratch. One receptor residue and one ligand
# residue cover both scan directions. Before 2026-10-05 the partner's complex
# Born radii were reused from the wild type (polar off by several kcal/mol),
# -eu=kCal/mol output was multiplied by 4.184 instead of divided, and the
# sentinel atoms still entered the SASA as probe-sized spheres (nonpolar off
# by up to 0.56 kcal/mol).
#
# -LibRoot needs AMBER99SB_aminoacids.rtp and ffAMBER99SBnonbonded.dat (the
# repository's lib_ff/ may be empty). About 15 s with three TITAN runs.
# -Alpb repeats the check with ALPB.

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$Retained = @('N', 'H', 'HN', 'CA', 'HA', 'HA1', 'HA2', 'C', 'O', 'OXT', 'CB',
              'HB', 'HB1', 'HB2', 'HB3')
$Mutations = @(
    @{ Chain = 'A'; Resnr = 59; File = 'cas_receptor.dat' },
    @{ Chain = 'D'; Resnr = 39; File = 'cas_ligand.dat' }
)

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Get-PdbLine([string]$PqrLine) {
    $name = $PqrLine.Substring(12, 4).Trim()
    $element = ($name.ToCharArray() | Where-Object { [char]::IsLetter($_) } | Select-Object -First 1)
    return $PqrLine.Substring(0, 54) + '  1.00  0.00          ' + ([string]$element).PadLeft(2)
}

function Write-Complex([string]$Path, [string[]]$Lines) {
    @('MODEL        1') + $Lines + @('ENDMDL', 'END') |
        Set-Content -LiteralPath $Path -Encoding ASCII
}

function Write-Index([string]$Path, [string[]]$Lines) {
    $groups = [ordered]@{}
    foreach ($title in @('System', 'barnase', 'barstar')) {
        $groups[$title] = New-Object System.Collections.Generic.List[int]
    }
    for ($i = 0; $i -lt $Lines.Count; ++$i) {
        $groups['System'].Add($i + 1)
        if ($Lines[$i].Substring(21, 1) -eq 'A') { $groups['barnase'].Add($i + 1) }
        else { $groups['barstar'].Add($i + 1) }
    }
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($title in $groups.Keys) {
        $out.Add("[ $title ]")
        $ids = $groups[$title].ToArray()
        for ($k = 0; $k -lt $ids.Count; $k += 15) {
            $last = [Math]::Min($k + 14, $ids.Count - 1)
            $out.Add((($ids[$k..$last] | ForEach-Object { '{0,5}' -f $_ }) -join ' '))
        }
    }
    $out | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Invoke-Gbsa([string]$CaseDir, [string]$Pdb, [string]$Index, [string[]]$Extra) {
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
    $arg = Join-Path $CaseDir 'arg.in'
    @(
        "-f=$(Convert-ToTitanPath $Pdb)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        "-n=$(Convert-ToTitanPath $Index)",
        '-g=1-2',
        '-mol-1-type=protein',
        '-mol-2-type=protein',
        "-path-a=$lib",
        "-path-n=$lib",
        '-nt=2',
        '-eu=kCal/mol',
        '-gbsa=true'
    ) + @(if ($Alpb) { '-alpb' }) + $Extra | Set-Content -LiteralPath $arg -Encoding ASCII
    Push-Location $CaseDir
    $oldErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        $output = & $TitanExe -a $arg 2>&1
        $code = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $oldErrorActionPreference
        Pop-Location
    }
    if ($code -ne 0) {
        throw "TITAN exited with code $code in $CaseDir. Output:`n$($output -join "`n")"
    }
}

function Read-Snapshot([string]$CaseDir) {
    $rows = Get-Content -LiteralPath (Join-Path $CaseDir 'traj_0_gbsa_snapshots.dat') |
        Where-Object { $_ -match '^\s*[0-9]' }
    if (@($rows).Count -ne 1) { throw "${CaseDir}: expected one GBSA snapshot" }
    $v = @($rows)[0].Split('|') | ForEach-Object { [double]::Parse($_.Trim(), [Globalization.CultureInfo]::InvariantCulture) }
    return @{ vdw = $v[1]; elec = $v[2]; polar = $v[3]; nonpolar = $v[4] }
}

foreach ($path in @($TitanExe, $LibRoot, $ComplexPqr)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing input: $path" }
}
foreach ($file in @('AMBER99SB_aminoacids.rtp', 'ffAMBER99SBnonbonded.dat')) {
    if (-not (Test-Path -LiteralPath (Join-Path $LibRoot $file))) {
        throw "-LibRoot $LibRoot lacks $file; pass a populated force-field library"
    }
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$complex = @(Get-Content -LiteralPath $ComplexPqr | Where-Object { $_ -match '^ATOM' } |
    ForEach-Object { Get-PdbLine $_ })

$wtDir = Join-Path $WorkRoot 'wild_type'
New-Item -ItemType Directory -Path $wtDir -Force | Out-Null
Write-Complex (Join-Path $wtDir 'complex.pdb') $complex
Write-Index (Join-Path $wtDir 'index.ndx') $complex
Invoke-Gbsa $wtDir (Join-Path $wtDir 'complex.pdb') (Join-Path $wtDir 'index.ndx') @('-cas=true')
$wt = Read-Snapshot $wtDir
if (-not (Select-String -Path (Join-Path $wtDir 'cas_summary.dat') -Pattern '^energy_unit: kcal/mol' -Quiet)) {
    throw 'cas_summary.dat does not report energy_unit: kcal/mol'
}

$failures = @()
foreach ($m in $Mutations) {
    $tag = "$($m.Chain)$($m.Resnr)"
    $caseDir = Join-Path $WorkRoot "deleted_$tag"
    New-Item -ItemType Directory -Path $caseDir -Force | Out-Null
    $kept = @($complex | Where-Object {
        -not ($_.Substring(21, 1) -eq $m.Chain -and [int]$_.Substring(22, 4) -eq $m.Resnr -and
              $Retained -notcontains $_.Substring(12, 4).Trim().ToUpperInvariant())
    })
    if ($kept.Count -ge $complex.Count) { throw "${tag}: no side-chain atoms removed" }
    Write-Complex (Join-Path $caseDir 'complex.pdb') $kept
    Write-Index (Join-Path $caseDir 'index.ndx') $kept
    Invoke-Gbsa $caseDir (Join-Path $caseDir 'complex.pdb') (Join-Path $caseDir 'index.ndx') @()
    $mut = Read-Snapshot $caseDir

    $row = Import-Csv -LiteralPath (Join-Path $wtDir $m.File) |
        Where-Object { $_.chain -eq $m.Chain -and [int]$_.resnr -eq $m.Resnr }
    if (-not $row -or $row.mutation_status -ne 'valid_mutation') {
        throw "${tag}: no valid CAS row in $($m.File)"
    }
    $inv = [Globalization.CultureInfo]::InvariantCulture
    foreach ($term in @(
            @{ Name = 'vdw'; Column = 'ddG_vdw_mean' },
            @{ Name = 'elec'; Column = 'ddG_elec_mean' },
            @{ Name = 'polar'; Column = 'ddG_polar_mean' },
            @{ Name = 'nonpolar'; Column = 'ddG_nonpolar_mean' })) {
        $reference = $mut[$term.Name] - $wt[$term.Name]
        $cas = [double]::Parse($row.($term.Column), $inv)
        $line = '{0} {1,-8} CAS {2,12:F4}  deletion {3,12:F4}  diff {4,9:F4}' -f
            $tag, $term.Name, $cas, $reference, ($cas - $reference)
        Write-Host $line
        if ([Math]::Abs($cas - $reference) -gt $Tolerance) { $failures += $line }
    }
}

if ($failures.Count -gt 0) {
    throw "CAS differs from explicit deletion beyond $Tolerance kcal/mol:`n$($failures -join "`n")"
}
Write-Host 'CAS deletion reference smoke PASS'
