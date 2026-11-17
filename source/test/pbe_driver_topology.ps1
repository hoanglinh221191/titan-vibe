# The PBE smokes start the engine on a small driver PDB and take the PB source
# from -pbe-test-pqr. With -g=1&1 TITAN also evaluates the driver's own energy,
# and that has required a full topology passed as -common_topo since bb06ac01;
# without one the run stops with unsupported_topology_self_energy.
#
# Write-PbeDriverTopology writes a neutral topology next to the driver, one
# atom per ATOM/HETATM line with the PDB's own atom names, in the same form as
# driver.top in pbe_bie_resolution_ses_smoke.py. It returns the .top path.
function Write-PbeDriverTopology([string]$PdbPath) {
    $atomLines = @(Get-Content -LiteralPath $PdbPath | Where-Object {
        $_.StartsWith('ATOM  ') -or $_.StartsWith('HETATM')
    })
    if ($atomLines.Count -eq 0) {
        throw "Driver PDB has no atoms: $PdbPath"
    }
    $lines = @(
        '[ defaults ]',
        '1 2 yes 0.5 0.8333333333',
        '[ atomtypes ]',
        'T 12.0 0.0 A 0.3 0.1',
        '[ moleculetype ]',
        'DRIVER 3',
        '[ atoms ]'
    )
    for ($i = 0; $i -lt $atomLines.Count; ++$i) {
        # TITAN reads the atom name from columns 13-16 (parse_pdb.cpp).
        $row = $atomLines[$i].PadRight(26)
        $atomName = $row.Substring(12, 4).Trim()
        $resName = $row.Substring(17, 3).Trim()
        $resId = $row.Substring(22, 4).Trim()
        if ([string]::IsNullOrEmpty($resName)) { $resName = 'UNK' }
        if ([string]::IsNullOrEmpty($resId)) { $resId = '1' }
        $lines += "$($i + 1) T $resId $resName $atomName $($i + 1) 0.0 12.0"
    }
    $lines += @(
        '[ system ]',
        'Neutral driver; the PB source is the PQR',
        '[ molecules ]',
        'DRIVER 1'
    )
    $topPath = [IO.Path]::ChangeExtension($PdbPath, '.top')
    Set-Content -LiteralPath $topPath -Value $lines -Encoding ASCII
    return $topPath
}
