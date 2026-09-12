param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_bie_cartesian_full_fmm'),
    [int]$Threads = 1,
    [int]$CaseTimeoutSeconds = 900,
    [ValidateSet(1, 4)][int[]]$Orders = @(1, 4),
    [double[]]$SaltMolar = @(0.0, 0.15),
    [ValidateSet('centroid', 'node-patch')]
    [string[]]$Discretizations = @('centroid', 'node-patch')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Test-Finite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Assert-UnderWorkRoot([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    $root = [IO.Path]::GetFullPath($WorkRoot).TrimEnd('\') + '\'
    if (-not ($full + '\').StartsWith($root,
            [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside Cartesian full FMM smoke root: $full"
    }
}

function New-CleanDirectory([string]$Path) {
    Assert-UnderWorkRoot $Path
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing TITAN executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
$libRoot = Join-Path $RuntimeRepoRoot 'lib_ff'
if (-not (Test-Path -LiteralPath (Join-Path $libRoot 'ffAMBER99SBnonbonded.dat'))) {
    $libRoot = Join-Path $RuntimeRepoRoot 'example\1\lib_ff'
}
if (-not (Test-Path -LiteralPath (Join-Path $libRoot 'ffAMBER99SBnonbonded.dat'))) {
    throw "Missing force-field directory: $libRoot"
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$rows = New-Object Collections.Generic.List[object]
foreach ($order in $Orders) {
foreach ($salt in $SaltMolar) {
foreach ($discretization in $Discretizations) {
    if (-not (Test-Finite $salt) -or $salt -lt 0.0) {
        throw "Invalid salt concentration: $salt"
    }
    $openingTheta = if ($order -eq 1) { 0.05 } else { 0.55 }
    $saltTag = $salt.ToString('G15', $InvariantCulture).Replace('.', 'p')
    $discretizationTag = $discretization.Replace('-', '')
    $runDir = Join-Path $WorkRoot "order${order}_salt${saltTag}_${discretizationTag}"
    New-CleanDirectory $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    ) | Set-Content -LiteralPath $pdbPath -Encoding ASCII
    $pqrPath = Join-Path $runDir 'four_atom.pqr'
    @(
        'ATOM      1  C   GLY A   1      -6.000   0.000   0.000  1.0000  1.7000',
        'ATOM      2  O   GLY A   1      -2.000   0.000   0.000 -1.0000  1.5000',
        'ATOM      3  N   GLY A   2       2.000   0.000   0.000  1.0000  1.5500',
        'ATOM      4  O   GLY A   2       6.000   0.000   0.000 -1.0000  1.5000',
        'END'
    ) | Set-Content -LiteralPath $pqrPath -Encoding ASCII

    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $libRoot).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
    @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=bie',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-max-iter=400',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        "-pbe-salt=$($salt.ToString('G15', $InvariantCulture))",
        '-pbe-temp=298.15',
        '-pbe-probe=1.4',
        '-pbe-timing=detail',
        '-pbe-bie-fast-operator=p1',
        "-pbe-bie-discretization=$discretization",
        '-pbe-bie-fmm-engine=cartesian-full-experimental',
        "-pbe-bie-fmm-cartesian-order=$order",
        '-pbe-bie-surface-points=8',
        '-pbe-bie-tree-leaf-size=8',
        "-pbe-bie-opening-theta=$($openingTheta.ToString('G15', $InvariantCulture))",
        '-pbe-bie-near-integration=analytic_singular',
        '-pbe-bie-quad-buffer=off',
        '-pbe-bie-precond=diagonal',
        '-pbe-bie-gmres-restart=64',
        '-pbe-bie-matvec-audit=full-small',
        '-pbe-bie-matvec-audit-sample-count=512',
        "-nt=$Threads"
    ) | Set-Content -LiteralPath $argPath -Encoding ASCII

    $stdoutPath = Join-Path $runDir 'stdout.log'
    $stderrPath = Join-Path $runDir 'stderr.log'
    $psi = New-Object Diagnostics.ProcessStartInfo
    $psi.FileName = $TitanExe
    $psi.WorkingDirectory = $runDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = '-a "' + ($argPath -replace '"', '\"') + '"'
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $psi
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
        $process.Kill($true)
        throw "Cartesian full FMM order $order timed out."
    }
    $process.WaitForExit()
    $stdout = $stdoutTask.GetAwaiter().GetResult()
    $stderr = $stderrTask.GetAwaiter().GetResult()
    Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding UTF8
    Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding UTF8
    $exitCode = $process.ExitCode
    $process.Dispose()
    if ($exitCode -ne 0) {
        throw "Cartesian full FMM order $order failed with exit ${exitCode}:`n$stdout`n$stderr"
    }

    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    if ($energy.status -ne 'ok') {
        throw "Cartesian full FMM order $order status=$($energy.status)"
    }
    if ($profile.fmm_local_gmres_active -ne 'true') {
        throw "Cartesian full FMM order $order did not use the prepared GMRES operator."
    }
    if ($profile.fmm_m2l_mode -ne 'cartesian-full-experimental') {
        throw "Cartesian full FMM order $order effective mode=$($profile.fmm_m2l_mode)"
    }
    if ([int]$profile.fmm_m2l_expansion_order -ne $order) {
        throw "Cartesian full FMM requested order $order but reported $($profile.fmm_m2l_expansion_order)."
    }
    if ($profile.bie_cartesian_full_fmm_active -ne 'true' -or
        $profile.bie_cartesian_full_fmm_requested -ne 'cartesian-full-experimental' -or
        $profile.bie_cartesian_full_fmm_effective -ne 'cartesian-full-experimental') {
        throw "Cartesian full FMM order $order requested/effective telemetry is inconsistent."
    }
    $cacheEntries = [int64]$profile.bie_cartesian_full_fmm_cache_entries
    $cacheBytes = [int64]$profile.bie_cartesian_full_fmm_cache_bytes
    $planBytes = [int64]$profile.bie_cartesian_full_fmm_plan_bytes
    if ($cacheEntries -gt 0 -and ($cacheBytes -le 0 -or
            $cacheBytes -gt $planBytes)) {
        throw "Cartesian full FMM order $order cache byte telemetry is inconsistent."
    }
    foreach ($parityColumn in @('bie_cartesian_full_fmm_m2m_parity_error',
            'bie_cartesian_full_fmm_l2l_parity_error')) {
        $parity = Convert-ToDouble $profile.$parityColumn
        if (-not (Test-Finite $parity) -or $parity -gt 1e-11) {
            throw "Cartesian full FMM order $order parity gate failed for ${parityColumn}: $parity"
        }
    }
    $audit = Convert-ToDouble $profile.tree_vs_direct_matvec_rel
    if (-not (Test-Finite $audit) -or $audit -gt 1e-4) {
        throw "Cartesian full FMM order $order direct matvec audit failed: $audit"
    }
    $trueResidual = Convert-ToDouble $profile.true_residual_rel
    if (-not (Test-Finite $trueResidual) -or $trueResidual -gt 1e-4) {
        throw "Cartesian full FMM order $order true residual failed: $trueResidual"
    }
    $rows.Add([pscustomobject]@{
        order = $order
        salt_molar = $salt
        discretization = $discretization
        status = $energy.status
        panels = [int64]$profile.panels
        iterations = [int64]$profile.iterations
        true_residual = $trueResidual
        direct_matvec_relative_l2 = $audit
        energy_kj_mol = Convert-ToDouble $energy.E_rxn_kJmol
        plan_seconds = Convert-ToDouble $profile.fmm_local_prepare_seconds
        matvec_seconds = Convert-ToDouble $profile.fmm_local_gmres_apply_seconds_total
        run_dir = $runDir
    }) | Out-Null
}
}
}

$rows.ToArray() | Export-Csv -LiteralPath (Join-Path $WorkRoot 'matrix.csv') `
    -NoTypeInformation -Encoding ASCII
Write-Host "PASS BIE Cartesian full FMM orders $($Orders -join ', ') for salt $($SaltMolar -join ', ') and discretizations $($Discretizations -join ', '): $WorkRoot"
$rows.ToArray() | Format-Table -AutoSize
