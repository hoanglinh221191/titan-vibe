param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"

# Some CI/agent environments expose both `Path` and `PATH`. Windows
# Start-Process treats those names case-insensitively and otherwise throws
# before launching TITAN, so normalize the process environment first.
if ([Environment]::OSVersion.Platform -eq [PlatformID]::Win32NT) {
  $processPath = $env:Path
  Remove-Item Env:Path -ErrorAction SilentlyContinue
  $env:Path = $processPath
}

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\energy_overlap_validation_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$sourcePdb = Join-Path $repo "example\2\md1_short.pdb"
$validPdb = Join-Path $work "valid_groups_4_3.pdb"
$fixturePdb = Join-Path $work "overlapping_groups_4_3.pdb"
$validArgPath = Join-Path $work "arg_valid.in"
$selfArgPath = Join-Path $work "arg_valid_self.in"
$argPath = Join-Path $work "arg_overlap.in"
$stdoutPath = Join-Path $work "stdout.txt"
$stderrPath = Join-Path $work "stderr.txt"

$sourceLines = [System.IO.File]::ReadAllLines($sourcePdb)
$modelStart = [Array]::IndexOf($sourceLines, 'MODEL        1')
$modelEnd = [Array]::IndexOf($sourceLines, 'ENDMDL')
if ($modelStart -lt 0 -or $modelEnd -le $modelStart) {
  throw "Could not locate the source PDB model."
}
$header = @($sourceLines[0..($modelStart - 1)])
$referenceModel = @($sourceLines[$modelStart..$modelEnd])
$secondHeader = @($header | ForEach-Object { $_ })
for ($i = 0; $i -lt $secondHeader.Length; ++$i) {
  if ($secondHeader[$i] -match '^TITLE') {
    $secondHeader[$i] = 'TITLE     Title t= 100100.00000 step= 50050000'
  }
}
$overlapModel = @($referenceModel | ForEach-Object { $_ })
$overlapModel[0] = 'MODEL        2'
$validLines = @($header + $referenceModel + $secondHeader + $overlapModel)
[System.IO.File]::WriteAllLines(
    $validPdb, $validLines, [System.Text.Encoding]::ASCII)

$atom1 = -1
$atom3448 = -1
for ($i = 0; $i -lt $overlapModel.Length; ++$i) {
  if ($overlapModel[$i] -match '^ATOM\s+1\s') { $atom1 = $i }
  if ($overlapModel[$i] -match '^ATOM\s+3448\s') { $atom3448 = $i }
}
if ($atom1 -lt 0 -or $atom3448 -lt 0) {
  throw "Could not locate the cross-group fixture atoms."
}

# Group 3 contains atoms 1..3447 and group 4 starts at atom 3448. Copy only
# the fixed-width xyz fields so the two selected atoms overlap exactly.
$coordinates = $overlapModel[$atom1].Substring(30, 24)
$target = $overlapModel[$atom3448]
$overlapModel[$atom3448] =
    $target.Substring(0, 30) + $coordinates + $target.Substring(54)
$lines = @($header + $referenceModel + $secondHeader + $overlapModel)
[System.IO.File]::WriteAllLines(
    $fixturePdb, $lines, [System.Text.Encoding]::ASCII)

# example/2/md1_short.pdb holds only the three protein chains of
# example/2/topol.top (PROA, PROB, PROC: 6788 atoms), not its ions or OPC
# water, and OPC.itp has a [ virtual_sites3 ] section, which the nonbonded
# topology reader rejects (unsupported_topology_self_energy). Use a topology
# of the same force field and molecule files with the protein chains only.
$toppar = (Join-Path $repo "example\2\toppar") -replace '\\', '/'
$proteinTop = Join-Path $work "md1_short_protein.top"
[System.IO.File]::WriteAllLines($proteinTop, @(
    '; Protein chains of example/2/topol.top, matching md1_short.pdb',
    "#include `"$toppar/forcefield.itp`"",
    "#include `"$toppar/PROA.itp`"",
    "#include `"$toppar/PROB.itp`"",
    "#include `"$toppar/PROC.itp`"",
    '',
    '[ system ]',
    'md1_short protein chains A, B, C',
    '',
    '[ molecules ]',
    'PROA 1',
    'PROB 1',
    'PROC 1'
), [System.Text.Encoding]::ASCII)

function New-ArgumentText([string]$PdbPath) {
  return @"
-ff=AMBER14SB
-f=$PdbPath
-1pdb
-g=4-3
-nt=1
-n=$(Join-Path $repo "example\2\index.ndx")
-common_topo=$proteinTop
-path-a=$(Join-Path $repo "lib_ff")
-path-n=$(Join-Path $repo "lib_ff")
"@
}
$argText = New-ArgumentText $fixturePdb
[System.IO.File]::WriteAllText(
    $argPath, $argText, [System.Text.Encoding]::ASCII)
$validArgText = New-ArgumentText $validPdb
[System.IO.File]::WriteAllText(
    $validArgPath, $validArgText, [System.Text.Encoding]::ASCII)
$selfArgText = $validArgText.Replace("-g=4-3", "-g=4-4")
[System.IO.File]::WriteAllText(
    $selfArgPath, $selfArgText, [System.Text.Encoding]::ASCII)

$validStdout = Join-Path $work "valid_stdout.txt"
$validStderr = Join-Path $work "valid_stderr.txt"
Remove-Item -LiteralPath $validStdout, $validStderr -Force `
    -ErrorAction SilentlyContinue
$validTimer = [System.Diagnostics.Stopwatch]::StartNew()
$validProcess = Start-Process -FilePath $titanPath `
    -ArgumentList @("-a", $validArgPath) -WorkingDirectory $work `
    -NoNewWindow -Wait -PassThru -RedirectStandardOutput $validStdout `
    -RedirectStandardError $validStderr
$validTimer.Stop()
if ($validProcess.ExitCode -ne 0) {
  throw "TITAN rejected a valid non-overlapping control frame."
}
Copy-Item -LiteralPath (Join-Path $work "log.txt") `
    -Destination (Join-Path $work "valid_log.txt") -Force

$selfStdout = Join-Path $work "valid_self_stdout.txt"
$selfStderr = Join-Path $work "valid_self_stderr.txt"
Remove-Item -LiteralPath $selfStdout, $selfStderr -Force `
    -ErrorAction SilentlyContinue
$selfProcess = Start-Process -FilePath $titanPath `
    -ArgumentList @("-a", $selfArgPath) -WorkingDirectory $work `
    -NoNewWindow -Wait -PassThru -RedirectStandardOutput $selfStdout `
    -RedirectStandardError $selfStderr
if ($selfProcess.ExitCode -ne 0) {
  throw "TITAN rejected a topology-bound self-interaction selection."
}
if ((Get-Content -LiteralPath (Join-Path $work "cal.xvg") -Raw) -match '(?i)(^|\s)(inf|nan)(\s|$)') {
  throw "TITAN emitted nonfinite self energy."
}

Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue
$process = Start-Process -FilePath $titanPath -ArgumentList @("-a", $argPath) `
    -WorkingDirectory $work -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath

if ($process.ExitCode -eq 0) {
  throw "TITAN unexpectedly accepted overlapping energy atoms."
}
$combined = (Get-Content -LiteralPath $stdoutPath -Raw) + "`n" +
            (Get-Content -LiteralPath $stderrPath -Raw)
if ($combined -notmatch 'invalid_overlapping_atoms') {
  throw "TITAN did not emit invalid_overlapping_atoms."
}
if ($combined -notmatch 'mol1_atom=1, mol2_atom=1') {
  throw "TITAN overlap diagnostic did not identify the selected atom pair."
}

Write-Host ("PASS energy_overlap_validation_smoke valid_control_seconds={0:F3}" `
    -f $validTimer.Elapsed.TotalSeconds)
