param(
    [string]$Distribution = 'Ubuntu',
    [string]$Version = '3.4.1',
    [ValidateRange(1,64)][int]$BuildJobs = 4,
    [switch]$VerifyOnly
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$scriptPath = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'pbe_jpcb_build_apbs_openmp_wsl.sh'))
if (-not (Test-Path -LiteralPath $scriptPath)) { throw "Missing WSL build script: $scriptPath" }
if (-not $VerifyOnly) {
    $wslScript = (& wsl.exe -d $Distribution -- wslpath -a ($scriptPath -replace '\\','/')).Trim()
    if ([string]::IsNullOrWhiteSpace($wslScript)) { throw "Cannot translate path into WSL: $scriptPath" }
    & wsl.exe -d $Distribution -- bash $wslScript $Version $BuildJobs
    if ($LASTEXITCODE -ne 0) { throw 'APBS OpenMP WSL build failed.' }
}
$wslHome = (& wsl.exe -d $Distribution -- /bin/sh -c 'printf "%s" "$HOME"').Trim()
if ([string]::IsNullOrWhiteSpace($wslHome)) { throw "Cannot resolve WSL home directory for '$Distribution'." }
$exe = "$wslHome/.local/apbs-$Version-openmp/bin/apbs"
$manifest = "$wslHome/.local/apbs-$Version-openmp/apbs_openmp_build_manifest.txt"
& wsl.exe -d $Distribution -- test -x $exe
if ($LASTEXITCODE -ne 0) { throw "Missing APBS OpenMP executable: $exe" }
& wsl.exe -d $Distribution -- test -f $manifest
if ($LASTEXITCODE -ne 0) { throw "Missing APBS OpenMP build manifest: $manifest" }
Write-Host "APBS OpenMP WSL executable: $exe"
Write-Host "Build manifest: $manifest"
Write-Host "Use OMP_NUM_THREADS=4 through the JPCB APBS wrapper."
