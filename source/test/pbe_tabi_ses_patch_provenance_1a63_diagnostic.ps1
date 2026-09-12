[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Exe,
    [string]$Fixture =
        'test\pbe_inputs_public\1a63_full_charmm_ph7.pqr',
    [Parameter(Mandatory = $true)]
    [string]$WorkRoot,
    [ValidateRange(0.01, 10.0)]
    [double]$Spacing = 1.0,
    [ValidateRange(1, 256)]
    [int]$Threads = 8,
    [ValidateRange(1, 120)]
    [int]$TimeoutMinutes = 30
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot

function Resolve-RepoPath([string]$Path) {
    if ([System.IO.Path]::IsPathRooted($Path)) {
        return [System.IO.Path]::GetFullPath($Path)
    }
    return [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $Path))
}

$Exe = Resolve-RepoPath $Exe
$Fixture = Resolve-RepoPath $Fixture
$WorkRoot = Resolve-RepoPath $WorkRoot
if (-not (Test-Path -LiteralPath $Exe -PathType Leaf)) {
    throw "Diagnostic executable not found: $Exe"
}
if (-not (Test-Path -LiteralPath $Fixture -PathType Leaf)) {
    throw "Fixture not found: $Fixture"
}
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite existing diagnostic folder: $WorkRoot"
}
New-Item -ItemType Directory -Path $WorkRoot | Out-Null

$rawPath = Join-Path $WorkRoot 'provenance_1a63_raw.csv'
$stderrPath = Join-Path $WorkRoot 'provenance_1a63.stderr.txt'
$manifestPath = Join-Path $WorkRoot 'run_manifest.txt'
$spacingText = $Spacing.ToString(
    'R', [System.Globalization.CultureInfo]::InvariantCulture)
$arguments = '"' + $Fixture.Replace('"', '\"') + '" ' +
    $spacingText + ' ' + $Threads

@(
    "fixture=$Fixture"
    "fixture_sha256=$((Get-FileHash -Algorithm SHA256 -LiteralPath $Fixture).Hash)"
    "executable=$Exe"
    "executable_sha256=$((Get-FileHash -Algorithm SHA256 -LiteralPath $Exe).Hash)"
    "spacing=$spacingText"
    "probe_radius=1.4"
    "threads=$Threads"
    "timeout_minutes=$TimeoutMinutes"
    "arguments=$arguments"
    "started_utc=$([DateTime]::UtcNow.ToString('o'))"
) | Set-Content -LiteralPath $manifestPath -Encoding UTF8

$startInfo = [System.Diagnostics.ProcessStartInfo]::new()
$startInfo.FileName = $Exe
$startInfo.Arguments = $arguments
$startInfo.UseShellExecute = $false
$startInfo.RedirectStandardOutput = $true
$startInfo.RedirectStandardError = $true
$startInfo.CreateNoWindow = $true
$process = [System.Diagnostics.Process]::new()
$process.StartInfo = $startInfo
if (-not $process.Start()) {
    throw 'Failed to start SES patch-provenance 1A63 diagnostic.'
}
$stdoutTask = $process.StandardOutput.ReadToEndAsync()
$stderrTask = $process.StandardError.ReadToEndAsync()
$completed = $process.WaitForExit($TimeoutMinutes * 60 * 1000)
if (-not $completed) {
    try {
        $process.Kill($true)
    } catch {
        $process.Kill()
    }
    [void]$process.WaitForExit(10000)
}
$stdoutCompleted = $stdoutTask.Wait(10000)
$stderrCompleted = $stderrTask.Wait(10000)
$stdout = if ($stdoutCompleted) {
    $stdoutTask.GetAwaiter().GetResult()
} else {
    'stdout drain timed out'
}
$stderr = if ($stderrCompleted) {
    $stderrTask.GetAwaiter().GetResult()
} else {
    'stderr drain timed out'
}
$exitCode = if ($process.HasExited) { $process.ExitCode } else { -1 }
$process.Dispose()

Set-Content -LiteralPath $rawPath -Value $stdout -Encoding UTF8
Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding UTF8
Add-Content -LiteralPath $manifestPath -Encoding UTF8 -Value @(
    "finished_utc=$([DateTime]::UtcNow.ToString('o'))"
    "completed=$completed"
    "exit_code=$exitCode"
    "stdout=$rawPath"
    "stderr=$stderrPath"
)
if (-not $completed) {
    throw "1A63 diagnostic timed out; stdout=$rawPath stderr=$stderrPath"
}
if ($exitCode -ne 0) {
    throw "1A63 diagnostic failed; stdout=$rawPath stderr=$stderrPath"
}
$rows = @()
try {
    $rows = @(Import-Csv -LiteralPath $rawPath)
} catch {
    $rows = @()
}
if ($rows.Count -ne 1 -or
    $rows[0].diagnostic_status -ne 'pass') {
    throw "Malformed 1A63 diagnostic CSV; stdout=$rawPath stderr=$stderrPath"
}
Write-Host "csv=$rawPath"
Write-Host "stderr=$stderrPath"
Write-Host ($rows[0] | ConvertTo-Csv -NoTypeInformation |
    Select-Object -Skip 1)
