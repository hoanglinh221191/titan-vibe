[CmdletBinding()]
param(
    [string]$WorkRoot =
        'test\tmp\pbe_tabi_ses_patch_provenance_benchmark',
    [ValidateRange(1, 30)]
    [int]$TimeoutMinutes = 30
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$WorkRoot = [System.IO.Path]::GetFullPath(
    (Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$exe = Join-Path $WorkRoot 'titan_tabi_ses_patch_provenance_benchmark.exe'
$compileArgs = @(
    '-std=c++23', '-O3', '-march=native', '-DNDEBUG',
    '-Wall', '-Wextra',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_patch_provenance_benchmark.cpp'),
    '-o', $exe
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile SES patch-provenance benchmark.'
}

$stdoutPath = Join-Path $WorkRoot 'provenance_microbench.csv'
$stderrPath = Join-Path $WorkRoot 'provenance_microbench.stderr.txt'
$startInfo = [System.Diagnostics.ProcessStartInfo]::new()
$startInfo.FileName = $exe
$startInfo.UseShellExecute = $false
$startInfo.RedirectStandardOutput = $true
$startInfo.RedirectStandardError = $true
$startInfo.CreateNoWindow = $true
$process = [System.Diagnostics.Process]::new()
$process.StartInfo = $startInfo
if (-not $process.Start()) {
    throw 'Failed to start SES patch-provenance benchmark.'
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
Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding UTF8
Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding UTF8
if (-not $completed) {
    throw "SES patch-provenance benchmark timed out; stdout=$stdoutPath stderr=$stderrPath"
}
if ($exitCode -ne 0) {
    throw "SES patch-provenance benchmark failed; stdout=$stdoutPath stderr=$stderrPath"
}
$rows = @(Import-Csv -LiteralPath $stdoutPath)
if ($rows.Count -ne 1 -or $rows[0].status -ne 'pass') {
    throw "Malformed SES patch-provenance CSV; stdout=$stdoutPath stderr=$stderrPath"
}
Write-Host "csv=$stdoutPath"
Write-Host ($rows[0] | ConvertTo-Csv -NoTypeInformation |
    Select-Object -Skip 1)
