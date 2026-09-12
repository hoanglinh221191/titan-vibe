[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\pbe_tabi_ses_patch_row_worklist_gate',
    [ValidateRange(1, 30)]
    [int]$TimeoutMinutes = 30
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite existing row-worklist gate folder: $WorkRoot"
}
New-Item -ItemType Directory -Path $WorkRoot | Out-Null

$exe = Join-Path $WorkRoot 'titan_tabi_ses_patch_row_worklist_smoke.exe'
$compileArgs = @(
    '-std=c++23', '-O2', '-DNDEBUG', '-Wall', '-Wextra', '-Werror',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_mapping.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_row_worklist.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_patch_row_worklist_smoke.cpp'),
    '-o', $exe
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile SES patch-row-worklist smoke test.'
}

$stdoutPath = Join-Path $WorkRoot 'row_worklist_smoke.stdout.txt'
$stderrPath = Join-Path $WorkRoot 'row_worklist_smoke.stderr.txt'
$startInfo = [System.Diagnostics.ProcessStartInfo]::new()
$startInfo.FileName = $exe
$startInfo.UseShellExecute = $false
$startInfo.RedirectStandardOutput = $true
$startInfo.RedirectStandardError = $true
$startInfo.CreateNoWindow = $true
$process = [System.Diagnostics.Process]::new()
$process.StartInfo = $startInfo
if (-not $process.Start()) {
    throw 'Failed to start SES patch-row-worklist smoke test.'
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
[void]$stdoutTask.Wait(10000)
[void]$stderrTask.Wait(10000)
$stdout = $stdoutTask.GetAwaiter().GetResult()
$stderr = $stderrTask.GetAwaiter().GetResult()
$exitCode = if ($process.HasExited) { $process.ExitCode } else { -1 }
$process.Dispose()
Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding UTF8
Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding UTF8
if (-not $completed) {
    throw "SES patch-row-worklist smoke timed out; stdout=$stdoutPath stderr=$stderrPath"
}
if ($exitCode -ne 0) {
    throw "SES patch-row-worklist smoke failed; stdout=$stdoutPath stderr=$stderrPath"
}
if ($stdout -notmatch 'tabi_ses_patch_row_worklist_smoke=pass') {
    throw "Malformed SES patch-row-worklist output; stdout=$stdoutPath"
}
Write-Host $stdout.Trim()
Write-Host "stdout=$stdoutPath"
Write-Host "stderr=$stderrPath"
