[CmdletBinding()]
param(
    [string]$WorkRoot =
        'test\tmp\pbe_tabi_ses_patch_projection_gate',
    [ValidateRange(1, 30)]
    [int]$TimeoutMinutes = 30
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$WorkRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite existing projection gate folder: $WorkRoot"
}
New-Item -ItemType Directory -Path $WorkRoot | Out-Null

$exe = Join-Path $WorkRoot 'titan_tabi_ses_patch_projection_smoke.exe'
$sources = @(
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_mapping.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_projection.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_patch_projection_smoke.cpp')
)
$compileArgs = @(
    '-std=c++23', '-O2', '-DNDEBUG', '-Wall', '-Wextra', '-Werror',
    '-fno-fast-math', '-frounding-math', '-ffp-contract=off'
) + $sources + @(
    '-o', $exe
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile SES patch-projection smoke test.'
}

$stdoutPath = Join-Path $WorkRoot 'projection_smoke.stdout.txt'
$stderrPath = Join-Path $WorkRoot 'projection_smoke.stderr.txt'
$startInfo = [System.Diagnostics.ProcessStartInfo]::new()
$startInfo.FileName = $exe
$startInfo.UseShellExecute = $false
$startInfo.RedirectStandardOutput = $true
$startInfo.RedirectStandardError = $true
$startInfo.CreateNoWindow = $true
$process = [System.Diagnostics.Process]::new()
$process.StartInfo = $startInfo
if (-not $process.Start()) {
    throw 'Failed to start SES patch-projection smoke test.'
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
    throw "SES patch-projection smoke timed out; stdout=$stdoutPath stderr=$stderrPath"
}
if ($exitCode -ne 0) {
    throw "SES patch-projection smoke failed; stdout=$stdoutPath stderr=$stderrPath"
}
if ($stdout -notmatch 'tabi_ses_patch_projection_smoke=pass') {
    throw "Malformed SES patch-projection output; stdout=$stdoutPath"
}

$fastExe = Join-Path $WorkRoot 'titan_tabi_ses_patch_projection_fast_math_smoke.exe'
$fastCompileArgs = @(
    '-std=c++23', '-O2', '-DNDEBUG', '-Wall', '-Wextra', '-Werror',
    '-ffast-math', '-Wno-nan-infinity-disabled'
) + $sources + @(
    '-o', $fastExe
)
& $Compiler @fastCompileArgs
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile fast-math SES patch-projection contract test.'
}
$fastStdout = & $fastExe '--rounding-contract-only'
if ($LASTEXITCODE -ne 0 -or
    ($fastStdout -join "`n") -notmatch 'c2b_rounding_contract_only=pass fast_math=1') {
    throw 'Fast-math SES patch-projection rounding contract did not fail closed.'
}
Write-Host $stdout.Trim()
Write-Host (($fastStdout -join "`n").Trim())
Write-Host "stdout=$stdoutPath"
Write-Host "stderr=$stderrPath"
