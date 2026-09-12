[CmdletBinding()]
param(
    [string]$Fixture = 'test\pbe_inputs_public\1a63_full_charmm_ph7.pqr',
    [ValidateSet('current','patch-ray-experimental','both')]
    [string]$Engine = 'both',
    [ValidateRange(0.25,4.0)][double]$Spacing = 0.75,
    [ValidateRange(1,256)][int]$Threads = 4,
    [ValidateRange(0,100000000)][int]$MaxAtoms = 0,
    [ValidateRange(1,120)][int]$TimeoutMinutes = 30,
    [string]$WorkRoot = 'test\tmp\pbe_tabi_ses_patch_ray_benchmark',
    [switch]$ProjectionDiagnostics
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Fixture = [IO.Path]::GetFullPath((Join-Path $RepoRoot $Fixture))
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (-not (Test-Path -LiteralPath $Fixture)) { throw "Missing PQR: $Fixture" }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$Exe = Join-Path $WorkRoot 'titan_tabi_ses_patch_ray_benchmark.exe'
$sources = @(
    'src\tabi\tabi_analytic_ses.cpp',
    'src\tabi\tabi_ses_quartic.cpp',
    'src\tabi\tabi_ses_patch_intersections.cpp',
    'src\tabi\tabi_ses_patch_provenance.cpp',
    'src\tabi\tabi_ses_patch_ray.cpp',
    'src\tabi\tabi_ses_bilevel_grid.cpp',
    'src\tabi\tabi_calderon_aux_surface.cpp',
    'src\tabi\tabi_calderon_trace_space.cpp',
    'src\titan_executor.cpp',
    'test\titan_tabi_ses_patch_ray_benchmark.cpp'
) | ForEach-Object { Join-Path $RepoRoot $_ }
$compileArguments = @(
    '-std=c++23', '-O3', '-march=native', '-Wall', '-Wextra'
)
if ($ProjectionDiagnostics) {
    $compileArguments += '-DTITAN_TABI_SES_TEST_HOOKS'
}
$compileArguments += ("-I" + (Join-Path $RepoRoot 'src'))
$compileArguments += $sources
$compileArguments += @('-o', $Exe)
& $Compiler @compileArguments
if ($LASTEXITCODE -ne 0) { throw 'Failed to compile patch-ray benchmark.' }

$engines = if ($Engine -eq 'both') {
    @('current','patch-ray-experimental')
} else { @($Engine) }
$previousDiagnosticEnvironment = $env:TITAN_TABI_SES_DEBUG_PROJECTION
try {
    if ($ProjectionDiagnostics) {
        $env:TITAN_TABI_SES_DEBUG_PROJECTION = '1'
    } else {
        Remove-Item Env:\TITAN_TABI_SES_DEBUG_PROJECTION -ErrorAction SilentlyContinue
    }
    $rows = foreach ($selected in $engines) {
        $stderrPath = Join-Path $WorkRoot ("surface_benchmark_{0}_stderr.log" -f $selected)
        $rawPath = Join-Path $WorkRoot ("surface_benchmark_{0}_raw_stdout.csv" -f $selected)
        $parsedPath = Join-Path $WorkRoot ("surface_benchmark_{0}_parsed.csv" -f $selected)
        $processArguments = @(
            ('"{0}"' -f $Fixture),
            $selected,
            $Spacing.ToString([Globalization.CultureInfo]::InvariantCulture),
            $Threads.ToString([Globalization.CultureInfo]::InvariantCulture),
            $MaxAtoms.ToString([Globalization.CultureInfo]::InvariantCulture)
        )
        $startInfo = [Diagnostics.ProcessStartInfo]::new()
        $startInfo.FileName = $Exe
        $startInfo.Arguments = $processArguments -join ' '
        $startInfo.UseShellExecute = $false
        $startInfo.RedirectStandardOutput = $true
        $startInfo.RedirectStandardError = $true
        $startInfo.CreateNoWindow = $true
        $process = [Diagnostics.Process]::new()
        $process.StartInfo = $startInfo
        if (-not $process.Start()) {
            throw "Failed to start patch-ray benchmark for $selected."
        }
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        $timeoutMilliseconds = [int]($TimeoutMinutes * 60 * 1000)
        $timedOut = -not $process.WaitForExit($timeoutMilliseconds)
        if ($timedOut) {
            try {
                if (-not $process.HasExited) {
                    $process.Kill()
                }
            } catch {
                # The child can exit between the timed wait and Kill().  The
                # blocking wait below remains safe once either event wins.
            }
            $process.WaitForExit()
        }
        $stdoutText = $stdoutTask.GetAwaiter().GetResult()
        $stderrText = $stderrTask.GetAwaiter().GetResult()
        $exitCode = $process.ExitCode
        $stdoutText | Set-Content -LiteralPath $rawPath -Encoding UTF8 -NoNewline
        $stderrText | Set-Content -LiteralPath $stderrPath -Encoding UTF8 -NoNewline
        $process.Dispose()
        if ($timedOut) {
            throw "Patch-ray benchmark timed out after $TimeoutMinutes minute(s) for $selected (raw=$rawPath, stderr=$stderrPath)."
        }
        $output = @(Get-Content -LiteralPath $rawPath)
        $parsed = @()
        if ($output.Count -ge 2) {
            try {
                $parsed = @($output | ConvertFrom-Csv)
                $parsed | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $parsedPath
            } catch {
                # Keep the raw stream and original exit code authoritative.  A
                # malformed diagnostic row must not erase the failing evidence.
                $parsed = @()
            }
        }
        if ($exitCode -ne 0) {
            throw "Patch-ray benchmark failed for $selected (exit=$exitCode, raw=$rawPath, parsed=$parsedPath, stderr=$stderrPath)."
        }
        if ($parsed.Count -eq 0) {
            throw "Patch-ray benchmark returned exit 0 without a parseable CSV row for $selected (raw=$rawPath, stderr=$stderrPath)."
        }
        $parsed
    }
} finally {
    if ($null -eq $previousDiagnosticEnvironment) {
        Remove-Item Env:\TITAN_TABI_SES_DEBUG_PROJECTION -ErrorAction SilentlyContinue
    } else {
        $env:TITAN_TABI_SES_DEBUG_PROJECTION = $previousDiagnosticEnvironment
    }
}
$csv = Join-Path $WorkRoot 'surface_benchmark.csv'
$rows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $csv
$rows | Format-Table engine,status,atoms,threads,panels,build_ms,repair_ms,projection_ms,orientation_certificate_ms,self_intersection_audit_ms,odd_rows
Write-Host "csv=$csv"
