[CmdletBinding()]
param(
    [string]$WorkRoot = 'test\tmp\pbe_tabi_ses_row_worklist_real_benchmark',
    [ValidateRange(5, 7)][int]$Samples = 5,
    [ValidateRange(1, 30)][int]$TimeoutMinutes = 30,
    [ValidateRange(1, 1800)][int]$AuditBudgetSeconds = 900,
    [ValidateSet('oneA63', 'huge8192', 'auto32768')][string]$MaximumStage = 'oneA63',
    [switch]$CompileOnly
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
if (Test-Path -LiteralPath $WorkRoot) {
    throw "Refusing to overwrite existing D2d benchmark folder: $WorkRoot"
}
New-Item -ItemType Directory -Path $WorkRoot | Out-Null

$oneA63 = [IO.Path]::GetFullPath((Join-Path $RepoRoot 'test\pbe_inputs_public\1a63_full_charmm_ph7.pqr'))
$huge = [IO.Path]::GetFullPath((Join-Path $RepoRoot 'example\pbe_huge\sphere_6co_full.pdb'))
$expectedOneA63 = '5058896FBA0C2BFDE704AFE2BA23D74E2AAFCEA1F1A95FAC67CE643E9264203B'
$expectedHuge = 'FA7B83D8D9A0B6CB6A6DBF07A3F93004B76951B51E8562184139060D3752D034'
foreach ($inputPath in @($oneA63, $huge)) {
    if (-not (Test-Path -LiteralPath $inputPath)) {
        throw "Missing D2d input: $inputPath"
    }
}
$oneA63Hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $oneA63).Hash.ToUpperInvariant()
$hugeHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $huge).Hash.ToUpperInvariant()
$hashLockError = $null
if ($oneA63Hash -ne $expectedOneA63) {
    $hashLockError = "1A63 SHA-256 mismatch: expected=$expectedOneA63 actual=$oneA63Hash"
} elseif ($hugeHash -ne $expectedHuge) {
    $hashLockError = "Huge PDB SHA-256 mismatch: expected=$expectedHuge actual=$hugeHash"
}

function Invoke-CapturedProcess {
    param(
        [Parameter(Mandatory)][string]$FileName,
        [Parameter(Mandatory)][string[]]$Arguments,
        [Parameter(Mandatory)][string]$StdoutPath,
        [Parameter(Mandatory)][string]$StderrPath,
        [Parameter(Mandatory)][int]$TimeoutMilliseconds
    )
    $startInfo = [Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = $FileName
    $startInfo.UseShellExecute = $false
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.CreateNoWindow = $true
    $quote = [char]34
    $slash = [char]92
    $quotedArguments = foreach ($argument in $Arguments) {
        $builder = [Text.StringBuilder]::new()
        [void]$builder.Append($quote)
        $backslashes = 0
        foreach ($character in $argument.ToCharArray()) {
            if ($character -eq $slash) {
                $backslashes++
                continue
            }
            if ($character -eq $quote) {
                [void]$builder.Append(($slash.ToString() * (2 * $backslashes + 1)))
                [void]$builder.Append($quote)
            } else {
                [void]$builder.Append(($slash.ToString() * $backslashes))
                [void]$builder.Append($character)
            }
            $backslashes = 0
        }
        [void]$builder.Append(($slash.ToString() * (2 * $backslashes)))
        [void]$builder.Append($quote)
        $builder.ToString()
    }
    $startInfo.Arguments = $quotedArguments -join ' '
    $stdoutStream = [IO.File]::Open($StdoutPath, [IO.FileMode]::Create,
        [IO.FileAccess]::Write, [IO.FileShare]::Read)
    $stderrStream = [IO.File]::Open($StderrPath, [IO.FileMode]::Create,
        [IO.FileAccess]::Write, [IO.FileShare]::Read)
    $process = [Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    $timer = [Diagnostics.Stopwatch]::StartNew()
    if (-not $process.Start()) {
        $stdoutStream.Dispose()
        $stderrStream.Dispose()
        throw "Failed to start process: $FileName"
    }
    $stdoutTask = $process.StandardOutput.BaseStream.CopyToAsync($stdoutStream)
    $stderrTask = $process.StandardError.BaseStream.CopyToAsync($stderrStream)
    $deadline = [DateTime]::UtcNow.AddMilliseconds($TimeoutMilliseconds)
    $peakWorkingSet = [int64]0
    $completed = $false
    while (-not $process.HasExited) {
        try {
            $process.Refresh()
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [int64]$process.WorkingSet64)
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [int64]$process.PeakWorkingSet64)
        } catch { }
        try {
            $snapshot = Get-Process -Id $process.Id -ErrorAction Stop
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [int64]$snapshot.WorkingSet64)
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [int64]$snapshot.PeakWorkingSet64)
        } catch { }
        $remainingMilliseconds = ($deadline - [DateTime]::UtcNow).TotalMilliseconds
        if ($remainingMilliseconds -le 0.0) {
            break
        }
        $pollMilliseconds = [Math]::Min(150, [Math]::Max(1, [int][Math]::Ceiling($remainingMilliseconds)))
        if ($process.WaitForExit($pollMilliseconds)) {
            $completed = $true
            break
        }
    }
    if ($process.HasExited) {
        $completed = $true
    }
    if (-not $completed) {
        try {
            $process.Kill($true)
        } catch {
            try { $process.Kill() } catch { }
        }
        [void]$process.WaitForExit(10000)
    }
    try { [void]$stdoutTask.Wait(10000) } catch { }
    try { [void]$stderrTask.Wait(10000) } catch { }
    $timer.Stop()
    try { $stdoutStream.Flush() } catch { }
    try { $stderrStream.Flush() } catch { }
    $stdoutStream.Dispose()
    $stderrStream.Dispose()
    $exitCode = if ($process.HasExited) { $process.ExitCode } else { -1 }
    $process.Dispose()
    [pscustomobject]@{
        ExitCode = $exitCode
        TimedOut = -not $completed
        WallMs = $timer.Elapsed.TotalMilliseconds
        PeakWorkingSetBytes = $peakWorkingSet
    }
}

$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$Exe = Join-Path $WorkRoot 'titan_tabi_ses_row_worklist_real_benchmark.exe'
$compileStdout = Join-Path $WorkRoot 'compile.stdout.txt'
$compileStderr = Join-Path $WorkRoot 'compile.stderr.txt'
$compilerVersionStdout = Join-Path $WorkRoot 'compiler_version.stdout.txt'
$compilerVersionStderr = Join-Path $WorkRoot 'compiler_version.stderr.txt'
$SourceInputs = @(
    Get-ChildItem -LiteralPath (Join-Path $RepoRoot 'src\tabi') -File |
        Where-Object { $_.Name -like 'tabi_ses_*' -or $_.Name -like 'tabi_analytic_ses*' } |
        Select-Object -ExpandProperty FullName
    (Join-Path $RepoRoot 'src\titan_executor.cpp')
    (Join-Path $RepoRoot 'src\titan_executor.h')
    (Join-Path $RepoRoot 'test\titan_tabi_ses_row_worklist_real_benchmark.cpp')
    (Join-Path $RepoRoot 'test\pbe_tabi_ses_row_worklist_real_benchmark.ps1')
) | Sort-Object -Unique
$compileArguments = @(
    '-std=c++23', '-O3', '-march=native', '-DNDEBUG', '-Wall', '-Wextra', '-Werror',
    ('-I' + (Join-Path $RepoRoot 'src')),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_mapping.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_row_worklist.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_torus_interval_bounds.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_torus_domain_bounds.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_ses_tiled_row_worklist.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    (Join-Path $RepoRoot 'test\titan_tabi_ses_row_worklist_real_benchmark.cpp'),
    '-o', $Exe
)
$fatalError = $hashLockError
$compile = $null
if ($null -eq $fatalError) {
    try {
        $compilerVersion = Invoke-CapturedProcess -FileName $Compiler -Arguments @('--version') `
            -StdoutPath $compilerVersionStdout -StderrPath $compilerVersionStderr -TimeoutMilliseconds 60000
        if ($compilerVersion.TimedOut -or $compilerVersion.ExitCode -ne 0) {
            $fatalError = "Unable to capture compiler version."
        } else {
            $compile = Invoke-CapturedProcess -FileName $Compiler -Arguments $compileArguments `
                -StdoutPath $compileStdout -StderrPath $compileStderr -TimeoutMilliseconds ($TimeoutMinutes * 60 * 1000)
            if ($compile.TimedOut -or $compile.ExitCode -ne 0) {
                $fatalError = "D2d compile failed or timed out; stdout=$compileStdout stderr=$compileStderr"
            }
        }
    } catch {
        $fatalError = $_.Exception.Message
        if (-not (Test-Path -LiteralPath $compileStdout)) {
            Set-Content -LiteralPath $compileStdout -Encoding UTF8 -Value 'compile did not start'
        }
        if (-not (Test-Path -LiteralPath $compileStderr)) {
            Set-Content -LiteralPath $compileStderr -Encoding UTF8 -Value $fatalError
        }
    }
} else {
    Set-Content -LiteralPath $compileStdout -Encoding UTF8 -Value 'compile skipped because input hash lock failed'
    Set-Content -LiteralPath $compileStderr -Encoding UTF8 -Value $fatalError
}

$rows = [Collections.Generic.List[object]]::new()
$events = [Collections.Generic.List[object]]::new()
function Invoke-D2dCase {
    param(
        [Parameter(Mandatory)][string]$CaseLabel,
        [Parameter(Mandatory)][string]$Fixture,
        [Parameter(Mandatory)][string]$InputHash,
        [Parameter(Mandatory)][int]$AtomCap,
        [Parameter(Mandatory)][string]$RadiiProvenance,
        [Parameter(Mandatory)][string]$ScopeLabel
    )
    $caseRoot = Join-Path $WorkRoot $CaseLabel
    New-Item -ItemType Directory -Path $caseRoot | Out-Null
    $rawPath = Join-Path $caseRoot 'benchmark.raw.csv'
    $stderrPath = Join-Path $caseRoot 'benchmark.stderr.txt'
    $parsedPath = Join-Path $caseRoot 'benchmark.parsed.csv'
    $arguments = @(
        $Fixture,
        $CaseLabel,
        '1.25',
        $AtomCap.ToString([Globalization.CultureInfo]::InvariantCulture),
        $RadiiProvenance,
        $InputHash,
        $Samples.ToString([Globalization.CultureInfo]::InvariantCulture),
        $AuditBudgetSeconds.ToString([Globalization.CultureInfo]::InvariantCulture)
    )
    $run = Invoke-CapturedProcess -FileName $Exe -Arguments $arguments `
        -StdoutPath $rawPath -StderrPath $stderrPath -TimeoutMilliseconds ($TimeoutMinutes * 60 * 1000)
    $parsed = $null
    if (-not $run.TimedOut -and (Test-Path -LiteralPath $rawPath)) {
        $lines = @(Get-Content -LiteralPath $rawPath)
        if ($lines.Count -ge 2) {
            try {
                $parsedRows = @($lines | ConvertFrom-Csv)
                if ($parsedRows.Count -eq 1) {
                    $parsed = $parsedRows[0]
                    $parsed | Add-Member -NotePropertyName runner_atom_cap -NotePropertyValue $AtomCap
                    $parsed | Add-Member -NotePropertyName runner_scope_label -NotePropertyValue $ScopeLabel
                    $parsed | Add-Member -NotePropertyName subprocess_wall_ms -NotePropertyValue $run.WallMs
                    $parsed | Add-Member -NotePropertyName peak_working_set_bytes -NotePropertyValue $run.PeakWorkingSetBytes
                    $parsed | Add-Member -NotePropertyName subprocess_exit_code -NotePropertyValue $run.ExitCode
                    $parsed | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $parsedPath
                    $rows.Add($parsed)
                }
            } catch {
                $parsed = $null
            }
        }
    }
    $events.Add([pscustomobject]@{
        case_label = $CaseLabel
        exit_code = $run.ExitCode
        timed_out = [int]$run.TimedOut
        subprocess_wall_ms = $run.WallMs
        peak_working_set_bytes = $run.PeakWorkingSetBytes
        parsed = [int]($null -ne $parsed)
        raw_path = $rawPath
        stderr_path = $stderrPath
    })
    [pscustomobject]@{ Row = $parsed; Run = $run }
}

$decision = 'not_started'
$scalingFactor = 5L
$projectedRssBytes = 'NA'
$projectedWallMs = 'NA'
$projectedPayloadBytes = 'NA'
$availableMemoryBytes = 0L
$headroomSource = 'unavailable_hard_guard_only'
if ($CompileOnly -and $null -eq $fatalError) {
    $decision = 'compile_only_pass'
}
if (-not $CompileOnly -and $null -eq $fatalError) {
    try {
        $oneA63Result = Invoke-D2dCase -CaseLabel '1a63_h1p25_full' -Fixture $oneA63 `
            -InputHash $oneA63Hash -AtomCap 0 -RadiiProvenance 'pqr_charmm_forcefield_radii' `
            -ScopeLabel 'full_1a63_charmm_broadphase_audit'
        $oneA63RssBytes = [int64]0
        $oneA63RssValid = $null -ne $oneA63Result.Row -and [int64]::TryParse(
            [string]$oneA63Result.Row.peak_working_set_bytes,
            [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture,
            [ref]$oneA63RssBytes) -and $oneA63RssBytes -gt 0L
        if ($null -eq $oneA63Result.Row -or $oneA63Result.Run.TimedOut -or
            $oneA63Result.Run.ExitCode -ne 0 -or $oneA63Result.Row.status -ne 'pass' -or
            -not $oneA63RssValid -or
            $oneA63Result.Row.audit_complete -ne '1' -or
            $oneA63Result.Row.false_negatives -ne '0' -or
            $oneA63Result.Row.d2b_d2c_csr_equal -ne '1' -or
            $oneA63Result.Row.d2c_nt_parity -ne '1' -or
            $oneA63Result.Row.d2c_reverse_parity -ne '1') {
            $fatalError = "1A63 D2d gate failed; raw=$($events[$events.Count - 1].raw_path) stderr=$($events[$events.Count - 1].stderr_path)"
        }
        if ($null -eq $fatalError -and $MaximumStage -eq 'oneA63') {
            $decision = 'stopped_after_1a63_requested'
        }
        if ($null -eq $fatalError -and $MaximumStage -ne 'oneA63') {
            $huge8192Result = Invoke-D2dCase -CaseLabel 'huge_prefix_8192_h1p25' -Fixture $huge `
                -InputHash $hugeHash -AtomCap 8192 -RadiiProvenance 'geometry_only_proxy_radii' `
                -ScopeLabel 'prefix_scaling_geometry_not_full_capsid'
            if ($null -ne $huge8192Result.Row -and
                $huge8192Result.Row.status -eq 'false_negative' -and
                $huge8192Result.Row.audit_complete -eq '1') {
                $fatalError = 'Huge 8192 full removed-candidate audit found a false negative.'
                $decision = 'stopped_after_8192_false_negative'
            } elseif ($null -eq $huge8192Result.Row -or $huge8192Result.Run.TimedOut -or
                $huge8192Result.Run.ExitCode -ne 0 -or $huge8192Result.Row.status -ne 'pass' -or
                $huge8192Result.Row.audit_complete -ne '1') {
                $fatalError = 'Huge 8192 stage failed, timed out, or returned a partial audit.'
            } else {
                $huge8192RssBytes = [int64]0
                $huge8192RssValid = [int64]::TryParse(
                    [string]$huge8192Result.Row.peak_working_set_bytes,
                    [Globalization.NumberStyles]::Integer,
                    [Globalization.CultureInfo]::InvariantCulture,
                    [ref]$huge8192RssBytes) -and $huge8192RssBytes -gt 0L
                if (-not $huge8192RssValid) {
                    $fatalError = 'Huge 8192 RSS is missing, nonnumeric, or nonpositive; refusing to open 32768.'
                    $decision = 'stopped_after_8192_missing_rss'
                }
                try {
                    $availableMemoryBytes = [int64](Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory * 1024L
                    $headroomSource = 'win32_operating_system_free_physical_memory'
                } catch {
                    $availableMemoryBytes = 0L
                    $headroomSource = 'unavailable_hard_guard_only'
                }
                $projectedWallMs = [double]$huge8192Result.Row.subprocess_wall_ms * $scalingFactor
                $projectedPayloadBytes = [int64]$huge8192Result.Row.d2c_accounted_peak_bytes * $scalingFactor
                $safeToOpen = $false
                if ($huge8192RssValid) {
                    $projectedRssBytes = $huge8192RssBytes * $scalingFactor
                    $headroomSafe = $availableMemoryBytes -eq 0L -or
                        $projectedRssBytes -lt [int64](0.70 * $availableMemoryBytes)
                    $safeToOpen = $projectedRssBytes -lt 6442450944L -and
                        $projectedWallMs -lt 1200000.0 -and
                        $projectedPayloadBytes -lt 2147483648L -and $headroomSafe
                }
                if ($null -ne $fatalError) {
                    $decision = 'stopped_after_8192_missing_rss'
                } elseif ($MaximumStage -eq 'huge8192') {
                    $decision = 'stopped_after_8192_requested'
                } elseif ($safeToOpen) {
                    $decision = 'opened_32768_after_projected_safety_guard'
                    $huge32768Result = Invoke-D2dCase -CaseLabel 'huge_prefix_32768_h1p25' -Fixture $huge `
                        -InputHash $hugeHash -AtomCap 32768 -RadiiProvenance 'geometry_only_proxy_radii' `
                        -ScopeLabel 'prefix_scaling_geometry_not_full_capsid'
                    if ($null -eq $huge32768Result.Row -or $huge32768Result.Run.TimedOut -or
                        $huge32768Result.Run.ExitCode -ne 0 -or $huge32768Result.Row.status -ne 'pass') {
                        $decision = '32768_stopped_failure_partial_or_timeout'
                    } else {
                        $decision = '32768_completed_full_removed_audit'
                    }
                } else {
                    $decision = 'stopped_after_8192_projected_safety_guard'
                }
            }
        }
    } catch {
        $fatalError = $_.Exception.Message
    }
}

$summaryPath = Join-Path $WorkRoot 'benchmark_summary.csv'
$eventsPath = Join-Path $WorkRoot 'run_events.csv'
if ($rows.Count -ne 0) {
    $rows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $summaryPath
} else {
    Set-Content -LiteralPath $summaryPath -Encoding UTF8 -Value 'no_parseable_benchmark_rows'
}
if ($events.Count -ne 0) {
    $events | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $eventsPath
} else {
    Set-Content -LiteralPath $eventsPath -Encoding UTF8 -Value 'no_benchmark_subprocess_started'
}
Set-Content -LiteralPath (Join-Path $WorkRoot 'staging_decision.txt') -Encoding UTF8 -Value @(
    "decision=$decision",
    "maximum_stage=$MaximumStage",
    "fatal_error=$fatalError",
    'huge_scope=prefix_scaling_geometry_not_full_capsid',
    'huge_radii=geometry_only_proxy_radii',
    'full_794880_atoms_run=not_run',
    "scaling_projection_factor=$scalingFactor",
    "projected_32768_rss_bytes=$projectedRssBytes",
    "projected_32768_wall_ms=$projectedWallMs",
    "projected_32768_accounted_payload_bytes=$projectedPayloadBytes",
    "available_memory_bytes_at_decision=$availableMemoryBytes",
    "headroom_source=$headroomSource",
    'projected_rss_hard_guard_bytes=6442450944',
    'projected_rss_available_memory_fraction=0.70',
    'projected_wall_guard_ms=1200000',
    'projected_accounted_payload_guard_bytes=2147483648',
    'benchmark_interpretation=same_working_tree_not_same_commit_ci'
)

$gitHead = (& git rev-parse HEAD 2>$null | Select-Object -First 1)
$gitStatus = @(& git status --porcelain=v1 2>$null)
Set-Content -LiteralPath (Join-Path $WorkRoot 'git_status_porcelain.txt') -Encoding UTF8 -Value $gitStatus
$sourceRows = [Collections.Generic.List[object]]::new()
foreach ($sourcePath in ($SourceInputs | Sort-Object)) {
    if (-not (Test-Path -LiteralPath $sourcePath)) { continue }
    $relative = $sourcePath.Substring($RepoRoot.Length).TrimStart([char]92, [char]47)
    $sourceRows.Add([pscustomobject]@{
        record_type = 'source'
        path = $relative
        sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $sourcePath).Hash
        value = ''
    })
}
$compilerVersionValue = if (Test-Path -LiteralPath $compilerVersionStdout) {
    ((Get-Content -LiteralPath $compilerVersionStdout | Select-Object -First 1) -join '')
} else { 'unavailable' }
$sourceRows.Add([pscustomobject]@{ record_type = 'metadata'; path = 'git_head'; sha256 = ''; value = $gitHead })
$sourceRows.Add([pscustomobject]@{ record_type = 'metadata'; path = 'git_dirty'; sha256 = ''; value = [int]($gitStatus.Count -ne 0) })
$sourceRows.Add([pscustomobject]@{ record_type = 'metadata'; path = 'compiler_path'; sha256 = ''; value = $Compiler })
$sourceRows.Add([pscustomobject]@{ record_type = 'metadata'; path = 'compiler_version'; sha256 = ''; value = $compilerVersionValue })
$oneA63Lock = if ($oneA63Hash -eq $expectedOneA63) { 'verified' } else { 'mismatch' }
$hugeLock = if ($hugeHash -eq $expectedHuge) { 'verified' } else { 'mismatch' }
$sourceRows.Add([pscustomobject]@{ record_type = 'input'; path = '1a63_full_charmm_ph7.pqr'; sha256 = $oneA63Hash; value = "expected=$expectedOneA63;lock_status=$oneA63Lock" })
$sourceRows.Add([pscustomobject]@{ record_type = 'input'; path = 'sphere_6co_full.pdb'; sha256 = $hugeHash; value = "expected=$expectedHuge;lock_status=$hugeLock;radii=geometry_only_proxy_radii" })
$compileFlags = '-std=c++23 -O3 -march=native -DNDEBUG -Wall -Wextra -Werror'
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'spacing'; sha256 = ''; value = '1.25' })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'samples'; sha256 = ''; value = $Samples })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'audit_budget_seconds'; sha256 = ''; value = $AuditBudgetSeconds })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'subprocess_timeout_minutes'; sha256 = ''; value = $TimeoutMinutes })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'compile_flags'; sha256 = ''; value = $compileFlags })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'scaling_projection_factor'; sha256 = ''; value = $scalingFactor })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'projected_rss_hard_guard_bytes'; sha256 = ''; value = '6442450944' })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'projected_wall_guard_ms'; sha256 = ''; value = '1200000' })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'projected_payload_guard_bytes'; sha256 = ''; value = '2147483648' })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'compile_only'; sha256 = ''; value = [int][bool]$CompileOnly })
$sourceRows.Add([pscustomobject]@{ record_type = 'config'; path = 'maximum_stage'; sha256 = ''; value = $MaximumStage })
$sourceRows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (Join-Path $WorkRoot 'source_manifest.csv')

$hashRows = [Collections.Generic.List[object]]::new()
foreach ($path in (Get-ChildItem -LiteralPath $WorkRoot -File -Recurse | Sort-Object FullName)) {
    if ($path.Name -eq 'artifact_hashes.csv') { continue }
    $hashRows.Add([pscustomobject]@{
        relative_path = $path.FullName.Substring($WorkRoot.Length).TrimStart([char]92, [char]47)
        sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $path.FullName).Hash
        bytes = $path.Length
    })
}
$hashRows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (Join-Path $WorkRoot 'artifact_hashes.csv')
$rows | Format-Table case_label,atoms,status,patch_build_ms,d1_candidates,d2b_candidates,d2c_candidates,removed_total,audit_complete,peak_working_set_bytes
Write-Host "decision=$decision"
Write-Host "summary=$(Join-Path $WorkRoot 'benchmark_summary.csv')"
Write-Host "events=$(Join-Path $WorkRoot 'run_events.csv')"
Write-Host "hashes=$(Join-Path $WorkRoot 'artifact_hashes.csv')"
if ($null -ne $fatalError) {
    throw $fatalError
}
