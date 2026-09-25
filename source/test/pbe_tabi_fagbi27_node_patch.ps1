param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi_fagbi27_node_patch',
    [string]$FixtureRoot = 'test\pbe_inputs_fagbi27',
    [string]$OutputCsv = '',
    [string[]]$Fixtures = @(),
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 3600,
    [int]$MaxIterations = 720,
    [int]$Restart = 128,
    [double]$Tolerance = 1.0e-4,
    [double]$AuditTolerance = 2.0e-2,
    [double]$SesGridSpacing = 1.25,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunFixtures = @(),
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$Inv = [Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function Format-Number([double]$Value) {
    return $Value.ToString('G15', $Inv)
}

function Convert-ToDouble($Value, [double]$Fallback = [double]::NaN) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
        return $Fallback
    }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float, $Inv, [ref]$parsed)) {
        return $parsed
    }
    return $Fallback
}

function Convert-ToInt64($Value, [int64]$Fallback = 0) {
    $parsed = 0L
    if ([int64]::TryParse([string]$Value, [ref]$parsed)) { return $parsed }
    return $Fallback
}

function Test-Finite([double]$Value) {
    return -not [double]::IsNaN($Value) -and
        -not [double]::IsInfinity($Value)
}

function Get-ProfileValue($Profile, [string]$Name, $Default = '') {
    if ($null -eq $Profile) { return $Default }
    $property = $Profile.PSObject.Properties[$Name]
    if ($null -eq $property) { return $Default }
    return $property.Value
}

function Normalize-List([string[]]$Values) {
    return @($Values | ForEach-Object { @([string]$_ -split ',') } |
        ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

function Get-Sha256([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Get-PqrAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
            $_ -match '^(ATOM|HETATM)'
        }).Count
}

function Get-PreconditionerPath([string]$Effective) {
    if ($Effective -match 'regularized') { return 'regularized_ilu0' }
    if ($Effective -match 'ilu0') { return 'ilu0' }
    if ($Effective -match '^diagonal') { return 'diagonal_fallback' }
    if ([string]::IsNullOrWhiteSpace($Effective)) { return 'missing' }
    return $Effective
}

function Invoke-NodePatchTitanCase(
    [string]$Exe,
    [string]$RunDir,
    [string]$ArgPath,
    [int]$TimeoutSeconds
) {
    $psi = [Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $Exe
    $psi.WorkingDirectory = $RunDir
    $psi.Arguments = '-a "' + ($ArgPath -replace '"', '\"') + '"'
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $process = [Diagnostics.Process]::new()
    $process.StartInfo = $psi
    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    $started = $false
    $peakBytes = 0L
    $timedOut = $false
    try {
        [void]$process.Start()
        $started = $true
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        $deadline = if ($TimeoutSeconds -gt 0) {
            [DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
        } else {
            [DateTime]::MaxValue
        }
        while (-not $process.HasExited -and [DateTime]::UtcNow -lt $deadline) {
            try {
                $process.Refresh()
                $peakBytes = [math]::Max($peakBytes, [int64]$process.WorkingSet64)
                $peakBytes = [math]::Max($peakBytes, [int64]$process.PeakWorkingSet64)
            } catch {}
            Start-Sleep -Milliseconds 100
        }
        $timedOut = -not $process.HasExited
        if ($timedOut) {
            try { $process.Kill() } catch {}
            $process.WaitForExit()
        }
        try {
            $process.Refresh()
            $peakBytes = [math]::Max($peakBytes, [int64]$process.PeakWorkingSet64)
        } catch {}
        $stdoutTask.Wait()
        $stderrTask.Wait()
        Set-Content -LiteralPath (Join-Path $RunDir 'stdout.log') `
            -Value $stdoutTask.Result -Encoding UTF8
        Set-Content -LiteralPath (Join-Path $RunDir 'stderr.log') `
            -Value $stderrTask.Result -Encoding UTF8
        return [pscustomobject]@{
            exit_code = if ($timedOut) { 124 } else { $process.ExitCode }
            timed_out = $timedOut
            wall_seconds = $stopwatch.Elapsed.TotalSeconds
            peak_rss_mb = $peakBytes / 1MB
        }
    } finally {
        if ($started) {
            try {
                if (-not $process.HasExited) {
                    $process.Kill()
                    [void]$process.WaitForExit(5000)
                }
            } catch {}
        }
        $stopwatch.Stop()
        $process.Dispose()
    }
}

if ($Threads -lt 1) { throw 'Threads must be positive.' }
if ($Tolerance -le 0.0) { throw 'Tolerance must be positive.' }
if ($AuditTolerance -le 0.0) { throw 'AuditTolerance must be positive.' }

if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) {
    $FixtureRoot = Join-Path $RepoRoot $FixtureRoot
}
if ([string]::IsNullOrWhiteSpace($OutputCsv)) {
    $OutputCsv = Join-Path $WorkRoot 'fagbi27_node_patch_results.csv'
} elseif (-not [IO.Path]::IsPathRooted($OutputCsv)) {
    $OutputCsv = Join-Path $RepoRoot $OutputCsv
}

$legacyMatrices = @(
    (Join-Path $RepoRoot 'test\tmp\pbe_tabi_fagbi27_analytic_ses_apbs\fagbi27_analytic_ses_apbs_matrix.csv'),
    (Join-Path $RepoRoot 'test\tmp\pbe_tabi_fagbi27_recovery\fagbi27_recovery_matrix.csv')
)
$resolvedOutput = [IO.Path]::GetFullPath($OutputCsv)
if ($legacyMatrices | Where-Object { [IO.Path]::GetFullPath($_) -eq $resolvedOutput }) {
    throw 'OutputCsv points to a protected historical FAGBI-27 result file.'
}

$candidateExe = if ([IO.Path]::IsPathRooted($TitanExe)) {
    $TitanExe
} else {
    Join-Path $RepoRoot $TitanExe
}
$TitanExe = (Resolve-Path -LiteralPath $candidateExe).Path
$helpText = (& $TitanExe -h 2>&1 | Out-String)
if ($helpText -notmatch [regex]::Escape('-pbe-tabi-discretization=centroid|node-patch')) {
    throw "The selected titan.exe does not expose the accelerated node-patch option. Build the node-patch main worktree or pass its executable with -TitanExe. Executable: $TitanExe"
}

$manifestPath = Join-Path $FixtureRoot 'fagbi27_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "Missing FAGBI-27 fixture manifest: $manifestPath"
}
$manifestRows = @(Import-Csv -LiteralPath $manifestPath | Sort-Object {
        Convert-ToInt64 $_.paper_index
    })

function Resolve-Fagbi27Pqr($ManifestRow) {
    $pdbId = ([string]$ManifestRow.pdb_id).ToLowerInvariant()
    $pqrFile = [string]$ManifestRow.pqr_file
    $candidates = @(
        [pscustomobject]@{
            source = 'manifest'
            path = Join-Path $FixtureRoot $pqrFile
        },
        [pscustomobject]@{
            source = 'recovery_v2_clean'
            path = Join-Path $RepoRoot `
                "test\tmp\pbe_tabi_fagbi27_recovery_v2\${pdbId}_clean_fixture\${pdbId}_titan\$pqrFile"
        },
        [pscustomobject]@{
            source = 'recovery'
            path = Join-Path $RepoRoot `
                "test\tmp\pbe_tabi_fagbi27_recovery\${pdbId}_titan\$pqrFile"
        }
    )
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate.path -PathType Leaf) {
            return [pscustomobject]@{
                source = $candidate.source
                path = (Resolve-Path -LiteralPath $candidate.path).Path
                sha256 = Get-Sha256 $candidate.path
                atom_count = Get-PqrAtomCount $candidate.path
            }
        }
    }
    return $null
}

$FixtureAvailability = @{}
foreach ($manifestRow in $manifestRows) {
    if ([string]::IsNullOrWhiteSpace([string]$manifestRow.pdb_id)) { continue }
    $resolved = Resolve-Fagbi27Pqr $manifestRow
    if ($null -ne $resolved) {
        $FixtureAvailability[([string]$manifestRow.pdb_id).ToUpperInvariant()] = `
            $resolved
    }
}
$availableManifestRows = @($manifestRows | Where-Object {
        -not [string]::IsNullOrWhiteSpace([string]$_.pdb_id) -and
        $FixtureAvailability.ContainsKey(
            ([string]$_.pdb_id).ToUpperInvariant())
    })
if ($availableManifestRows.Count -lt 27 -and -not $Quick -and
    $Fixtures.Count -eq 0) {
    $missing = @($manifestRows | Where-Object {
            -not [string]::IsNullOrWhiteSpace([string]$_.pdb_id) -and
            -not $FixtureAvailability.ContainsKey(
                ([string]$_.pdb_id).ToUpperInvariant())
        } | ForEach-Object { $_.pdb_id })
    throw "Only $($availableManifestRows.Count)/27 FAGBI fixtures can be resolved. Missing: $($missing -join ', ')"
}

$Fixtures = @(Normalize-List $Fixtures | ForEach-Object { $_.ToUpperInvariant() })
$RerunFixtures = @(Normalize-List $RerunFixtures |
    ForEach-Object { $_.ToUpperInvariant() })
if ($Quick) { $Fixtures = @([string]$availableManifestRows[0].pdb_id) }
if ($Fixtures.Count -gt 0) {
    $selectedRows = @($availableManifestRows | Where-Object {
            ([string]$_.pdb_id).ToUpperInvariant() -in $Fixtures
        })
    $missingIds = @($Fixtures | Where-Object {
            $_ -notin @($selectedRows | ForEach-Object {
                    ([string]$_.pdb_id).ToUpperInvariant()
                })
        })
    if ($missingIds.Count -gt 0) {
        throw "Unknown or unavailable FAGBI-27 fixtures: $($missingIds -join ', ')"
    }
} else {
    $selectedRows = $availableManifestRows
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
New-Item -ItemType Directory -Path (Split-Path -Parent $OutputCsv) -Force | Out-Null

$TitanSha = Get-Sha256 $TitanExe
$Protocol = @(
    'fagbi27-node-patch-v1',
    'analytic-ses',
    "h=$(Format-Number $SesGridSpacing)",
    'analytical-fmm-m2l',
    'theta=0.45',
    'near=analytic_singular',
    'precond=fmm-nearfield-ilu0',
    'ilu-cap=64',
    'ilu-regularization=adaptive',
    "tol=$(Format-Number $Tolerance)",
    "restart=$Restart",
    "maxiter=$MaxIterations",
    'audit=sampled512',
    "threads=$Threads"
) -join '|'

$Checkpoint = @{}
if (Test-Path -LiteralPath $OutputCsv) {
    foreach ($row in @(Import-Csv -LiteralPath $OutputCsv)) {
        $Checkpoint[[string]$row.case_id] = $row
    }
}

function Save-Results {
    @($Checkpoint.Values | Sort-Object {
            Convert-ToInt64 $_.paper_index
        }) | Export-Csv -LiteralPath $OutputCsv -NoTypeInformation
}

function New-ResultRow($Fixture, $FixtureInfo, [string]$CaseId,
    [string]$PqrSha, [string]$RunDir) {
    return [pscustomobject][ordered]@{
        case_id = $CaseId
        paper_index = $Fixture.paper_index
        pdb_id = ([string]$Fixture.pdb_id).ToUpperInvariant()
        protocol = $Protocol
        executable_sha256 = $TitanSha
        pqr_sha256 = $PqrSha
        status = 'planned'
        failure = ''
        exit_code = ''
        timed_out = ''
        fixture_source = if ($null -ne $FixtureInfo) {
            $FixtureInfo.source
        } else { '' }
        fixture_pqr_path = if ($null -ne $FixtureInfo) {
            $FixtureInfo.path
        } else { '' }
        atom_count = if ($null -ne $FixtureInfo) {
            $FixtureInfo.atom_count
        } else { $Fixture.atom_count }
        surface_generator = ''
        discretization = ''
        panels = ''
        mesh_vertices = ''
        mesh_triangles = ''
        unknown_count = ''
        energy_kj_mol = ''
        energy_kcal_mol = ''
        iterations = ''
        true_residual = ''
        operator_audit_rows = ''
        operator_audit_rms = ''
        process_wall_seconds = ''
        solver_total_seconds = ''
        surface_seconds = ''
        preconditioner_setup_seconds = ''
        gmres_seconds = ''
        matvec_seconds = ''
        peak_rss_mb = ''
        preconditioner_requested = 'fmm-nearfield-ilu0'
        preconditioner_effective = ''
        preconditioner_path = ''
        ilu_regularization_mode = ''
        ilu_regularization_selected = ''
        ilu_regularization_attempts = ''
        ilu_regularization_rejected_attempts = ''
        ilu_regularization_shifted_rows = ''
        ilu_quality_status = ''
        ilu_factor_status = ''
        near_kernel_policy = ''
        fmm_mode = ''
        run_dir = $RunDir
    }
}

function Should-Run($OldRow, [string]$PqrSha, [string]$PdbId) {
    if ($Force -or $PdbId -in $RerunFixtures) { return $true }
    if ($null -eq $OldRow) { return $true }
    $sameContract = [string]$OldRow.protocol -eq $Protocol -and
        [string]$OldRow.executable_sha256 -eq $TitanSha -and
        [string]$OldRow.pqr_sha256 -eq $PqrSha
    if (-not $sameContract) { return $true }
    if ([string]$OldRow.status -in @('planned', 'running', 'incomplete')) {
        return $true
    }
    if ($RerunFailed -and [string]$OldRow.status -ne 'ok') { return $true }
    return $false
}

foreach ($fixture in $selectedRows) {
    $pdbId = ([string]$fixture.pdb_id).ToUpperInvariant()
    $caseId = $pdbId.ToLowerInvariant() + '_node_patch'
    $fixtureInfo = $FixtureAvailability[$pdbId]
    $pqrPath = if ($null -ne $fixtureInfo) { $fixtureInfo.path } else { '' }
    $runDir = Join-Path $WorkRoot $caseId
    if ([string]::IsNullOrWhiteSpace($pqrPath) -or
        -not (Test-Path -LiteralPath $pqrPath)) {
        $row = New-ResultRow $fixture $fixtureInfo $caseId '' $runDir
        $row.status = 'fixture_gap'
        $row.failure = "missing PQR: $pqrPath"
        $Checkpoint[$caseId] = $row
        Save-Results
        continue
    }
    $pqrSha = $fixtureInfo.sha256
    if ($fixtureInfo.source -eq 'manifest' -and
        -not [string]::IsNullOrWhiteSpace([string]$fixture.pqr_sha256) -and
        $pqrSha -ne ([string]$fixture.pqr_sha256).ToLowerInvariant()) {
        $row = New-ResultRow $fixture $fixtureInfo $caseId $pqrSha $runDir
        $row.status = 'fixture_gap'
        $row.failure = 'PQR checksum drift'
        $Checkpoint[$caseId] = $row
        Save-Results
        continue
    }
    $oldRow = if ($Checkpoint.ContainsKey($caseId)) {
        $Checkpoint[$caseId]
    } else { $null }
    if (-not (Should-Run $oldRow $pqrSha $pdbId)) {
        Write-Host "resume: $caseId"
        continue
    }

    if (Test-Path -LiteralPath $runDir) {
        Remove-Item -LiteralPath $runDir -Recurse -Force
    }
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $localPqr = Join-Path $runDir (Split-Path -Leaf $pqrPath)
    Copy-Item -LiteralPath $pqrPath -Destination $localPqr
    $row = New-ResultRow $fixture $fixtureInfo $caseId $pqrSha $runDir
    $row.status = 'running'
    $Checkpoint[$caseId] = $row
    Save-Results

    $argPath = Write-FmmCliArgFile -RepoRoot $RepoRoot -RunDir $runDir `
        -PqrPath $localPqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -FmmLocal off `
        -SurfacePoints 8 -Threads $Threads -MaxIter $MaxIterations `
        -Restart $Restart -Tolerance $Tolerance -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -ProbeRadius 1.4 -OpeningTheta 0.45 `
        -AuditMode sampled -AuditSampleCount 512 `
        -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -PreconditionerIluRegularization adaptive `
        -PreconditionerIluRegularizationMax 0.01 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -PrimarySurface analytic-ses -SesGridSpacing $SesGridSpacing `
        -DgSolv $true
    Add-Content -LiteralPath $argPath -Encoding ASCII `
        -Value '-pbe-tabi-discretization=node-patch'

    Write-Host "run: $caseId (node-patch, timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-NodePatchTitanCase -Exe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -TimeoutSeconds $CaseTimeoutSeconds
    $row.exit_code = $run.exit_code
    $row.timed_out = $run.timed_out
    $row.process_wall_seconds = $run.wall_seconds
    $row.peak_rss_mb = $run.peak_rss_mb
    if ($run.timed_out) {
        $row.status = 'timeout'
        $row.failure = 'titan_timeout'
    } elseif ($run.exit_code -ne 0) {
        $row.status = 'failed'
        $row.failure = 'titan_exit_nonzero'
    } else {
        try {
            $profile = @(Import-Csv -LiteralPath `
                    (Join-Path $runDir 'pbe_tabi_profile.csv'))[-1]
            $energy = @(Import-Csv -LiteralPath `
                    (Join-Path $runDir 'pbe_energy.csv'))[-1]
            $energyKj = Convert-ToDouble $energy.E_rxn_kJmol
            $trueResidual = Convert-ToDouble `
                (Get-ProfileValue $profile 'true_residual_rel')
            $auditRms = Convert-ToDouble `
                (Get-ProfileValue $profile 'tree_vs_direct_matvec_rel')
            $row.surface_generator = Get-ProfileValue $profile 'surface_generator'
            $row.discretization = Get-ProfileValue $profile 'discretization'
            $row.panels = Get-ProfileValue $profile 'panels'
            $row.mesh_vertices = Get-ProfileValue $profile 'mesh_vertices'
            $row.mesh_triangles = Get-ProfileValue $profile 'mesh_triangles'
            $row.unknown_count = Get-ProfileValue $profile 'unknown_count'
            $row.energy_kj_mol = $energyKj
            $row.energy_kcal_mol = $energyKj / 4.184
            $row.iterations = Get-ProfileValue $profile 'iterations'
            $row.true_residual = $trueResidual
            $row.operator_audit_rows = Get-ProfileValue $profile 'matvec_audit_rows'
            $row.operator_audit_rms = $auditRms
            $row.solver_total_seconds = Get-ProfileValue $profile 'total_seconds'
            $row.surface_seconds = Get-ProfileValue $profile 'surface_build_seconds'
            $row.preconditioner_setup_seconds = `
                Get-ProfileValue $profile 'preconditioner_setup_seconds'
            $row.gmres_seconds = Get-ProfileValue $profile 'gmres_seconds'
            $row.matvec_seconds = Get-ProfileValue $profile 'matvec_seconds'
            $row.preconditioner_effective = `
                Get-ProfileValue $profile 'preconditioner_effective'
            $row.preconditioner_path = `
                Get-PreconditionerPath $row.preconditioner_effective
            $row.ilu_regularization_mode = `
                Get-ProfileValue $profile 'preconditioner_ilu_regularization'
            $row.ilu_regularization_selected = Get-ProfileValue $profile `
                'preconditioner_ilu_regularization_selected'
            $row.ilu_regularization_attempts = Get-ProfileValue $profile `
                'preconditioner_ilu_regularization_attempts'
            $row.ilu_regularization_rejected_attempts = `
                Get-ProfileValue $profile `
                    'preconditioner_ilu_regularization_rejected_attempts'
            $row.ilu_regularization_shifted_rows = Get-ProfileValue $profile `
                'preconditioner_ilu_regularization_shifted_rows'
            $row.ilu_quality_status = `
                Get-ProfileValue $profile 'preconditioner_quality_status'
            $row.ilu_factor_status = `
                Get-ProfileValue $profile 'preconditioner_factor_status'
            $row.near_kernel_policy = `
                Get-ProfileValue $profile 'near_kernel_policy'
            $row.fmm_mode = Get-ProfileValue $profile 'fmm_m2l_mode'

            if ([string]$profile.status -ne 'ok') {
                $row.status = 'solver_gap'
                $row.failure = "profile_status=$($profile.status)"
            } elseif ($row.discretization -ne 'node_patch') {
                $row.status = 'contract_gap'
                $row.failure = "discretization=$($row.discretization)"
            } elseif (-not (Test-Finite $energyKj)) {
                $row.status = 'accuracy_gap'
                $row.failure = 'non_finite_energy'
            } elseif (-not (Test-Finite $trueResidual) -or
                $trueResidual -gt $Tolerance) {
                $row.status = 'accuracy_gap'
                $row.failure = "true_residual=$trueResidual"
            } elseif (-not (Test-Finite $auditRms) -or
                $auditRms -gt $AuditTolerance) {
                $row.status = 'accuracy_gap'
                $row.failure = "operator_audit_rms=$auditRms"
            } else {
                $row.status = 'ok'
                $row.failure = ''
            }
        } catch {
            $row.status = 'schema_gap'
            $row.failure = $_.Exception.Message
        }
    }
    $Checkpoint[$caseId] = $row
    Save-Results
}

$selectedIds = @($selectedRows | ForEach-Object {
        ([string]$_.pdb_id).ToUpperInvariant()
    })
$selectedResults = @($Checkpoint.Values | Where-Object {
        ([string]$_.pdb_id).ToUpperInvariant() -in $selectedIds
    })
$okRows = @($selectedResults | Where-Object { [string]$_.status -eq 'ok' })
$timeouts = @($selectedResults | Where-Object { [string]$_.status -eq 'timeout' })
$gaps = @($selectedResults | Where-Object {
        [string]$_.status -notin @('ok', 'timeout')
    })
$decision = if ($selectedResults.Count -lt $selectedRows.Count) {
    'fagbi27_node_patch_incomplete'
} elseif ($timeouts.Count -gt 0) {
    'fagbi27_node_patch_runtime_limited'
} elseif ($gaps.Count -gt 0) {
    'fagbi27_node_patch_gap'
} elseif ($okRows.Count -eq $selectedRows.Count) {
    if ($Quick) { 'fagbi27_node_patch_quick_pass' }
    else { 'fagbi27_node_patch_pass' }
} else {
    'fagbi27_node_patch_incomplete'
}

$regularized = @($okRows | Where-Object {
        [string]$_.preconditioner_path -eq 'regularized_ilu0'
    }).Count
$plainIlu = @($okRows | Where-Object {
        [string]$_.preconditioner_path -eq 'ilu0'
    }).Count
$diagonalFallback = @($okRows | Where-Object {
        [string]$_.preconditioner_path -eq 'diagonal_fallback'
    }).Count

Write-Host "Decision: $decision"
Write-Host "Completed: $($okRows.Count)/$($selectedRows.Count)"
Write-Host "Preconditioner paths: ILU0=$plainIlu regularized_ILU0=$regularized diagonal_fallback=$diagonalFallback"
Write-Host "Results: $OutputCsv"
Write-Host 'Resume: rerun unchanged; add -RerunFailed for timeout/failed rows, -RerunFixtures <PDB IDs>, or -Force for all rows.'

if (-not $AllowFailedDecision -and
    $decision -notin @('fagbi27_node_patch_pass',
        'fagbi27_node_patch_quick_pass')) {
    throw "FAGBI-27 node-patch gate did not pass: $decision"
}
