param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_cli_config_parser')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

function Invoke-PbctopoParserCase(
    [string]$Name,
    [string[]]$CaseArgs,
    [string]$ExpectedMessage
) {
    $stdout = Join-Path $WorkRoot "$Name.stdout.log"
    $stderr = Join-Path $WorkRoot "$Name.stderr.log"
    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = $TitanExe
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.Arguments = (($CaseArgs | ForEach-Object {
        '"' + $_.Replace('"', '\"') + '"'
    }) -join ' ')
    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    if (-not $process.Start()) {
        throw "$Name failed to start TITAN"
    }
    $stdoutText = $process.StandardOutput.ReadToEnd()
    $stderrText = $process.StandardError.ReadToEnd()
    $process.WaitForExit()
    $exitCode = $process.ExitCode
    [System.IO.File]::WriteAllText($stdout, $stdoutText)
    [System.IO.File]::WriteAllText($stderr, $stderrText)
    $log = $stdoutText + $stderrText
    $process.Dispose()
    if ($exitCode -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if ($log -notmatch $ExpectedMessage) {
        throw "$Name error was unclear: $log"
    }
}

Invoke-PbctopoParserCase 'bad_bond_cutoff' `
    @('-pbctopo', '-pbctopo-bond-cutoff', '-1') `
    'pbctopo-bond-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_clash_cutoff' `
    @('-pbctopo', '-pbctopo-clash-cutoff', '0') `
    'pbctopo-clash-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_contact_cutoff' `
    @('-pbctopo', '-pbctopo-contact-cutoff', '-0.5') `
    'pbctopo-contact-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_continuity_cutoff' `
    @('-pbctopo', '-pbctopo-continuity-cutoff', '0') `
    'pbctopo-continuity-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_exact_limit' `
    @('-pbctopo', '-pbctopo-exact-component-limit', '0') `
    'pbctopo-exact-component-limit.*> 0'
Invoke-PbctopoParserCase 'bad_exact_limit_too_high' `
    @('-pbctopo', '-pbctopo-exact-component-limit', '17') `
    'pbctopo-exact-component-limit.*<= 16'
Invoke-PbctopoParserCase 'bad_beam_width' `
    @('-pbctopo', '-pbctopo-beam-width', '0') `
    'pbctopo-beam-width.*> 0'
Invoke-PbctopoParserCase 'bad_beam_width_too_high' `
    @('-pbctopo', '-pbctopo-beam-width', '4097') `
    'pbctopo-beam-width.*<= 4096'
Invoke-PbctopoParserCase 'bad_anchor_backend' `
    @('-pbctopo', '-pbctopo-anchor-backend', 'bad') `
    'pbctopo-anchor-backend.*auto.*prim.*spatial-all'
Invoke-PbctopoParserCase 'bad_pbctopo_mode' `
    @('-pbctopo', '-pbctopo-mode', 'bad') `
    'pbctopo-mode.*vibe.*llps-local'
Invoke-PbctopoParserCase 'bad_pbctopo_search_policy' `
    @('-pbctopo', '-pbctopo-search-policy', 'bad') `
    'pbctopo-search-policy.*fast-certified.*audit-exact'
Invoke-PbctopoParserCase 'bad_audit_max_states' `
    @('-pbctopo', '-pbctopo-audit-max-states', '0') `
    'pbctopo-audit-max-states.*1\.\.1000000000'
Invoke-PbctopoParserCase 'bad_audit_max_seconds' `
    @('-pbctopo', '-pbctopo-audit-max-seconds', '0') `
    'pbctopo-audit-max-seconds.*\(0, 86400\]'
Invoke-PbctopoParserCase 'bad_audit_max_memory' `
    @('-pbctopo', '-pbctopo-audit-max-memory-mb', '0') `
    'pbctopo-audit-max-memory-mb.*1\.\.1048576'
Invoke-PbctopoParserCase 'bad_anchor_metrics' `
    @('-pbctopo', '-pbctopo-anchor-metrics', 'bad') `
    'pbctopo-anchor-metrics.*full.*mst-only.*no-order.*no-contact'
Invoke-PbctopoParserCase 'bad_temporal_mode' `
    @('-pbctopo', '-pbctopo-temporal', 'bad') `
    'pbctopo-temporal.*greedy.*viterbi'
Invoke-PbctopoParserCase 'bad_temporal_k' `
    @('-pbctopo', '-pbctopo-temporal-k', '33') `
    'pbctopo-temporal-k.*1\.\.32'
Invoke-PbctopoParserCase 'bad_temporal_penalty' `
    @('-pbctopo', '-pbctopo-temporal-penalty', '-1') `
    'pbctopo-temporal-penalty.*>= 0'
Invoke-PbctopoParserCase 'bad_vdw_scale' `
    @('-pbctopo', '-pbctopo-vdw-clash-scale', '1.2') `
    'pbctopo-vdw-clash-scale.*\(0, 1\]'
Invoke-PbctopoParserCase 'missing_pbctopo_mode' `
    @('-pbctopo', '-pbctopo-mode') `
    'pbctopo-mode.*requires .*value'
Invoke-PbctopoParserCase 'missing_pbctopo_search_policy' `
    @('-pbctopo', '-pbctopo-search-policy') `
    'pbctopo-search-policy.*requires .*value'
Invoke-PbctopoParserCase 'missing_bond_cutoff' `
    @('-pbctopo', '-pbctopo-bond-cutoff') `
    'pbctopo-bond-cutoff.*requires .*value'
Invoke-PbctopoParserCase 'missing_required_contacts' `
    @('-pbctopo', '-pbctopo-required-contacts') `
    'pbctopo-required-contacts.*requires .*value'
Invoke-PbctopoParserCase 'required_contacts_llps_conflict' `
    @('-pbctopo-llps-local', '-pbctopo-required-contacts', 'contacts.dat') `
    'pbctopo-required-contacts.*only valid in vibe mode'
Invoke-PbctopoParserCase 'audit_exact_llps_conflict' `
    @('-pbctopo-llps-local', '-pbctopo-search-policy', 'audit-exact') `
    'audit-exact.*only valid in vibe mode'

Write-Host 'PASS pbctopo_cli_config_parser_smoke'
