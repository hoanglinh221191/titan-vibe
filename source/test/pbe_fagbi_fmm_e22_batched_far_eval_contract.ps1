param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Require-Text {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][string]$Pattern,
        [Parameter(Mandatory = $true)][string]$Message
    )
    $full = Join-Path $RepoRoot $Path
    if (-not (Test-Path $full)) {
        throw "missing required file: $Path"
    }
    $content = Get-Content -Raw $full
    if ($content -notmatch $Pattern) {
        throw $Message
    }
}

Require-Text -Path "src/tabi/tabi_fmm_kifmm_real_panel_cloud.h" `
    -Pattern "full_apply_policy" `
    -Message "KifMM real-panel config must expose full_apply_policy."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.h" `
    -Pattern "kifmm_full_apply_policy" `
    -Message "Full-panel matvec result must report kifmm_full_apply_policy."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.h" `
    -Pattern "kifmm_target_tasks" `
    -Message "Full-panel matvec result must report target-node task count."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.h" `
    -Pattern "kifmm_channel_build_seconds" `
    -Message "Full-panel matvec result must split unknown-dependent channel build time."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.h" `
    -Pattern "TitanExecutor" `
    -Message "Prepared full-panel apply must expose an executor-aware overload."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.cpp" `
    -Pattern "target-node-parallel" `
    -Message "Full-panel matvec implementation must include target-node-parallel policy."

Require-Text -Path "src/tabi/tabi_fmm_kifmm_full_panel_matvec.cpp" `
    -Pattern "noncontiguous_target_node_worklist" `
    -Message "Target-node parallel apply must reject noncontiguous target-node worklists before parallel writes."

Require-Text -Path "src/arguments.h" `
    -Pattern "tabi_kifmm_full_apply" `
    -Message "CLI argument state must include tabi_kifmm_full_apply."

Require-Text -Path "src/arguments.cpp" `
    -Pattern "pbe-tabi-kifmm-full-apply" `
    -Message "Arguments parser must accept -pbe-tabi-kifmm-full-apply."

Require-Text -Path "test/titan_fmm_kifmm_full_panel_matvec_smoke.cpp" `
    -Pattern "e22-profile" `
    -Message "Full-panel smoke driver must expose --e22-profile."

Require-Text -Path "test/pbe_fagbi_fmm_e22_batched_far_eval_gate.ps1" `
    -Pattern "fmm_e22_batched_far_eval" `
    -Message "E22 gate harness is missing."

Require-Text -Path "test/pbe_fagbi_fmm_e22_batched_far_eval_gate.ps1" `
    -Pattern "Assert-E22CanonicalPair" `
    -Message "E22 gate must reject stale duplicate canonical KIFMM source/header prefixes."

Require-Text -Path "test/pbe_fagbi_fmm_e22_batched_far_eval_gate.ps1" `
    -Pattern "fmm_e22_target_node_parallel_apply_quick_pass" `
    -Message "E22 closeout decision must use the narrower target-node-parallel quick-pass label."

Require-Text -Path "TITAN_FMM_ENGINE_E22_BATCHED_FAR_EVAL_REPORT.txt" `
    -Pattern "FMM-E22" `
    -Message "E22 root report is missing."

Require-Text -Path "TITAN_FMM_ENGINE_E22_BATCHED_FAR_EVAL_REPORT.txt" `
    -Pattern "worklist_hash must match" `
    -Message "E22 report must state the serial/parallel worklist-hash invariant."

Require-Text -Path "TITAN_FMM_ENGINE_E22_BATCHED_FAR_EVAL_REPORT.txt" `
    -Pattern "policy_hash is allowed to differ" `
    -Message "E22 report must explain why policy hash differs between serial and target-node-parallel rows."

Require-Text -Path "TITAN_FMM_ENGINE_E22_BATCHED_FAR_EVAL_REPORT.txt" `
    -Pattern "sublinear but useful" `
    -Message "E22 report must avoid overclaiming near-linear scaling."

Write-Host "FMM-E22 contract ok"
