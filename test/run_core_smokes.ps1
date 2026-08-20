param(
    [string]$Compiler = 'g++'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$tests = @(
    'mol_topology_smoke.ps1',
    'pbctopo_lattice_math_smoke.ps1',
    'pbctopo_lattice_stability_smoke.ps1',
    'pbctopo_lattice_basis_invariance_smoke.ps1',
    'pbctopo_candidate_order_smoke.ps1',
    'pbctopo_candidate_scoring_smoke.ps1',
    'pbctopo_component_metric_smoke.ps1',
    'pbctopo_pair_graph_math_smoke.ps1',
    'pbctopo_sparse_pair_shift_table_smoke.ps1',
    'pbctopo_component_search_header_smoke.ps1',
    'pbctopo_component_search_optimality_smoke.ps1',
    'pbctopo_adaptive_image_shell_smoke.ps1',
    'pbctopo_density_cut_planes_smoke.ps1',
    'pbctopo_kinematic_seed_smoke.ps1',
    'pbctopo_virtual_discrete_search_smoke.ps1',
    'pbctopo_anchor_mst_exact_smoke.ps1',
    'pbctopo_certification_independence_smoke.ps1',
    'pbctopo_required_interface_contacts_smoke.ps1',
    'pbctopo_identifiability_smoke.ps1',
    'pbctopo_image_assignment_smoke.ps1',
    'pbctopo_assignment_set_report_smoke.ps1',
    'pbctopo_observation_schema_smoke.ps1',
    'pbctopo_evidence_semantics_smoke.ps1',
    'pbctopo_complete_uncertainty_band_policy_smoke.ps1',
    'pbctopo_temporal_viterbi_smoke.ps1',
    'pbctopo_temporal_assignment_replay_smoke.ps1',
    'pbctopo_reference_strain_smoke.ps1',
    'pbctopo_mathematical_objectives_smoke.ps1',
    'pbctopo_metadata_local_graph_smoke.ps1',
    'pbctopo_warm_state_smoke.ps1',
    'pbctopo_logging_smoke.ps1'
)

foreach ($test in $tests) {
    $path = Join-Path $PSScriptRoot $test
    Write-Host "==> $test"
    $command = Get-Command $path
    if ($command.Parameters.ContainsKey('Compiler')) {
        & $path -Compiler $Compiler
    } else {
        & $path
    }
    if ($LASTEXITCODE -ne 0) {
        throw "$test failed with exit code $LASTEXITCODE"
    }
}

Write-Host "PASS all $($tests.Count) TITAN-VIBE core smoke tests"
