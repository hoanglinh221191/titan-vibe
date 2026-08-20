param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_logging_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_logging_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_logging_smoke.exe'

@'
#include "pbctopo/logging.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main() {
  std::filesystem::remove("pbctopo_debug.log");

  titan_pbctopo::reset_log_file();
  titan_pbctopo::append_log_line("first chunk");

  titan_pbctopo::reset_log_file();
  titan_pbctopo::append_log_line("second chunk");

  std::ifstream in("pbctopo_debug.log");
  std::ostringstream buffer;
  buffer << in.rdbuf();
  const std::string contents = buffer.str();

  if (contents.find("first chunk") == std::string::npos ||
      contents.find("second chunk") == std::string::npos) {
    std::cerr << "FAIL pbctopo_logging_smoke: reset_log_file truncated chunk diagnostics\n";
    std::cerr << contents;
    return 1;
  }

  std::cout << "PASS pbctopo_logging_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\pbctopo\logging.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile pbctopo logging smoke'
}

Push-Location $WorkRoot
try {
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "pbctopo logging smoke failed with exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}

$solver = Join-Path $repo 'src\solv_pbc_pbctopo.cpp'
$solverText = Get-Content -LiteralPath $solver -Raw
$requiredTimingKeys = @(
    'pbctopo_kinematic_seed_candidates',
    'pbctopo_kinematic_seed_evaluated',
    'pbctopo_kinematic_seed_deduped',
    'pbctopo_kinematic_seed_best',
    'pbctopo_kinematic_seed_history_resets',
    'pbctopo_density_cut_planes_frames',
    'pbctopo_density_cut_planes_x',
    'pbctopo_density_cut_planes_y',
    'pbctopo_density_cut_planes_z',
    'pbctopo_density_cut_planes_dense_equiv_x',
    'pbctopo_density_cut_planes_dense_equiv_y',
    'pbctopo_density_cut_planes_dense_equiv_z',
    'pbctopo_density_cut_plane_collisions',
    'pbctopo_density_cut_plane_dense_fallbacks',
    'pbctopo_density_cut_plane_build',
    'pbctopo_virtual_discrete_frames',
    'pbctopo_virtual_axis_states_x',
    'pbctopo_virtual_axis_states_y',
    'pbctopo_virtual_axis_states_z',
    'pbctopo_virtual_state_space',
    'pbctopo_virtual_states_evaluated',
    'pbctopo_virtual_search_exhaustive_frames',
    'pbctopo_virtual_seam_enumeration_exhaustive_frames',
    'pbctopo_virtual_projection_scoring_exhaustive_frames',
    'pbctopo_virtual_full_atom_verification_exhaustive_frames',
    'pbctopo_virtual_search_shortlist_frames',
    'pbctopo_virtual_branch_and_bound_frames',
    'pbctopo_virtual_sweep_workspace_reuse_frames',
    'pbctopo_virtual_sweep_events_x',
    'pbctopo_virtual_sweep_events_y',
    'pbctopo_virtual_sweep_events_z',
    'pbctopo_virtual_sweep_workspace_peak_bytes',
    'pbctopo_verify_geometry_buffer_reuses',
    'pbctopo_frame_point_buffer_peak_bytes_estimate',
    'pbctopo_virtual_omitted_seam_lower_bound_min',
    'pbctopo_virtual_lower_bound_gap_min',
    'pbctopo_virtual_seam_sweep',
    'pbctopo_certification_attempts',
    'pbctopo_certification_successes',
    'pbctopo_certification_assignment_rejects',
    'pbctopo_certification_hard_rejects',
    'pbctopo_output_certification_band_precheck_rejects',
    'pbctopo_output_certification_band_rebuilds',
    'pbctopo_output_certification_band_offset_bytes_peak',
    'pbctopo_output_certification_band_point_bytes_avoided_peak',
    'pbctopo_heavy_search_frames',
    'pbctopo_heavy_search_all_atom_candidate_rejects',
    'pbctopo_heavy_search_top_candidate_rejects',
    'pbctopo_heavy_search_all_atom_ranking_recoveries',
    'pbctopo_heavy_search_all_atom_band_exhausted',
    'pbctopo_runtime_simd_forced_scalar',
    'pbctopo_certification_topology_edges_checked',
    'pbctopo_certification_topology_edge_residuals',
    'pbctopo_certification_topology_cycle_residuals',
    'pbctopo_warm_state_commits',
    'pbctopo_warm_state_resets',
    'pbctopo_warm_state_commit_rejects',
    'pbctopo_warm_projection_rejects',
    'pbctopo_warm_component_rejects',
    'pbctopo_weak_observation_output_frames',
    'pbctopo_weak_observation_warm_preserved',
    'pbctopo_contradicted_output_frames',
    'pbctopo_selected_evidence_audits',
    'pbctopo_selected_evidence_contacts_checked',
    'pbctopo_selected_evidence_contacts_lost',
    'pbctopo_eobs_construction_attempts',
    'pbctopo_eobs_construction_completed',
    'pbctopo_eobs_construction_failures',
    'pbctopo_eobs_cell_setup_failures',
    'pbctopo_eobs_atom_index_failures',
    'pbctopo_eobs_mic_query_failures',
    'pbctopo_selected_evidence_audit_s',
    'pbctopo_warm_layout_mismatches',
    'pbctopo_certification_atom_image_residuals',
    'pbctopo_certification_assignment_point_residuals',
    'pbctopo_certification_steric_pair_checks',
    'pbctopo_certification_steric_clashes',
    'pbctopo_certification_required_contacts_checked',
    'pbctopo_certification_required_contacts_lost',
    'pbctopo_required_interface_contacts_loaded',
    'pbctopo_certification_explicit_interface_contacts_checked',
    'pbctopo_certification_explicit_interface_contacts_lost',
    'component_required_contact_constraint_builds',
    'component_required_contact_constraint_pairs',
    'component_required_contact_relative_deltas_checked',
    'component_required_contact_relative_deltas_rejected',
    'component_required_contact_infeasible',
    'pbctopo_geometry_classifications',
    'pbctopo_compact_unit_box_certified_frames',
    'pbctopo_noncompact_certified_frames',
    'pbctopo_winding_rejected_candidates',
    'pbctopo_extended_winding_rejected_candidates',
    'pbctopo_metadata_inter_hard_edges',
    'pbctopo_metadata_inter_soft_edges',
    'pbctopo_owner_com_split_observations',
    'pbctopo_anchor_intra_owner_contact_pairs',
    'pbctopo_soft_observed_contact_pairs',
    'pbctopo_hard_required_owner_pairs',
    'pbctopo_diagnostic_only_candidates',
    'pbctopo_virtual_geometry_outlier_observations',
    'component_graph_outlier_observations',
    'component_graph_primary_attempts',
    'component_graph_primary_certified',
    'component_graph_primary_virtual_skips',
    'component_graph_primary_fallbacks',
    'component_graph_primary_unidentified',
    'component_graph_unidentified_output_rejects',
    'component_graph_spatial_shortlist_builds',
    'component_graph_spatial_shortlist_fallbacks',
    'component_graph_spatial_shortlist_pairs',
    'component_graph_spatial_shortlist_shift_scores_avoided',
    'component_image_shell_expansions',
    'component_image_shell_retry_successes',
    'component_image_shell_max_radius',
    'component_image_shell_max_shift_count',
    'component_large_mode_internal_offsets_preserved',
    'component_large_mode_internal_offset_nonzero',
    'component_search_exact_budget_skips',
    'component_audit_exact_incomplete',
    'component_audit_exact_state_limit',
    'component_audit_exact_wall_time_limit',
    'component_audit_exact_memory_limit',
    'component_search_exact_state_upper_bound',
    'component_search_filtered_state_upper_bound',
    'component_search_arc_consistency_passes',
    'component_search_arc_consistency_revisions',
    'component_search_arc_consistency_removals',
    'component_search_arc_consistency_failures',
    'component_search_exact_recoveries',
    'component_search_arc_consistency_time',
    'component_pair_table_relative_deltas',
    'component_pair_table_relative_reuses',
    'component_search_shift_descriptors',
    'component_search_shift_materialized',
    'component_search_shift_domain_max',
    'component_search_lazy_materialize_time',
    'component_search_outlier_observations',
    'component_large_mode_beam_final_rescore_all_frames',
    'component_search_beam_provisional_children',
    'component_search_beam_state_materializations',
    'component_search_beam_state_copies_avoided',
    'component_search_beam_parent_state_moves',
    'component_search_certification_survivors',
    'component_search_certification_band_retained',
    'component_graph_refine_anchor_spatial_expansions',
    'component_graph_refine_anchor_spatial_exact_certified',
    'component_graph_refine_anchor_spatial_uncertified',
    'component_graph_refine_anchor_spatial_max_k',
    'pbctopo_lattice_audits',
    'pbctopo_lattice_invalid_frames',
    'pbctopo_lattice_reduction_fallback_frames',
    'pbctopo_lattice_fatal_query_failures',
    'pbctopo_lattice_optional_query_failures',
    'pbctopo_lattice_condition_max',
    'temporal_diagnostic_only_frames',
    'temporal_proposed_changes',
    'temporal_certified_candidates',
    'temporal_assignment_deduped',
    'temporal_replay_frames',
    'temporal_replay_applied',
    'temporal_replay_skipped'
)
foreach ($key in $requiredTimingKeys) {
    if (-not $solverText.Contains($key)) {
        throw "Missing PBCTOPO kinematic timing key: $key"
    }
}
