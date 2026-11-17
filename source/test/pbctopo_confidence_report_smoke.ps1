param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_confidence_report_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
    $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$pdb = Join-Path $WorkRoot 'confidence_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO confidence report smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA     1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY     2      12.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA     1      50.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY     2      52.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      5  CA  DUM C   1      25.000  25.000  25.000  1.00  0.00      DUM  C
ENDMDL
TITLE     Title t= 1.1234567890123457 step= 1
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA     1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY     2      12.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA     1      50.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY     2      52.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      5  CA  DUM C   1      25.000  25.000  25.000  1.00  0.00      DUM  C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ PROT ]
1 2 3 4

[ DUMMY ]
5
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
MOLA 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
MOLB 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
DUM 3

[ atoms ]
1 CT 1 DUM CA 1 0.0 12.0

[ system ]
pbctopo confidence report smoke

[ molecules ]
MOLA 1
MOLB 1
DUM 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f=$pdb
-pbctopo=true
-pbctopo-bond-cutoff=1.10
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=64
-pbctopo-search-policy=audit-exact
-pbctopo-audit-max-states=1
-pbctopo-anchor-backend=auto
-g=0-1
-nt=1
-1pdb=true
-nocal=true
-n=$ndx
-common_topo=$top
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=confidence_report_out.pdb
"@
[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a $arg
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

$report = Join-Path $WorkRoot 'pbctopo_frame_report.csv'
if (-not (Test-Path -LiteralPath $report)) {
    throw "Missing frame report: $report"
}

function Assert-CsvFieldParity {
    param([Parameter(Mandatory = $true)][string]$Path)

    $lines = @(Get-Content -LiteralPath $Path)
    if ($lines.Count -lt 2) {
        throw "CSV has no data rows: $Path"
    }
    $expected = $lines[0].Split(',').Count
    for ($lineIndex = 1; $lineIndex -lt $lines.Count; ++$lineIndex) {
        if ([string]::IsNullOrWhiteSpace($lines[$lineIndex])) {
            continue
        }
        $actual = $lines[$lineIndex].Split(',').Count
        if ($actual -ne $expected) {
            throw "CSV field mismatch in $Path at line $($lineIndex + 1): expected $expected, found $actual"
        }
    }
}

Assert-CsvFieldParity -Path $report
$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -ne 2) {
    throw "Expected exactly 2 frame report rows, found $($rows.Count)"
}
if ($rows[1].time_ps -ne '1.1234567890123457') {
    throw "Frame-report time_ps did not round-trip: $($rows[1].time_ps)"
}

$requiredColumns = @(
    'schema_version',
    'producer_version',
    'coordinate_unit',
    'time_unit',
    'endianness',
    'identity_hash_algorithm',
    'binding_hash_algorithm',
    'payload_checksum_algorithm',
    'observation_scope',
    'temporal_policy_status',
    'frame',
    'time_ps',
    'trajectory_frame_index',
    'frame_binding_valid',
    'box_b00',
    'box_b01',
    'box_b02',
    'box_b10',
    'box_b11',
    'box_b12',
    'box_b20',
    'box_b21',
    'box_b22',
    'box_hash',
    'atom_selection_hash',
    'canonical_atom_universe_hash',
    'wrapped_coordinate_hash',
    'status',
    'certificate_scope',
    'certificate_graph_source',
    'hard_graph_source',
    'search_evidence_source',
    'hard_feasible',
    'evidence_state',
    'evidence_consistency_evaluated',
    'evidence_consistent',
    'evidence_graph_connected',
    'soft_observed_contact_pair_hash',
    'soft_observed_lost_pair_hash',
    'soft_observed_contact_hash_algorithm',
    'soft_observed_contact_hash_domain',
    'soft_observed_contact_cutoff_policy',
    'soft_observed_contact_audit_policy_version',
    'soft_observed_contact_max_pairs_per_unit_pair',
    'soft_observed_contact_selection_tie_tolerance_d2_A2',
    'soft_observed_contact_loss_tolerance_d2_A2',
    'soft_observed_contact_cutoff_A',
    'soft_observed_contact_construction_attempted',
    'soft_observed_contact_construction_complete',
    'soft_observed_contact_cell_setup_failures',
    'soft_observed_contact_atom_index_failures',
    'soft_observed_contact_mic_query_failures',
    'identified',
    'unique_within_search_domain',
    'uniqueness_search_exhaustive',
    'objective_uniqueness_known',
    'objective_unique_within_search_domain',
    'evidence_uniqueness_known',
    'evidence_unique_within_search_domain',
    'feasible_assignment_set_enumerated',
    'feasible_assignment_set_semantics',
    'feasible_assignments_within_domain',
    'evidence_component_count',
    'evidence_edge_count',
    'explicit_contact_evidence_edges',
    'scored_contact_evidence_edges',
    'evidence_connected_component_count',
    'unconstrained_relative_components',
    'relative_lattice_degrees_of_freedom',
    'equivalent_best_assignments_within_domain',
    'certificate_atoms_total',
    'hard_graph_atoms_covered',
    'hard_graph_atom_fraction',
    'topology_bond_atoms_covered',
    'topology_bond_atom_fraction',
    'metadata_only_atoms',
    'metadata_only_atom_fraction',
    'mode',
    'component_count',
    'soft_valid_candidate_count',
    'nonclashing_candidate_count',
    'unique_certified_assignment_count',
    'committed_candidate_count',
    'selected_score_vector',
    'pool_best_score_vector',
    'pool_second_score_vector',
    'best_score_vector',
    'second_score_vector',
    'score_gap_metric',
    'score_gap',
    'clash_margin_A',
    'continuity_margin_A',
    'continuity_margin_atom_a',
    'continuity_margin_atom_b',
    'anchor_max_edge_A',
    'split_owner_gap_scaled',
    'split_owner_pairs',
    'geometry_class',
    'current_basis_fundamental_parallelepiped_certified',
    'compact_unit_box_certified',
    'scaled_span_x',
    'scaled_span_y',
    'scaled_span_z',
    'fits_current_basis_fundamental_parallelepiped',
    'fits_one_fundamental_cell',
    'cartesian_extent_x_A',
    'cartesian_extent_y_A',
    'cartesian_extent_z_A',
    'topological_winding_rank',
    'topological_winding_b0_x',
    'topological_winding_b0_y',
    'topological_winding_b0_z',
    'topological_winding_b1_x',
    'topological_winding_b1_y',
    'topological_winding_b1_z',
    'topological_winding_b2_x',
    'topological_winding_b2_y',
    'topological_winding_b2_z',
    'hard_evaluated',
    'hard_certified',
    'hard_connectivity_edges_available',
    'topology_edges_available',
    'metadata_edges_available',
    'hard_topology_edges_checked',
    'hard_topology_edge_fraction',
    'hard_ambiguous_topology_edges',
    'hard_metadata_edges_checked',
    'hard_metadata_edge_fraction',
    'hard_ambiguous_metadata_edges',
    'hard_metadata_edge_image_residuals',
    'hard_metadata_edge_distance_residuals',
    'hard_metadata_edge_max_distance_error_A',
    'hard_edge_residuals',
    'hard_atom_image_residuals',
        'hard_atom_image_max_error',
        'hard_assignment_point_residuals',
        'hard_assignment_point_max_error',
    'hard_steric_audit_available',
    'hard_steric_scan_complete',
        'hard_steric_pairs_checked',
        'steric_radius_force_field_atoms',
        'steric_radius_explicit_zero_lj_atoms',
        'steric_zero_lj_policy',
        'steric_radius_inferred_atoms',
        'steric_radius_generic_fallback_atoms',
        'steric_radius_unavailable_atoms',
        'hard_min_steric_margin_A',
        'soft_observed_contacts_checked',
        'soft_observed_contacts_lost',
        'soft_observed_contact_lost_fraction',
    'hard_required_contacts_checked',
    'hard_topology_cycle_residuals',
    'hard_unsupported_winding_rank',
    'hard_clashes',
    'hard_metadata_internal_pairs_checked',
    'hard_metadata_internal_clashes',
    'hard_required_contacts_lost',
    'hard_explicit_interface_contacts_checked',
    'hard_explicit_interface_contacts_lost',
    'hard_explicit_interface_weight_checked',
    'hard_explicit_interface_weight_lost',
    'topology_epoch_id',
    'component_layout_signature',
    'assignment_signature',
    'assignment_component_count',
    'virtual_search_exhaustive',
    'virtual_seam_enumeration_exhaustive',
    'virtual_projection_scoring_exhaustive',
    'virtual_full_atom_verification_exhaustive',
    'output_certificate_search_exhaustive',
    'virtual_branch_and_bound_used',
    'virtual_state_space',
    'virtual_states_evaluated',
    'virtual_omitted_seam_lower_bound',
    'virtual_lower_bound_gap',
    'component_search_attempted',
    'component_search_completed_under_admissible_bounds',
    'component_search_incumbent_bound_pruning_used',
    'component_search_best_objective_proven',
    'component_search_equivalent_optima_complete',
    'component_feasible_assignment_set_enumerated',
    'component_feasible_assignments_within_domain',
    'component_soft_score_valid_assignments_within_domain',
    'component_bounded_domain_feasibility_known',
    'component_bounded_domain_feasible',
    'component_search_shell_radius',
    'selected_assignment_boundary_known',
    'selected_assignment_touches_shell_boundary',
    'domain_expansion_attempted',
    'domain_expansion_exhausted',
    'bounded_domain_only',
    'component_search_exhaustive',
    'beam_survivors_exactly_rescored',
    'component_offset_enumeration_exhaustive',
    'component_all_atom_feasibility_exhaustive',
    'component_certificate_search_exhaustive',
    'component_fixed_certificate_audit_evaluated',
    'component_fixed_certificate_input_valid',
    'component_fixed_certificate_radius_model_valid',
    'component_fixed_certificate_internal_steric_scan_complete',
    'component_fixed_certificate_internal_pairs_checked',
    'component_fixed_certificate_internal_clashes',
    'component_fixed_certificate_required_contact_scan_complete',
    'component_fixed_certificate_required_contacts_checked',
    'component_fixed_certificate_required_contacts_lost',
    'component_fixed_certificate_passed',
    'component_audit_exact_requested',
    'component_audit_exact_completed',
    'component_audit_exact_incomplete',
    'component_audit_exact_incomplete_reason',
    'component_audit_exact_states_visited',
    'component_audit_exact_wall_seconds',
    'component_audit_exact_memory_bytes_estimate',
    'component_audit_exact_memory_budget_semantics',
    'component_audit_exact_reference_comparison_available',
    'component_audit_exact_reference_assignment_match',
    'component_search_optimality_gap_known',
    'component_search_optimality_gap',
    'fallback_reason',
    'time_total_ms',
    'time_graph_ms',
    'time_virtual_ms',
    'time_component_ms',
    'soft_pair_cycle_support',
    'soft_pair_cycle_mismatch',
    'soft_pair_inconsistent_edges'
)

$header = (Get-Content -LiteralPath $report -First 1).Split(',')
foreach ($column in $requiredColumns) {
    if ($header -notcontains $column) {
        throw "Missing frame report column: $column"
    }
}

$certificateScopeIndex = [Array]::IndexOf($header, 'certificate_scope')
if ($certificateScopeIndex -lt 0 -or
    $certificateScopeIndex + 1 -ge $header.Count -or
    $header[$certificateScopeIndex + 1] -ne 'certificate_graph_source') {
    throw 'certificate_graph_source is not immediately after certificate_scope'
}

foreach ($row in $rows) {
    if ($row.schema_version -ne '8' -or
        [string]::IsNullOrWhiteSpace($row.producer_version) -or
        $row.coordinate_unit -ne 'angstrom' -or $row.time_unit -ne 'ps' -or
        $row.endianness -ne 'not_applicable_text' -or
        $row.identity_hash_algorithm -ne 'sha256-128' -or
        $row.binding_hash_algorithm -ne 'sha256-128' -or
        $row.payload_checksum_algorithm -ne 'sha256-256' -or
        $row.canonical_atom_universe_hash -notmatch '^[0-9a-fA-F]{32}$' -or
        $row.observation_scope -ne 'framewise_bounded_domain' -or
        $row.temporal_policy_status -ne 'production_framewise_greedy') {
        throw 'Frame report observation schema metadata is invalid'
    }
    if ($row.evidence_state -notmatch '^(not_evaluated|not_applicable|supported_connected|weak_disconnected|weak_ambiguous|contradicted)$' -or
        $row.soft_observed_contact_hash_algorithm -ne 'sha256-128' -or
        $row.soft_observed_contact_hash_domain -ne 'titan-pbctopo-eobs-pairs-v1' -or
        $row.soft_observed_contact_cutoff_policy -ne 'wrapped-mic-heavy-closest8-per-owner-fragment-pair-v1' -or
        $row.soft_observed_contact_audit_policy_version -ne '1' -or
        $row.soft_observed_contact_max_pairs_per_unit_pair -ne '8' -or
        $row.soft_observed_contact_construction_attempted -ne '1' -or
        $row.soft_observed_contact_construction_complete -ne '1' -or
        $row.soft_observed_contact_cell_setup_failures -ne '0' -or
        $row.soft_observed_contact_atom_index_failures -ne '0' -or
        $row.soft_observed_contact_mic_query_failures -ne '0') {
        throw 'Frame report evidence schema metadata is invalid'
    }
    if ($row.identified -ne $row.evidence_graph_connected) {
        throw 'Legacy identified alias disagrees with evidence_graph_connected'
    }
    if ($row.steric_zero_lj_policy -ne 'force_field_no_hard_core') {
        throw "Unexpected zero-LJ steric policy: $($row.steric_zero_lj_policy)"
    }
    $expectedFeasibleSemantics = if ($row.feasible_assignment_set_enumerated -eq '1') {
        'hard_constraints_within_search_domain'
    } else {
        'not_enumerated'
    }
    if ($row.feasible_assignment_set_semantics -ne $expectedFeasibleSemantics) {
        throw 'Frame feasible-set semantics is inconsistent'
    }
    if ($row.status -notmatch '^(certified|weak_observation|rescued|fallback|contradicted)$') {
        throw "Unexpected status: $($row.status)"
    }
    if ($row.fallback_reason -notmatch '^(none|no_seed|projection_invalid|full_verify_clash|continuity_reject|anchor_reject|split_owner_reject|component_graph_failed|beam_empty|not_better_than_incumbent|evidence_disconnected|evidence_underdetermined|evidence_contradicted|evidence_hypothesis_unavailable|certification_reject|fallback_pbci|lattice_query_failed|invalid_lattice)$') {
        throw "Unexpected fallback reason: $($row.fallback_reason)"
    }
    if ([double]$row.time_total_ms -lt 0.0) {
        throw "Negative total timing"
    }
    if ($row.component_audit_exact_requested -ne '1' -or
        $row.component_audit_exact_completed -ne '0' -or
        $row.component_audit_exact_incomplete -ne '1' -or
        $row.component_audit_exact_incomplete_reason -ne 'state_limit' -or
        [uint64]$row.component_audit_exact_states_visited -ne 1 -or
        $row.component_audit_exact_reference_comparison_available -ne '0') {
        throw 'Resource-limited audit-exact did not fail closed'
    }
    if ($row.component_audit_exact_memory_budget_semantics -ne
        'estimated_workspace') {
        throw 'Audit memory-budget semantics are not explicit'
    }
    if ($row.status -eq 'certified' -and
        $row.best_score_vector -match '(?:^|;)shift=na(?:;|$)') {
        throw "Certified candidate reported a non-finite shift norm"
    }
    if ($row.status -eq 'certified' -and $row.mode -ne 'llps-local') {
        if ($row.certificate_scope -ne 'image_intercomponent_steric') {
            throw "Unexpected certificate scope: $($row.certificate_scope)"
        }
        if ($row.certificate_graph_source -ne 'gromacs_topology') {
            throw "Unexpected certificate graph source: $($row.certificate_graph_source)"
        }
        if ($row.hard_evaluated -ne '1' -or $row.hard_certified -ne '1') {
            throw 'Certified frame is missing a successful hard-feasibility check'
        }
        $hardEdgesAvailable = [uint64]$row.hard_connectivity_edges_available
        $topologyEdgesAvailable = [uint64]$row.topology_edges_available
        $topologyEdgesChecked = [uint64]$row.hard_topology_edges_checked
        if ($hardEdgesAvailable -eq 0 -or $topologyEdgesAvailable -eq 0 -or
            $topologyEdgesChecked -gt $topologyEdgesAvailable) {
            throw 'Certified frame reported invalid hard-connectivity coverage'
        }
        $expectedTopologyFraction =
            [double]$topologyEdgesChecked / [double]$topologyEdgesAvailable
        if ([math]::Abs([double]$row.hard_topology_edge_fraction -
                        $expectedTopologyFraction) -gt 1.0e-9) {
            throw 'Certified frame reported an inconsistent topology-edge fraction'
        }
        if ([uint64]$row.metadata_edges_available -ne 0) {
            throw 'Topology-only fixture unexpectedly reported metadata edges'
        }
        if ([uint64]$row.certificate_atoms_total -ne 5 -or
            [uint64]$row.hard_graph_atoms_covered -ne 4 -or
            [uint64]$row.topology_bond_atoms_covered -ne 4 -or
            [uint64]$row.metadata_only_atoms -ne 0 -or
            [math]::Abs([double]$row.hard_graph_atom_fraction - 0.8) -gt 1.0e-9 -or
            [math]::Abs([double]$row.topology_bond_atom_fraction - 0.8) -gt 1.0e-9 -or
            [math]::Abs([double]$row.metadata_only_atom_fraction) -gt 1.0e-9) {
            throw 'Certified frame reported inconsistent atom-level graph coverage'
        }
        if ([uint64]$row.hard_unsupported_winding_rank -ne 0) {
            throw 'Certified frame reported unsupported topological winding'
        }
        if ([uint64]$row.assignment_signature -eq 0 -or
            [int]$row.assignment_component_count -le 0 -or
            [uint64]$row.topology_epoch_id -eq 0) {
            throw 'Certified frame is missing a canonical image assignment'
        }
    }
}

$assignmentReport = Join-Path $WorkRoot 'pbctopo_assignment_candidates.csv'
if (-not (Test-Path -LiteralPath $assignmentReport)) {
    throw "Missing assignment candidate report: $assignmentReport"
}
Assert-CsvFieldParity -Path $assignmentReport
$assignmentRows = @(Import-Csv -LiteralPath $assignmentReport)
if ($assignmentRows.Count -eq 0) {
    throw 'Assignment candidate report contains no certified assignments'
}
$assignmentHeader =
    (Get-Content -LiteralPath $assignmentReport -First 1).Split(',')
$requiredAssignmentColumns = @(
    'schema_version',
    'producer_version',
    'coordinate_unit',
    'time_unit',
    'endianness',
    'identity_hash_algorithm',
    'binding_hash_algorithm',
    'payload_checksum_algorithm',
    'observation_scope',
    'topology_epoch_id',
    'frame',
    'candidate_rank',
    'selected_framewise',
    'candidate_set_scope',
    'observation_class',
    'search_equivalence_known',
    'equivalent_to_search_best',
    'equivalent_under_output_order',
    'equivalent_set_complete',
    'search_equivalent_assignments_exported',
    'equivalent_best_assignments_within_domain',
    'objective_uniqueness_known',
    'objective_unique_within_search_domain',
    'evidence_uniqueness_known',
    'evidence_unique_within_search_domain',
    'feasible_assignment_set_enumerated',
    'feasible_assignment_set_semantics',
    'feasible_assignments_within_domain',
    'soft_score_valid_assignments_within_domain',
    'search_domain_assignments_exported',
    'credible_alternative_set_complete',
    'component_search_shell_radius',
    'search_domain_audit_known',
    'search_domain_evidence_complete',
    'selected_assignment_touches_shell_boundary',
    'domain_expansion_attempted',
    'domain_expansion_exhausted',
    'bounded_domain_only',
    'layout_signature',
    'assignment_signature',
    'component_assignment',
    'normalized_gram_g00',
    'normalized_gram_g22'
)
foreach ($column in $requiredAssignmentColumns) {
    if ($assignmentHeader -notcontains $column) {
        throw "Missing assignment candidate column: $column"
    }
}
foreach ($assignmentRow in $assignmentRows) {
    if ($assignmentRow.schema_version -ne '8' -or
        $assignmentRow.coordinate_unit -ne 'angstrom' -or
        $assignmentRow.time_unit -ne 'ps' -or
        [uint64]$assignmentRow.topology_epoch_id -eq 0) {
        throw 'Assignment artifact schema metadata is invalid'
    }
    if ($assignmentRow.hard_feasible -ne '1' -or
        $assignmentRow.candidate_set_scope -ne 'retained_certified_band') {
        throw 'Assignment artifact contains a non-certified or unlabeled row'
    }
    $expectedFeasibleSemantics = if ($assignmentRow.feasible_assignment_set_enumerated -eq '1') {
        'hard_constraints_within_search_domain'
    } else {
        'not_enumerated'
    }
    if ($assignmentRow.feasible_assignment_set_semantics -ne $expectedFeasibleSemantics) {
        throw 'Assignment feasible-set semantics is inconsistent'
    }
    if ($assignmentRow.component_assignment -notmatch
        '^\d+:-?\d+:-?\d+:-?\d+(\|\d+:-?\d+:-?\d+:-?\d+)*$') {
        throw "Invalid stable component assignment: $($assignmentRow.component_assignment)"
    }
    if ($assignmentRow.normalized_gram_g00 -eq 'na' -or
        $assignmentRow.normalized_gram_g22 -eq 'na') {
        throw 'Assignment artifact is missing the normalized lattice metric'
    }
    if ($assignmentRow.search_equivalence_known -eq '0' -and
        $assignmentRow.equivalent_to_search_best -ne 'na') {
        throw 'Unknown exact-search equivalence was emitted as a boolean'
    }
    if ($assignmentRow.search_equivalence_known -eq '1' -and
        $assignmentRow.equivalent_to_search_best -notin @('0', '1')) {
        throw 'Known exact-search equivalence is not a boolean'
    }
    if ($assignmentRow.equivalent_set_complete -eq '1' -and
        $assignmentRow.search_equivalent_assignments_exported -ne
            $assignmentRow.equivalent_best_assignments_within_domain) {
        throw 'Complete equivalent set count does not match exact search'
    }
    if ($assignmentRow.credible_alternative_set_complete -eq '1' -and
        $assignmentRow.search_domain_assignments_exported -ne
            $assignmentRow.feasible_assignments_within_domain) {
        throw 'Complete credible alternative set does not match exact search'
    }
    if ($assignmentRow.observation_class -eq 'strong_global_anchor') {
        throw 'Schema-v6 production output emitted unsupported global certification'
    }
    if ($assignmentRow.observation_class -eq 'strong_bounded_anchor' -and
        ($assignmentRow.selected_framewise -ne '1' -or
         $assignmentRow.search_equivalence_known -ne '1' -or
         $assignmentRow.equivalent_to_search_best -ne '1' -or
         $assignmentRow.objective_unique_within_search_domain -ne '1' -or
         $assignmentRow.evidence_unique_within_search_domain -ne '1' -or
         $assignmentRow.credible_alternative_set_complete -ne '1' -or
         $assignmentRow.search_domain_audit_known -ne '1' -or
         $assignmentRow.search_domain_evidence_complete -ne '1' -or
         $assignmentRow.selected_assignment_touches_shell_boundary -ne '0')) {
        throw 'Strong anchor lacks its objective/evidence proof contract'
    }
}
if (@($assignmentRows | Where-Object selected_framewise -eq '1').Count -eq 0) {
    throw 'Assignment artifact did not mark a framewise selected state'
}

$layoutReport = Join-Path $WorkRoot 'pbctopo_component_layout.csv'
if (-not (Test-Path -LiteralPath $layoutReport)) {
    throw "Missing component layout report: $layoutReport"
}
Assert-CsvFieldParity -Path $layoutReport
$layoutRows = @(Import-Csv -LiteralPath $layoutReport)
if ($layoutRows.Count -eq 0) {
    throw 'Component layout report contains no exact membership rows'
}
foreach ($layoutRow in $layoutRows) {
    if ($layoutRow.schema_version -ne '8' -or
        $layoutRow.coordinate_unit -ne 'angstrom' -or
        $layoutRow.canonical_atom_universe_hash -notmatch '^[0-9a-fA-F]{32}$' -or
        $layoutRow.hard_graph_source -ne 'explicit_topology' -or
        $layoutRow.exact_atom_membership -notmatch
            '^-?\d+:\d+(\|-?\d+:\d+)*$' -or
        [uint64]$layoutRow.topology_epoch_id -eq 0 -or
        [uint64]$layoutRow.component_id -eq 0) {
        throw 'Component layout row violates schema v7'
    }
}
$hardEdgeReport = Join-Path $WorkRoot 'pbctopo_hard_edges.csv'
if (-not (Test-Path -LiteralPath $hardEdgeReport)) {
    throw "Missing hard-edge report: $hardEdgeReport"
}
Assert-CsvFieldParity -Path $hardEdgeReport
$hardEdgeRows = @(Import-Csv -LiteralPath $hardEdgeReport)
if ($hardEdgeRows.Count -eq 0) {
    throw 'Hard-edge report contains no topology bonds'
}
foreach ($hardEdgeRow in $hardEdgeRows) {
    if ($hardEdgeRow.schema_version -ne '8' -or
        $hardEdgeRow.observation_scope -ne 'static_hard_edges' -or
        $hardEdgeRow.edge_kind -ne 'topology_bond' -or
        $hardEdgeRow.hard_graph_source -ne 'explicit_topology' -or
        [uint64]$hardEdgeRow.topology_epoch_id -eq 0) {
        throw 'Hard-edge row violates schema v7'
    }
}
$reportedEpochs = @($rows | Where-Object status -eq 'certified' |
    Select-Object -ExpandProperty topology_epoch_id -Unique)
$layoutEpochs = @($layoutRows |
    Select-Object -ExpandProperty topology_epoch_id -Unique)
foreach ($epoch in $reportedEpochs) {
    if ($layoutEpochs -notcontains $epoch) {
        throw "Frame report topology epoch is missing from layout sidecar: $epoch"
    }
}

$pythonPath = (Get-Command python -ErrorAction Stop).Source
$validator = Join-Path $repo 'tools\validate_vibe_observation.py'
& $pythonPath $validator $WorkRoot
if ($LASTEXITCODE -ne 0) {
    throw "VIBE observation validator failed: $LASTEXITCODE"
}

function Assert-ValidatorRejects {
    param(
        [string]$CaseName,
        [scriptblock]$Mutate
    )
    $tamperRoot = Join-Path $WorkRoot "validator_tamper_$CaseName"
    if (Test-Path -LiteralPath $tamperRoot) {
        Remove-Item -LiteralPath $tamperRoot -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $tamperRoot | Out-Null
    foreach ($artifact in @(
            'pbctopo_component_layout.csv',
            'pbctopo_hard_edges.csv',
            'pbctopo_frame_report.csv',
            'pbctopo_assignment_candidates.csv')) {
        Copy-Item -LiteralPath (Join-Path $WorkRoot $artifact) `
            -Destination (Join-Path $tamperRoot $artifact)
    }
    & $Mutate $tamperRoot
    $savedErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        & $pythonPath $validator $tamperRoot 2>&1 | Out-Null
        $validatorExitCode = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
    if ($validatorExitCode -eq 0) {
        throw "VIBE validator accepted tampered artifact: $CaseName"
    }
    $global:LASTEXITCODE = 0
}

Assert-ValidatorRejects 'component_identity' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_component_layout.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].component_id =
        ([uint64]$tampered[0].component_id + [uint64]1).ToString()
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'strong_anchor_contract' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_assignment_candidates.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].observation_class = 'strong_bounded_anchor'
    $tampered[0].selected_framewise = '0'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'selected_frame_signature' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].assignment_signature = '1'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'derived_candidate_set_contract' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_assignment_candidates.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].credible_alternative_set_complete =
        if ($tampered[0].credible_alternative_set_complete -eq '1') {
            '0'
        } else {
            '1'
        }
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'canonical_atom_universe' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].canonical_atom_universe_hash =
        '11111111111111111111111111111111'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'spatial_ambiguity_bound' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].spatial_lift_ambiguous_hard_edges = '999'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'spatial_cycle_bound' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].spatial_lift_cycle_residuals = '999'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'payload_hash_contract' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_assignment_candidates.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].payload_checksum_algorithm = 'sha256-128'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'unsupported_global_claim' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_assignment_candidates.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].domain_completeness = 'global_certified'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'fallback_stale_evidence_payload' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].status = 'fallback'
    $tampered[0].evidence_state = 'not_evaluated'
    $tampered[0].evidence_consistency_evaluated = '0'
    $tampered[0].evidence_consistent = '0'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Assert-ValidatorRejects 'local_unwrap_stale_evidence_payload' {
    param($tamperRoot)
    $path = Join-Path $tamperRoot 'pbctopo_frame_report.csv'
    $tampered = @(Import-Csv -LiteralPath $path)
    $tampered[0].status = 'local_unwrap'
    $tampered[0].evidence_state = 'not_applicable'
    $tampered[0].evidence_consistency_evaluated = '0'
    $tampered[0].evidence_consistent = '0'
    $tampered | Export-Csv -LiteralPath $path -NoTypeInformation `
        -Encoding ASCII
}

Write-Host 'PASS pbctopo_confidence_report_smoke'
