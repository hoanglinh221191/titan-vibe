param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_evidence_semantics_smoke')
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/evidence_audit.h"
#include "pbctopo/report.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace titan_pbctopo;

struct Point { double x, y, z; };

void expect(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL pbctopo_evidence_semantics_smoke: " << message << '\n';
    std::exit(1);
  }
}

int main() {
  expect(pbctopo_resolve_fallback_reason(
             false, PbctopoFailureReason::EvidenceHypothesisUnavailable) ==
             PbctopoFailureReason::EvidenceHypothesisUnavailable,
         "specific evidence failure was overwritten by PBCI fallback");
  expect(pbctopo_resolve_fallback_reason(false,
                                         PbctopoFailureReason::None) ==
             PbctopoFailureReason::FallbackPbci,
         "generic PBCI fallback reason was not selected");
  expect(!pbctopo_frame_status_has_selected_assignment(
             PbctopoFrameStatus::Fallback) &&
             !pbctopo_frame_status_has_selected_assignment(
                 PbctopoFrameStatus::LocalUnwrap) &&
             pbctopo_frame_status_has_selected_assignment(
                 PbctopoFrameStatus::WeakObservation),
         "fallback selected-assignment status contract is invalid");
  const std::array<PbctopoCanonicalAtomPair, 2> ordered{{{7, 11}, {20, 30}}};
  const auto hash = hash_pbctopo_canonical_atom_pairs(ordered);
  expect(hash.has_value() && !hash->empty(), "pair hash was unavailable");
  const std::array<PbctopoCanonicalAtomPair, 2> changed{{{7, 11}, {20, 31}}};
  const auto changed_hash = hash_pbctopo_canonical_atom_pairs(changed);
  expect(changed_hash && *changed_hash != *hash,
         "pair identity change did not change hash");
  const std::array<PbctopoCanonicalAtomPair, 2> reversed_order{{{20, 30}, {7, 11}}};
  expect(!hash_pbctopo_canonical_atom_pairs(reversed_order),
         "noncanonical pair order was accepted");

  std::array<PbctopoAtom, 2> atoms{};
  atoms[0].global_atom_id = 11;
  atoms[1].global_atom_id = 7;
  const std::array<PbctopoAnchorContactEdge, 1> edges{{{0, 1, 16.0}}};
  const std::array<Point, 2> supported_points{{{0, 0, 0}, {3, 0, 0}}};
  const auto supported = audit_pbctopo_selected_evidence(
      std::span<const PbctopoAtom>(atoms), std::span<const Point>(supported_points),
      std::span<const PbctopoAnchorContactEdge>(edges), true, true);
  expect(supported && supported->state == PbctopoEvidenceState::SupportedConnected &&
             supported->contacts_lost == 0,
         "supported connected state was not produced");
  const auto weak = audit_pbctopo_selected_evidence(
      std::span<const PbctopoAtom>(atoms), std::span<const Point>(supported_points),
      std::span<const PbctopoAnchorContactEdge>{}, true, false);
  expect(weak && weak->state == PbctopoEvidenceState::WeakDisconnected,
         "empty evidence was treated as contradiction");
  const std::array<Point, 2> contradicted_points{{{0, 0, 0}, {13, 0, 0}}};
  const auto contradicted = audit_pbctopo_selected_evidence(
      std::span<const PbctopoAtom>(atoms),
      std::span<const Point>(contradicted_points),
      std::span<const PbctopoAnchorContactEdge>(edges), true, true);
  expect(contradicted &&
             contradicted->state == PbctopoEvidenceState::Contradicted &&
             contradicted->evidence_graph_connected &&
             contradicted->contacts_lost == 1,
         "graph-connected contradiction was not detected by direct distance");
  expect(contradicted->observed_pair_hash == supported->observed_pair_hash,
         "pair orientation/global-id canonicalization was unstable");
  std::array<PbctopoAtom, 4> order_atoms{};
  order_atoms[0].global_atom_id = 11;
  order_atoms[1].global_atom_id = 7;
  order_atoms[2].global_atom_id = 30;
  order_atoms[3].global_atom_id = 20;
  const std::array<Point, 4> order_points{{
      {0, 0, 0}, {3, 0, 0}, {0, 0, 0}, {13, 0, 0}}};
  const std::array<PbctopoAnchorContactEdge, 2> forward_edges{{
      {0, 1, 16.0}, {2, 3, 16.0}}};
  const std::array<PbctopoAnchorContactEdge, 2> reversed_edges{{
      {3, 2, 16.0}, {1, 0, 16.0}}};
  const auto forward = audit_pbctopo_selected_evidence(
      std::span<const PbctopoAtom>(order_atoms),
      std::span<const Point>(order_points),
      std::span<const PbctopoAnchorContactEdge>(forward_edges), true, true);
  const auto reversed = audit_pbctopo_selected_evidence(
      std::span<const PbctopoAtom>(order_atoms),
      std::span<const Point>(order_points),
      std::span<const PbctopoAnchorContactEdge>(reversed_edges), true, true);
  expect(forward && reversed &&
             forward->contacts_checked == reversed->contacts_checked &&
             forward->contacts_lost == reversed->contacts_lost &&
             forward->observed_pair_hash == reversed->observed_pair_hash &&
             forward->lost_pair_hash == reversed->lost_pair_hash,
         "edge input order/orientation changed audit identity");
  atoms[1].global_atom_id = atoms[0].global_atom_id;
  expect(!audit_pbctopo_selected_evidence(
              std::span<const PbctopoAtom>(atoms),
              std::span<const Point>(supported_points),
              std::span<const PbctopoAnchorContactEdge>(edges), true, true),
         "duplicate global atom identities were not rejected");

  std::array<PbctopoAtom, 4> hypothesis_atoms{};
  for (std::size_t i = 0; i < hypothesis_atoms.size(); ++i)
    hypothesis_atoms[i].global_atom_id = 100 + i;
  const std::array<std::size_t, 4> atom_components{{0, 1, 0, 1}};
  const std::array<Int3, 4> lift_images{};
  PbctopoImageAssignment identity;
  identity.component_offsets.resize(2);
  identity.component_ids = {10, 20};
  identity.valid = true;
  const std::array<PbctopoAnchorContactEdge, 2> hypothesis_edges{{
      {0, 1, 16.0, {1, 0, 0}, 1, true},
      {2, 3, 16.0, {0, 0, 0}, 1, true}}};
  const auto hypotheses = build_pbctopo_contact_hypotheses(
      hypothesis_atoms, hypothesis_edges, atom_components, lift_images,
      identity);
  expect(hypotheses.complete && hypotheses.hypotheses.size() == 2 &&
             hypotheses.component_relation_count == 1 &&
             hypotheses.ambiguous_component_relation_count == 1,
         "incompatible shift hypotheses were not kept as alternatives");

  PbctopoImageAssignment selected = identity;
  selected.component_offsets[1] = {1, 0, 0};
  const std::array<Point, 4> selected_points{{
      {0, 0, 0}, {1, 0, 0}, {0, 3, 0}, {10, 3, 0}}};
  const auto ambiguous = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>(hypothesis_atoms),
      std::span<const Point>(selected_points), hypotheses, selected, true);
  expect(ambiguous && ambiguous->state == PbctopoEvidenceState::WeakAmbiguous &&
             ambiguous->contacts_checked == 2 &&
             ambiguous->selected_compatible_contact_count == 1 &&
             ambiguous->contacts_lost == 0 &&
             !ambiguous->selected_hypothesis_hash.empty(),
         "incompatible alternative was incorrectly counted as a lost contact");
  expect(ambiguous->relation_audits.size() == 1 &&
             ambiguous->relation_audits[0].component_a == 0 &&
             ambiguous->relation_audits[0].component_b == 1 &&
             ambiguous->relation_audits[0].component_id_a == 10 &&
             ambiguous->relation_audits[0].component_id_b == 20 &&
             ambiguous->relation_audits[0].observed_hypothesis_relation &&
             ambiguous->relation_audits[0].support.compatibility_match &&
             ambiguous->relation_audits[0].support.active_at_threshold &&
             ambiguous->relation_audits[0].support.support_primary == 1.0 &&
             ambiguous->relation_audits[0]
                     .support.closest_shift_linf_distance == 0 &&
             ambiguous->relation_audits[0].compatible_contact_count == 1 &&
             ambiguous->relation_audits[0]
                     .alternative_hypothesis_contact_count == 1,
         "compatible relation support record was not provenance complete");
  selected.component_offsets[1] = {2, 0, 0};
  const auto no_support = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>(hypothesis_atoms),
      std::span<const Point>(selected_points), hypotheses, selected, true);
  expect(no_support && no_support->state == PbctopoEvidenceState::WeakAmbiguous &&
             no_support->compatible_hypothesis_count == 0 &&
             no_support->no_supported_contact_relation_count == 1 &&
             no_support->selected_hypothesis_hash.empty(),
         "no-support complement was serialized as a selected hypothesis");
  expect(no_support->relation_audits.size() == 1 &&
             !no_support->relation_audits[0].support.compatibility_match &&
             !no_support->relation_audits[0].support.active_at_threshold &&
             no_support->relation_audits[0].support.support_primary == 0.0 &&
             no_support->relation_audits[0]
                     .support.closest_shift_linf_distance_defined &&
             no_support->relation_audits[0]
                     .support.closest_shift_linf_distance == 1 &&
             no_support->relation_audits[0].support.reason ==
                 PbctopoEvidenceSupportReason::NoCompatibleHypothesis &&
             no_support->evidence_edge_count == 0,
         "inactive relation support did not replay the disconnected graph");

  const std::array<std::pair<std::size_t, std::size_t>, 1>
      explicit_relation{{{1, 0}}};
  const auto explicit_no_support = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>(hypothesis_atoms),
      std::span<const Point>(selected_points), hypotheses, selected, true,
      explicit_relation);
  expect(explicit_no_support && explicit_no_support->relation_audits.size() == 1 &&
             explicit_no_support->relation_audits[0]
                 .support.explicit_relation &&
             explicit_no_support->relation_audits[0]
                 .support.active_at_threshold &&
             explicit_no_support->relation_audits[0]
                     .support.support_primary == 1.0 &&
             explicit_no_support->relation_audits[0].support.reason ==
                 PbctopoEvidenceSupportReason::ExplicitRelation &&
             explicit_no_support->evidence_edge_count == 1 &&
             explicit_no_support->explicit_contact_evidence_edges == 1 &&
             explicit_no_support->scored_contact_evidence_edges == 0 &&
             explicit_no_support->evidence_graph_connected,
         "explicit relation did not override an incompatible hypothesis");

  PbctopoContactHypothesisSet empty_hypotheses;
  empty_hypotheses.complete = true;
  const auto explicit_only = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>(hypothesis_atoms),
      std::span<const Point>(selected_points), empty_hypotheses, selected, true,
      explicit_relation);
  expect(explicit_only && explicit_only->relation_audits.size() == 1 &&
             !explicit_only->relation_audits[0]
                  .observed_hypothesis_relation &&
             explicit_only->relation_audits[0].support.explicit_relation &&
             explicit_only->relation_audits[0].support.active_at_threshold &&
             !explicit_only->relation_audits[0]
                  .support.closest_shift_linf_distance_defined &&
             explicit_only->evidence_edge_count == 1 &&
             explicit_only->evidence_graph_connected,
         "explicit-only relation was omitted from the support audit");

  PbctopoContactHypothesisSet chain_hypotheses;
  chain_hypotheses.hypotheses = {
      {0, 1, 10, 20, {0, 0, 0}, {}},
      {1, 2, 20, 30, {1, 0, 0}, {}}};
  chain_hypotheses.component_relation_count = 2;
  chain_hypotheses.complete = true;
  PbctopoImageAssignment chain_selected;
  chain_selected.component_ids = {10, 20, 30};
  chain_selected.component_offsets = {{0, 0, 0}, {0, 0, 0}, {2, 0, 0}};
  chain_selected.valid = true;
  const auto weak_chain = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>{}, std::span<const Point>{},
      chain_hypotheses, chain_selected, true);
  expect(weak_chain && weak_chain->relation_audits.size() == 2 &&
             weak_chain->relation_audits[0].support.active_at_threshold &&
             !weak_chain->relation_audits[1].support.active_at_threshold &&
             weak_chain->evidence_edge_count == 1 &&
             weak_chain->scored_contact_evidence_edges == 1 &&
             weak_chain->evidence_connected_component_count == 2 &&
             weak_chain->relative_lattice_degrees_of_freedom == 3 &&
             !weak_chain->evidence_graph_connected,
         "relation reduction did not reproduce a weak three-component chain");
  const std::array<std::pair<std::size_t, std::size_t>, 1>
      chain_explicit{{{2, 1}}};
  const auto connected_chain = audit_pbctopo_selected_hypotheses(
      std::span<const PbctopoAtom>{}, std::span<const Point>{},
      chain_hypotheses, chain_selected, true, chain_explicit);
  expect(connected_chain && connected_chain->relation_audits.size() == 2 &&
             connected_chain->relation_audits[0].support.active_at_threshold &&
             connected_chain->relation_audits[1].support.active_at_threshold &&
             connected_chain->evidence_edge_count == 2 &&
             connected_chain->explicit_contact_evidence_edges == 1 &&
             connected_chain->scored_contact_evidence_edges == 1 &&
             connected_chain->evidence_connected_component_count == 1 &&
             connected_chain->relative_lattice_degrees_of_freedom == 0 &&
             connected_chain->evidence_graph_connected,
         "explicit support did not close the three-component evidence chain");

  auto duplicate_ids = selected;
  duplicate_ids.component_ids = {10, 10};
  expect(!audit_pbctopo_selected_hypotheses(
              std::span<const PbctopoAtom>(hypothesis_atoms),
              std::span<const Point>(selected_points), hypotheses,
              duplicate_ids, true),
         "duplicate stable component IDs were accepted by the support audit");

  const Int3 n_ab{1, 0, 0};
  Int3 n_ba{};
  expect(checked_int3_negate(n_ab, n_ba) && n_ba.x == -1,
         "pair reversal did not negate the MIC image");
  const double b_minus_a = -9.0;
  const double a_minus_b = 9.0;
  expect(b_minus_a + 10.0 * static_cast<double>(n_ab.x) ==
             -(a_minus_b + 10.0 * static_cast<double>(n_ba.x)),
         "reversed MIC triplet does not reproduce the same Cartesian vector");
  std::cout << "PASS pbctopo_evidence_semantics_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src src\pbctopo\evidence_audit.cpp -o $exe
    if ($LASTEXITCODE -ne 0) { throw 'Failed to compile evidence semantics smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw "Evidence semantics smoke failed: $LASTEXITCODE" }
}
finally { Pop-Location }
