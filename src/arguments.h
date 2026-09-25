#ifndef TITAN_VIBE_ARGUMENTS_H
#define TITAN_VIBE_ARGUMENTS_H

#include <cstddef>
#include <string>

// Minimal host-configuration surface used by the published VIBE adapter.
// The full TITAN command-line parser is intentionally outside this release.
struct PbcArgs {
  int chain_ref{-1};
  std::string pbctopo_mode{"vibe"};
  std::string pbctopo_search_policy{"fast-certified"};
  int pbctopo_audit_max_states{10000000};
  double pbctopo_audit_max_seconds{60.0};
  int pbctopo_audit_max_memory_mb{1024};
  std::string pbctopo_anchor_backend{"auto"};
  std::string pbctopo_anchor_metrics{"full"};
  std::string pbctopo_temporal{"greedy"};
  int pbctopo_temporal_k{4};
  double pbctopo_temporal_penalty{2.0};
  double pbctopo_vdw_clash_scale{0.55};
  double pbctopo_bond_cutoff{0.0};
  bool pbctopo_bond_cutoff_set{false};
  double pbctopo_clash_cutoff{0.0};
  bool pbctopo_clash_cutoff_set{false};
  double pbctopo_contact_cutoff{0.0};
  bool pbctopo_contact_cutoff_set{false};
  std::string pbctopo_required_contacts_file{};
  double pbctopo_continuity_cutoff{0.0};
  bool pbctopo_continuity_cutoff_set{false};
  int pbctopo_exact_component_limit{0};
  bool pbctopo_exact_component_limit_set{false};
  int pbctopo_beam_width{0};
  bool pbctopo_beam_width_set{false};
  bool pbctopo_evidence_export{false};
  // Host fields present at the bound TITAN source commit; defaults unchanged.
  int pbctopo_component_shell_radius{0};
  bool pbctopo_component_shell_radius_set{false};
  bool pbctopo_h9_orbit{false};
};

struct arguments {
  PbcArgs pbc;
  std::size_t nt{0};
  std::string force_f{"AMBER99SB"};
};

#endif
