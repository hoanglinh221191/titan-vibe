#include "tabi/tabi_fmm_kifmm_cartesian_bridge.h"

#include <iostream>

int main() {
  const auto result = titan::tabi::fmm_kifmm::run_cartesian_kifmm_bridge_smoke_suite();
  bool ok = result.decision == "fmm_e26_synthetic_single_box_hybrid_bridge_pass" ||
            result.decision == "fmm_e26_q8_required_q6_gap";
  for (const auto &row : result.rows) {
    std::cout << "E26_ROW"
              << " row_id=" << row.row_id
              << " execution_status=" << row.execution_status
              << " accuracy_status=" << row.accuracy_status
              << " accuracy_failure=" << row.accuracy_failure
              << " aggregate_rms=" << row.aggregate_relative_rms
              << " k1_rms=" << row.k1_relative_rms
              << " k2_rms=" << row.k2_relative_rms
              << " k3_rms=" << row.k3_relative_rms
              << " k4_rms=" << row.k4_relative_rms
              << " max_channel_name=" << row.max_channel_name
              << " max_channel_rms=" << row.max_channel_relative_rms
              << " m2e_residual=" << row.m2e_max_relative_residual
              << " e2l_coeff_rms=" << row.e2l_coefficient_relative_rms
              << " l2p_gradient_method=" << row.l2p_gradient_method
              << " l2p_gradient_rms=" << row.l2p_target_gradient_relative_rms
              << " surface_order=" << row.surface_order
              << " cartesian_order=" << row.cartesian_order
              << " policy_hash=" << row.policy_hash << '\n';
    ok = ok && row.execution_status == "pass";
    ok = ok && row.e2l_coefficient_relative_rms < 0.0;
    if (row.row_id == "scalar_e2l_laplace_q8")
      ok = ok && row.l2p_gradient_method == "analytic_cartesian_local";
  }
  std::cout << "E26_DECISION"
            << " decision=" << result.decision
            << " rows=" << result.rows.size()
            << " required_rows=" << result.required_rows << '\n';
  return ok ? 0 : 2;
}
