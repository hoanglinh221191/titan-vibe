#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_weak_form.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::MixedTraceCalderonConfig;
using TitanPBE::Tabi::calderon::assemble_jcp_mixed_trace_pairing_calibration;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

} // namespace

int main() {
  try {
    const std::vector<Atom> sphere = {{0.0, 0.0, 0.0, 0.0, 1.5}};
    const auto auxiliary = build_implicit_vdw_auxiliary_surface(sphere, 1.00);
    require(auxiliary.status == "pass", "auxiliary sphere must be watertight");

    MixedTraceCalderonConfig config;
    config.eps_ratio = 80.0;
    config.kappa = 0.125;
    config.adjacent_pair_subdivision_depth = 4;
    const auto result = assemble_jcp_mixed_trace_pairing_calibration(
        auxiliary.mesh, auxiliary.trace_spaces, config);
    std::cout << "calderon_c5c1_mixed_trace_diagnostic status=" << result.status
              << " p1=" << result.p1_dofs << " dual_p0=" << result.dual_p0_dofs
              << " patches=" << result.dual_p0_subtriangles
              << " partition=" << result.dual_p0_partition_relative
              << " pivot=" << result.pairing_pivot_ratio
              << " inverse=" << result.pairing_inverse_relative
              << " k1_sym=" << result.k1_symmetry_relative
              << " k23_negative_eps="
              << result.k2_k3_negative_epsilon_closure_relative
              << " k23_positive_eps="
              << result.k2_k3_positive_epsilon_closure_relative
              << " k23_nonclosure="
              << result.k2_k3_separated_nonclosure_relative
              << " k4_separated_sym="
              << result.k4_separated_pair_symmetry_relative
              << " zero=" << result.zero_screening_kernel_absolute << '\n';
    require(result.p1_dofs == result.dual_p0_dofs,
            "dual P0 trace must have one cell per P1 vertex");
    require(result.dual_p0_subtriangles == 6 * auxiliary.mesh.triangles.size(),
            "each primal triangle must contribute six dual subtriangles");
    require(result.dual_p0_partition_relative < 1.0e-12,
            "dual P0 cells must partition the closed surface");
    require(result.pairing_inverse_relative < 1.0e-11,
            "mixed P1/dual-P0 Riesz pairing must be invertible");
    require(result.status == "pass",
            "mixed trace calibration must pass the JCP kernel-parity gate");
    std::cout << "calderon_c5c1_mixed_trace_smoke=pass k23_relation="
              << result.derivative_block_relation
              << " k4_status=" << result.screened_hypersingular_status << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c1_mixed_trace_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
