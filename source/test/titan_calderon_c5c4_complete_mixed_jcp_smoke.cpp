#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_weak_form.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::CompleteMixedJcpConfig;
using TitanPBE::Tabi::calderon::assemble_complete_mixed_jcp_calibration;
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
    const auto coarse_auxiliary =
        build_implicit_vdw_auxiliary_surface(sphere, 1.00);
    const auto fine_auxiliary =
        build_implicit_vdw_auxiliary_surface(sphere, 0.75);
    require(coarse_auxiliary.status == "pass" && fine_auxiliary.status == "pass",
            "auxiliary spheres must be watertight");

    CompleteMixedJcpConfig config;
    config.eps_ratio = 80.0;
    config.kappa = 0.125;
    config.adjacent_pair_subdivision_depth = 2;
    config.direct_action_relative_tolerance = 0.025;
    const auto coarse = assemble_complete_mixed_jcp_calibration(
        coarse_auxiliary.mesh, coarse_auxiliary.trace_spaces, config);
    const auto fine = assemble_complete_mixed_jcp_calibration(
        fine_auxiliary.mesh, fine_auxiliary.trace_spaces, config);
    std::cout << "calderon_c5c4_complete_mixed_diagnostic coarse_status="
              << coarse.status << " coarse_p1=" << coarse.p1_dofs
              << " coarse_system=" << coarse.strong_system_relative
              << " fine_status=" << fine.status << " fine_p1=" << fine.p1_dofs
              << " fine_dual=" << fine.dual_p0_dofs
              << " fine_patches=" << fine.dual_p0_subtriangles
              << " fine_pairing_pivot=" << fine.pairing_pivot_ratio
              << " fine_pairing_inverse=" << fine.pairing_inverse_relative
              << " fine_k4=" << fine.direct_k4_weak_relative
              << " fine_interaction=" << fine.strong_interaction_relative
              << " fine_system=" << fine.strong_system_relative
              << " fine_probe=" << fine.deterministic_probe_relative
              << " lhs_phi=" << fine.lhs_phi
              << " lhs_phin=" << fine.lhs_phin
              << " weak_sign=" << fine.weak_form_sign << '\n';
    require(coarse.status == "pass" && fine.status == "pass",
            "complete mixed JCP calibration must pass on both sphere levels");
    require(fine.p1_dofs > coarse.p1_dofs,
            "fine auxiliary mesh must add trace degrees of freedom");
    require(fine.dual_p0_dofs == fine.p1_dofs,
            "dual-P0 and P1 spaces must remain square");
    require(fine.pairing_pivot_ratio > 1.0e-12 &&
                fine.pairing_inverse_relative <= 1.0e-11,
            "mixed Riesz pairing must remain well conditioned");
    require(fine.strong_system_relative <= 0.025 &&
                fine.deterministic_probe_relative <= 0.025,
            "full JCP strong action must match the dense direct oracle");
    require(fine.strong_system_relative <= 1.05 * coarse.strong_system_relative,
            "complete mixed JCP parity must not worsen under h-refinement");
    std::cout << "calderon_c5c4_complete_mixed_smoke=pass system_rel="
              << fine.strong_system_relative << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c4_complete_mixed_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
