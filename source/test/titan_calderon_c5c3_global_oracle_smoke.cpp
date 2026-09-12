#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_weak_form.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::ScreenedHypersingularConfig;
using TitanPBE::Tabi::calderon::assemble_jcp_screened_hypersingular_weak_form;
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

    ScreenedHypersingularConfig config;
    config.kappa = 0.125;
    config.adjacent_pair_subdivision_depth = 2;
    const auto coarse = assemble_jcp_screened_hypersingular_weak_form(
        coarse_auxiliary.mesh, coarse_auxiliary.trace_spaces, config);
    const auto fine = assemble_jcp_screened_hypersingular_weak_form(
        fine_auxiliary.mesh, fine_auxiliary.trace_spaces, config);
    std::cout << "calderon_c5c3_global_oracle_diagnostic coarse_status="
              << coarse.status << " coarse_p1=" << coarse.p1_dofs
              << " coarse_global=" << coarse.global_direct_weak_relative
              << " fine_status=" << fine.status << " fine_p1=" << fine.p1_dofs
              << " fine_direct_pre_sym=" << fine.direct_pre_symmetry_relative
              << " fine_direct_sym=" << fine.direct_symmetry_relative
              << " fine_global_negative=" << fine.global_direct_negative_sign_relative
              << " fine_global_positive=" << fine.global_direct_positive_sign_relative
              << " fine_global=" << fine.global_direct_weak_relative
              << " weak_sign=" << fine.weak_form_sign << '\n';
    require(coarse.status == "pass" && fine.status == "pass",
            "screened weak assembly must pass C5C.2");
    require(fine.p1_dofs > coarse.p1_dofs,
            "fine auxiliary mesh must add trace degrees of freedom");
    require(fine.direct_symmetry_relative < 1.0e-12,
            "global direct K4 matrix must be symmetric after Galerkin symmetrization");
    require(fine.global_direct_positive_sign_relative <
                fine.global_direct_negative_sign_relative,
            "positive Maue sign must outperform the opposite K4 sign");
    require(fine.global_direct_weak_relative <= 0.25,
            "screened weak K4 must match the global direct action gate");
    require(fine.global_direct_weak_relative <=
                1.05 * coarse.global_direct_weak_relative,
            "global direct/weak discrepancy must not worsen under h-refinement");
    std::cout << "calderon_c5c3_global_oracle_smoke=pass global_rel="
              << fine.global_direct_weak_relative << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c3_global_oracle_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
