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
    const auto auxiliary = build_implicit_vdw_auxiliary_surface(sphere, 1.00);
    require(auxiliary.status == "pass", "auxiliary sphere must be watertight");

    ScreenedHypersingularConfig config;
    config.kappa = 0.125;
    config.adjacent_pair_subdivision_depth = 2;
    const auto result = assemble_jcp_screened_hypersingular_weak_form(
        auxiliary.mesh, auxiliary.trace_spaces, config);
    std::cout << "calderon_c5c2_screened_diagnostic status=" << result.status
              << " p1=" << result.p1_dofs
              << " triangles=" << result.triangles
              << " separated_evaluations=" << result.separated_source_evaluations
              << " weak_pre_sym=" << result.weak_pre_symmetry_relative
              << " weak_sym=" << result.weak_symmetry_relative
              << " negative_rel=" << result.negative_sign_separated_relative
              << " positive_rel=" << result.positive_sign_separated_relative
              << " selected_rel=" << result.selected_sign_separated_relative
              << " sign=" << result.weak_form_sign << '\n';
    require(result.p1_dofs > 0 && result.triangles > 0,
            "screened weak form needs a closed P1 mesh");
    require(result.separated_source_evaluations > 0,
            "screened weak form needs non-adjacent calibration pairs");
    require(result.status == "pass",
            "screened weak K4 form must pass the direct separated-pair gate");
    std::cout << "calderon_c5c2_screened_smoke=pass weak_sym="
              << result.weak_symmetry_relative
              << " selected_rel=" << result.selected_sign_separated_relative
              << " sign=" << result.weak_form_sign << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c2_screened_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
