#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_transfer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::calderon::PrimaryAuxiliaryTransferConfig;
using TitanPBE::Tabi::calderon::apply_auxiliary_to_primary_adjoint_trace_transfer;
using TitanPBE::Tabi::calderon::apply_primary_to_auxiliary_trace_transfer;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;
using TitanPBE::Tabi::calderon::build_primary_auxiliary_trace_transfer;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::vector<Panel> shifted_panels(std::vector<Panel> panels, double shift) {
  for (Panel &panel : panels) {
    panel.x += shift * panel.nx;
    panel.y += shift * panel.ny;
    panel.z += shift * panel.nz;
  }
  return panels;
}

} // namespace

int main() {
  try {
    const std::vector<Atom> sphere = {{0.0, 0.0, 0.0, 0.0, 1.5}};
    const auto auxiliary = build_implicit_vdw_auxiliary_surface(sphere, 0.75);
    require(auxiliary.status == "pass", "auxiliary sphere must be watertight");
    require(!auxiliary.mesh.panels.empty(),
            "auxiliary sphere must provide panel geometry");

    PrimaryAuxiliaryTransferConfig config;
    config.nearest_primary_panels = 4;
    config.max_primary_to_auxiliary_distance = 0.85;
    config.max_auxiliary_to_primary_distance = 0.85;
    config.min_normal_dot = 0.0;
    config.max_relative_area_mismatch = 1.0e-12;
    const auto matched = build_primary_auxiliary_trace_transfer(
        auxiliary.mesh.panels, auxiliary, config);
    const auto repeat = build_primary_auxiliary_trace_transfer(
        auxiliary.mesh.panels, auxiliary, config);
    std::cout << "calderon_c5c5_transfer_diagnostic matched_status="
              << matched.status << " primary_panels=" << matched.primary_panels
              << " auxiliary_vertices=" << matched.auxiliary_vertices
              << " p2a_max=" << matched.primary_to_auxiliary_max_distance
              << " a2p_max=" << matched.auxiliary_to_primary_max_distance
              << " area_rel=" << matched.relative_area_mismatch
              << " forward_constant=" << matched.forward_constant_relative
              << " adjoint=" << matched.adjoint_relative
              << " reverse_constant="
              << matched.adjoint_reverse_constant_relative << '\n';
    require(matched.status == "pass" && repeat.status == "pass",
            "matched primary and auxiliary surfaces must pass transfer gate");
    require(matched.primary_to_auxiliary_p1.size() == matched.auxiliary_vertices &&
                matched.auxiliary_to_primary_adjoint.size() == matched.primary_panels,
            "transfer stencils must cover both trace spaces");
    require(matched.forward_constant_relative <= 1.0e-13 &&
                matched.adjoint_relative <= 1.0e-13,
            "forward constants and Riesz adjointness must be preserved");
    require(matched.primary_to_auxiliary_p1 == repeat.primary_to_auxiliary_p1 &&
                matched.auxiliary_to_primary_adjoint ==
                    repeat.auxiliary_to_primary_adjoint,
            "transfer construction must be deterministic");

    std::vector<double> primary_values(matched.primary_panels, 0.0);
    for (std::size_t i = 0; i < primary_values.size(); ++i) {
      primary_values[i] = std::sin(0.25 * static_cast<double>(i + 1));
    }
    const std::vector<double> auxiliary_values =
        apply_primary_to_auxiliary_trace_transfer(matched, primary_values);
    const std::vector<double> primary_adjoint =
        apply_auxiliary_to_primary_adjoint_trace_transfer(matched, auxiliary_values);
    require(auxiliary_values.size() == matched.auxiliary_vertices &&
                primary_adjoint.size() == matched.primary_panels,
            "transfer application must preserve expected vector dimensions");

    PrimaryAuxiliaryTransferConfig mismatch_config = config;
    mismatch_config.max_primary_to_auxiliary_distance = 0.25;
    mismatch_config.max_auxiliary_to_primary_distance = 0.25;
    const auto shifted = build_primary_auxiliary_trace_transfer(
        shifted_panels(auxiliary.mesh.panels, 1.0), auxiliary, mismatch_config);
    require(shifted.status == "geometry_coverage_gap",
            "large primary/auxiliary displacement must fail geometry coverage");
    std::cout << "calderon_c5c5_primary_aux_transfer_smoke=pass status="
              << matched.status << " shifted_status=" << shifted.status << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c5_primary_aux_transfer_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
