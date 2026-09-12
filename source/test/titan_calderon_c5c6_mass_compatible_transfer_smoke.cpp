#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_transfer.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::calderon::MassCompatibleTransferConfig;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;
using TitanPBE::Tabi::calderon::build_mass_compatible_primary_auxiliary_trace_transfer;

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

    MassCompatibleTransferConfig config;
    config.geometry.nearest_primary_panels = 8;
    config.geometry.max_primary_to_auxiliary_distance = 0.85;
    config.geometry.max_auxiliary_to_primary_distance = 0.85;
    config.geometry.min_normal_dot = 0.0;
    config.geometry.max_relative_area_mismatch = 1.0e-12;
    config.max_balance_iterations = 8000;
    config.balance_relative_tolerance = 1.0e-10;
    const auto matched = build_mass_compatible_primary_auxiliary_trace_transfer(
        auxiliary.mesh.panels, auxiliary, config);
    const auto repeat = build_mass_compatible_primary_auxiliary_trace_transfer(
        auxiliary.mesh.panels, auxiliary, config);
    std::cout << "calderon_c5c6_mass_transfer_diagnostic matched_status="
              << matched.status << " geometry_status=" << matched.geometry_status
              << " primary_panels=" << matched.primary_panels
              << " auxiliary_vertices=" << matched.auxiliary_vertices
              << " support_edges=" << matched.support_edges
              << " support_components=" << matched.support_components
              << " unsupported_before="
              << matched.unsupported_primary_panels_before_repair
              << " repair_edges=" << matched.support_repair_edges
              << " balance_iterations=" << matched.balance_iterations
              << " balance_relative=" << matched.balance_relative
              << " primary_mass_scale=" << matched.primary_mass_global_scale
              << " component_scale_min=" << matched.component_mass_scale_min
              << " component_scale_max=" << matched.component_mass_scale_max
              << " forward_constant=" << matched.forward_constant_relative
              << " reverse_constant=" << matched.reverse_constant_relative
              << " adjoint=" << matched.adjoint_relative << '\n';
    require(matched.status == "pass" && repeat.status == "pass",
            "matched mass-compatible transfer must pass");
    require(matched.balance_relative <= 1.0e-10 &&
                matched.forward_constant_relative <= 1.0e-9 &&
                matched.reverse_constant_relative <= 1.0e-9 &&
                matched.adjoint_relative <= 1.0e-9,
            "balanced transfer must preserve both constants and adjointness");
    require(matched.primary_to_auxiliary_p1 == repeat.primary_to_auxiliary_p1 &&
                matched.auxiliary_to_primary_balanced_adjoint ==
                    repeat.auxiliary_to_primary_balanced_adjoint,
            "mass-compatible transfer must be deterministic");

    MassCompatibleTransferConfig shifted_config = config;
    shifted_config.geometry.max_primary_to_auxiliary_distance = 0.25;
    shifted_config.geometry.max_auxiliary_to_primary_distance = 0.25;
    const auto shifted = build_mass_compatible_primary_auxiliary_trace_transfer(
        shifted_panels(auxiliary.mesh.panels, 1.0), auxiliary, shifted_config);
    require(shifted.status == "real_primary_geometry_gap",
            "mass balancing must not hide a geometry mismatch");
    require(shifted.forward_constant_relative <= 1.0e-9 &&
                shifted.reverse_constant_relative <= 1.0e-9,
            "geometry rejection must remain independent of transfer algebra");
    std::cout << "calderon_c5c6_mass_compatible_transfer_smoke=pass status="
              << matched.status << " shifted_status=" << shifted.status << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c6_mass_compatible_transfer_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
