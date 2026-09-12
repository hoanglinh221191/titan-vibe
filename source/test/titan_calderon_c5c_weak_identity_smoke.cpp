#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_weak_form.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::WeakCalderonConfig;
using TitanPBE::Tabi::calderon::assemble_laplace_weak_calderon_identity;
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

    WeakCalderonConfig config;
    config.adjacent_pair_subdivision_depth = 2;
    const auto coarse = assemble_laplace_weak_calderon_identity(
        coarse_auxiliary.mesh, coarse_auxiliary.trace_spaces, config);
    const auto fine = assemble_laplace_weak_calderon_identity(
        fine_auxiliary.mesh, fine_auxiliary.trace_spaces, config);
    std::cout << "calderon_c5c_weak_identity_diagnostic coarse_status="
              << coarse.status << " coarse_p1=" << coarse.p1_dofs
              << " coarse_strong=" << coarse.strong_identity_relative
              << " fine_status=" << fine.status << " fine_p1=" << fine.p1_dofs
              << " fine_strong=" << fine.strong_identity_relative << '\n';
    require(fine.status != "invalid_trace_space" && fine.status != "invalid_mesh" &&
                fine.status != "invalid_triangle_geometry" &&
                fine.status != "assembly_failure" && fine.status != "nonfinite",
            "weak Calderon assembly must be finite");
    require(fine.p1_dofs > coarse.p1_dofs && fine.triangles > coarse.triangles,
            "finer auxiliary mesh must add trace degrees of freedom");
    require(fine.single_layer_symmetry_relative < 1.0e-12,
            "single-layer weak matrix must be symmetric");
    require(fine.hypersingular_symmetry_relative < 1.0e-12,
            "Maue hypersingular weak matrix must be symmetric");
    require(fine.hypersingular_constant_relative < 1.0e-10,
            "weak hypersingular operator must annihilate constants");
    require(coarse.status == "pass" && fine.status == "pass",
            "mass-scaled Calderon identity must pass the synthetic-mesh gate");
    require(fine.strong_identity_relative < coarse.strong_identity_relative,
            "h-refinement must reduce the strong-form identity discrepancy");
    std::cout << "calderon_c5c_weak_identity_smoke=pass coarse_status="
              << coarse.status << " fine_status=" << fine.status
              << " fine_strong_identity=" << fine.strong_identity_relative
              << " fine_weak_identity_mass_scaled="
              << fine.weak_identity_mass_scaled_relative
              << " fine_weak_identity_raw="
              << fine.weak_identity_raw_frobenius_relative << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c_weak_identity_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
