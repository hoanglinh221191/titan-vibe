#include "tabi/tabi_calderon_aux_surface.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::AuxiliarySurfaceResult;
using TitanPBE::Tabi::calderon::build_implicit_ses_surface;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::size_t connected_components(const AuxiliarySurfaceResult &surface) {
  std::vector<std::size_t> parent(surface.mesh.vertices.size());
  std::iota(parent.begin(), parent.end(), 0U);
  const auto find = [&](std::size_t value, const auto &self) -> std::size_t {
    if (parent[value] != value)
      parent[value] = self(parent[value], self);
    return parent[value];
  };
  const auto join = [&](std::size_t a, std::size_t b) {
    a = find(a, find);
    b = find(b, find);
    if (a != b)
      parent[b] = a;
  };
  for (const auto &triangle : surface.mesh.triangles) {
    join(triangle.v0, triangle.v1);
    join(triangle.v1, triangle.v2);
  }
  std::size_t count = 0;
  for (std::size_t index = 0; index < parent.size(); ++index) {
    count += find(index, find) == index ? 1U : 0U;
  }
  return count;
}

} // namespace

int main() {
  try {
    const std::vector<Atom> sphere = {{0.0, 0.0, 0.0, 0.0, 1.5}};
    const AuxiliarySurfaceResult one =
        build_implicit_vdw_auxiliary_surface(sphere, 0.30);
    require(one.status == "pass", "single sphere must be watertight");
    require(one.trace_spaces.quality.boundary_edge_count == 0,
            "single sphere has no boundary edges");
    const double exact_area = 4.0 * std::acos(-1.0) * 1.5 * 1.5;
    require(std::abs(one.trace_spaces.quality.total_area - exact_area) / exact_area <
                0.15,
            "single sphere area must be reasonable at this grid spacing");

    const std::vector<Atom> overlap = {{-0.9, 0.0, 0.0, 0.0, 1.5},
                                       {0.9, 0.0, 0.0, 0.0, 1.5}};
    const AuxiliarySurfaceResult two =
        build_implicit_vdw_auxiliary_surface(overlap, 0.30);
    require(two.status == "pass", "overlapping spheres must be watertight");
    require(two.trace_spaces.quality.nonmanifold_edge_count == 0,
            "overlapping spheres must be manifold");

    const AuxiliarySurfaceResult ses_sphere =
        build_implicit_ses_surface(sphere, 0.25, 1.4);
    require(ses_sphere.status == "pass", "single-sphere SES must be watertight");
    require(ses_sphere.trace_spaces.quality.boundary_edge_count == 0,
            "single-sphere SES has no boundary edges");
    require(std::abs(ses_sphere.trace_spaces.quality.total_area - exact_area) /
                    exact_area <
                0.20,
            "closing must preserve a convex sphere within grid tolerance");

    const std::vector<Atom> narrow_gap = {
        {-1.3, 0.0, 0.0, 0.0, 1.0}, {1.3, 0.0, 0.0, 0.0, 1.0}};
    const AuxiliarySurfaceResult gap_vdw =
        build_implicit_vdw_auxiliary_surface(narrow_gap, 0.25);
    const AuxiliarySurfaceResult gap_ses =
        build_implicit_ses_surface(narrow_gap, 0.25, 1.0);
    require(gap_vdw.status == "pass" && gap_ses.status == "pass",
            "narrow-gap VDW and SES meshes must be watertight");
    require(connected_components(gap_vdw) == 2,
            "separated VDW spheres must remain separate");
    require(connected_components(gap_ses) == 1,
            "probe closing must bridge a narrow solvent-excluded gap");

    const AuxiliarySurfaceResult repeat =
        build_implicit_vdw_auxiliary_surface(overlap, 0.30);
    require(repeat.mesh.triangles.size() == two.mesh.triangles.size() &&
                repeat.mesh.vertices.size() == two.mesh.vertices.size(),
            "auxiliary mesh generation must be deterministic");
    std::cout << "calderon_c5_aux_surface_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5_aux_surface_smoke=fail reason=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
