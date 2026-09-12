#include "tabi/tabi_calderon_trace_space.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

using TitanPBE::Tabi::SurfaceMesh;
using TitanPBE::Tabi::Triangle;
using TitanPBE::Tabi::Vertex;
using TitanPBE::Tabi::calderon::TraceSpaces;
using TitanPBE::Tabi::calderon::build_closed_surface_trace_spaces;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

SurfaceMesh oriented_tetrahedron() {
  SurfaceMesh mesh;
  mesh.vertices = {
      Vertex{0.0, 0.0, 0.0}, Vertex{1.0, 0.0, 0.0},
      Vertex{0.0, 1.0, 0.0}, Vertex{0.0, 0.0, 1.0}};
  const double half = 0.5;
  const double sloped = std::sqrt(3.0) / 2.0;
  mesh.triangles = {
      Triangle{0, 2, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, half},
      Triangle{0, 1, 3, 0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, half},
      Triangle{0, 3, 2, 0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, half},
      Triangle{1, 2, 3, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, sloped}};
  return mesh;
}

double sum(const std::vector<double> &values) {
  double total = 0.0;
  for (const double value : values) {
    total += value;
  }
  return total;
}

} // namespace

int main() {
  try {
    const SurfaceMesh tetrahedron = oriented_tetrahedron();
    const TraceSpaces closed = build_closed_surface_trace_spaces(tetrahedron);
    require(closed.quality.valid, "closed tetrahedron must be accepted");
    require(closed.quality.edge_count == 6, "tetrahedron must have six edges");
    require(closed.quality.boundary_edge_count == 0, "closed mesh has no boundary");
    require(closed.quality.orientation_mismatch_count == 0,
            "closed mesh must have consistent orientation");
    require(closed.quality.topology_invariants_valid,
            "closed mesh topology invariants must be valid");
    require(closed.quality.connected_component_count == 1,
            "tetrahedron must have one component");
    require(closed.quality.euler_characteristic == 2,
            "tetrahedron Euler characteristic must be two");
    require(closed.quality.genus_sum == 0,
            "tetrahedron genus must be zero");
    require(closed.p0_triangle_mass.size() == 4, "P0 mass count");
    require(closed.p1_vertex_lumped_mass.size() == 4, "P1 mass count");
    require(std::abs(sum(closed.p0_triangle_mass) - closed.quality.total_area) <
                1.0e-14,
            "P0 mass must integrate the surface area");
    require(std::abs(sum(closed.p1_vertex_lumped_mass) - closed.quality.total_area) <
                1.0e-14,
            "lumped P1 mass must integrate the surface area");

    SurfaceMesh disconnected = tetrahedron;
    const std::size_t offset = disconnected.vertices.size();
    for (Vertex vertex : tetrahedron.vertices) {
      vertex.x += 3.0;
      disconnected.vertices.push_back(vertex);
    }
    for (Triangle triangle : tetrahedron.triangles) {
      triangle.v0 += offset;
      triangle.v1 += offset;
      triangle.v2 += offset;
      disconnected.triangles.push_back(triangle);
    }
    const TraceSpaces two_components =
        build_closed_surface_trace_spaces(disconnected);
    require(two_components.quality.valid,
            "two closed components must be accepted");
    require(two_components.quality.connected_component_count == 2,
            "disconnected tetrahedra must have two components");
    require(two_components.quality.euler_characteristic == 4,
            "two spheres must have Euler characteristic four");
    require(two_components.quality.genus_sum == 0,
            "two tetrahedra have zero total genus");

    SurfaceMesh open = tetrahedron;
    open.triangles.pop_back();
    const TraceSpaces open_spaces = build_closed_surface_trace_spaces(open);
    require(!open_spaces.quality.valid && open_spaces.quality.boundary_edge_count != 0,
            "decimated mesh must be rejected as open");

    SurfaceMesh reversed = tetrahedron;
    std::swap(reversed.triangles[0].v1, reversed.triangles[0].v2);
    const TraceSpaces reversed_spaces = build_closed_surface_trace_spaces(reversed);
    require(!reversed_spaces.quality.valid &&
                reversed_spaces.quality.orientation_mismatch_count != 0,
            "reversed triangle must be rejected");

    SurfaceMesh orphan = tetrahedron;
    orphan.vertices.push_back(Vertex{2.0, 2.0, 2.0});
    const TraceSpaces orphan_spaces = build_closed_surface_trace_spaces(orphan);
    require(!orphan_spaces.quality.valid &&
                orphan_spaces.quality.unreferenced_vertex_count == 1,
            "orphan vertex must be rejected");

    std::cout << "calderon_c5_trace_space_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5_trace_space_smoke=fail reason=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
