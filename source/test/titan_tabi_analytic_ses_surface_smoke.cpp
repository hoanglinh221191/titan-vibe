#include "tabi/tabi_calderon_aux_surface.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::SurfaceMesh;
using TitanPBE::Tabi::calderon::audit_surface_self_intersections;
using TitanPBE::Tabi::calderon::build_analytic_ses_surface;
#ifdef TITAN_TABI_SES_TEST_HOOKS
using TitanPBE::Tabi::calderon::testing::
    collapse_single_manifold_edge_for_test;
using TitanPBE::Tabi::calderon::testing::
    commit_projection_proposals_for_test;
#endif

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

struct MeshStats {
  double min_area{std::numeric_limits<double>::infinity()};
  double max_edge_ratio{0.0};
  double enclosed_volume{0.0};
};

struct OrientationStats {
  std::size_t inverted_faces{0};
  double minimum_cosine{1.0};
  std::size_t first_inverted{static_cast<std::size_t>(-1)};
};

OrientationStats vertex_normal_orientation_stats(
    const TitanPBE::Tabi::SurfaceMesh &mesh) {
  OrientationStats stats;
  for (std::size_t triangle_index = 0U;
       triangle_index < mesh.triangles.size(); ++triangle_index) {
    const auto &triangle = mesh.triangles[triangle_index];
    const auto &a = mesh.vertices[triangle.v0];
    const auto &b = mesh.vertices[triangle.v1];
    const auto &c = mesh.vertices[triangle.v2];
    const double abx = b.x - a.x;
    const double aby = b.y - a.y;
    const double abz = b.z - a.z;
    const double acx = c.x - a.x;
    const double acy = c.y - a.y;
    const double acz = c.z - a.z;
    const double fnx = aby * acz - abz * acy;
    const double fny = abz * acx - abx * acz;
    const double fnz = abx * acy - aby * acx;
    const double anx = a.nx + b.nx + c.nx;
    const double any = a.ny + b.ny + c.ny;
    const double anz = a.nz + b.nz + c.nz;
    const double denominator =
        std::hypot(fnx, fny, fnz) * std::hypot(anx, any, anz);
    if (!(denominator > 1.0e-30)) {
      ++stats.inverted_faces;
      stats.minimum_cosine = -1.0;
      continue;
    }
    const double cosine =
        std::clamp((fnx * anx + fny * any + fnz * anz) / denominator,
                   -1.0, 1.0);
    stats.minimum_cosine = std::min(stats.minimum_cosine, cosine);
    if (!(cosine > 0.0)) {
      if (stats.first_inverted == static_cast<std::size_t>(-1))
        stats.first_inverted = triangle_index;
      ++stats.inverted_faces;
    }
  }
  return stats;
}

MeshStats mesh_stats(const TitanPBE::Tabi::SurfaceMesh &mesh) {
  MeshStats stats;
  double signed_volume_six = 0.0;
  for (const auto &triangle : mesh.triangles) {
    stats.min_area = std::min(stats.min_area, triangle.area);
    const auto &a = mesh.vertices[triangle.v0];
    const auto &b = mesh.vertices[triangle.v1];
    const auto &c = mesh.vertices[triangle.v2];
    const auto length = [](const auto &p, const auto &q) {
      return std::hypot(p.x - q.x, p.y - q.y, p.z - q.z);
    };
    const double ab = length(a, b);
    const double bc = length(b, c);
    const double ca = length(c, a);
    const double shortest = std::min({ab, bc, ca});
    if (shortest > 1.0e-14)
      stats.max_edge_ratio =
          std::max(stats.max_edge_ratio, std::max({ab, bc, ca}) / shortest);
    signed_volume_six +=
        a.x * (b.y * c.z - b.z * c.y) +
        a.y * (b.z * c.x - b.x * c.z) +
        a.z * (b.x * c.y - b.y * c.x);
  }
  stats.enclosed_volume = std::abs(signed_volume_six) / 6.0;
  return stats;
}

void check_self_intersection_audit() {
  SurfaceMesh crossing;
  crossing.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {0.5, 0.5, -1.0}, {0.5, 0.5, 1.0}, {1.5, 0.5, 0.0}};
  crossing.triangles = {{0U, 1U, 2U}, {3U, 4U, 5U}};
  const auto intersection =
      audit_surface_self_intersections(crossing, 0.25);
  require(intersection.count == 1U && intersection.first == 0U &&
              intersection.second == 1U,
          "self-intersection audit missed a proper triangle crossing");

  crossing.vertices[3].x += 10.0;
  crossing.vertices[4].x += 10.0;
  crossing.vertices[5].x += 10.0;
  require(audit_surface_self_intersections(crossing, 0.25).count == 0U,
          "self-intersection audit reported disjoint triangles");

  SurfaceMesh shared_vertex_crossing;
  shared_vertex_crossing.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {0.5, 0.5, -1.0}, {0.5, 0.5, 1.0}};
  shared_vertex_crossing.triangles = {{0U, 1U, 2U}, {0U, 3U, 4U}};
  require(audit_surface_self_intersections(shared_vertex_crossing, 0.25)
              .count == 1U,
          "self-intersection audit missed a crossing beyond a shared vertex");

  SurfaceMesh shared_vertex_valid;
  shared_vertex_valid.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {-1.0, 0.0, 1.0}, {0.0, -1.0, 1.0}};
  shared_vertex_valid.triangles = {{0U, 1U, 2U}, {0U, 3U, 4U}};
  require(audit_surface_self_intersections(shared_vertex_valid, 0.25)
              .count == 0U,
          "self-intersection audit rejected a valid shared vertex");

  SurfaceMesh shared_edge_valid;
  shared_edge_valid.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {0.0, -2.0, 0.0}};
  shared_edge_valid.triangles = {{0U, 1U, 2U}, {1U, 0U, 3U}};
  require(audit_surface_self_intersections(shared_edge_valid, 0.25)
              .count == 0U,
          "self-intersection audit rejected a valid coplanar shared edge");

  SurfaceMesh shared_edge_overlap;
  shared_edge_overlap.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {0.5, 0.5, 0.0}};
  shared_edge_overlap.triangles = {{0U, 1U, 2U}, {1U, 0U, 3U}};
  require(audit_surface_self_intersections(shared_edge_overlap, 0.25)
              .count == 1U,
          "self-intersection audit missed a coplanar shared-edge overlap");

  SurfaceMesh duplicate;
  duplicate.vertices = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}};
  duplicate.triangles = {{0U, 1U, 2U}, {0U, 2U, 1U}};
  require(audit_surface_self_intersections(duplicate, 0.25).count == 1U,
          "self-intersection audit missed a duplicate triangle");
}

#ifdef TITAN_TABI_SES_TEST_HOOKS
SurfaceMesh collapse_orientation_fixture(double keep_x) {
  SurfaceMesh mesh;
  mesh.vertices = {
      {keep_x, 0.0, 0.0}, // kept endpoint
      {1.0, 0.0, 0.0},   // removed endpoint
      {0.0, 0.0, 1.0},   // first link vertex
      {0.0, 0.0, -1.0},  // second link vertex
      {0.0, 1.0, 0.0},   // surviving face
      {0.0, -1.0, 0.0}};
  mesh.triangles = {
      {0U, 1U, 2U}, {1U, 0U, 3U}, {1U, 4U, 5U}};
  mesh.panels.resize(mesh.triangles.size());
  return mesh;
}

void check_collapse_orientation_guard() {
  const SurfaceMesh reversing = collapse_orientation_fixture(-1.0);
  require(!collapse_single_manifold_edge_for_test(
              reversing, 0U, 1U, false),
          "orientation-reversing manifold collapse was accepted");

  const SurfaceMesh safe = collapse_orientation_fixture(0.75);
  require(collapse_single_manifold_edge_for_test(safe, 0U, 1U, false),
          "orientation-preserving manifold collapse was rejected");
}

SurfaceMesh projection_transaction_fixture() {
  SurfaceMesh mesh;
  mesh.vertices = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}};
  mesh.triangles = {{0U, 1U, 2U}};
  mesh.panels.resize(1U);
  return mesh;
}

void check_coupled_projection_transaction() {
  SurfaceMesh coupled = projection_transaction_fixture();
  std::vector<TitanPBE::Tabi::Vertex> coupled_proposals = coupled.vertices;
  // Moving either endpoint alone reverses the xy-oriented face:
  // A'=(0.8,0.8) gives signed area 1-A'x-A'y < 0 and B'x=-0.2 < 0.
  // Moving both together gives signed area 0.76 > 0, so the deterministic
  // coupled-component transaction must accept the pair atomically.
  coupled_proposals[0].x = 0.8;
  coupled_proposals[0].y = 0.8;
  coupled_proposals[1].x = -0.2;
  coupled_proposals[1].y = 2.0;
  require(commit_projection_proposals_for_test(coupled, coupled_proposals),
          "orientation-preserving coupled projection was rejected");
  require(std::abs(coupled.vertices[0].x - 0.8) < 1.0e-15 &&
              std::abs(coupled.vertices[1].y - 2.0) < 1.0e-15,
          "coupled projection did not commit both endpoints");

  SurfaceMesh reversing = projection_transaction_fixture();
  std::vector<TitanPBE::Tabi::Vertex> reversing_proposals =
      reversing.vertices;
  reversing_proposals[0].x = 0.8;
  reversing_proposals[0].y = 0.8;
  reversing_proposals[1].x = -0.2;
  reversing_proposals[1].y = 0.0;
  require(
      !commit_projection_proposals_for_test(reversing, reversing_proposals),
      "orientation-reversing coupled projection was accepted");
  require(std::abs(reversing.vertices[0].x) < 1.0e-15 &&
              std::abs(reversing.vertices[1].x - 1.0) < 1.0e-15,
          "rejected coupled projection mutated an endpoint");
}
#endif

} // namespace

int main() {
  try {
    check_self_intersection_audit();
#ifdef TITAN_TABI_SES_TEST_HOOKS
    check_collapse_orientation_guard();
    check_coupled_projection_transaction();
#endif
    const std::vector<Atom> sphere = {{0.0, 0.0, 0.0, 0.0, 1.5}};
    const auto coarse = build_analytic_ses_surface(sphere, 0.40, 1.4);
    const auto one = build_analytic_ses_surface(sphere, 0.25, 1.4);
    require(one.status == "pass", "analytic sphere must be closed");
    require(one.trace_spaces.quality.boundary_edge_count == 0,
            "analytic sphere has boundary edges");
    const double exact_area = 4.0 * std::acos(-1.0) * 1.5 * 1.5;
    require(std::abs(one.trace_spaces.quality.total_area - exact_area) /
                    exact_area <
                0.08,
            "analytic sphere area must converge to VDW sphere area");
    require(std::abs(one.trace_spaces.quality.total_area -
                     coarse.trace_spaces.quality.total_area) /
                    exact_area <
                0.03,
            "sphere area must remain stable under grid refinement");

    const std::vector<Atom> pair = {{-1.47, 0.0, 0.0, 0.0, 1.0},
                                    {1.46, 0.0, 0.0, 0.0, 1.0}};
    const auto two = build_analytic_ses_surface(pair, 0.20, 1.0);
    require(two.status == "pass", "analytic two-atom SES must be closed");
    require(two.trace_spaces.quality.nonmanifold_edge_count == 0,
            "analytic two-atom SES must be manifold");
    const std::vector<Atom> translated_pair = {
        {pair[0].x + 1.0e6, pair[0].y - 17.0, pair[0].z + 29.0, 0.0,
         pair[0].radius},
        {pair[1].x + 1.0e6, pair[1].y - 17.0, pair[1].z + 29.0, 0.0,
         pair[1].radius}};
    const auto translated_two =
        build_analytic_ses_surface(translated_pair, 0.20, 1.0);
    const double current_area_relative =
        std::abs(translated_two.trace_spaces.quality.total_area -
                 two.trace_spaces.quality.total_area) /
        two.trace_spaces.quality.total_area;
    const double current_panel_relative =
        std::abs(static_cast<double>(translated_two.mesh.panels.size()) -
                 static_cast<double>(two.mesh.panels.size())) /
        static_cast<double>(two.mesh.panels.size());
    if (translated_two.status != "pass" || translated_two.odd_ray_rows != 0U ||
        current_area_relative >= 1.0e-4 || current_panel_relative >= 0.01) {
      std::cerr << "current_translation_mismatch status="
                << translated_two.status << " panels="
                << two.mesh.panels.size() << '/'
                << translated_two.mesh.panels.size() << " area="
                << two.trace_spaces.quality.total_area << '/'
                << translated_two.trace_spaces.quality.total_area
                << " odd_rows=" << two.odd_ray_rows << '/'
                << translated_two.odd_ray_rows << " area_rel="
                << current_area_relative << " panel_rel="
                << current_panel_relative << '\n';
    }
    require(translated_two.status == "pass" &&
                translated_two.odd_ray_rows == 0U &&
                current_area_relative < 1.0e-4 &&
                current_panel_relative < 0.01,
            "current local torus root path lost translation stability");

    const double root3 = std::sqrt(3.0);
    const std::vector<Atom> triple = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                      {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                      {0.0, root3, 0.0, 0.0, 1.0}};
    const auto three = build_analytic_ses_surface(triple, 0.20, 1.0);
    require(three.status == "pass", "analytic three-atom SES must be closed");
    require(three.trace_spaces.quality.orientation_mismatch_count == 0,
            "analytic reentrant patches must orient consistently");
    const MeshStats stats = mesh_stats(three.mesh);
    require(stats.min_area > 1.0e-8,
            "analytic remeshing retained a micro-area sliver");
    require(stats.max_edge_ratio < 200.0,
            "analytic remeshing retained an extreme edge ratio");

    const auto patch_ray_sphere = build_analytic_ses_surface(
        sphere, 0.30, 1.4, nullptr, 4000000, true,
        "patch-ray-experimental");
    if (patch_ray_sphere.status != "pass")
      std::cerr << "patch_ray_sphere_status=" << patch_ray_sphere.status
                << " audit=" << patch_ray_sphere.analytic_audit_status
                << " vertex_misses="
                << patch_ray_sphere.analytic_vertex_projection_misses
                << " rejected_projection_rounds="
                << patch_ray_sphere.rejected_projection_rounds
                << " flips=" << patch_ray_sphere.flipped_edges
                << " projection_rms="
                << patch_ray_sphere.analytic_projection_distance_rms
                << " residual_rms="
                << patch_ray_sphere.analytic_patch_residual_rms
                << " normal_rms="
                << patch_ray_sphere.analytic_normal_angle_rms_deg
                << " normal_max="
                << patch_ray_sphere.analytic_normal_angle_max_deg
                << " odd_rows=" << patch_ray_sphere.odd_ray_rows << '\n';
    require(patch_ray_sphere.status == "pass",
            "patch-ray analytic sphere must be closed");
    require(patch_ray_sphere.analytic_ses_engine ==
                "patch-ray-experimental",
            "patch-ray engine selection was not retained");
    require(patch_ray_sphere.patch_ray_events > 0,
            "patch-ray engine emitted no compact events");
    require(patch_ray_sphere.topology_retained_bytes > 0,
            "patch-ray engine emitted no compressed topology");
    require(std::abs(patch_ray_sphere.trace_spaces.quality.total_area -
                     exact_area) /
                    exact_area <
                0.10,
            "patch-ray sphere area must converge to VDW sphere area");
    const auto patch_ray_fine_sphere = build_analytic_ses_surface(
        sphere, 0.10, 1.4, nullptr, 4000000, true,
        "patch-ray-experimental");
    const OrientationStats patch_ray_sphere_orientation =
        vertex_normal_orientation_stats(patch_ray_sphere.mesh);
    const OrientationStats patch_ray_fine_sphere_orientation =
        vertex_normal_orientation_stats(patch_ray_fine_sphere.mesh);
    if (patch_ray_sphere_orientation.inverted_faces != 0U ||
        patch_ray_fine_sphere_orientation.inverted_faces != 0U)
      std::cerr << "patch_ray_inverted_faces coarse="
                << patch_ray_sphere_orientation.inverted_faces
                << " fine="
                << patch_ray_fine_sphere_orientation.inverted_faces
                << " coarse_min_cos="
                << patch_ray_sphere_orientation.minimum_cosine
                << " fine_min_cos="
                << patch_ray_fine_sphere_orientation.minimum_cosine
                << " fine_first="
                << patch_ray_fine_sphere_orientation.first_inverted << '\n';
    if (patch_ray_fine_sphere_orientation.first_inverted !=
        static_cast<std::size_t>(-1)) {
      const auto &triangle = patch_ray_fine_sphere.mesh.triangles[
          patch_ray_fine_sphere_orientation.first_inverted];
      for (const std::size_t vertex_index :
           {triangle.v0, triangle.v1, triangle.v2}) {
        const auto &vertex = patch_ray_fine_sphere.mesh.vertices[vertex_index];
        std::cerr << " inverted_vertex id=" << vertex_index << " p=("
                  << vertex.x << ',' << vertex.y << ',' << vertex.z
                  << ") n=(" << vertex.nx << ',' << vertex.ny << ','
                  << vertex.nz << ")\n";
      }
    }
    const MeshStats patch_ray_fine_stats =
        mesh_stats(patch_ray_fine_sphere.mesh);
    const double exact_volume =
        (4.0 / 3.0) * std::acos(-1.0) * 1.5 * 1.5 * 1.5;
    const double coarse_patch_area_error =
        std::abs(patch_ray_sphere.trace_spaces.quality.total_area - exact_area);
    const double fine_patch_area_error = std::abs(
        patch_ray_fine_sphere.trace_spaces.quality.total_area - exact_area);
    if (!(patch_ray_fine_sphere.status == "pass" &&
          fine_patch_area_error / exact_area < 0.005 &&
          fine_patch_area_error < coarse_patch_area_error &&
          std::abs(patch_ray_fine_stats.enclosed_volume - exact_volume) /
                  exact_volume <
              0.01))
      std::cerr << "patch_ray_fine_sphere_status="
                << patch_ray_fine_sphere.status
                << " area_rel=" << fine_patch_area_error / exact_area
                << " volume_rel="
                << std::abs(patch_ray_fine_stats.enclosed_volume -
                            exact_volume) /
                       exact_volume
                << " vertex_misses="
                << patch_ray_fine_sphere.analytic_vertex_projection_misses
                << " flips=" << patch_ray_fine_sphere.flipped_edges
                << " normal_max="
                << patch_ray_fine_sphere.analytic_normal_angle_max_deg
                << '\n';
    require(patch_ray_fine_sphere.status == "pass" &&
                fine_patch_area_error / exact_area < 0.005 &&
                fine_patch_area_error < coarse_patch_area_error &&
                std::abs(patch_ray_fine_stats.enclosed_volume - exact_volume) /
                        exact_volume <
                    0.01,
            "patch-ray sphere area/volume must converge under refinement");
    require((patch_ray_sphere.topology_source == "three_axis_ray_parity" ||
             patch_ray_sphere.topology_source ==
                 "three_axis_ray_parity_symbolic_perturbation" ||
             patch_ray_sphere.topology_source ==
                 "three_axis_ray_parity_morphology_fallback") &&
                patch_ray_sphere.parity_fallback_vertices <
                    patch_ray_sphere.grid_vertices,
            "isolated sphere must use global parity or compact local repair");
    require(patch_ray_sphere.morphology_fallback_vertices == 0U &&
                patch_ray_fine_sphere.morphology_fallback_vertices == 0U,
            "strict patch-ray SES retained a morphology topology fallback");
    require(patch_ray_sphere.trace_spaces.quality.topology_invariants_valid &&
                patch_ray_sphere.trace_spaces.quality.connected_component_count ==
                    1 &&
                patch_ray_sphere.trace_spaces.quality.euler_characteristic ==
                    2 &&
                patch_ray_sphere.trace_spaces.quality.genus_sum == 0,
            "patch-ray sphere topology invariants must describe one sphere");
    require(patch_ray_sphere.analytic_audit_status == "pass" &&
                patch_ray_sphere.analytic_audit_samples > 0 &&
                patch_ray_sphere.analytic_audit_hits > 0 &&
                patch_ray_sphere.analytic_vertex_projection_misses == 0U &&
                patch_ray_sphere.analytic_vertex_residual_max <= 1.0e-8 &&
                patch_ray_sphere.projected_vertices >=
                    patch_ray_sphere.mesh.vertices.size(),
            "patch-ray sphere analytic projection audit must pass");
    require(patch_ray_sphere.analytic_orientation_status == "pass" &&
                patch_ray_fine_sphere.analytic_orientation_status == "pass" &&
                patch_ray_sphere.analytic_inverted_faces == 0U &&
                patch_ray_fine_sphere.analytic_inverted_faces == 0U &&
                patch_ray_sphere.analytic_self_intersections == 0U &&
                patch_ray_fine_sphere.analytic_self_intersections == 0U &&
                patch_ray_sphere.analytic_orientation_certified_faces ==
                    patch_ray_sphere.mesh.triangles.size() &&
                patch_ray_fine_sphere.analytic_orientation_certified_faces ==
                    patch_ray_fine_sphere.mesh.triangles.size(),
            "full-face analytic orientation certificate must cover the mesh");
    require(
        patch_ray_sphere.analytic_normal_quality_status == "pass" &&
            patch_ray_fine_sphere.analytic_normal_quality_status == "pass" &&
            patch_ray_sphere.analytic_orientation_smooth_faces +
                    patch_ray_sphere.analytic_orientation_seam_faces ==
                patch_ray_sphere.mesh.triangles.size() &&
            patch_ray_fine_sphere.analytic_orientation_smooth_faces +
                    patch_ray_fine_sphere.analytic_orientation_seam_faces ==
                patch_ray_fine_sphere.mesh.triangles.size() &&
            std::isfinite(
                patch_ray_sphere.analytic_normal_area_weighted_rms_deg) &&
            std::isfinite(patch_ray_fine_sphere.analytic_normal_p99_deg),
        "full-face seam-aware analytic normal quality gate must pass");
    const std::vector<Atom> shifted_sphere = {
        {1.0e3, -1.0e3, 0.5e3, 0.0, 1.5}};
    const auto patch_ray_shifted = build_analytic_ses_surface(
        shifted_sphere, 0.30, 1.4, nullptr, 4000000, true,
        "patch-ray-experimental");
    require(patch_ray_shifted.status == "pass",
            "translated patch-ray sphere must be closed");
    if (patch_ray_shifted.mesh.panels.size() !=
            patch_ray_sphere.mesh.panels.size() ||
        std::abs(patch_ray_shifted.trace_spaces.quality.total_area -
                 patch_ray_sphere.trace_spaces.quality.total_area) >=
            1.0e-7 * exact_area) {
      double max_vertex_delta = 0.0;
      if (patch_ray_shifted.mesh.vertices.size() ==
          patch_ray_sphere.mesh.vertices.size()) {
        for (std::size_t index = 0;
             index < patch_ray_sphere.mesh.vertices.size(); ++index) {
          const auto &a = patch_ray_sphere.mesh.vertices[index];
          const auto &b = patch_ray_shifted.mesh.vertices[index];
          max_vertex_delta = std::max(
              max_vertex_delta,
              std::hypot((b.x - 1.0e3) - a.x, (b.y + 1.0e3) - a.y,
                         (b.z - 0.5e3) - a.z));
        }
      }
      std::cerr << "translation_mismatch base_panels="
                << patch_ray_sphere.mesh.panels.size()
                << " shifted_panels=" << patch_ray_shifted.mesh.panels.size()
                << " base_area="
                << patch_ray_sphere.trace_spaces.quality.total_area
                << " shifted_area="
                << patch_ray_shifted.trace_spaces.quality.total_area
                << " base_vertices=" << patch_ray_sphere.mesh.vertices.size()
                << " shifted_vertices="
                << patch_ray_shifted.mesh.vertices.size()
                << " max_vertex_delta=" << max_vertex_delta << '\n';
    }
    require(patch_ray_shifted.mesh.panels.size() ==
                patch_ray_sphere.mesh.panels.size() &&
                std::abs(patch_ray_shifted.trace_spaces.quality.total_area -
                         patch_ray_sphere.trace_spaces.quality.total_area) <
                    1.0e-7 * exact_area,
            "patch-ray mesh is not translation invariant");

    const auto patch_ray_pair = build_analytic_ses_surface(
        pair, 0.25, 1.0, nullptr, 4000000, true,
        "patch-ray-experimental");
    if (patch_ray_pair.status != "pass")
      std::cerr << "patch_ray_pair_status=" << patch_ray_pair.status
                << " audit=" << patch_ray_pair.analytic_audit_status
                << " misses=" << patch_ray_pair.analytic_audit_misses
                << " residual_rms="
                << patch_ray_pair.analytic_patch_residual_rms
                << " normal_rms="
                << patch_ray_pair.analytic_normal_angle_rms_deg
                << " odd_rows=" << patch_ray_pair.odd_ray_rows << '\n';
    require(patch_ray_pair.status == "pass",
            "patch-ray two-atom SES must be closed");
    require(patch_ray_pair.trace_spaces.quality.nonmanifold_edge_count == 0,
            "patch-ray two-atom SES must be manifold");
    require(patch_ray_pair.trace_spaces.quality.topology_invariants_valid,
            "patch-ray pair topology invariants must be valid");

    const double exact_grid_spacing = 0.125;
    const double next_grid_spacing =
        std::nextafter(exact_grid_spacing, 1.0);
    const auto patch_ray_exact_grid = build_analytic_ses_surface(
        triple, exact_grid_spacing, 1.0, nullptr, 4000000, true,
        "patch-ray-experimental");
    const auto patch_ray_next_grid = build_analytic_ses_surface(
        triple, next_grid_spacing, 1.0, nullptr, 4000000, true,
        "patch-ray-experimental");
    const OrientationStats exact_grid_vertex_orientation =
        vertex_normal_orientation_stats(patch_ray_exact_grid.mesh);
    if (patch_ray_exact_grid.status != "pass" ||
        patch_ray_next_grid.status != "pass")
      std::cerr << "patch_ray_nextafter_grid exact_status="
                << patch_ray_exact_grid.status << " exact_boundary="
                << patch_ray_exact_grid.trace_spaces.quality.boundary_edge_count
                << " exact_triangles=" << patch_ray_exact_grid.mesh.triangles.size()
                << " exact_inverted="
                << patch_ray_exact_grid.analytic_inverted_faces
                << " exact_uncertified_orientation="
                << patch_ray_exact_grid.analytic_orientation_uncertified_faces
                << " exact_orientation_queries="
                << patch_ray_exact_grid.analytic_orientation_centroid_queries
                << " exact_min_cos="
                << patch_ray_exact_grid.analytic_orientation_min_cosine
                << " exact_vertex_inverted="
                << exact_grid_vertex_orientation.inverted_faces
                << " exact_vertex_first="
                << exact_grid_vertex_orientation.first_inverted
                << " next_status=" << patch_ray_next_grid.status
                << " next_boundary="
                << patch_ray_next_grid.trace_spaces.quality.boundary_edge_count
                << " next_triangles=" << patch_ray_next_grid.mesh.triangles.size()
                << " next_inverted="
                << patch_ray_next_grid.analytic_inverted_faces
                << '\n';
    if (exact_grid_vertex_orientation.first_inverted !=
        static_cast<std::size_t>(-1)) {
      const auto &triangle = patch_ray_exact_grid.mesh.triangles[
          exact_grid_vertex_orientation.first_inverted];
      for (const std::size_t vertex_index :
           {triangle.v0, triangle.v1, triangle.v2}) {
        const auto &vertex = patch_ray_exact_grid.mesh.vertices[vertex_index];
        std::cerr << " exact_inverted_vertex id=" << vertex_index << " p=("
                  << vertex.x << ',' << vertex.y << ',' << vertex.z
                  << ") n=(" << vertex.nx << ',' << vertex.ny << ','
                  << vertex.nz << ")\n";
      }
    }
    require(patch_ray_exact_grid.status == "pass" &&
                patch_ray_next_grid.status == "pass" &&
                patch_ray_exact_grid.trace_spaces.quality.boundary_edge_count ==
                    0U &&
                patch_ray_next_grid.trace_spaces.quality.boundary_edge_count ==
                    0U,
            "patch-ray exact/nextafter grid degeneracy must remain closed");
    const double nextafter_panel_relative =
        std::abs(static_cast<double>(patch_ray_exact_grid.mesh.panels.size()) -
                 static_cast<double>(patch_ray_next_grid.mesh.panels.size())) /
        static_cast<double>(patch_ray_exact_grid.mesh.panels.size());
    const double nextafter_area_relative =
        std::abs(patch_ray_exact_grid.trace_spaces.quality.total_area -
                 patch_ray_next_grid.trace_spaces.quality.total_area) /
        patch_ray_exact_grid.trace_spaces.quality.total_area;
    require(nextafter_panel_relative <= 1.0e-3 &&
                nextafter_area_relative <= 1.0e-8 &&
                patch_ray_exact_grid.morphology_fallback_vertices == 0U &&
                patch_ray_next_grid.morphology_fallback_vertices == 0U &&
                patch_ray_exact_grid.analytic_inverted_faces == 0U &&
                patch_ray_next_grid.analytic_inverted_faces == 0U,
            "one-ULP grid perturbation changed certified SES geometry");
    require(patch_ray_sphere_orientation.inverted_faces == 0U &&
                patch_ray_fine_sphere_orientation.inverted_faces == 0U,
            "patch-ray sphere contains an analytically inverted face");

    std::cout << "tabi_analytic_ses_surface_smoke=pass sphere_panels="
              << one.mesh.panels.size()
              << " pair_panels=" << two.mesh.panels.size()
              << " triple_panels=" << three.mesh.panels.size()
              << " patch_ray_sphere_panels="
              << patch_ray_sphere.mesh.panels.size()
              << " patch_ray_fine_sphere_panels="
              << patch_ray_fine_sphere.mesh.panels.size()
              << " patch_ray_fine_area_rel="
              << fine_patch_area_error / exact_area
              << " patch_ray_fine_volume_rel="
              << std::abs(patch_ray_fine_stats.enclosed_volume - exact_volume) /
                     exact_volume
              << " patch_ray_fine_morphology_fallback="
              << patch_ray_fine_sphere.morphology_fallback_vertices
              << " patch_ray_pair_panels=" << patch_ray_pair.mesh.panels.size()
              << " patch_ray_events=" << patch_ray_pair.patch_ray_events
              << " min_area=" << stats.min_area
              << " max_edge_ratio=" << stats.max_edge_ratio << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_analytic_ses_surface_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
