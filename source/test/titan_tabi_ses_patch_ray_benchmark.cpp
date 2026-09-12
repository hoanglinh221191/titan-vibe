#include "tabi/tabi_calderon_aux_surface.h"
#include "titan_executor.h"

#include <bit>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::calderon::build_analytic_ses_surface;

std::string trimmed(std::string value) {
  const auto first = value.find_first_not_of(' ');
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(' ');
  return value.substr(first, last - first + 1U);
}

double pdb_radius(std::string element) {
  element = trimmed(std::move(element));
  if (element == "H")
    return 1.20;
  if (element == "N")
    return 1.55;
  if (element == "O")
    return 1.52;
  if (element == "S" || element == "P")
    return 1.80;
  return 1.70;
}

std::vector<Atom> read_atoms(const std::string &path,
                             std::size_t max_atoms) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open PQR: " + path);
  std::vector<Atom> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM")))
      continue;
    Atom atom;
    const bool fixed_pdb = line.size() >= 54U &&
                           !path.ends_with(".pqr") &&
                           !path.ends_with(".PQR");
    if (fixed_pdb) {
      atom.x = std::stod(line.substr(30U, 8U));
      atom.y = std::stod(line.substr(38U, 8U));
      atom.z = std::stod(line.substr(46U, 8U));
      atom.q = 0.0;
      atom.radius = pdb_radius(line.size() >= 78U ? line.substr(76U, 2U)
                                                  : std::string{});
      atoms.push_back(atom);
      if (max_atoms != 0U && atoms.size() >= max_atoms)
        break;
      continue;
    }
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while (stream >> field)
      fields.push_back(field);
    if (fields.size() < 9U)
      continue;
    const std::size_t n = fields.size();
    atom.x = std::stod(fields[n - 5U]);
    atom.y = std::stod(fields[n - 4U]);
    atom.z = std::stod(fields[n - 3U]);
    atom.q = std::stod(fields[n - 2U]);
    atom.radius = std::stod(fields[n - 1U]);
    atoms.push_back(atom);
    if (max_atoms != 0U && atoms.size() >= max_atoms)
      break;
  }
  if (atoms.empty())
    throw std::runtime_error("PQR contains no atoms: " + path);
  return atoms;
}

std::uint64_t mesh_checksum(const TitanPBE::Tabi::SurfaceMesh &mesh) noexcept {
  std::uint64_t hash = 14695981039346656037ULL;
  const auto mix_word = [&](std::uint64_t word) {
    for (unsigned byte = 0U; byte < 8U; ++byte) {
      hash ^= (word >> (8U * byte)) & 0xffULL;
      hash *= 1099511628211ULL;
    }
  };
  const auto mix_double = [&](double value) {
    mix_word(std::bit_cast<std::uint64_t>(value));
  };
  mix_word(static_cast<std::uint64_t>(mesh.vertices.size()));
  mix_word(static_cast<std::uint64_t>(mesh.triangles.size()));
  for (const auto &vertex : mesh.vertices) {
    mix_double(vertex.x);
    mix_double(vertex.y);
    mix_double(vertex.z);
    mix_double(vertex.nx);
    mix_double(vertex.ny);
    mix_double(vertex.nz);
    mix_word(static_cast<std::uint64_t>(vertex.atom_index));
  }
  for (const auto &triangle : mesh.triangles) {
    mix_word(static_cast<std::uint64_t>(triangle.v0));
    mix_word(static_cast<std::uint64_t>(triangle.v1));
    mix_word(static_cast<std::uint64_t>(triangle.v2));
    mix_word(static_cast<std::uint64_t>(triangle.atom_index));
    mix_double(triangle.cx);
    mix_double(triangle.cy);
    mix_double(triangle.cz);
    mix_double(triangle.nx);
    mix_double(triangle.ny);
    mix_double(triangle.nz);
    mix_double(triangle.area);
  }
  return hash;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2 || argc > 6) {
      std::cerr << "usage: benchmark <input.pqr|pdb> [engine] [spacing] "
                   "[threads] [max_atoms]\n";
      return EXIT_FAILURE;
    }
    const std::string engine =
        argc >= 3 ? argv[2] : "patch-ray-experimental";
    const double spacing = argc >= 4 ? std::stod(argv[3]) : 0.75;
    const std::size_t threads =
        argc >= 5 ? static_cast<std::size_t>(std::stoull(argv[4])) : 4U;
    const std::size_t max_atoms =
        argc >= 6 ? static_cast<std::size_t>(std::stoull(argv[5])) : 0U;
    const std::vector<Atom> atoms = read_atoms(argv[1], max_atoms);
    TitanExecutor executor(threads);
    const auto result = build_analytic_ses_surface(
        atoms, spacing, 1.4, &executor, 100000000U, true, engine);
    std::cout << std::setprecision(12)
              << "engine,topology_source,status,atoms,threads,spacing,grid_vertices,panels,"
                 "build_ms,patch_ms,ray_ms,topology_ms,boundary_magnitude_ms,"
                 "polygonize_ms,repair_ms,current_pre_smooth_repair_ms,"
                 "current_component_filter_ms,current_smooth_analytic_slivers_ms,"
                 "current_post_smooth_repair_ms,"
                 "current_final_orient_trace_space_ms,"
                 "projection_ms,local_repair_ms,orientation_certificate_ms,"
                 "self_intersection_audit_ms,sampled_analytic_audit_ms,"
                 "initial_collapse_ms,stable_sliver_ms,post_projection_collapse_ms,"
                 "local_sliver_ms,full_projection_ms,seam_weld_ms,"
                 "certificate_cleanup_ms,orientation_repair_ms,wide_repair_ms,"
                 "projection_index_bytes,repair_scratch_bytes,"
                 "projected_vertices,collapsed_edges,flipped_edges,"
                 "projection_queries,projection_hits,rejected_projection_vertices,"
                 "rejected_projection_rounds,"
                 "projection_orientation_reversals,"
                 "coupled_projection_components_attempted,"
                 "coupled_projection_components_accepted,"
                 "coupled_projection_vertices_rescued,"
                 "coupled_projection_max_component_size,"
                 "coupled_projection_unresolved_vertices,"
                 "orientation_repair_rounds,"
                 "orientation_repair_faces,wide_repair_attempts,"
                 "topology_certificate_accepts,topology_fallback_audits,"
                 "welded_vertices,"
                 "odd_rows,degenerate_rows_repaired,uncertified_intersection_queries,"
                 "uncertified_rows_repaired,unresolved_uncertified_rows,"
                 "events,event_duplicates,event_bytes,event_peak_temp_bytes,topology_bytes,"
                 "mixed_blocks,active_cells,active_cell_bytes,geometric_certificate,"
                 "parity_confident,"
                 "parity_local_repair,morphology_fallback,parity_fallback,"
                 "surface_components,euler_characteristic,genus_sum,"
                 "filled_cavity_regions,filled_cavity_vertices,"
                 "poor_triangles_before,poor_triangles_after,"
                 "triangle_quality_min_before,triangle_quality_min_after,"
                 "max_edge_ratio_before,max_edge_ratio_after,"
                 "analytic_audit_status,analytic_audit_samples,"
                 "analytic_audit_hits,analytic_audit_misses,"
                 "projection_distance_rms,projection_distance_max,"
                 "patch_residual_rms,patch_residual_max,"
                 "normal_angle_rms_deg,normal_angle_max_deg,"
                 "sampled_hausdorff_proxy,analytic_vertex_projection_misses,"
                 "analytic_vertex_residual_max,analytic_vertex_certificate,"
                 "analytic_orientation_status,"
                 "analytic_orientation_certified_faces,analytic_orientation_centroid_queries,"
                 "analytic_orientation_uncertified_faces,analytic_inverted_faces,"
                 "analytic_self_intersections,analytic_self_intersection_first,"
                 "analytic_self_intersection_second,analytic_orientation_min_cosine,"
                 "analytic_normal_quality_status,analytic_orientation_smooth_faces,"
                 "analytic_orientation_seam_faces,"
                 "analytic_normal_cosine_below_025_faces,"
                 "analytic_normal_cosine_below_050_faces,"
                 "analytic_normal_area_weighted_rms_deg,analytic_normal_p99_deg,"
                 "analytic_smooth_normal_p99_deg,analytic_seam_normal_p99_deg,"
                 "analytic_normal_area_below_cosine_025_fraction,"
                 "analytic_normal_area_below_cosine_050_fraction,"
                 "area,mesh_checksum\n"
              << engine << ',' << result.topology_source << ','
              << result.status << ',' << atoms.size() << ','
              << threads << ',' << spacing << ',' << result.grid_vertices
              << ',' << result.mesh.panels.size() << ',' << result.build_ms
              << ',' << result.patch_build_ms << ',' << result.ray_build_ms
              << ',' << result.topology_ms << ','
              << result.boundary_magnitude_ms << ',' << result.polygonize_ms
              << ',' << result.repair_ms << ','
              << result.current_pre_smooth_repair_ms << ','
              << result.current_component_filter_ms << ','
              << result.current_smooth_analytic_slivers_ms << ','
              << result.current_post_smooth_repair_ms << ','
              << result.current_final_orient_trace_space_ms << ','
              << result.projection_ms << ','
              << result.local_repair_ms << ','
              << result.orientation_certificate_ms << ','
              << result.self_intersection_audit_ms << ','
              << result.sampled_analytic_audit_ms << ','
              << result.initial_collapse_ms << ','
              << result.stable_sliver_ms << ','
              << result.post_projection_collapse_ms << ','
              << result.local_sliver_ms << ','
              << result.full_projection_ms << ',' << result.seam_weld_ms
              << ',' << result.certificate_cleanup_ms << ','
              << result.orientation_repair_ms << ','
              << result.wide_repair_ms << ','
              << result.projection_index_bytes << ','
              << result.repair_scratch_bytes << ','
              << result.projected_vertices << ',' << result.collapsed_edges
              << ',' << result.flipped_edges << ','
              << result.projection_queries << ',' << result.projection_hits
              << ',' << result.rejected_projection_vertices
              << ',' << result.rejected_projection_rounds << ','
              << result.projection_orientation_reversals << ','
              << result.coupled_projection_components_attempted << ','
              << result.coupled_projection_components_accepted << ','
              << result.coupled_projection_vertices_rescued << ','
              << result.coupled_projection_max_component_size << ','
              << result.coupled_projection_unresolved_vertices << ','
              << result.orientation_repair_rounds << ','
              << result.orientation_repair_faces << ','
              << result.wide_repair_attempts << ','
              << result.topology_certificate_accepts << ','
              << result.topology_fallback_audits << ','
              << result.welded_vertices << ','
              << result.odd_ray_rows << ','
              << result.patch_ray_degenerate_rows_repaired << ','
              << result.patch_ray_uncertified_intersection_queries << ','
              << result.patch_ray_uncertified_rows_repaired << ','
              << result.patch_ray_unresolved_uncertified_rows << ','
              << result.patch_ray_events << ',' << result.patch_ray_duplicates
              << ',' << result.patch_ray_retained_bytes << ','
              << result.patch_ray_peak_temporary_bytes << ','
              << result.topology_retained_bytes << ','
              << result.topology_mixed_blocks << ',' << result.active_cells
              << ',' << result.active_cell_bytes << ','
              << result.geometric_certificate_vertices << ','
              << result.parity_confident_vertices << ','
              << result.parity_local_repair_vertices << ','
              << result.morphology_fallback_vertices << ','
              << result.parity_fallback_vertices << ','
              << result.trace_spaces.quality.connected_component_count << ','
              << result.trace_spaces.quality.euler_characteristic << ','
              << result.trace_spaces.quality.genus_sum << ','
              << result.filled_cavity_regions << ','
              << result.filled_cavity_vertices << ','
              << result.poor_triangles_before << ','
              << result.poor_triangles_after << ','
              << result.triangle_quality_min_before << ','
              << result.triangle_quality_min_after << ','
              << result.max_edge_ratio_before << ','
              << result.max_edge_ratio_after << ','
              << result.analytic_audit_status << ','
              << result.analytic_audit_samples << ','
              << result.analytic_audit_hits << ','
              << result.analytic_audit_misses << ','
              << result.analytic_projection_distance_rms << ','
              << result.analytic_projection_distance_max << ','
              << result.analytic_patch_residual_rms << ','
              << result.analytic_patch_residual_max << ','
              << result.analytic_normal_angle_rms_deg << ','
              << result.analytic_normal_angle_max_deg << ','
              << result.sampled_hausdorff_proxy << ','
              << result.analytic_vertex_projection_misses << ','
              << result.analytic_vertex_residual_max << ','
              << (result.analytic_vertex_certificate ? 1 : 0) << ','
              << result.analytic_orientation_status << ','
              << result.analytic_orientation_certified_faces << ','
              << result.analytic_orientation_centroid_queries << ','
              << result.analytic_orientation_uncertified_faces << ','
              << result.analytic_inverted_faces << ','
              << result.analytic_self_intersections << ','
              << result.analytic_self_intersection_first << ','
              << result.analytic_self_intersection_second << ','
              << result.analytic_orientation_min_cosine << ','
              << result.analytic_normal_quality_status << ','
              << result.analytic_orientation_smooth_faces << ','
              << result.analytic_orientation_seam_faces << ','
              << result.analytic_normal_cosine_below_025_faces << ','
              << result.analytic_normal_cosine_below_050_faces << ','
              << result.analytic_normal_area_weighted_rms_deg << ','
              << result.analytic_normal_p99_deg << ','
              << result.analytic_smooth_normal_p99_deg << ','
              << result.analytic_seam_normal_p99_deg << ','
              << result.analytic_normal_area_below_cosine_025_fraction << ','
              << result.analytic_normal_area_below_cosine_050_fraction << ','
              << result.trace_spaces.quality.total_area << ','
              << mesh_checksum(result.mesh) << '\n';
    if (std::getenv("TITAN_SES_DEBUG_SELF_INTERSECTION") != nullptr &&
        result.analytic_self_intersections != 0U &&
        result.analytic_self_intersection_first < result.mesh.triangles.size() &&
        result.analytic_self_intersection_second < result.mesh.triangles.size()) {
      for (const std::size_t triangle_index :
           {result.analytic_self_intersection_first,
            result.analytic_self_intersection_second}) {
        const auto &triangle = result.mesh.triangles[triangle_index];
        std::cerr << "self_intersection_triangle=" << triangle_index;
        for (const std::size_t vertex_index :
             {triangle.v0, triangle.v1, triangle.v2}) {
          const auto &vertex = result.mesh.vertices[vertex_index];
          std::cerr << " v" << vertex_index << "=(" << vertex.x << ','
                    << vertex.y << ',' << vertex.z << ')';
        }
        std::cerr << '\n';
      }
    }
    return result.status == "pass" ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_ray_benchmark=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
