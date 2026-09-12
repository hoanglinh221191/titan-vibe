#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_intersections.h"
#include "tabi/tabi_ses_patch_provenance.h"
#include "tabi/tabi_ses_patch_ray.h"
#include "titan_executor.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

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

// Kept byte-for-byte equivalent in semantics to the established patch-ray
// benchmark parser so this diagnostic uses the same CHARMM PQR atoms.
std::vector<Atom> read_atoms(const std::string &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open PQR: " + path);
  std::vector<Atom> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM")))
      continue;
    Atom atom;
    const bool fixed_pdb = line.size() >= 54U && !path.ends_with(".pqr") &&
                           !path.ends_with(".PQR");
    if (fixed_pdb) {
      atom.x = std::stod(line.substr(30U, 8U));
      atom.y = std::stod(line.substr(38U, 8U));
      atom.z = std::stod(line.substr(46U, 8U));
      atom.q = 0.0;
      atom.radius =
          pdb_radius(line.size() >= 78U ? line.substr(76U, 2U) : std::string{});
      atoms.push_back(atom);
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
  }
  if (atoms.empty())
    throw std::runtime_error("PQR contains no atoms: " + path);
  return atoms;
}

std::uint64_t grid_checksum(const PatchRayGrid &grid) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(grid.raw_events);
  mix(grid.duplicate_events);
  mix(grid.uncertified_intersection_queries);
  mix(grid.uncertified_rows_repaired);
  mix(grid.unresolved_uncertified_rows);
  mix(grid.degenerate_rows_repaired);
  for (const AxisRayCrossings &axis : grid.axes) {
    mix(static_cast<std::uint64_t>(axis.axis));
    mix(axis.u_count);
    mix(axis.v_count);
    mix(axis.odd_rows);
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const double coordinate : axis.coordinates)
      mix(std::bit_cast<std::uint64_t>(coordinate));
  }
  return hash;
}

void require_same_grid(const PatchRayGrid &off, const PatchRayGrid &on) {
  require(off.status == on.status && off.raw_events == on.raw_events &&
              off.duplicate_events == on.duplicate_events &&
              off.uncertified_intersection_queries ==
                  on.uncertified_intersection_queries &&
              off.uncertified_rows_repaired == on.uncertified_rows_repaired &&
              off.unresolved_uncertified_rows ==
                  on.unresolved_uncertified_rows &&
              off.degenerate_rows_repaired == on.degenerate_rows_repaired,
          "provenance changed patch-ray scalar contract");
  for (std::size_t axis = 0U; axis < off.axes.size(); ++axis) {
    require(off.axes[axis].axis == on.axes[axis].axis &&
                off.axes[axis].u_count == on.axes[axis].u_count &&
                off.axes[axis].v_count == on.axes[axis].v_count &&
                off.axes[axis].offsets == on.axes[axis].offsets &&
                off.axes[axis].coordinates == on.axes[axis].coordinates &&
                off.axes[axis].odd_rows == on.axes[axis].odd_rows,
            "provenance changed patch-ray CSR contract");
  }
  require(grid_checksum(off) == grid_checksum(on),
          "provenance changed patch-ray checksum");
}

std::size_t
exact_provenance_payload_bytes(const PatchRayProvenance &provenance) {
  std::size_t bytes = 0U;
  for (const AxisPatchProvenance &axis : provenance.axes) {
    std::size_t axis_bytes = 0U;
    require(axis.retained_payload_bytes(axis_bytes),
            "provenance payload byte count overflow");
    require(bytes <= std::numeric_limits<std::size_t>::max() - axis_bytes,
            "provenance payload byte sum overflow");
    bytes += axis_bytes;
  }
  return bytes;
}

struct GridSpec {
  Vec3 lower;
  Vec3 upper;
  std::size_t nx{0U};
  std::size_t ny{0U};
  std::size_t nz{0U};
  std::size_t vertices{0U};
};

// This intentionally mirrors build_analytic_ses_surface: expanded atom
// bounds, a 2*h margin, and ceil(range/h)+1 vertices on each axis.
GridSpec make_surface_grid_spec(const PatchComplex &patches,
                                double probe_radius, double grid_spacing) {
  Vec3 lower{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
  Vec3 upper{-std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
  for (const Atom &atom : patches.atoms) {
    const double radius = atom.radius + probe_radius;
    lower.x = std::min(lower.x, atom.x - radius);
    lower.y = std::min(lower.y, atom.y - radius);
    lower.z = std::min(lower.z, atom.z - radius);
    upper.x = std::max(upper.x, atom.x + radius);
    upper.y = std::max(upper.y, atom.y + radius);
    upper.z = std::max(upper.z, atom.z + radius);
  }
  const double margin = 2.0 * grid_spacing;
  lower.x -= margin;
  lower.y -= margin;
  lower.z -= margin;
  upper.x += margin;
  upper.y += margin;
  upper.z += margin;
  const auto axis_count = [&](double lo, double hi) {
    return static_cast<std::size_t>(std::ceil((hi - lo) / grid_spacing)) + 1U;
  };
  const std::size_t nx = axis_count(lower.x, upper.x);
  const std::size_t ny = axis_count(lower.y, upper.y);
  const std::size_t nz = axis_count(lower.z, upper.z);
  constexpr std::size_t cap = 100000000U;
  require(nx >= 2U && ny >= 2U && nz >= 2U && nx <= cap / ny &&
              nx * ny <= cap / nz,
          "surface grid exceeds the diagnostic cap");
  return {lower, upper, nx, ny, nz, nx * ny * nz};
}

class CsvRow {
public:
  template <class Value> void add(std::string key, const Value &value) {
    std::ostringstream stream;
    stream << std::setprecision(17);
    if constexpr (std::is_same_v<std::remove_cvref_t<Value>, bool>)
      stream << static_cast<int>(value);
    else
      stream << value;
    fields_.emplace_back(std::move(key), stream.str());
  }

  void print() const {
    for (std::size_t index = 0U; index < fields_.size(); ++index) {
      if (index != 0U)
        std::cout << ',';
      std::cout << fields_[index].first;
    }
    std::cout << '\n';
    for (std::size_t index = 0U; index < fields_.size(); ++index) {
      if (index != 0U)
        std::cout << ',';
      std::cout << fields_[index].second;
    }
    std::cout << '\n';
  }

private:
  std::vector<std::pair<std::string, std::string>> fields_;
};

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4) {
      std::cerr << "usage: diagnostic <fixture.pqr> <spacing> <threads>\n";
      return EXIT_FAILURE;
    }
    constexpr double probe_radius = 1.4;
    const double spacing = std::stod(argv[2]);
    const std::size_t threads = static_cast<std::size_t>(std::stoull(argv[3]));
    require(std::isfinite(spacing) && spacing > 0.0 && threads > 0U,
            "invalid spacing or thread count");

    const std::vector<Atom> atoms = read_atoms(argv[1]);
    TitanExecutor executor(threads);

    const auto patch_begin = std::chrono::steady_clock::now();
    const PatchComplex complex =
        build_patch_complex(atoms, probe_radius, &executor, false);
    const double patch_build_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - patch_begin)
            .count();
    require(complex.status == "pass", "analytic patch build failed");
    const GridSpec grid =
        make_surface_grid_spec(complex, probe_radius, spacing);

    const PatchRayGrid off = build_patch_ray_grid(
        complex, grid.lower, spacing, grid.nx, grid.ny, grid.nz, &executor);
    PatchRayBuildOptions provenance_options;
    provenance_options.retain_provenance = true;
    const PatchRayGrid on =
        build_patch_ray_grid(complex, grid.lower, spacing, grid.nx, grid.ny,
                             grid.nz, &executor, provenance_options);
    require_same_grid(off, on);
    require(!off.provenance.enabled && off.provenance.retained_bytes == 0U &&
                on.provenance.enabled,
            "provenance off/on state is inconsistent");
    const std::size_t provenance_payload =
        exact_provenance_payload_bytes(on.provenance);
    require(provenance_payload == on.provenance.retained_bytes,
            "provenance retained payload is not exact");
    require(on.retained_bytes >= off.retained_bytes &&
                on.retained_bytes - off.retained_bytes == provenance_payload,
            "grid retained-byte delta does not equal provenance payload");

    const auto index_begin = std::chrono::steady_clock::now();
    const LocalPatchProjectionIndex projection_index =
        build_local_patch_projection_index(
            complex, LocalPatchProjectionIndexBuildOptions{true});
    const double projection_index_wall_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - index_begin)
            .count();
    const auto oracle_begin = std::chrono::steady_clock::now();
    const AnalyticPatchAdjacencyOracle oracle =
        build_analytic_patch_adjacency_oracle(complex, projection_index);
    const double oracle_wall_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - oracle_begin)
            .count();
    const PatchRayProvenanceAudit audit =
        audit_patch_ray_provenance(on, oracle);
    const AnalyticPatchSeamCertificate &seams = oracle.seam_certificate();

    const auto empty_count = [](const std::vector<std::uint8_t> &flags) {
      return static_cast<std::size_t>(
          std::count(flags.begin(), flags.end(), std::uint8_t{1}));
    };
    CsvRow row;
    row.add("diagnostic_status", "pass");
    row.add("atoms", atoms.size());
    row.add("threads", threads);
    row.add("spacing", spacing);
    row.add("probe_radius", probe_radius);
    row.add("lower_x", grid.lower.x);
    row.add("lower_y", grid.lower.y);
    row.add("lower_z", grid.lower.z);
    row.add("upper_x", grid.upper.x);
    row.add("upper_y", grid.upper.y);
    row.add("upper_z", grid.upper.z);
    row.add("nx", grid.nx);
    row.add("ny", grid.ny);
    row.add("nz", grid.nz);
    row.add("grid_vertices", grid.vertices);
    row.add("patch_status", complex.status);
    row.add("contacts", complex.contacts.size());
    row.add("tori", complex.tori.size());
    row.add("reentrant", complex.reentrant.size());
    row.add("patch_build_ms", patch_build_ms);
    row.add("off_status", off.status);
    row.add("on_status", on.status);
    row.add("exact_grid_contract", true);
    row.add("off_checksum", grid_checksum(off));
    row.add("on_checksum", grid_checksum(on));
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const std::string prefix = axis == 0U ? "x_" : (axis == 1U ? "y_" : "z_");
      row.add(prefix + "rows", on.axes[axis].row_count());
      row.add(prefix + "crossings", on.axes[axis].coordinates.size());
      row.add(prefix + "odd_rows", on.axes[axis].odd_rows);
    }
    row.add("raw_events", on.raw_events);
    row.add("duplicate_events", on.duplicate_events);
    row.add("uncertified_intersection_queries",
            on.uncertified_intersection_queries);
    row.add("uncertified_rows_repaired", on.uncertified_rows_repaired);
    row.add("unresolved_uncertified_rows", on.unresolved_uncertified_rows);
    row.add("degenerate_rows_repaired", on.degenerate_rows_repaired);
    row.add("off_ray_build_ms", off.build_ms);
    row.add("on_ray_build_ms", on.build_ms);
    row.add("on_over_off_ray_ratio", on.build_ms / off.build_ms);
    row.add("off_retained_bytes", off.retained_bytes);
    row.add("on_retained_bytes", on.retained_bytes);
    row.add("off_peak_temporary_bytes", off.peak_temporary_bytes);
    row.add("on_peak_temporary_bytes", on.peak_temporary_bytes);
    row.add("peak_temporary_ratio",
            static_cast<double>(on.peak_temporary_bytes) /
                static_cast<double>(off.peak_temporary_bytes));
    row.add("provenance_status", on.provenance.status);
    row.add("provenance_certified", on.provenance.certified);
    row.add("provenance_retained_bytes", provenance_payload);
    row.add("provenance_peak_temporary_bytes",
            on.provenance.peak_temporary_bytes);
    row.add("provenance_build_ms", on.provenance.build_ms);
    row.add("unique_crossings", on.provenance.unique_crossings);
    row.add("seam_groups", on.provenance.seam_groups);
    row.add("refs_total", on.provenance.refs_total);
    row.add("max_refs_per_crossing", on.provenance.max_refs_per_crossing);
    row.add("roots_total", on.provenance.roots_total);
    row.add("repeated_root_groups", on.provenance.repeated_root_groups);
    row.add("same_primitive_near_pair_unclassified",
            on.provenance.same_primitive_near_pair_unclassified);
    row.add("coordinate_merge_tolerance",
            on.provenance.coordinate_merge_tolerance);
    row.add("near_pair_tolerance", on.provenance.near_pair_tolerance);
    row.add("projection_index_wall_ms", projection_index_wall_ms);
    row.add("projection_index_retained_bytes",
            projection_index.retained_bytes());
    row.add("projection_stamp_bytes_per_worker",
            projection_index.stamp_bytes_per_worker());
    row.add("bounded_domains_certified",
            projection_index.bounded_domains_certified);
    row.add("bvh_invariants_certified",
            projection_index.bvh_invariants_certified);
    row.add("active_pair_cache_build_ms",
            projection_index.active_pair_cache_build_ms);
    row.add("active_pair_cache_bytes",
            projection_index.active_pair_cache_bytes());
    row.add("torus_reentrant_audit_ms",
            projection_index.torus_reentrant_audit_ms);
    row.add("torus_reentrant_pairs_audited",
            projection_index.torus_reentrant_pairs_audited);
    row.add("empty_contact_domains",
            empty_count(projection_index.contact_domain_empty));
    row.add("empty_reentrant_domains",
            empty_count(projection_index.reentrant_domain_empty));
    row.add("oracle_status", oracle.status());
    row.add("oracle_certified", oracle.certified());
    row.add("oracle_wall_ms", oracle_wall_ms);
    row.add("oracle_build_ms", oracle.build_ms());
    row.add("oracle_retained_bytes", oracle.retained_payload_bytes());
    row.add("seam_structure_certified", seams.structure_certified);
    row.add("all_seams_certified", seams.all_seams_certified);
    row.add("certified_torus_contact_seams",
            seams.certified_torus_contact_seams);
    row.add("certified_torus_reentrant_seams",
            seams.certified_torus_reentrant_seams);
    row.add("uncertified_seams", seams.uncertified_seams);
    row.add("seam_certificate_retained_bytes", seams.retained_payload_bytes());
    row.add("audit_status", audit.status);
    row.add("audit_ms", audit.audit_ms);
    row.add("relation_audit_complete", audit.relation_audit_complete);
    row.add("crossing_location_certified", audit.crossing_location_certified);
    row.add("active_cell_mapping_certified",
            audit.active_cell_mapping_certified);
    row.add("active_cell_mapping_status", audit.active_cell_mapping_status);
    row.add("mesh_mapping_certified", audit.mesh_mapping_certified);
    row.add("mesh_mapping_status", audit.mesh_mapping_status);
    row.add("relation_same_surface", audit.relation_counts[0]);
    row.add("relation_curve_adjacent", audit.relation_counts[1]);
    row.add("relation_point_only", audit.relation_counts[2]);
    row.add("relation_separated", audit.relation_counts[3]);
    row.add("relation_uncertified", audit.relation_counts[4]);
    row.add("seam_groups_audited", audit.seam_groups_audited);
    row.add("curve_adjacent_groups", audit.curve_adjacent_groups);
    row.add("point_only_groups", audit.point_only_groups);
    row.add("incompatible_crossing_groups", audit.incompatible_crossing_groups);
    row.add("uncertified_crossing_groups", audit.uncertified_crossing_groups);
    row.add("same_primitive_merged_groups", audit.same_primitive_merged_groups);
    row.add("same_primitive_adjacent_near_pairs",
            audit.same_primitive_adjacent_near_pairs);
    row.add("audit_same_primitive_near_pair_unclassified",
            audit.same_primitive_near_pair_unclassified);
    row.add("adjacent_near_pairs", audit.adjacent_near_pairs);
    row.add("near_relation_same_surface", audit.near_pair_relation_counts[0]);
    row.add("near_relation_curve_adjacent", audit.near_pair_relation_counts[1]);
    row.add("near_relation_point_only", audit.near_pair_relation_counts[2]);
    row.add("near_relation_separated", audit.near_pair_relation_counts[3]);
    row.add("near_relation_uncertified", audit.near_pair_relation_counts[4]);
    row.add("incompatible_near_pairs", audit.incompatible_near_pairs);
    row.add("uncertified_near_pairs", audit.uncertified_near_pairs);
    row.add("incompatible_active_cells", audit.incompatible_active_cells);
    row.add("uncertified_active_cells", audit.uncertified_active_cells);
    row.add("incompatible_mesh_edges", audit.incompatible_mesh_edges);
    row.add("uncertified_mesh_edges", audit.uncertified_mesh_edges);
    row.add("incompatible_mesh_faces", audit.incompatible_mesh_faces);
    row.add("uncertified_mesh_faces", audit.uncertified_mesh_faces);
    row.add("surface_area_status", "not_computed_audit_only");
    row.add("surface_volume_status", "not_computed_audit_only");
    row.print();
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_provenance_1a63_diagnostic=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
