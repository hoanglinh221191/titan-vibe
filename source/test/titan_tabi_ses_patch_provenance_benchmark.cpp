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
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

std::vector<Atom> make_tiled_triples(std::size_t tile_x,
                                     std::size_t tile_y,
                                     std::size_t tile_z) {
  const double root3 = std::sqrt(3.0);
  const std::vector<Atom> base = {
      {-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
      {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
      {0.0, root3, 0.0, 0.0, 1.0}};
  std::vector<Atom> atoms;
  atoms.reserve(3U * tile_x * tile_y * tile_z);
  constexpr double separation = 9.0;
  for (std::size_t iz = 0U; iz < tile_z; ++iz) {
    for (std::size_t iy = 0U; iy < tile_y; ++iy) {
      for (std::size_t ix = 0U; ix < tile_x; ++ix) {
        for (Atom atom : base) {
          atom.x += separation * static_cast<double>(ix);
          atom.y += separation * static_cast<double>(iy);
          atom.z += separation * static_cast<double>(iz);
          atoms.push_back(atom);
        }
      }
    }
  }
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
  require(off.status == on.status &&
              off.raw_events == on.raw_events &&
              off.duplicate_events == on.duplicate_events &&
              off.uncertified_intersection_queries ==
                  on.uncertified_intersection_queries &&
              off.uncertified_rows_repaired ==
                  on.uncertified_rows_repaired &&
              off.unresolved_uncertified_rows ==
                  on.unresolved_uncertified_rows &&
              off.degenerate_rows_repaired ==
                  on.degenerate_rows_repaired,
          "provenance changed patch-ray scalar contract");
  for (std::size_t axis = 0U; axis < off.axes.size(); ++axis) {
    require(off.axes[axis].axis == on.axes[axis].axis &&
                off.axes[axis].u_count == on.axes[axis].u_count &&
                off.axes[axis].v_count == on.axes[axis].v_count &&
                off.axes[axis].offsets == on.axes[axis].offsets &&
                off.axes[axis].coordinates ==
                    on.axes[axis].coordinates &&
                off.axes[axis].odd_rows == on.axes[axis].odd_rows,
            "provenance changed patch-ray CSR contract");
  }
  require(grid_checksum(off) == grid_checksum(on),
          "provenance changed patch-ray checksum");
}

double median(std::vector<double> values) {
  require(!values.empty(), "empty benchmark timing sample");
  std::sort(values.begin(), values.end());
  return values[values.size() / 2U];
}

std::string samples_string(const std::vector<double> &values) {
  std::string result;
  for (std::size_t index = 0U; index < values.size(); ++index) {
    if (index != 0U)
      result += ';';
    result += std::to_string(values[index]);
  }
  return result;
}

std::size_t exact_provenance_payload_bytes(
    const PatchRayProvenance &provenance) {
  std::size_t bytes = 0U;
  for (const AxisPatchProvenance &axis : provenance.axes) {
    std::size_t axis_bytes = 0U;
    require(axis.retained_payload_bytes(axis_bytes),
            "provenance payload byte count overflow");
    require(bytes <=
                std::numeric_limits<std::size_t>::max() - axis_bytes,
            "provenance payload byte sum overflow");
    bytes += axis_bytes;
  }
  return bytes;
}

struct GridSpec {
  Vec3 lower;
  std::size_t nx{0U};
  std::size_t ny{0U};
  std::size_t nz{0U};
};

GridSpec make_grid_spec(const std::vector<Atom> &atoms,
                        double probe_radius, double spacing) {
  Vec3 lower{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
  Vec3 upper{-std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
  for (const Atom &atom : atoms) {
    const double extent = atom.radius + probe_radius + 2.0 * spacing;
    lower.x = std::min(lower.x, atom.x - extent);
    lower.y = std::min(lower.y, atom.y - extent);
    lower.z = std::min(lower.z, atom.z - extent);
    upper.x = std::max(upper.x, atom.x + extent);
    upper.y = std::max(upper.y, atom.y + extent);
    upper.z = std::max(upper.z, atom.z + extent);
  }
  const auto count = [spacing](double low, double high) {
    return static_cast<std::size_t>(
               std::ceil((high - low) / spacing)) +
           1U;
  };
  return {lower, count(lower.x, upper.x),
          count(lower.y, upper.y), count(lower.z, upper.z)};
}

} // namespace

int main() {
  try {
    constexpr double probe_radius = 1.0;
    constexpr double spacing = 0.25;
    constexpr std::size_t repeats = 7U;
    const std::vector<Atom> atoms =
        make_tiled_triples(3U, 2U, 2U);
    const PatchComplex complex =
        build_patch_complex(atoms, probe_radius, nullptr, false);
    require(complex.status == "pass" && !complex.tori.empty() &&
                !complex.reentrant.empty(),
            "tiled molecular patch fixture failed");

    const GridSpec grid_spec =
        make_grid_spec(atoms, probe_radius, spacing);
    const Vec3 lower = grid_spec.lower;
    const std::size_t nx = grid_spec.nx;
    const std::size_t ny = grid_spec.ny;
    const std::size_t nz = grid_spec.nz;

    TitanExecutor executor(8);
    PatchRayBuildOptions provenance_options;
    provenance_options.retain_provenance = true;

    // Warm both paths twice before recording alternating samples.
    PatchRayGrid warm_off;
    PatchRayGrid warm_on;
    for (std::size_t warmup = 0U; warmup < 2U; ++warmup) {
      warm_off = build_patch_ray_grid(
          complex, lower, spacing, nx, ny, nz, &executor);
      warm_on = build_patch_ray_grid(
          complex, lower, spacing, nx, ny, nz, &executor,
          provenance_options);
      require_same_grid(warm_off, warm_on);
      require(warm_off.status == "pass" &&
                  !warm_off.provenance.enabled &&
                  warm_off.provenance.retained_bytes == 0U &&
                  warm_on.provenance.certified,
              "warm patch-ray build was not certified");
    }

    std::vector<double> off_ms;
    std::vector<double> on_ms;
    off_ms.reserve(repeats);
    on_ms.reserve(repeats);
    PatchRayGrid retained_off;
    PatchRayGrid retained_on;
    for (std::size_t repeat = 0U; repeat < repeats; ++repeat) {
      if ((repeat & 1U) == 0U) {
        retained_off = build_patch_ray_grid(
            complex, lower, spacing, nx, ny, nz, &executor);
        retained_on = build_patch_ray_grid(
            complex, lower, spacing, nx, ny, nz, &executor,
            provenance_options);
      } else {
        retained_on = build_patch_ray_grid(
            complex, lower, spacing, nx, ny, nz, &executor,
            provenance_options);
        retained_off = build_patch_ray_grid(
            complex, lower, spacing, nx, ny, nz, &executor);
      }
      require_same_grid(retained_off, retained_on);
      require(retained_off.status == "pass" &&
                  !retained_off.provenance.enabled &&
                  !retained_off.provenance.certified &&
                  retained_off.provenance.retained_bytes == 0U &&
                  retained_on.provenance.certified &&
                  retained_on.provenance.status == "pass",
              "timed patch-ray build was not certified");
      const std::size_t sample_payload =
          exact_provenance_payload_bytes(retained_on.provenance);
      require(sample_payload ==
                  retained_on.provenance.retained_bytes &&
                  retained_on.retained_bytes >=
                      retained_off.retained_bytes &&
                  retained_on.retained_bytes -
                          retained_off.retained_bytes ==
                      sample_payload &&
                  retained_off.peak_temporary_bytes > 0U &&
                  static_cast<double>(
                      retained_on.peak_temporary_bytes) /
                          static_cast<double>(
                              retained_off.peak_temporary_bytes) <=
                      1.60,
              "timed sample violated provenance memory invariants");
      off_ms.push_back(retained_off.build_ms);
      on_ms.push_back(retained_on.build_ms);
    }

    const auto index_begin = std::chrono::steady_clock::now();
    const LocalPatchProjectionIndex projection_index =
        build_local_patch_projection_index(
            complex, LocalPatchProjectionIndexBuildOptions{true});
    const double projection_index_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - index_begin)
            .count();
    const AnalyticPatchAdjacencyOracle oracle =
        build_analytic_patch_adjacency_oracle(complex,
                                              projection_index);
    require(oracle.certified(),
            "provenance benchmark oracle was not certified");
    const PatchRayProvenanceAudit audit =
        audit_patch_ray_provenance(retained_on, oracle);

    const std::vector<Atom> half_atoms =
        make_tiled_triples(3U, 2U, 1U);
    const PatchComplex half_complex =
        build_patch_complex(half_atoms, probe_radius, nullptr, false);
    const GridSpec half_spec =
        make_grid_spec(half_atoms, probe_radius, spacing);
    const PatchRayGrid half_on = build_patch_ray_grid(
        half_complex, half_spec.lower, spacing, half_spec.nx,
        half_spec.ny, half_spec.nz, &executor,
        provenance_options);
    require(half_complex.status == "pass" &&
                half_on.status == "pass" &&
                half_on.provenance.certified,
            "half-size provenance memory fixture failed");
    const std::size_t half_payload =
        exact_provenance_payload_bytes(half_on.provenance);
    require(half_payload == half_on.provenance.retained_bytes,
            "half-size provenance retained bytes are not exact");

    const std::size_t exact_payload =
        exact_provenance_payload_bytes(retained_on.provenance);
    require(exact_payload == retained_on.provenance.retained_bytes,
            "reported provenance retained bytes are not exact");
    require(!retained_off.provenance.enabled &&
                retained_off.provenance.retained_bytes == 0U &&
                retained_on.retained_bytes >= retained_off.retained_bytes &&
                retained_on.retained_bytes -
                        retained_off.retained_bytes ==
                    exact_payload,
            "on/off retained-byte delta does not equal provenance payload");
    require(retained_off.peak_temporary_bytes > 0U,
            "off-path peak temporary bytes were not recorded");

    const double off_median_ms = median(off_ms);
    const double on_median_ms = median(on_ms);
    const double slowdown = on_median_ms / off_median_ms;
    const double peak_ratio =
        static_cast<double>(retained_on.peak_temporary_bytes) /
        static_cast<double>(retained_off.peak_temporary_bytes);
    const std::size_t half_units =
        half_on.provenance.unique_crossings +
        half_on.provenance.refs_total +
        half_on.provenance.seam_groups;
    const std::size_t full_units =
        retained_on.provenance.unique_crossings +
        retained_on.provenance.refs_total +
        retained_on.provenance.seam_groups;
    require(half_units > 0U && full_units > half_units,
            "provenance linear-memory normalization is empty");
    const double half_bytes_per_unit =
        static_cast<double>(half_payload) /
        static_cast<double>(half_units);
    const double full_bytes_per_unit =
        static_cast<double>(exact_payload) /
        static_cast<double>(full_units);
    const double normalized_memory_ratio =
        std::max(half_bytes_per_unit, full_bytes_per_unit) /
        std::min(half_bytes_per_unit, full_bytes_per_unit);

    std::cout << std::setprecision(12)
              << "status,atoms,contacts,tori,reentrant,nx,ny,nz,"
                 "crossings,seam_groups,refs_total,checksum,"
                 "off_median_ms,on_median_ms,slowdown_ratio,"
                 "off_min_ms,off_max_ms,on_min_ms,on_max_ms,"
                 "off_samples_ms,on_samples_ms,"
                 "off_peak_bytes,on_peak_bytes,peak_ratio,"
                 "off_retained_bytes,on_retained_bytes,"
                 "off_provenance_retained_bytes,"
                 "provenance_retained_bytes,projection_index_ms,"
                 "oracle_build_ms,oracle_retained_bytes,audit_ms,"
                 "half_crossings,half_refs,half_groups,"
                 "half_provenance_bytes,half_bytes_per_unit,"
                 "full_bytes_per_unit,normalized_memory_ratio,"
                 "audit_status,relation_complete,"
                 "crossing_location_certified,"
                 "active_cell_mapping_certified,"
                 "mesh_mapping_certified\n"
              << "pass," << atoms.size() << ','
              << complex.contacts.size() << ',' << complex.tori.size()
              << ',' << complex.reentrant.size() << ',' << nx << ','
              << ny << ',' << nz << ','
              << retained_on.provenance.unique_crossings << ','
              << retained_on.provenance.seam_groups << ','
              << retained_on.provenance.refs_total << ','
              << grid_checksum(retained_on) << ',' << off_median_ms
              << ',' << on_median_ms << ',' << slowdown << ','
              << *std::min_element(off_ms.begin(), off_ms.end()) << ','
              << *std::max_element(off_ms.begin(), off_ms.end()) << ','
              << *std::min_element(on_ms.begin(), on_ms.end()) << ','
              << *std::max_element(on_ms.begin(), on_ms.end()) << ','
              << samples_string(off_ms) << ','
              << samples_string(on_ms) << ','
              << retained_off.peak_temporary_bytes << ','
              << retained_on.peak_temporary_bytes << ',' << peak_ratio
              << ',' << retained_off.retained_bytes << ','
              << retained_on.retained_bytes << ','
              << retained_off.provenance.retained_bytes << ','
              << exact_payload
              << ',' << projection_index_ms << ',' << oracle.build_ms()
              << ',' << oracle.retained_payload_bytes() << ','
              << audit.audit_ms << ','
              << half_on.provenance.unique_crossings << ','
              << half_on.provenance.refs_total << ','
              << half_on.provenance.seam_groups << ','
              << half_payload << ',' << half_bytes_per_unit << ','
              << full_bytes_per_unit << ','
              << normalized_memory_ratio << ',' << audit.status << ','
              << static_cast<int>(audit.relation_audit_complete) << ','
              << static_cast<int>(audit.crossing_location_certified)
              << ','
              << static_cast<int>(
                     audit.active_cell_mapping_certified)
              << ','
              << static_cast<int>(audit.mesh_mapping_certified) << '\n';

    require(slowdown <= 1.20,
            "provenance ray-build slowdown exceeds 20 percent");
    require(peak_ratio <= 1.60,
            "provenance peak temporary ratio exceeds 1.6");
    require(normalized_memory_ratio <= 2.0,
            "provenance retained memory is not empirically linear");
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_provenance_benchmark=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
