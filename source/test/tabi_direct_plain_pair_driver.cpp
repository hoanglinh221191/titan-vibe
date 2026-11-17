// Check of accumulate_jcp_near_ordered() against the per-pair dispatcher it
// used before the plain-centroid fast path (TABI-02).
//
// The reference row is evaluate_jcp_near_pair_kernel() for every source in
// order, accumulated as apply_kernel() does, followed by the end-of-row stats
// assignments: exactly the loop accumulate_jcp_near_ordered() had. On flat
// triangles, sphere patches, node patches, curved node patches, point panels
// and a non-finite panel, for coupled/uncoupled rows, kappa 0/>0, every near
// integration mode and every quad-buffer mode, the near-field counters and
// any exception must be identical, and the row sums must agree to 1e-13 of
// the sum of absolute terms. With TITAN_FORCE_SCALAR=1 (no AVX2 plain-pair
// kernel) every row is also identical bit for bit.
//
// Build: link with the engine objects of the same profile (all but Main.o).
// Usage: tabi_direct_plain_pair_driver [points_per_atom]
#include "tabi/tabi_near.h"
#include "tabi/tabi_node_patch.h"
#include "tabi/tabi_node_patch_curved.h"
#include "titan_executor.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace TitanPBE::Tabi;

namespace {

struct RowResult {
  bool threw{false};
  std::string error;
  JcpNearAccum acc;
  JcpNearAccum scale; // sum of |terms| per row, reference only
  TreecodeStats stats;
};

RowResult reference_row(const std::vector<Panel> &panels, std::size_t i,
                        const std::vector<double> &x, double eps_ratio,
                        double kappa, bool coupled,
                        const std::vector<std::size_t> &order,
                        const std::string &near, const std::string &buffer) {
  RowResult out;
  const std::size_t n = panels.size();
  const Panel &target = panels[i];
  if (curved_dual_cell(target)) {
    // Rows served from a prepared curved pair-block cache never reach the
    // pair loop; that branch is unchanged and is not what this checks.
    const auto &cache = *target.node_patch_geometry->curved;
    if (!cache.pair_blocks.empty() && n == cache.cells.size() &&
        target.node_patch_index == i) {
      out.error = "cached_row";
      return out;
    }
  }
  try {
    for (std::size_t pos = 0; pos < n; ++pos) {
      const std::size_t j = order[pos];
      const JcpNearPairKernel k = evaluate_jcp_near_pair_kernel(
          panels[i], i, panels[j], j, eps_ratio, kappa, coupled, near, buffer,
          &out.stats);
      if (!k.valid)
        continue;
      // apply_kernel(), verbatim.
      const double phi_j = x[j];
      if (!coupled) {
        out.acc.phi += k.k2 * phi_j;
        out.scale.phi += std::fabs(k.k2 * phi_j);
        continue;
      }
      const double phin_j = x[n + j];
      out.acc.phi += k.k1 * phin_j + k.k2 * phi_j;
      out.acc.phin += k.k3 * phin_j + k.k4 * phi_j;
      out.scale.phi += std::fabs(k.k1 * phin_j) + std::fabs(k.k2 * phi_j);
      out.scale.phin += std::fabs(k.k3 * phin_j) + std::fabs(k.k4 * phi_j);
    }
  } catch (const std::exception &e) {
    out.threw = true;
    out.error = e.what();
  }
  return out;
}

RowResult fast_row(const std::vector<Panel> &panels, std::size_t i,
                   const std::vector<double> &x, double eps_ratio,
                   double kappa, bool coupled,
                   const std::vector<std::size_t> &order,
                   const std::string &near, const std::string &buffer,
                   const NearSourceCache *cache) {
  RowResult out;
  try {
    out.acc = accumulate_jcp_near_ordered(panels[i], i, panels, x, eps_ratio,
                                          kappa, coupled, order, 0,
                                          panels.size(), near, buffer,
                                          &out.stats, cache);
  } catch (const std::exception &e) {
    out.threw = true;
    out.error = e.what();
  }
  return out;
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool same_stats(const TreecodeStats &a, const TreecodeStats &b) {
  return a.near_self_pairs == b.near_self_pairs &&
         a.near_one_ring_pairs == b.near_one_ring_pairs &&
         a.near_duffy_pairs == b.near_duffy_pairs &&
         a.near_centroid_pairs == b.near_centroid_pairs &&
         a.near_fallbacks == b.near_fallbacks &&
         a.near_geometric_buffer_pairs == b.near_geometric_buffer_pairs &&
         a.near_geometric_buffer_adaptive_pairs ==
             b.near_geometric_buffer_adaptive_pairs &&
         a.near_geometric_buffer_gauss_pairs ==
             b.near_geometric_buffer_gauss_pairs &&
         same_bits(a.duffy_subdivisions_mean, b.duffy_subdivisions_mean) &&
         a.duffy_subdivisions_max == b.duffy_subdivisions_max &&
         a.near_k1_accepted == b.near_k1_accepted &&
         a.near_k1_rejected == b.near_k1_rejected &&
         a.near_k2_accepted == b.near_k2_accepted &&
         a.near_k2_rejected == b.near_k2_rejected &&
         a.near_k3_accepted == b.near_k3_accepted &&
         a.near_k3_rejected == b.near_k3_rejected &&
         a.near_k4_accepted == b.near_k4_accepted &&
         a.near_k4_rejected == b.near_k4_rejected &&
         a.near_integration == b.near_integration &&
         a.quad_buffer == b.quad_buffer &&
         a.near_kernel_policy == b.near_kernel_policy;
}

struct Tally {
  std::size_t rows{0};
  std::size_t rows_threw{0};
  std::size_t rows_cached{0};
  std::size_t plain_pairs{0};
  std::size_t rows_bitwise{0};
  double max_error{0.0};
  std::size_t failures{0};
};

void compare_surface(const std::string &label,
                     const std::vector<Panel> &panels, Tally &tally) {
  const std::size_t n = panels.size();
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::mt19937_64 rng(20261004);
  std::uniform_real_distribution<double> unit(-1.0, 1.0);
  std::vector<double> x(2 * n);
  for (double &v : x)
    v = unit(rng);

  const char *nears[] = {"centroid", "analytic_singular", "duffy_adaptive"};
  const char *buffers[] = {"off", "gauss7", "adaptive"};
  const double kappas[] = {0.0, 0.1235};
  const double eps_ratio = 78.5 / 2.0;
  const std::size_t stride = n > 160 ? n / 160 : 1;
  std::size_t node_patch = 0, triangles = 0, quadrature = 0;
  for (const Panel &p : panels) {
    node_patch += p.node_patch_geometry != nullptr;
    triangles += p.has_triangle_geometry;
    quadrature += p.sphere_patch_quadrature != nullptr;
  }
  double ref_seconds = 0.0, fast_seconds = 0.0, aos_seconds = 0.0;
  const NearSourceCache source_cache = build_near_source_cache(panels);
  const auto seconds_since = [](auto start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
        .count();
  };
  for (const bool coupled : {true, false})
    for (const double kappa : kappas)
      for (const char *near : nears)
        for (const char *buffer : buffers) {
          std::size_t failures = 0;
          for (std::size_t i = 0; i < n; i += stride) {
            auto start = std::chrono::steady_clock::now();
            const RowResult ref = reference_row(panels, i, x, eps_ratio,
                                                kappa, coupled, order, near,
                                                buffer);
            ref_seconds += seconds_since(start);
            if (ref.error == "cached_row") {
              ++tally.rows_cached;
              continue;
            }
            // Both the per-Panel path and the source-cache path the direct
            // operators use.
            for (const NearSourceCache *cache : {(const NearSourceCache *)nullptr,
                                                 &source_cache}) {
            start = std::chrono::steady_clock::now();
            const RowResult got = fast_row(panels, i, x, eps_ratio, kappa,
                                           coupled, order, near, buffer, cache);
            (cache ? fast_seconds : aos_seconds) += seconds_since(start);
            ++tally.rows;
            bool ok = ref.threw == got.threw;
            if (ok && ref.threw) {
              ok = ref.error == got.error;
              ++tally.rows_threw;
            } else if (ok) {
              const double err_phi =
                  std::fabs(ref.acc.phi - got.acc.phi) /
                  std::max(ref.scale.phi, 1e-300);
              const double err_phin =
                  std::fabs(ref.acc.phin - got.acc.phin) /
                  std::max(ref.scale.phin, 1e-300);
              const double err = std::max(ref.scale.phi > 0 ? err_phi : 0.0,
                                          ref.scale.phin > 0 ? err_phin : 0.0);
              tally.max_error = std::max(tally.max_error, err);
              if (same_bits(ref.acc.phi, got.acc.phi) &&
                  same_bits(ref.acc.phin, got.acc.phin))
                ++tally.rows_bitwise;
              ok = err <= 1e-13 && same_stats(ref.stats, got.stats);
              tally.plain_pairs += got.stats.near_centroid_pairs;
            }
            if (!ok) {
              if (failures < 3)
                std::cerr << "MISMATCH " << label << " coupled=" << coupled
                          << " kappa=" << kappa << " near=" << near
                          << " buffer=" << buffer << " row=" << i
                          << " ref=(" << ref.acc.phi << "," << ref.acc.phin
                          << (ref.threw ? " threw: " + ref.error : "")
                          << ") got=(" << got.acc.phi << "," << got.acc.phin
                          << (got.threw ? " threw: " + got.error : "")
                          << ") centroid " << ref.stats.near_centroid_pairs
                          << "/" << got.stats.near_centroid_pairs
                          << (cache ? " (source cache)" : " (panels)") << "\n";
              ++failures;
            }
            }
          }
          tally.failures += failures;
        }
  std::cout << label << ": panels=" << n << " node_patch=" << node_patch
            << " triangle=" << triangles << " sphere_quadrature=" << quadrature
            << " reference_s=" << ref_seconds << " panels_s=" << aos_seconds
            << " cache_s=" << fast_seconds << "\n";
}

std::vector<Atom> cluster() {
  // Three overlapping spheres and one detached sphere: buried seams, close
  // but non-adjacent panels for the geometric buffer, and far pairs.
  return {{0.0, 0.0, 0.0, -0.5, 1.9},
          {1.6, 0.4, 0.0, 0.3, 1.7},
          {0.7, 1.5, 0.6, 0.2, 1.5},
          {6.5, -1.0, 0.8, 0.0, 1.6}};
}

} // namespace

int main(int argc, char **argv) {
  const int points = argc > 1 ? std::stoi(argv[1]) : 96;
  TitanExecutor executor(4);
  Options o;
  o.primary_surface = "sphere-patch";
  o.surface_mode = "vdw";
  o.surface_only = true;
  o.sphere_mesh_mode = "icosphere";
  o.surface_point_count = points;
  const std::vector<Panel> flat = build_surface(cluster(), o, nullptr,
                                                &executor);
  // Curved sphere-patch quadrature is attached only to full (not
  // surface-only) builds with a non-centroid near rule.
  Options curved_options = o;
  curved_options.surface_only = false;
  curved_options.near_integration = "analytic_singular";
  const std::vector<Panel> sphere =
      build_surface(cluster(), curved_options, nullptr, &executor);

  std::vector<Panel> flat_qc = flat;
  for (std::size_t i = 0; i < flat_qc.size(); ++i) {
    Panel &p = flat_qc[i];
    p.has_source_quadrature_center = true;
    p.source_x = p.x + 1.0e-3 * static_cast<double>(i % 7);
    p.source_y = p.y - 2.0e-3 * static_cast<double>(i % 5);
    p.source_z = p.z + 1.5e-3 * static_cast<double>(i % 3);
  }

  std::vector<Panel> points_only = flat;
  for (Panel &p : points_only) {
    p.has_triangle_geometry = false;
    p.v0_index = p.v1_index = p.v2_index = static_cast<std::size_t>(-1);
  }

  const std::vector<Panel> node = build_node_patch_panels(flat);
  // Curved node patches need a closed surface: use one isolated sphere.
  std::vector<Panel> single = build_surface({{0.0, 0.0, 0.0, 0.4, 1.7}}, o,
                                            nullptr, &executor);
  std::vector<Panel> curved = build_node_patch_panels(single);
  bool have_curved = false;
  try {
    have_curved = prepare_curved_node_patches(single, curved, 78.5, 0.1235,
                                              0, executor);
  } catch (const std::exception &e) {
    std::cout << "curved_node_patch: preparation failed: " << e.what()
              << "\n";
  }

  std::vector<Panel> poisoned = flat;
  poisoned[poisoned.size() / 2].x = std::numeric_limits<double>::quiet_NaN();

  Tally tally;
  compare_surface("sphere_patch", sphere, tally);
  compare_surface("flat_triangle", flat, tally);
  compare_surface("flat_quadrature_center", flat_qc, tally);
  compare_surface("point_panel", points_only, tally);
  compare_surface("node_patch", node, tally);
  if (have_curved)
    compare_surface("curved_node_patch", curved, tally);
  else
    std::cout << "curved_node_patch: not prepared, skipped\n";
  compare_surface("non_finite_panel", poisoned, tally);

  std::cout << "rows=" << tally.rows << " rows_threw_identically="
            << tally.rows_threw << " rows_cached_skipped=" << tally.rows_cached
            << " centroid_pairs=" << tally.plain_pairs
            << " rows_bitwise=" << tally.rows_bitwise
            << " max_normalised_error=" << tally.max_error
            << " failures=" << tally.failures << "\n";
  if (tally.failures != 0) {
    std::cout << "FAIL tabi_direct_plain_pair_driver\n";
    return 1;
  }
  std::cout << "PASS tabi_direct_plain_pair_driver\n";
  return 0;
}
