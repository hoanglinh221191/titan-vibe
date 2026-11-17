// Check of the AVX2 treecode leaf kernels against the scalar ones.
//
// direct_jcp_ordered_sum() with "avx2-gather" (the "auto" default on AVX2
// CPUs) and direct_jcp_ordered_tile_sum() with "tile-avx2" are compared with
// "scalar" / "tile-scalar" on random panels, index lists of every length
// mod 4, the target in its own list, coincident sources (r2 <= 1e-28),
// tile spans with unaligned starts and exclusion lists, coupled and Poisson
// rows, kappa from 0 to values where exp(-kappa r) takes the libm branch.
// Near-interaction counts must be identical and each row must agree to 1e-13
// of its piece scale: the summed magnitudes of the pieces each term is built
// from, before they cancel (at kappa = 0 the screened and unscreened mixed
// terms cancel exactly, so a row's sum of |terms| is rounding noise there).
// A NaN panel must give NaN on both paths.
// With TITAN_FORCE_SCALAR=1, "auto" must resolve to "scalar".
//
// Build: link with the engine objects of the same profile (all but Main.o).
// Usage: tabi_leaf_jcp_simd_driver [panels] [targets]
#include "tabi/tabi_treecode_internal.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::detail;

namespace {

std::vector<Panel> make_panels(std::size_t n, std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> box(-12.0, 12.0);
  std::normal_distribution<double> gauss(0.0, 1.0);
  std::uniform_real_distribution<double> area(0.02, 0.6);
  std::vector<Panel> panels(n);
  for (std::size_t i = 0; i < n; ++i) {
    Panel &p = panels[i];
    p.x = box(rng);
    p.y = box(rng);
    p.z = box(rng);
    double nx = gauss(rng), ny = gauss(rng), nz = gauss(rng);
    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    p.nx = nx / norm;
    p.ny = ny / norm;
    p.nz = nz / norm;
    p.area = area(rng);
    // A few panels carry a separate source quadrature centre.
    if (i % 17 == 5) {
      p.has_source_quadrature_center = true;
      p.source_x = p.x + 0.01;
      p.source_y = p.y - 0.02;
      p.source_z = p.z + 0.015;
    }
  }
  // Coincident pairs: source centre equal to another panel's target point.
  for (std::size_t i = 1; i + 1 < n; i += 97) {
    panels[i + 1].x = panels[i].x;
    panels[i + 1].y = panels[i].y;
    panels[i + 1].z = panels[i].z;
    panels[i + 1].has_source_quadrature_center = false;
  }
  return panels;
}

struct Row {
  JcpAccum acc;
  std::size_t count{0};
};

Row run_list(const std::vector<Panel> &panels, std::size_t target,
             const TreePanelCache &pc, const TreeUnknownCache &uc,
             double eps_ratio, double kappa, bool coupled,
             const std::vector<std::size_t> &list, const std::string &kernel) {
  Row row;
  row.acc = direct_jcp_ordered_sum(panels[target], target, pc, uc, eps_ratio,
                                   kappa, coupled, list, 0, list.size(),
                                   kernel, row.count);
  return row;
}

JcpAccum piece_scale_list(const Panel &target, std::size_t target_index,
                          const TreePanelCache &pc, const TreeUnknownCache &uc,
                          double eps_ratio, double kappa, bool coupled,
                          const std::vector<std::size_t> &list) {
  const double inv_four_pi = 1.0 / (4.0 * 3.14159265358979323846);
  const double kk = std::max(0.0, kappa);
  const double inv_eps = 1.0 / std::max(1.0e-12, eps_ratio);
  JcpAccum s;
  for (std::size_t j : list) {
    if (j == target_index)
      continue;
    const double ux = target.x - pc.x[j];
    const double uy = target.y - pc.y[j];
    const double uz = target.z - pc.z[j];
    const double r2 = ux * ux + uy * uy + uz * uz;
    if (r2 <= 1.0e-28)
      continue;
    const double r = std::sqrt(r2);
    const double ir = 1.0 / r;
    const double ir3 = ir * ir * ir;
    const double sp = std::abs(pc.nx[j] * ux + pc.ny[j] * uy + pc.nz[j] * uz);
    const double tp =
        std::abs(target.nx * ux + target.ny * uy + target.nz * uz);
    const double a = pc.area[j];
    const double phi = std::abs(uc.phi[j]);
    const double phin = std::abs(uc.phin[j]);
    const double d0s = inv_four_pi * sp * ir3;
    if (!coupled) {
      s.phi += a * std::abs(eps_ratio - 1.0) * d0s * phi;
      continue;
    }
    const double e = std::exp(-kk * r);
    const double opk = 1.0 + kk * r;
    const double g0 = inv_four_pi * ir;
    s.phi += a * ((g0 + e * g0) * phin + (eps_ratio * e * opk * d0s + d0s) * phi);
    const double d0t = inv_four_pi * tp * ir3;
    const double nxu = tp * ir;
    const double nyu = sp * ir;
    const double nxny = std::abs(target.nx * pc.nx[j] + target.ny * pc.ny[j] +
                                 target.nz * pc.nz[j]);
    const double m0 = inv_four_pi * ir3 * (nxny + 3.0 * nxu * nyu);
    const double fsec =
        e * inv_four_pi * (kk * kk * ir + 2.0 * kk * ir * ir + 2.0 * ir3);
    const double fpr = e * opk * inv_four_pi * ir3;
    const double mk = fsec * nxu * nyu + fpr * (nxny + nxu * nyu);
    s.phin += a * ((d0t + inv_eps * e * opk * d0t) * phin + (mk + m0) * phi);
  }
  return s;
}

bool close(double a, double b, double scale, double &worst) {
  if (std::isnan(a) || std::isnan(b))
    return std::isnan(a) && std::isnan(b);
  const double err = std::abs(a - b) / std::max(scale, 1.0e-300);
  if (scale > 0.0)
    worst = std::max(worst, err);
  return std::abs(a - b) <= 1.0e-13 * scale;
}

} // namespace

int main(int argc, char **argv) {
  const std::size_t n = argc > 1 ? std::stoul(argv[1]) : 3000;
  const std::size_t targets = argc > 2 ? std::stoul(argv[2]) : 600;
  const char *forced = std::getenv("TITAN_FORCE_SCALAR");
  const bool force_scalar = forced != nullptr && std::string(forced) == "1";
  const std::string resolved = resolve_leaf_panel_kernel("auto");
  std::cout << "auto resolves to " << resolved << '\n';
  if (force_scalar && resolved != "scalar") {
    std::cout << "FAIL: TITAN_FORCE_SCALAR=1 but auto -> " << resolved << '\n';
    return 1;
  }
  const std::string vec = force_scalar ? "scalar" : "avx2-gather";
  const std::string vec_tile = force_scalar ? "tile-scalar" : "tile-avx2";
  if (!force_scalar && resolved != "avx2-gather") {
    std::cout << "FAIL: AVX2 expected for auto on this CPU\n";
    return 1;
  }

  std::mt19937_64 rng(20261011);
  std::vector<Panel> panels = make_panels(n, rng);
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::shuffle(order.begin(), order.end(), rng);
  std::uniform_real_distribution<double> unit(-1.0, 1.0);
  std::vector<double> unknowns(2 * n);
  for (double &u : unknowns)
    u = unit(rng);

  std::size_t rows = 0, failures = 0;
  double worst_phi = 0.0, worst_phin = 0.0;
  const double kappas[] = {0.0, 0.125, 1.0, 40.0};
  const double eps_ratios[] = {80.0 / 4.0, 1.0 / 20.0};
  for (int pass = 0; pass < 2; ++pass) {
    std::vector<Panel> work = panels;
    if (pass == 1)
      work[n / 3].y = std::numeric_limits<double>::quiet_NaN();
    const TreePanelCache pc = build_tree_panel_cache(work, order);
    for (bool coupled : {true, false}) {
      const TreeUnknownCache uc =
          build_tree_unknown_cache(work, unknowns, coupled, order);
      for (double eps_ratio : eps_ratios) {
        for (double kappa : kappas) {
          std::mt19937_64 list_rng(static_cast<unsigned long long>(
              1000 * pass + 100 * coupled + 10 * kappa + eps_ratio));
          std::uniform_int_distribution<std::size_t> pick(0, n - 1);
          std::uniform_int_distribution<std::size_t> length(0, 71);
          for (std::size_t r = 0; r < targets; ++r) {
            const std::size_t i = pick(list_rng);
            std::vector<std::size_t> list(length(list_rng));
            for (std::size_t &j : list)
              j = pick(list_rng);
            if (!list.empty() && r % 3 == 0)
              list[r % list.size()] = i; // self pair
            if (list.size() > 2 && r % 5 == 0 && i + 1 < n)
              list[1] = i + 1; // possibly coincident with the target
            if (pass == 1 && list.size() > 3 && r % 7 == 0)
              list[3] = n / 3; // NaN source
            const Row a = run_list(work, i, pc, uc, eps_ratio, kappa,
                                   coupled, list, "scalar");
            const Row b =
                run_list(work, i, pc, uc, eps_ratio, kappa, coupled, list, vec);
            const JcpAccum scale = piece_scale_list(
                work[i], i, pc, uc, eps_ratio, kappa, coupled, list);
            ++rows;
            const bool ok =
                a.count == b.count &&
                close(a.acc.phi, b.acc.phi, scale.phi, worst_phi) &&
                close(a.acc.phin, b.acc.phin, scale.phin, worst_phin) &&
                (!force_scalar ||
                 (std::memcmp(&a.acc, &b.acc, sizeof(JcpAccum)) == 0));
            if (!ok && ++failures <= 10) {
              std::cout << "FAIL list pass=" << pass << " coupled=" << coupled
                        << " eps=" << eps_ratio << " kappa=" << kappa
                        << " i=" << i << " len=" << list.size()
                        << " count " << a.count << '/' << b.count << " phi "
                        << a.acc.phi << '/' << b.acc.phi << " phin "
                        << a.acc.phin << '/' << b.acc.phin << '\n';
            }
          }
          // Tile spans in tree order, unaligned starts, with exclusions.
          std::uniform_int_distribution<std::size_t> span_len(1, 70);
          for (std::size_t r = 0; r < targets / 4; ++r) {
            const std::size_t i = pick(list_rng);
            TreeNode node;
            node.panel_count = std::min(span_len(list_rng), n);
            node.first_panel =
                std::uniform_int_distribution<std::size_t>(
                    0, n - node.panel_count)(list_rng);
            std::vector<std::size_t> exclusions;
            if (r % 2 == 0) {
              exclusions.push_back(order[node.first_panel]);
              exclusions.push_back(
                  order[node.first_panel + node.panel_count / 2]);
            }
            std::size_t ca = 0, cb = 0;
            const JcpAccum ta = direct_jcp_ordered_tile_sum(
                work[i], i, pc, uc, eps_ratio, kappa, coupled, node,
                "tile-scalar", ca, &exclusions, 0, exclusions.size());
            const JcpAccum tb = direct_jcp_ordered_tile_sum(
                work[i], i, pc, uc, eps_ratio, kappa, coupled, node, vec_tile,
                cb, &exclusions, 0, exclusions.size());
            std::vector<std::size_t> span(order.begin() + node.first_panel,
                                          order.begin() + node.first_panel +
                                              node.panel_count);
            std::vector<std::size_t> kept;
            for (std::size_t j : span)
              if (std::find(exclusions.begin(), exclusions.end(), j) ==
                  exclusions.end())
                kept.push_back(j);
            const JcpAccum scale = piece_scale_list(
                work[i], i, pc, uc, eps_ratio, kappa, coupled, kept);
            ++rows;
            const bool ok = ca == cb &&
                            close(ta.phi, tb.phi, scale.phi, worst_phi) &&
                            close(ta.phin, tb.phin, scale.phin, worst_phin);
            if (!ok && ++failures <= 10) {
              std::cout << "FAIL tile pass=" << pass << " coupled=" << coupled
                        << " kappa=" << kappa << " i=" << i << " first="
                        << node.first_panel << " len=" << node.panel_count
                        << " count " << ca << '/' << cb << " phi " << ta.phi
                        << '/' << tb.phi << " phin " << ta.phin << '/'
                        << tb.phin << '\n';
            }
          }
        }
      }
    }
  }

  // Throughput of the index-list kernel on long lists (one thread).
  const TreePanelCache pc = build_tree_panel_cache(panels, order);
  const TreeUnknownCache uc = build_tree_unknown_cache(panels, unknowns, true,
                                                       order);
  std::vector<std::size_t> list(256);
  std::mt19937_64 bench_rng(7);
  std::uniform_int_distribution<std::size_t> pick(0, n - 1);
  for (std::size_t &j : list)
    j = pick(bench_rng);
  for (bool coupled : {true, false}) {
    for (const std::string &kernel : {std::string("scalar"), vec}) {
      double sink = 0.0;
      std::size_t count = 0;
      const auto t0 = std::chrono::steady_clock::now();
      for (int rep = 0; rep < 4000; ++rep) {
        const std::size_t i = static_cast<std::size_t>(rep) % n;
        const JcpAccum acc = direct_jcp_ordered_sum(
            panels[i], i, pc, uc, 20.0, 0.125, coupled, list, 0, list.size(),
            kernel, count);
        sink += acc.phi + acc.phin;
      }
      const double ns = std::chrono::duration<double, std::nano>(
                            std::chrono::steady_clock::now() - t0)
                            .count() /
                        (4000.0 * static_cast<double>(list.size()));
      std::cout << "bench coupled=" << coupled << ' ' << kernel << ' ' << ns
                << " ns/pair (sink " << sink << ")\n";
    }
  }

  std::cout << "rows=" << rows << " failures=" << failures
            << " worst_phi=" << worst_phi << " worst_phin=" << worst_phin
            << '\n';
  std::cout << (failures == 0 ? "PASS" : "FAIL") << '\n';
  return failures == 0 ? 0 : 1;
}
