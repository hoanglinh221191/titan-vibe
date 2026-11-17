// The FMM near-field branches must apply one operator, and that operator must
// converge to the direct one.  A small PN sphere (an icosphere of 162
// vertices, its sphere-patch parents dropped so the cells are PN lifts) on
// one mesh and one density, against one direct reference (the audit's rows,
// also built pair by pair as the four kernel blocks K1..K4):
//
//  1. Parity.  Under a leaf-1 tree every leaf holds one cell, so the tree
//     puts one-ring and curved near-zone pairs in the far field, where only a
//     correction (curved minus one-point) gives them their curved value.  The
//     p1 packed and scalar near loops, and the p4 packed and scalar
//     fallbacks, must agree to rounding (the scalar p1 loop was 4% off before
//     it took the correction).
//  2. Convergence.  On the same tree, tightening theta must shrink the error
//     of p1 and p4 in every block, p4 faster: a truncation error, not a
//     wrong moment, derivative sign or correction (which would not shrink).
//  3. All pairs direct (the root a leaf): every branch must reproduce the
//     direct reference to rounding.
//
// Each block is reported against its own size, ||dKx|| / ||Kx||, and against
// the cancellation-free size ||dKx|| / || |K||x| ||, so a small or cancelling
// block is neither hidden by a large one nor inflated by its own
// cancellation.
//
// The scalar p1 loop runs where AVX2 is withheld (TITAN_FORCE_SCALAR=1, read
// once per process), so the program runs itself again for it.
#include "tabi/tabi_fmm_m2l_fallback.h"
#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_near.h"
#include "tabi/tabi_node_patch.h"
#include "tabi/tabi_node_patch_curved.h"
#include "titan_executor.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
using namespace TitanPBE::Tabi;
namespace fmm = titan::tabi::fmm_m2l;
namespace {
constexpr double kEps = 80.0, kKappa = 0.125, kZone = 8.0;
struct Config {
  int leaf;
  double theta;
};
// The leaf-1 sweep (theta falling), then the root a leaf.
const std::vector<Config> kConfigs{
    {1, 0.6}, {1, 0.45}, {1, 0.3}, {1 << 20, 0.45}};
double seconds_since(std::chrono::steady_clock::time_point begin) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin)
      .count();
}
using Vec = std::vector<double>;
// The three densities: (phi, dphi), (phi, 0) and (0, dphi).  An output is
// the interaction part of the operator (the audit's sums), phi rows then
// dphi rows: from (phi, 0) the phi rows are K2 phi and the dphi rows K4 phi;
// from (0, dphi) they are K1 dphi and K3 dphi.
struct Case {
  std::vector<Panel> panels;
  std::array<Vec, 3> phi, phin;
};
using Outputs = std::array<Vec, 3>;
Case build(TitanExecutor &ex) {
  Options o;
  o.primary_surface = "sphere-patch";
  o.surface_mode = "vdw";
  o.sphere_mesh_mode = "icosphere";
  o.surface_point_count = 300;
  const std::vector<Atom> atoms{{0, 0, 0, 1, 3.0}};
  auto triangles = build_surface(atoms, o, nullptr, &ex);
  for (auto &t : triangles)
    t.sphere_patch_quadrature.reset();
  Case c;
  c.panels = build_node_patch_panels(triangles);
  // The FMM route: no eager pair cache, curved cells within kZone radii.
  if (!prepare_curved_node_patches(triangles, c.panels, kEps, kKappa,
                                   64u * 1024 * 1024, ex, false, 4, kZone))
    throw std::runtime_error("no curved geometry");
  if (c.panels.front().node_patch_geometry->curved->representation !=
      "curved_pn")
    throw std::runtime_error("not a PN surface");
  const std::size_t n = c.panels.size();
  for (auto &v : c.phi)
    v.assign(n, 0.0);
  for (auto &v : c.phin)
    v.assign(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    c.phi[0][j] = c.phi[1][j] = std::sin(0.137 * j) + 0.2;
    c.phin[0][j] = c.phin[2][j] = std::cos(0.071 * j) - 0.3;
  }
  return c;
}
fmm::FmmM2lAnalyticContext context(const Case &c, const Config &k,
                                   TitanExecutor &ex) {
  fmm::PrototypeConfig config;
  config.kernel = fmm::KernelKind::Yukawa;
  config.kappa = kKappa;
  config.max_leaf_size = k.leaf;
  config.expansion_order = 4;
  config.separation_theta = k.theta;
  config.threads = static_cast<int>(ex.capacity());
  std::vector<fmm::PrototypePanel> prototype;
  for (const Panel &p : c.panels) {
    fmm::PrototypePanel q;
    q.x = p.x; q.y = p.y; q.z = p.z;
    q.nx = p.nx; q.ny = p.ny; q.nz = p.nz;
    q.area = p.area;
    q.source_x = p.source_x; q.source_y = p.source_y; q.source_z = p.source_z;
    q.source_nx = p.source_nx; q.source_ny = p.source_ny;
    q.source_nz = p.source_nz;
    q.has_source_quadrature_center = p.has_source_quadrature_center;
    prototype.push_back(q);
  }
  return fmm::FmmM2lAnalyticContext(config, std::move(prototype), kEps, kKappa);
}
Vec joined(Vec a, const Vec &b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}
struct Run {
  Outputs out;
  std::size_t far_pairs{0}, refined_cells{0}, fallback_pairs{0};
  std::string branch;
};
Run p1(const Case &c, const Config &k, TitanExecutor &ex) {
  const auto ctx = context(c, k, ex);
  Run run;
  for (int d = 0; d < 3; ++d) {
    Vec phi, phin;
    const auto r = ctx.apply_into(c.phi[d], c.phin[d], ex, phi, phin, true,
                                  "analytic_singular", "off", &c.panels);
    if (!r.finite_outputs)
      throw std::runtime_error("p1 apply not finite");
    run.out[d] = joined(phi, phin);
    run.far_pairs = r.far_pairs;
    run.branch = r.near_kernel_policy;
  }
  return run;
}
Run p4(const Case &c, const Config &k, TitanExecutor &ex,
       const std::string &fallback) {
  auto ctx = context(c, k, ex);
  fmm::JcpLocalP4OperatorOptions o;
  o.near_integration = "analytic_singular";
  o.refine_mode = "operator-tail-static";
  o.refine_tolerance = 0.002;
  o.cache_mode = "green-jet-compressed";
  o.cache_budget_mb = 1024;
  o.fallback_mode = fallback;
  o.plan_builder = "metric-cache";
  o.refine_probe_count = 4;
  (void)ctx.prepare_local_p4_operator(c.panels, ex, o);
  Run run;
  for (int d = 0; d < 3; ++d) {
    const auto r =
        ctx.apply_prepared_local_p4_operator(c.phi[d], c.phin[d], ex, true);
    if (!r.finite_outputs || r.status != "ok")
      throw std::runtime_error("p4 apply (" + fallback + "): " + r.status);
    run.out[d] = joined(r.phi, r.phin);
    run.far_pairs = r.far_pairs;
    run.refined_cells = r.refined_cells;
    run.fallback_pairs = r.direct_fallback_pairs;
    run.branch = r.fallback_effective;
  }
  return run;
}
// The direct reference, pair by pair: K[b][t * n + j] for b = K1..K4.
struct Direct {
  std::array<Vec, 4> K;
  Outputs out;
  double row_check{0.0};
};
Direct direct(const Case &c) {
  const std::size_t n = c.panels.size();
  Direct d;
  for (auto &k : d.K)
    k.assign(n * n, 0.0);
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  Vec unit(2 * n, 0.0);
  for (std::size_t t = 0; t < n; ++t)
    for (std::size_t j = 0; j < n; ++j) {
      unit[j] = 1.0;
      const JcpNearAccum a = accumulate_jcp_near_ordered(
          c.panels[t], t, c.panels, unit, kEps, kKappa, true, order, j, 1,
          "analytic_singular", "off");
      unit[j] = 0.0;
      unit[n + j] = 1.0;
      const JcpNearAccum b = accumulate_jcp_near_ordered(
          c.panels[t], t, c.panels, unit, kEps, kKappa, true, order, j, 1,
          "analytic_singular", "off");
      unit[n + j] = 0.0;
      d.K[1][t * n + j] = a.phi;
      d.K[3][t * n + j] = a.phin;
      d.K[0][t * n + j] = b.phi;
      d.K[2][t * n + j] = b.phin;
    }
  for (int s = 0; s < 3; ++s) {
    d.out[s].assign(2 * n, 0.0);
    for (std::size_t t = 0; t < n; ++t)
      for (std::size_t j = 0; j < n; ++j) {
        d.out[s][t] += d.K[1][t * n + j] * c.phi[s][j] +
                       d.K[0][t * n + j] * c.phin[s][j];
        d.out[s][n + t] += d.K[3][t * n + j] * c.phi[s][j] +
                           d.K[2][t * n + j] * c.phin[s][j];
      }
  }
  // The pairwise blocks against the audit's whole rows.
  Vec u = joined(c.phi[0], c.phin[0]);
  double diff = 0.0, size = 0.0;
  for (std::size_t t = 0; t < n; ++t) {
    const JcpNearAccum a = accumulate_jcp_near_ordered(
        c.panels[t], t, c.panels, u, kEps, kKappa, true, order, 0, n,
        "analytic_singular", "off");
    diff = std::max({diff, std::abs(a.phi - d.out[0][t]),
                     std::abs(a.phin - d.out[0][n + t])});
    size = std::max({size, std::abs(a.phi), std::abs(a.phin)});
  }
  d.row_check = diff / size;
  return d;
}
// Errors of a run: [0] the whole matvec, [1..4] the blocks K1..K4; each
// relative to the block (rel) and to its cancellation-free size (abs).
struct Error {
  std::array<double, 5> rel{}, abs{};
};
Error error(const Case &c, const Direct &d, const Outputs &out) {
  const std::size_t n = c.panels.size();
  Error e;
  // Block b: density s (1 = phi, 2 = dphi), rows offset (0 = phi, n = dphi).
  const std::array<int, 4> density{2, 1, 2, 1};
  const std::array<std::size_t, 4> rows{0, 0, n, n};
  std::array<Vec, 5> magnitude;
  for (auto &m : magnitude)
    m.assign(2 * n, 0.0);
  for (int b = 0; b < 4; ++b) {
    const Vec &x = density[b] == 1 ? c.phi[1] : c.phin[2];
    double dd = 0.0, ss = 0.0, mm = 0.0;
    for (std::size_t t = 0; t < n; ++t) {
      double m = 0.0;
      for (std::size_t j = 0; j < n; ++j)
        m += std::abs(d.K[b][t * n + j] * x[j]);
      const double ref = d.out[density[b]][rows[b] + t];
      const double delta = out[density[b]][rows[b] + t] - ref;
      dd += delta * delta;
      ss += ref * ref;
      mm += m * m;
      magnitude[0][rows[b] + t] += m;
    }
    e.rel[b + 1] = std::sqrt(dd / std::max(ss, 1e-300));
    e.abs[b + 1] = std::sqrt(dd / std::max(mm, 1e-300));
  }
  double dd = 0.0, ss = 0.0, mm = 0.0;
  for (std::size_t i = 0; i < 2 * n; ++i) {
    const double delta = out[0][i] - d.out[0][i];
    dd += delta * delta;
    ss += d.out[0][i] * d.out[0][i];
    mm += magnitude[0][i] * magnitude[0][i];
  }
  e.rel[0] = std::sqrt(dd / std::max(ss, 1e-300));
  e.abs[0] = std::sqrt(dd / std::max(mm, 1e-300));
  return e;
}
double parity(const Outputs &a, const Outputs &b) {
  double worst = 0.0;
  for (int s = 0; s < 3; ++s) {
    double dd = 0.0, ss = 0.0;
    for (std::size_t i = 0; i < a[s].size(); ++i) {
      dd += (a[s][i] - b[s][i]) * (a[s][i] - b[s][i]);
      ss += b[s][i] * b[s][i];
    }
    worst = std::max(worst, std::sqrt(dd / std::max(ss, 1e-300)));
  }
  return worst;
}
void save(const std::string &path, const std::vector<Outputs> &runs) {
  std::ofstream f(path, std::ios::binary);
  for (const auto &run : runs)
    for (const auto &v : run) {
      const std::size_t n = v.size();
      f.write(reinterpret_cast<const char *>(&n), sizeof n);
      f.write(reinterpret_cast<const char *>(v.data()),
              static_cast<std::streamsize>(n * sizeof(double)));
    }
  if (!f)
    throw std::runtime_error("cannot write " + path);
}
std::vector<Outputs> load(const std::string &path, std::size_t count) {
  std::ifstream f(path, std::ios::binary);
  std::vector<Outputs> runs(count);
  for (auto &run : runs)
    for (auto &v : run) {
      std::size_t n = 0;
      f.read(reinterpret_cast<char *>(&n), sizeof n);
      v.resize(n);
      f.read(reinterpret_cast<char *>(v.data()),
             static_cast<std::streamsize>(n * sizeof(double)));
    }
  if (!f)
    throw std::runtime_error("cannot read " + path);
  return runs;
}
// The far-classified pairs the corrections cover, from the near-field log.
std::size_t far_classified(const std::string &log) {
  std::ifstream f(log);
  std::string line;
  std::size_t count = 0;
  const std::string head = "TABI node-patch: ",
                    tail = " one-ring or curved near-zone pair(s) the tree puts";
  while (std::getline(f, line)) {
    const auto at = line.find(tail), from = line.find(head);
    if (at != std::string::npos && from != std::string::npos)
      count = std::max<std::size_t>(
          count, std::stoull(line.substr(from + head.size(),
                                         at - from - head.size())));
  }
  return count;
}
void print(const char *name, const Error &e) {
  std::cout << "    " << std::left << std::setw(10) << name << std::right;
  const char *label[5] = {"matvec", "K1", "K2", "K3", "K4"};
  for (int b = 0; b < 5; ++b)
    std::cout << " " << label[b] << " " << std::setw(9) << e.rel[b] << "/"
              << std::setw(9) << e.abs[b];
  std::cout << "\n";
}
} // namespace
int main(int argc, char **argv) {
  try {
    TitanExecutor ex(4);
    if (argc == 3 && std::string(argv[1]) == "--scalar-p1") {
      // The child: AVX2 withheld, so the p1 near loop is the scalar one.
      const Case c = build(ex);
      std::vector<Outputs> runs;
      std::string branch;
      for (const Config &k : kConfigs) {
        Run run = p1(c, k, ex);
        branch = run.branch;
        runs.push_back(std::move(run.out));
      }
      save(argv[2], runs);
      std::cout << "CHILD p1 near policy " << branch << "\n";
      return 0;
    }
    const std::string self = argv[0];
    const std::string vec = self + ".scalar_p1.bin", log = self + ".scalar_p1.log";
    std::remove(vec.c_str());
#ifdef _WIN32
    _putenv_s("TITAN_FORCE_SCALAR", "1");
#else
    setenv("TITAN_FORCE_SCALAR", "1", 1);
#endif
    std::string command = "\"" + self + "\" --scalar-p1 \"" + vec + "\" > \"" +
                          log + "\" 2>&1";
#ifdef _WIN32
    command = "\"" + command + "\""; // cmd.exe strips one pair of quotes
#endif
    const int child = std::system(command.c_str());
    // The packed branches below need AVX2, and runtime_has_avx2() reads the
    // variable once, on first use (nothing here has used it yet), so it is
    // removed again before this process applies anything.
#ifdef _WIN32
    _putenv_s("TITAN_FORCE_SCALAR", "");
#else
    unsetenv("TITAN_FORCE_SCALAR");
#endif
    if (child != 0)
      throw std::runtime_error("scalar p1 child failed (see " + log + ")");
    const std::vector<Outputs> scalar_p1 = load(vec, kConfigs.size());
    const std::size_t far_scalar = far_classified(log);
    // Without AVX2 the packed branches do not exist here; the scalar ones
    // are still checked against the direct operator.
    const bool avx2 = fmm::packed_jcp_fallback_avx2_available();
    const std::string packed_mode = avx2 ? "packed-avx2" : "scalar";
    if (!avx2)
      std::cout << "NOTE no AVX2: the packed branches are not exercised\n";

    auto clock = std::chrono::steady_clock::now();
    const Case c = build(ex);
    const std::size_t n = c.panels.size();
    const Direct d = direct(c);
    const double direct_seconds = seconds_since(clock);
    std::cout << std::setprecision(3)
              << "BRANCHES cells=" << n << " zone=" << kZone
              << " eps=" << kEps << " kappa=" << kKappa
              << "; pairwise K1..K4 reproduce the audit rows to "
              << d.row_check << " (surface and reference " << direct_seconds
              << " s)\n"
              << "  error rel/abs: rel = ||dKx||/||Kx||, abs = ||dKx||/|| |K||x| ||\n";
    bool ok = d.row_check < 1e-13 && far_scalar > 0;
    std::vector<Error> sweep_p1, sweep_p4;
    for (std::size_t i = 0; i < kConfigs.size(); ++i) {
      const Config &k = kConfigs[i];
      clock = std::chrono::steady_clock::now();
      const Run packed1 = p1(c, k, ex);
      const Run packed4 = p4(c, k, ex, packed_mode);
      const Run scalar4 = p4(c, k, ex, "scalar");
      const double parity1 = parity(scalar_p1[i], packed1.out);
      const double parity4 = parity(scalar4.out, packed4.out);
      const Error e1 = error(c, d, packed1.out), e4 = error(c, d, packed4.out);
      const Error s1 = error(c, d, scalar_p1[i]), s4 = error(c, d, scalar4.out);
      std::cout << "  leaf=" << (k.leaf > int(n) ? std::string("all")
                                                  : std::to_string(k.leaf))
                << " theta=" << k.theta << ": far pairs p1 " << packed1.far_pairs
                << ", p4 " << packed4.far_pairs << " (refined cells "
                << packed4.refined_cells << ", direct fallback pairs "
                << packed4.fallback_pairs << "); scalar vs packed p1 "
                << parity1 << ", p4 " << parity4 << " (" << seconds_since(clock)
                << " s)\n";
      print("p1 packed", e1);
      print("p1 scalar", s1);
      print("p4 packed", e4);
      print("p4 scalar", s4);
      ok = ok && parity1 < 1e-10 && parity4 < 1e-10 &&
           packed4.branch == packed_mode && scalar4.branch == "scalar";
      if (k.leaf > int(n)) {
        // Every pair direct: the branches are the direct operator.
        ok = ok && packed1.far_pairs == 0;
        for (const Error *e : {&e1, &s1, &e4, &s4})
          for (double r : e->rel)
            ok = ok && r < 1e-12;
      } else if (k.leaf == 1) {
        sweep_p1.push_back(e1);
        sweep_p4.push_back(e4);
      }
    }
    // Theta falls along the sweep: each error must fall with it, p4 below
    // p1.  (The single-cell orders, h^2 for p1 in every block, are
    // tabi_fmm_p1_degree_probe.cpp.)
    for (std::size_t i = 1; i < sweep_p1.size(); ++i)
      for (int b = 0; b < 5; ++b)
        ok = ok && sweep_p1[i].abs[b] < sweep_p1[i - 1].abs[b] &&
             sweep_p4[i].abs[b] < sweep_p4[i - 1].abs[b] &&
             sweep_p4[i].abs[b] < sweep_p1[i].abs[b];
    std::cout << "  far-classified pairs corrected (leaf 1): " << far_scalar
              << "\n"
              << (ok ? "PASS" : "FAIL")
              << " FMM near-branch parity, theta convergence, all-direct\n";
    return ok ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << "\n";
    return 2;
  }
}
