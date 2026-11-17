// KIN-14 oracle: histogram transfer entropy (TE) without the 0.5 pseudocount.
//
// TE(X->Y) = H(Y_t, Y_t-tau) + H(Y_t-tau, X_t-tau) - H(Y_t-tau)
//            - H(Y_t, Y_t-tau, X_t-tau), each entropy a plug-in estimate with
// the Miller-Madow correction (m - 1)/(2N) (m occupied cells, N samples).
//
// Bound for independent series. Under X independent of Y the plug-in TE is a
// conditional mutual information and 2 N TE ~ chi^2(df), df = (B-1)^2 B (the
// G-test of conditional independence, all cells occupied). Miller-Madow
// subtracts df/(2N) (its four corrections sum to -B (B-1)^2 / (2N)), so the
// corrected TE has mean ~0 and standard deviation sigma0 = sqrt(2 df)/(2N).
//
//   helper (new helper titan_info::transfer_entropy only): AR(1) series with
//     uniform marginals (Gaussian AR(1) through the normal CDF, so every bin
//     is equally populated and the chi^2 limit applies), N = 20000, B = 8,
//     tau = 5: sigma0 = 7.0e-4.
//       independent (phi_x 0.9, phi_y 0.5): |TE| <= 5 sigma0 both ways
//       coupled Y_t += 0.3 X_t-1:  TE(X->Y) - TE(Y->X) > 20 sigma0
//     (the pseudocount estimator gives TE(X->Y) = +0.0095 +- 0.0006 for the
//     independent pair, 13 sigma0: a 10-run simulation, same settings)
//   pipeline (DCCMNetwork::calculate_allostery_information_optimized, both
//     builds): four nodes moving along x, N = 2000 frames, default bins (8),
//     tau = 5. A: AR(1) phi 0.98, B: white noise, both independent of
//     everything; C: AR(1) phi 0.9 driving D (D_t = 0.3 D_t-1 + e + 0.8 C_t-1).
//     network_te_edges.csv (|net TE| > 0.05) must hold exactly the edge
//     C -> D. With the pseudocount the independent pair A,B has a net TE of
//     +0.11 +- 0.01 (simulation) and is written as a spurious edge A -> B.
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "aligned_allocator.h"
#include "arguments.h"
#include "Multi_MOL.h"
#include "titan_executor.h"
#include "titan_sparse_graph.h"
#if __has_include("titan_info_entropy.h")
#include "titan_info_entropy.h"
#define KIN14_HAVE_HELPER 1
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// the pipeline part fills the trajectory cache directly
#define private public
#include "dccm_network.h"
#undef private

namespace {

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

std::vector<double> ar1(std::size_t n, double phi, std::mt19937_64 &rng,
                        const std::vector<double> *drive = nullptr,
                        double coupling = 0.0) {
  std::normal_distribution<double> normal(0.0, 1.0);
  std::vector<double> x(n, 0.0);
  for (std::size_t t = 1; t < n; ++t)
    x[t] = phi * x[t - 1] + normal(rng) +
           (drive ? coupling * (*drive)[t - 1] : 0.0);
  return x;
}

#ifdef KIN14_HAVE_HELPER
std::vector<std::uint8_t> bins_of(const std::vector<double> &x, int B) {
  double lo = x[0], hi = x[0];
  for (double v : x) {
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  const double range = std::max(hi - lo, 1e-15);
  std::vector<std::uint8_t> b(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    b[i] = static_cast<std::uint8_t>(
        std::clamp(static_cast<int>((x[i] - lo) / range * B), 0, B - 1));
  return b;
}

std::vector<double> to_uniform(const std::vector<double> &x) {
  double m = 0.0, s = 0.0;
  for (double v : x)
    m += v;
  m /= x.size();
  for (double v : x)
    s += (v - m) * (v - m);
  s = std::sqrt(s / x.size());
  std::vector<double> u(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    u[i] = 0.5 * std::erfc(-(x[i] - m) / s / std::sqrt(2.0));
  return u;
}

void helper_part() {
  constexpr std::size_t n = 20000, B = 8, tau = 5;
  const double df = (B - 1.0) * (B - 1.0) * B;
  const double sigma0 = std::sqrt(2.0 * df) / (2.0 * (n - tau));
  std::mt19937_64 rng(4242u);
  const auto x = ar1(n, 0.9, rng);
  const auto y = ar1(n, 0.5, rng);
  const auto bx = bins_of(to_uniform(x), B), by = bins_of(to_uniform(y), B);
  const double te_xy = titan_info::transfer_entropy(bx.data(), by.data(), n, B, tau);
  const double te_yx = titan_info::transfer_entropy(by.data(), bx.data(), n, B, tau);
  std::cout << std::setprecision(6) << std::fixed
            << "case helper independent (N=20000, B=8, tau=5, sigma0="
            << sigma0 << "): TE(x->y) = " << te_xy << ", TE(y->x) = " << te_yx
            << " nats\n";
  check(std::abs(te_xy) <= 5 * sigma0 && std::abs(te_yx) <= 5 * sigma0,
        "independent: |TE| <= 5 sigma0 both ways");

  const auto y2 = ar1(n, 0.5, rng, &x, 0.3);
  const auto by2 = bins_of(to_uniform(y2), B);
  const double c_xy = titan_info::transfer_entropy(bx.data(), by2.data(), n, B, tau);
  const double c_yx = titan_info::transfer_entropy(by2.data(), bx.data(), n, B, tau);
  std::cout << "case helper coupled: TE(x->y) = " << c_xy << ", TE(y->x) = "
            << c_yx << " nats\n";
  check(c_xy - c_yx > 20 * sigma0, "coupled: TE(x->y) - TE(y->x) > 20 sigma0");
}
#endif

void pipeline_part() {
  constexpr std::size_t n = 2000;
  std::mt19937_64 rng(99u);
  const auto a = ar1(n, 0.98, rng);
  const auto b = ar1(n, 0.0, rng);
  const auto c = ar1(n, 0.9, rng);
  const auto d = ar1(n, 0.3, rng, &c, 0.8);
  const std::vector<double> *series[4] = {&a, &b, &c, &d};

  DCCMNetwork net;
  net.ca_indices = {0, 1, 2, 3};
  net.ca_resnrs = {1, 2, 3, 4};
  net.ca_is_lig = {0, 0, 0, 0};
  net.ca_chain = {"A", "A", "A", "A"};
  net.ca_resname = {"ALA", "GLY", "SER", "THR"};
  net.store_full_traj_ = true;
  for (std::size_t f = 0; f < n; ++f) {
    DCCMNetwork::TrajectoryRow xr(4), yr(4), zr(4);
    for (std::size_t k = 0; k < 4; ++k) {
      xr[k] = 10.0 * k + 0.3 * (*series[k])[f];
      yr[k] = 1.0 + 0.5 * k;
      zr[k] = -2.0;
    }
    net.traj_x.push_back(xr);
    net.traj_y.push_back(yr);
    net.traj_z.push_back(zr);
  }
  net.total_frames = n;
  DCCMNetwork::AllosteryConfig cfg; // tau 5, default bins, histogram
  TitanExecutor executor(2);
  net.calculate_allostery_information_optimized(executor, cfg, "gc.csv",
                                                "mi.csv", "te.csv");
  std::ifstream in("te.csv");
  std::string line;
  std::getline(in, line);
  std::vector<std::tuple<std::string, std::string, double>> edges;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    std::string s, t, v;
    std::getline(row, s, ',');
    std::getline(row, t, ',');
    std::getline(row, v, ',');
    edges.emplace_back(s, t, std::stod(v));
  }
  std::cout << "case pipeline (N=2000, default bins, tau=5): " << edges.size()
            << " TE edge(s)";
  for (const auto &[s, t, v] : edges)
    std::cout << "  " << s << "->" << t << " " << v;
  std::cout << '\n';
  bool spurious = false, driven = false;
  for (const auto &[s, t, v] : edges) {
    if (s == "SER3_A" && t == "THR4_A")
      driven = true;
    else
      spurious = true;
  }
  check(driven, "pipeline: edge C -> D (SER3 -> THR4) found");
  check(!spurious, "pipeline: no edge between independent nodes");
}

} // namespace

int main() {
#ifdef KIN14_HAVE_HELPER
  helper_part();
#else
  check(false, "no Miller-Madow transfer-entropy helper (titan_info_entropy.h)");
#endif
  pipeline_part();
  if (failures) {
    std::cout << "FAIL KIN-14 transfer entropy: " << failures
              << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-14 transfer entropy: TE ~ 0 for independent series "
               "(chi^2 bound), X->Y for a coupled pair, no spurious edge\n";
  return 0;
}
