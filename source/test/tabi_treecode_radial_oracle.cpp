// TABI-03 oracle: the radial-derivative Taylor coefficients that the treecode
// far field now uses, against an independent long-double reference, with the
// truncated power-series path they replaced measured on the same cases.
//
// 1. Coefficients. The reference expands the four JCP kernel channels about a
//    node centre with its own long-double multivariate series algebra
//    (composition through f^(k)(s0)/k! t^k, not the recurrences of the
//    production series code). For random targets, normals, node centres and
//    moments, orders 2-5, kappa 0..5, eps 1..80, coupled or not, each channel
//    of taylor_jcp_node_channel_sum (radial) and of
//    taylor_jcp_node_channel_sum_series is compared with it, normalised by the
//    sum of absolute moment terms.
// 2. Expansion. Point-source clusters of radius rho seen from rho/theta: the
//    Taylor far field from exact long-double moments against the direct sum of
//    the point kernels. Radial and series must leave the same truncation error,
//    and that error must fall with the order.
// 3. taylor_jcp_node_sum equals the sum of its channels.
// 4. The fast far-field kernel of the treecode apply (radial_jcp_far_kernel),
//    row by row (phi = K1+K2, phin = K3+K4), against the same reference and
//    against taylor_jcp_node_sum; the timing block compares the two as the
//    apply's far loop calls them.
// 5. The four-node kernel (radial_jcp_far_kernel4): every lane bit-identical
//    to the one-node kernel, degenerate lanes included; and the scalar
//    fallback kernel (radial_fast::node_sum, used without vector types) equal
//    to the one-node kernel to rounding.
//
// Build: link with the engine objects of the same profile (all but Main.o).
#include "tabi/tabi_radial_taylor.h"
#include "tabi/tabi_treecode_internal.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::detail;

namespace {

using LD = long double;
constexpr int kMaxDegree = 5;
constexpr LD kPiL = 3.141592653589793238462643383279502884L;

// Truncated series in the source displacement h = (hx, hy, hz).
struct Series {
  int degree{kMaxDegree};
  LD c[kMaxDegree + 1][kMaxDegree + 1][kMaxDegree + 1]{};
};

Series constant(LD value, int degree) {
  Series s;
  s.degree = degree;
  s.c[0][0][0] = value;
  return s;
}

Series operator+(const Series &a, const Series &b) {
  Series s = a;
  for (int i = 0; i <= a.degree; ++i)
    for (int j = 0; i + j <= a.degree; ++j)
      for (int k = 0; i + j + k <= a.degree; ++k)
        s.c[i][j][k] += b.c[i][j][k];
  return s;
}

Series operator*(LD x, const Series &a) {
  Series s = a;
  for (auto &plane : s.c)
    for (auto &row : plane)
      for (LD &v : row)
        v *= x;
  return s;
}

Series operator-(const Series &a, const Series &b) { return a + (-1.0L) * b; }

Series operator*(const Series &a, const Series &b) {
  Series s;
  s.degree = a.degree;
  const int d = a.degree;
  for (int i = 0; i <= d; ++i)
    for (int j = 0; i + j <= d; ++j)
      for (int k = 0; i + j + k <= d; ++k) {
        LD sum = 0;
        for (int i1 = 0; i1 <= i; ++i1)
          for (int j1 = 0; j1 <= j; ++j1)
            for (int k1 = 0; k1 <= k; ++k1)
              sum += a.c[i1][j1][k1] * b.c[i - i1][j - j1][k - k1];
        s.c[i][j][k] = sum;
      }
  return s;
}

// f(s0 + t) = sum_k f_k t^k with f_k = f^(k)(s0) / k!, t = s - s0.
template <typename Coefficient> Series compose(const Series &s, Coefficient f) {
  const LD s0 = s.c[0][0][0];
  Series t = s;
  t.c[0][0][0] = 0;
  Series power = constant(1, s.degree);
  Series out = constant(f(0, s0), s.degree);
  for (int k = 1; k <= s.degree; ++k) {
    power = power * t;
    out = out + f(k, s0) * power;
  }
  return out;
}

LD inv_factorial(int k) {
  LD f = 1;
  for (int i = 2; i <= k; ++i)
    f *= i;
  return 1 / f;
}

Series sqrt_series(const Series &s) {
  return compose(s, [](int k, LD s0) {
    LD binomial = 1; // binom(1/2, k)
    for (int i = 0; i < k; ++i)
      binomial *= (0.5L - i) / (i + 1);
    return binomial * std::pow(s0, 0.5L - k);
  });
}

Series reciprocal_series(const Series &s) {
  return compose(s, [](int k, LD s0) {
    return (k % 2 ? -1.0L : 1.0L) * std::pow(s0, -1.0L - k);
  });
}

// exp(s) - 1 without cancellation in the constant term.
Series expm1_series(const Series &s) {
  return compose(s, [](int k, LD s0) {
    return k == 0 ? std::expm1(s0) : std::exp(s0) * inv_factorial(k);
  });
}

struct Vec {
  LD x, y, z;
};

// Kernel channels of the JCP rows as functions of the source s = c + h, for a
// target at t with unit normal n: the Laplace kernels and the screened-minus-
// Laplace differences. The differences are formed from exp(-kr) - 1 and
// (1 + kr) exp(-kr) - 1 directly, never by subtracting two kernels, so they
// keep full relative accuracy when kr is small.
struct KernelSeries {
  Series green;                   // G0
  Series target_normal;           // dG0/dn_t
  std::array<Series, 3> source;   // source-normal component c of G0
  std::array<Series, 3> mixed;    // mixed second derivative of G0
  Series d_green;                 // Gk - G0, and so on
  Series d_target_normal;
  std::array<Series, 3> d_source;
  std::array<Series, 3> d_mixed;
};

KernelSeries kernel_series(const Vec &t, const Vec &n, const Vec &centre,
                           LD kappa, int degree) {
  std::array<Series, 3> u;
  const LD d[3] = {t.x - centre.x, t.y - centre.y, t.z - centre.z};
  for (int i = 0; i < 3; ++i) {
    u[i] = constant(d[i], degree);
    if (degree >= 1) {
      if (i == 0) u[i].c[1][0][0] = -1;
      if (i == 1) u[i].c[0][1][0] = -1;
      if (i == 2) u[i].c[0][0][1] = -1;
    }
  }
  const Series r2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  const Series r = sqrt_series(r2);
  const Series ir = reciprocal_series(r);
  const Series ir2 = ir * ir;
  const Series ir3 = ir2 * ir;
  const Series em1 = expm1_series((-kappa) * r);   // exp(-kr) - 1
  const Series e = constant(1, degree) + em1;
  // (1 + kr) exp(-kr) - 1 = (1 + kr) (exp(-kr) - 1) + kr
  const Series one_kr = constant(1, degree) + kappa * r;
  const Series g1 = one_kr * em1 + kappa * r;
  const LD q = 1 / (4 * kPiL);
  const Series n_dot_u = n.x * u[0] + n.y * u[1] + n.z * u[2];
  const Series nxu = n_dot_u * ir;
  const LD nc[3] = {n.x, n.y, n.z};
  KernelSeries k;
  k.green = q * ir;
  k.d_green = q * (em1 * ir);
  k.target_normal = (-q) * (n_dot_u * ir3);
  k.d_target_normal = (-q) * (g1 * n_dot_u * ir3);
  // fprime = -q (1 + kr) e / r^2; fsecond = q e (k^2/r + 2k/r^2 + 2/r^3)
  const Series fprime0 = (-q) * ir2;
  const Series d_fprime = (-q) * (g1 * ir2);
  const Series fsecond0 = (2 * q) * ir3;
  const Series d_fsecond =
      q * (em1 * (2.0L * ir3) +
           e * ((kappa * kappa) * ir + (2 * kappa) * ir2));
  for (int c = 0; c < 3; ++c) {
    k.source[c] = q * (u[c] * ir3);
    k.d_source[c] = q * (g1 * u[c] * ir3);
    const Series nyu = u[c] * ir;
    const Series metric = constant(nc[c], degree) - nxu * nyu;
    k.mixed[c] = (-1.0L) * (fsecond0 * nxu * nyu + fprime0 * ir * metric);
    k.d_mixed[c] = (-1.0L) * (d_fsecond * nxu * nyu + d_fprime * ir * metric);
  }
  return k;
}

struct ChannelRef {
  LD value[4]{};
  LD scale[4]{};
};

LD coefficient(const Series &s, const MomentIndex &a) {
  return s.c[a.ax][a.ay][a.az];
}

ChannelRef reference_channels(const Vec &t, const Vec &n, const Vec &centre,
                              const NodeMoments &scalar,
                              const std::array<NodeMoments, 3> &normal,
                              LD eps, LD kappa, bool coupled, int order) {
  const KernelSeries k = kernel_series(t, n, centre, kappa, order);
  const auto &indices = moment_indices(order);
  ChannelRef out;
  const auto add = [&](std::size_t channel, LD term) {
    out.value[channel] += term;
    out.scale[channel] += std::fabs(term);
  };
  const LD inv_eps = 1 / std::max<LD>(1e-12L, eps);
  for (std::size_t m = 0; m < indices.size(); ++m) {
    const MomentIndex &a = indices[m];
    if (coupled) {
      // G0 - Gk
      add(kJcpK1PhiChannel, -scalar.values[m] * coefficient(k.d_green, a));
      // T0 - Tk/eps = (1 - 1/eps) T0 - (Tk - T0)/eps
      add(kJcpK3PhinChannel,
          scalar.values[m] * ((1 - inv_eps) * coefficient(k.target_normal, a) -
                              inv_eps * coefficient(k.d_target_normal, a)));
    }
    for (int c = 0; c < 3; ++c) {
      const LD d0 = coefficient(k.source[c], a);
      if (!coupled) {
        add(kJcpK2PhiNormalChannel, normal[c].values[m] * (eps - 1) * d0);
        continue;
      }
      // eps Dk - D0 = (eps - 1) D0 + eps (Dk - D0)
      add(kJcpK2PhiNormalChannel,
          normal[c].values[m] *
              ((eps - 1) * d0 + eps * coefficient(k.d_source[c], a)));
      add(kJcpK4MixedChannel,
          normal[c].values[m] * coefficient(k.d_mixed[c], a));
    }
  }
  return out;
}

double channel_value(const JcpChannelAccum &a, std::size_t channel) {
  return channel == kJcpK1PhiChannel || channel == kJcpK2PhiNormalChannel
             ? a.channel[channel].phi
             : a.channel[channel].phin;
}

Vec random_unit(std::mt19937_64 &rng) {
  std::normal_distribution<double> g;
  double x = g(rng), y = g(rng), z = g(rng);
  const double s = std::sqrt(x * x + y * y + z * z);
  return {x / s, y / s, z / s};
}

// Point kernels for the direct sum, from the same closed forms in long double.
struct PointChannels {
  LD v[4]{};
};

PointChannels point_channels(const Vec &t, const Vec &n, const Vec &s, LD q,
                             const LD w[3], LD eps, LD kappa, bool coupled) {
  PointChannels out;
  const LD u[3] = {t.x - s.x, t.y - s.y, t.z - s.z};
  const LD r2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2], r = std::sqrt(r2);
  const LD qq = 1 / (4 * kPiL);
  const LD n_dot_u = n.x * u[0] + n.y * u[1] + n.z * u[2];
  const LD nc[3] = {n.x, n.y, n.z};
  const auto green = [&](LD k) { return qq * std::exp(-k * r) / r; };
  const auto tn = [&](LD k) {
    return -qq * std::exp(-k * r) * (1 + k * r) * n_dot_u / (r2 * r);
  };
  const auto src = [&](LD k, int c) {
    return qq * std::exp(-k * r) * (1 + k * r) * u[c] / (r2 * r);
  };
  const auto mixed = [&](LD k, int c) {
    const LD e = std::exp(-k * r);
    const LD fprime = -qq * e * (1 + k * r) / r2;
    const LD fsecond = qq * e * (k * k / r + 2 * k / r2 + 2 / (r2 * r));
    const LD nxu = n_dot_u / r, nyu = u[c] / r;
    return -(fsecond * nxu * nyu + fprime / r * (nc[c] - nxu * nyu));
  };
  const LD inv_eps = 1 / std::max<LD>(1e-12L, eps);
  if (coupled) {
    out.v[kJcpK1PhiChannel] = q * (green(0) - green(kappa));
    out.v[kJcpK3PhinChannel] = q * (tn(0) - inv_eps * tn(kappa));
  }
  for (int c = 0; c < 3; ++c) {
    if (!coupled) {
      out.v[kJcpK2PhiNormalChannel] += w[c] * (eps - 1) * src(0, c);
      continue;
    }
    out.v[kJcpK2PhiNormalChannel] += w[c] * (eps * src(kappa, c) - src(0, c));
    out.v[kJcpK4MixedChannel] += w[c] * (mixed(kappa, c) - mixed(0, c));
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const bool skip_timing =
      argc > 1 && std::strcmp(argv[1], "--no-timing") == 0;
  std::mt19937_64 rng(20261005);
  std::uniform_real_distribution<double> unit(-1.0, 1.0);
  int failures = 0;

  // ---- 1. coefficients against the long-double reference -----------------
  constexpr std::size_t kKappaCount = 8;
  const double kappas[kKappaCount] = {0.0,    1e-9, 1e-4, 0.05,
                                      0.1235, 0.5,  2.0,  5.0};
  const double epss[] = {1.0, 4.0, 39.25, 80.0};
  double radial_max = 0, series_max = 0, fast_max = 0, fast_vs_radial_max = 0;
  double radial_max_screened_small = 0, series_max_screened_small = 0;
  std::array<double, kKappaCount> radial_by_kappa{}, series_by_kappa{};
  std::size_t cases = 0;
  for (int order = 2; order <= 5; ++order)
    for (std::size_t kappa_index = 0; kappa_index < kKappaCount; ++kappa_index)
      for (const double eps : epss)
        for (const bool coupled : {false, true})
          for (int sample = 0; sample < 10; ++sample) {
            const double kappa = kappas[kappa_index];
            const double distance = 1.5 + 28.5 * (0.5 + 0.5 * unit(rng));
            const Vec dir = random_unit(rng);
            const Vec centre{40 * unit(rng), 40 * unit(rng), 40 * unit(rng)};
            Panel target;
            target.x = static_cast<double>(centre.x + distance * dir.x);
            target.y = static_cast<double>(centre.y + distance * dir.y);
            target.z = static_cast<double>(centre.z + distance * dir.z);
            const Vec n = random_unit(rng);
            target.nx = static_cast<double>(n.x);
            target.ny = static_cast<double>(n.y);
            target.nz = static_cast<double>(n.z);
            TreeNode node;
            node.cx = static_cast<double>(centre.x);
            node.cy = static_cast<double>(centre.y);
            node.cz = static_cast<double>(centre.z);
            const double rho = distance * (0.1 + 0.6 * (0.5 + 0.5 * unit(rng)));
            NodeMoments scalar;
            std::array<NodeMoments, 3> normal{};
            const auto &indices = moment_indices(order);
            const std::size_t count = indices.size();
            for (std::size_t m = 0; m < count; ++m) {
              const auto &a = indices[m];
              const double scale = std::pow(rho, a.ax + a.ay + a.az);
              scalar.values[m] = scale * unit(rng);
              for (auto &component : normal)
                component.values[m] = scale * unit(rng);
            }
            // The reference uses the double-rounded target and centre so the
            // comparison measures the expansion, not input rounding.
            const Vec tL{target.x, target.y, target.z};
            const Vec nL{target.nx, target.ny, target.nz};
            const Vec cL{node.cx, node.cy, node.cz};
            const ChannelRef ref = reference_channels(
                tL, nL, cL, scalar, normal, eps, kappa, coupled, order);
            const JcpChannelAccum radial = taylor_jcp_node_channel_sum(
                target, node, scalar, normal, eps, kappa, coupled, order,
                "analytic");
            const JcpChannelAccum series = taylor_jcp_node_channel_sum_series(
                target, node, scalar, normal, eps, kappa, coupled, order,
                "analytic");
            for (std::size_t ch = 0; ch < 4; ++ch) {
              if (ref.scale[ch] == 0)
                continue;
              const double er = static_cast<double>(
                  std::fabs(channel_value(radial, ch) - ref.value[ch]) /
                  ref.scale[ch]);
              const double es = static_cast<double>(
                  std::fabs(channel_value(series, ch) - ref.value[ch]) /
                  ref.scale[ch]);
              radial_max = std::max(radial_max, er);
              series_max = std::max(series_max, es);
              radial_by_kappa[kappa_index] =
                  std::max(radial_by_kappa[kappa_index], er);
              series_by_kappa[kappa_index] =
                  std::max(series_by_kappa[kappa_index], es);
              // Below kr = 0.25 the screened-minus-Laplace channels are a
              // tiny fraction of either kernel (about (kr)^2/2), and their
              // own relative accuracy is what is measured here.
              const bool small_screening = kappa > 0 && kappa * distance < 0.25;
              if (small_screening) {
                radial_max_screened_small = std::max(radial_max_screened_small, er);
                series_max_screened_small = std::max(series_max_screened_small, es);
              }
              const double gate = small_screening ? 1e-9 : 1e-12;
              if (!std::isfinite(er) || er > gate) {
                if (failures < 5)
                  std::printf("COEFF FAIL order=%d kappa=%g eps=%g coupled=%d "
                              "channel=%zu radial_err=%.3e series_err=%.3e\n",
                              order, kappa, eps, coupled ? 1 : 0, ch, er, es);
                ++failures;
              }
            }
            // 3. the plain sum equals its channels.
            const JcpAccum sum = taylor_jcp_node_sum(
                target, node, scalar, normal, eps, kappa, coupled, order,
                "analytic");
            const JcpAccum from_channels = sum_jcp_channels(radial);
            const double sum_scale = static_cast<double>(
                ref.scale[0] + ref.scale[1] + ref.scale[2] + ref.scale[3]);
            const double sum_err =
                std::max(std::fabs(sum.phi - from_channels.phi),
                         std::fabs(sum.phin - from_channels.phin)) /
                std::max(sum_scale, 1e-300);
            if (!(sum_err <= 1e-14)) {
              if (failures < 5)
                std::printf("SUM FAIL order=%d kappa=%g err=%.3e\n", order,
                            kappa, sum_err);
              ++failures;
            }
            // 4. the fast far-field kernel, per row, against the reference.
            const JcpAccum quick = radial_jcp_far_kernel(order, coupled)(
                target, node, scalar.values.data(), normal[0].values.data(),
                normal[1].values.data(), normal[2].values.data(), eps, kappa);
            const LD rows[2][3] = {
                {ref.value[kJcpK1PhiChannel] + ref.value[kJcpK2PhiNormalChannel],
                 ref.scale[kJcpK1PhiChannel] + ref.scale[kJcpK2PhiNormalChannel],
                 static_cast<LD>(quick.phi)},
                {ref.value[kJcpK3PhinChannel] + ref.value[kJcpK4MixedChannel],
                 ref.scale[kJcpK3PhinChannel] + ref.scale[kJcpK4MixedChannel],
                 static_cast<LD>(quick.phin)}};
            const double plain[2] = {sum.phi, sum.phin};
            for (int row = 0; row < 2; ++row) {
              if (rows[row][1] == 0)
                continue;
              const double ef = static_cast<double>(
                  std::fabs(rows[row][2] - rows[row][0]) / rows[row][1]);
              const double eq = static_cast<double>(
                  std::fabs(rows[row][2] - plain[row]) / rows[row][1]);
              fast_max = std::max(fast_max, ef);
              fast_vs_radial_max = std::max(fast_vs_radial_max, eq);
              const bool small_screening = kappa > 0 && kappa * distance < 0.25;
              if (!std::isfinite(ef) || ef > (small_screening ? 1e-9 : 1e-12) ||
                  !(eq <= 1e-13)) {
                if (failures < 5)
                  std::printf("FAST FAIL order=%d kappa=%g eps=%g coupled=%d "
                              "row=%d ref_err=%.3e radial_gap=%.3e\n",
                              order, kappa, eps, coupled ? 1 : 0, row, ef, eq);
                ++failures;
              }
            }
            ++cases;
          }
  std::printf("coefficients: cases=%zu max normalised error radial=%.3e "
              "series=%.3e; screened-small (0<kappa*r<0.25) radial=%.3e "
              "series=%.3e\n",
              cases, radial_max, series_max, radial_max_screened_small,
              series_max_screened_small);
  for (std::size_t k = 0; k < kKappaCount; ++k)
    std::printf("  kappa=%-7g radial max=%.3e  series max=%.3e\n", kappas[k],
                radial_by_kappa[k], series_by_kappa[k]);
  std::printf("fast far kernel: max normalised row error %.3e, max gap to "
              "taylor_jcp_node_sum %.3e\n",
              fast_max, fast_vs_radial_max);

  // ---- 4. four-node lanes against the one-node kernel ---------------------
  // Each lane must carry exactly the one-node result, degenerate lanes
  // included, so the apply's answer does not depend on how the far list is
  // cut into fours.
  {
    std::size_t lanes = 0, differing = 0;
    double lane_gap = 0, scalar_gap = 0;
    const RadialJcpFarKernel scalar_kernels[2][4] = {
        {&radial_fast::node_sum<2, false>, &radial_fast::node_sum<3, false>,
         &radial_fast::node_sum<4, false>, &radial_fast::node_sum<5, false>},
        {&radial_fast::node_sum<2, true>, &radial_fast::node_sum<3, true>,
         &radial_fast::node_sum<4, true>, &radial_fast::node_sum<5, true>}};
    for (int order = 2; order <= 5; ++order)
      for (const double kappa : kappas)
        for (const bool coupled : {false, true})
          for (int sample = 0; sample < 20; ++sample) {
            const RadialJcpFarKernel one = radial_jcp_far_kernel(order, coupled);
            const RadialJcpFarKernel4 four =
                radial_jcp_far_kernel4(order, coupled);
            if (four == nullptr) {
              std::printf("FOUR FAIL: no four-node kernel for order %d\n",
                          order);
              ++failures;
              continue;
            }
            Panel target;
            target.x = 30 * unit(rng);
            target.y = 30 * unit(rng);
            target.z = 30 * unit(rng);
            const Vec n = random_unit(rng);
            target.nx = static_cast<double>(n.x);
            target.ny = static_cast<double>(n.y);
            target.nz = static_cast<double>(n.z);
            TreeNode nodes[4];
            NodeMoments density[4][4];
            const TreeNode *node_ptrs[4];
            const double *density_ptrs[16];
            for (int lane = 0; lane < 4; ++lane) {
              const double distance = 1.5 + 28.5 * (0.5 + 0.5 * unit(rng));
              const Vec dir = random_unit(rng);
              nodes[lane].cx = target.x - distance * static_cast<double>(dir.x);
              nodes[lane].cy = target.y - distance * static_cast<double>(dir.y);
              nodes[lane].cz = target.z - distance * static_cast<double>(dir.z);
              if (sample == 0 && lane == 2) { // degenerate: node on the target
                nodes[lane].cx = target.x;
                nodes[lane].cy = target.y;
                nodes[lane].cz = target.z;
              }
              const double rho = distance * 0.5;
              const auto &indices = moment_indices(order);
              for (std::size_t m = 0; m < indices.size(); ++m) {
                const auto &a = indices[m];
                const double scale = std::pow(rho, a.ax + a.ay + a.az);
                for (int kind = 0; kind < 4; ++kind)
                  density[lane][kind].values[m] = scale * unit(rng);
              }
              node_ptrs[lane] = &nodes[lane];
              for (int kind = 0; kind < 4; ++kind)
                density_ptrs[4 * lane + kind] = density[lane][kind].values.data();
            }
            const double eps = epss[sample % 4];
            JcpAccum got[4];
            four(target, node_ptrs, density_ptrs, eps, kappa, got);
            for (int lane = 0; lane < 4; ++lane) {
              const JcpAccum want = one(
                  target, nodes[lane], density[lane][0].values.data(),
                  density[lane][1].values.data(), density[lane][2].values.data(),
                  density[lane][3].values.data(), eps, kappa);
              ++lanes;
              const double scale =
                  std::max({std::fabs(want.phi), std::fabs(want.phin), 1e-300});
              if (std::memcmp(&got[lane].phi, &want.phi, sizeof(double)) != 0 ||
                  std::memcmp(&got[lane].phin, &want.phin, sizeof(double)) != 0) {
                ++differing;
                lane_gap = std::max(
                    lane_gap, std::max(std::fabs(got[lane].phi - want.phi),
                                       std::fabs(got[lane].phin - want.phin)) /
                                  scale);
              }
              const JcpAccum scalar = scalar_kernels[coupled][order - 2](
                  target, nodes[lane], density[lane][0].values.data(),
                  density[lane][1].values.data(), density[lane][2].values.data(),
                  density[lane][3].values.data(), eps, kappa);
              scalar_gap = std::max(
                  scalar_gap, std::max(std::fabs(scalar.phi - want.phi),
                                       std::fabs(scalar.phin - want.phin)) /
                                  scale);
            }
          }
    std::printf("four-node lanes: %zu lanes, %zu differ from the one-node "
                "kernel (max relative gap %.3e); scalar fallback max relative "
                "gap %.3e\n",
                lanes, differing, lane_gap, scalar_gap);
    if (differing != 0) {
      std::printf("FOUR FAIL: lanes are not bit-identical to the one-node "
                  "kernel\n");
      ++failures;
    }
    if (!(scalar_gap <= 1e-12)) {
      std::printf("FOUR FAIL: the scalar fallback kernel is not the one-node "
                  "kernel to rounding\n");
      ++failures;
    }
  }

  // ---- 2. expansion against the direct sum --------------------------------
  // The truncation error of one configuration need not fall at every order
  // (an even-order term can vanish by symmetry), so the order check is on the
  // geometric mean over the configurations at each opening angle.
  double worst_gap_ratio = 0;
  constexpr double thetas[3] = {0.3, 0.5, 0.7};
  std::array<std::array<double, 6>, 3> log_error_sum{};
  std::array<int, 3> configurations{};
  for (const double kappa : {0.0, 0.1235, 1.0})
    for (std::size_t theta_index = 0; theta_index < 3; ++theta_index)
      for (const bool coupled : {false, true})
        for (int sample = 0; sample < 4; ++sample) {
          const double theta = thetas[theta_index];
          ++configurations[theta_index];
          const double rho = 2.0;
          const Vec centre{5 * unit(rng), 5 * unit(rng), 5 * unit(rng)};
          std::vector<Vec> sources;
          std::vector<LD> q;
          std::vector<std::array<LD, 3>> w;
          for (int s = 0; s < 60; ++s) {
            Vec p;
            do {
              p = {unit(rng), unit(rng), unit(rng)};
            } while (p.x * p.x + p.y * p.y + p.z * p.z > 1);
            sources.push_back({centre.x + rho * p.x, centre.y + rho * p.y,
                               centre.z + rho * p.z});
            q.push_back(unit(rng));
            const Vec nn = random_unit(rng);
            const LD phi = unit(rng);
            w.push_back({phi * nn.x, phi * nn.y, phi * nn.z});
          }
          const Vec dir = random_unit(rng);
          const double distance = rho / theta;
          Panel target;
          target.x = static_cast<double>(centre.x + distance * dir.x);
          target.y = static_cast<double>(centre.y + distance * dir.y);
          target.z = static_cast<double>(centre.z + distance * dir.z);
          const Vec n = random_unit(rng);
          target.nx = static_cast<double>(n.x);
          target.ny = static_cast<double>(n.y);
          target.nz = static_cast<double>(n.z);
          const Vec tL{target.x, target.y, target.z};
          const Vec nL{target.nx, target.ny, target.nz};
          TreeNode node;
          node.cx = static_cast<double>(centre.x);
          node.cy = static_cast<double>(centre.y);
          node.cz = static_cast<double>(centre.z);
          const double eps = 39.25;
          LD direct[4]{}, direct_scale = 0;
          for (std::size_t s = 0; s < sources.size(); ++s) {
            const PointChannels pc = point_channels(
                tL, nL, sources[s], q[s], w[s].data(), eps, kappa, coupled);
            for (int ch = 0; ch < 4; ++ch) {
              direct[ch] += pc.v[ch];
              direct_scale += std::fabs(pc.v[ch]);
            }
          }
          for (int order = 2; order <= 5; ++order) {
            const auto &indices = moment_indices(order);
            NodeMoments scalar;
            std::array<NodeMoments, 3> normal{};
            for (std::size_t m = 0; m < indices.size(); ++m) {
              const auto &a = indices[m];
              LD sq = 0, sn[3] = {0, 0, 0};
              for (std::size_t s = 0; s < sources.size(); ++s) {
                const LD h = std::pow(sources[s].x - node.cx, (LD)a.ax) *
                             std::pow(sources[s].y - node.cy, (LD)a.ay) *
                             std::pow(sources[s].z - node.cz, (LD)a.az);
                sq += q[s] * h;
                for (int c = 0; c < 3; ++c)
                  sn[c] += w[s][c] * h;
              }
              scalar.values[m] = static_cast<double>(sq);
              for (int c = 0; c < 3; ++c)
                normal[c].values[m] = static_cast<double>(sn[c]);
            }
            const JcpChannelAccum radial = taylor_jcp_node_channel_sum(
                target, node, scalar, normal, eps, kappa, coupled, order,
                "analytic");
            const JcpChannelAccum series = taylor_jcp_node_channel_sum_series(
                target, node, scalar, normal, eps, kappa, coupled, order,
                "analytic");
            double truncation = 0, gap = 0;
            for (std::size_t ch = 0; ch < 4; ++ch) {
              truncation = std::max(
                  truncation,
                  static_cast<double>(std::fabs(channel_value(radial, ch) -
                                                direct[ch]) /
                                      direct_scale));
              gap = std::max(gap, std::fabs(channel_value(radial, ch) -
                                            channel_value(series, ch)) /
                                      static_cast<double>(direct_scale));
            }
            worst_gap_ratio =
                std::max(worst_gap_ratio, gap / std::max(truncation, 1e-300));
            log_error_sum[theta_index][order] +=
                std::log(std::max(truncation, 1e-300));
            if (gap > 1e-9 * truncation + 1e-15) {
              if (failures < 5)
                std::printf("EXPANSION FAIL kappa=%g theta=%g order=%d "
                            "truncation=%.3e radial-series=%.3e\n",
                            kappa, theta, order, truncation, gap);
              ++failures;
            }
          }
        }
  std::printf("expansion: worst |radial-series|/truncation=%.3e\n",
              worst_gap_ratio);
  for (std::size_t t = 0; t < 3; ++t) {
    std::printf("  theta=%.1f geometric-mean truncation by order:", thetas[t]);
    for (int order = 2; order <= 5; ++order) {
      const double mean =
          std::exp(log_error_sum[t][order] / configurations[t]);
      std::printf(" p%d=%.3e", order, mean);
      if (order > 2 && t < 2 &&
          !(log_error_sum[t][order] < log_error_sum[t][order - 1])) {
        std::printf(" (FAIL: did not fall)");
        ++failures;
      }
    }
    std::printf("\n");
  }

  // ---- timing (information only) -------------------------------------------
  // Old and new alternate in blocks, order flipped per round; the median of
  // the per-round ratios is what is reported. --no-timing (the CI gate) skips
  // it.
  if (!skip_timing) {
    constexpr int kPairs = 4000;
    std::vector<Panel> targets(kPairs);
    std::vector<TreeNode> nodes(kPairs);
    for (int i = 0; i < kPairs; ++i) {
      const Vec dir = random_unit(rng), n = random_unit(rng);
      const double distance = 3 + 20 * (0.5 + 0.5 * unit(rng));
      nodes[i].cx = 10 * unit(rng);
      nodes[i].cy = 10 * unit(rng);
      nodes[i].cz = 10 * unit(rng);
      targets[i].x = nodes[i].cx + distance * static_cast<double>(dir.x);
      targets[i].y = nodes[i].cy + distance * static_cast<double>(dir.y);
      targets[i].z = nodes[i].cz + distance * static_cast<double>(dir.z);
      targets[i].nx = static_cast<double>(n.x);
      targets[i].ny = static_cast<double>(n.y);
      targets[i].nz = static_cast<double>(n.z);
    }
    for (const int order : {2, 4, 5}) {
      NodeMoments scalar;
      std::array<NodeMoments, 3> normal{};
      for (std::size_t m = 0; m < moment_indices(order).size(); ++m) {
        scalar.values[m] = unit(rng);
        for (auto &c : normal)
          c.values[m] = unit(rng);
      }
      std::vector<double> ratios;
      double sink = 0, best_radial = 1e300, best_series = 1e300;
      for (int round = 0; round < 7; ++round) {
        double seconds[2] = {0, 0};
        for (const bool radial : {round % 2 == 0, round % 2 != 0}) {
          const auto start = std::chrono::steady_clock::now();
          for (int i = 0; i < kPairs; ++i) {
            const JcpAccum a =
                radial ? taylor_jcp_node_sum(targets[i], nodes[i], scalar,
                                             normal, 39.25, 0.1235, true,
                                             order, "analytic")
                       : taylor_jcp_node_sum_series(
                             targets[i], nodes[i], scalar, normal, 39.25,
                             0.1235, true, order, "analytic");
            sink += a.phi + a.phin;
          }
          seconds[radial] = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - start)
                                .count();
        }
        best_series = std::min(best_series, seconds[0]);
        best_radial = std::min(best_radial, seconds[1]);
        ratios.push_back(seconds[0] / seconds[1]);
      }
      std::sort(ratios.begin(), ratios.end());
      std::printf("timing order=%d coupled: series %.2f us/pair, radial %.2f "
                  "us/pair, median ratio %.1fx (checksum %.3g)\n",
                  order, 1e6 * best_series / kPairs,
                  1e6 * best_radial / kPairs, ratios[ratios.size() / 2], sink);
    }
    // The apply's far loop: taylor_jcp_node_sum with the moments copied, as
    // the loop did, against the fast kernel reading them in place.
    for (const int order : {2, 3, 4, 5})
      for (const bool coupled : {true, false}) {
        if (!coupled && order != 4)
          continue;
        std::vector<NodeMoments> scalar(kPairs);
        std::array<std::vector<NodeMoments>, 3> normal;
        for (auto &c : normal)
          c.resize(kPairs);
        for (int i = 0; i < kPairs; ++i)
          for (std::size_t m = 0; m < moment_indices(order).size(); ++m) {
            scalar[i].values[m] = unit(rng);
            for (auto &c : normal)
              c[i].values[m] = unit(rng);
          }
        const std::string mode = "analytic";
        const JcpFarEvaluator evaluator =
            make_jcp_far_evaluator(order, coupled, mode);
        std::vector<double> ratios;
        double sink[2] = {0, 0}, best[2] = {1e300, 1e300};
        for (int round = 0; round < 9; ++round) {
          double seconds[2] = {0, 0};
          for (const bool fast : {round % 2 == 0, round % 2 != 0}) {
            const auto start = std::chrono::steady_clock::now();
            for (int repeat = 0; repeat < 4; ++repeat)
              for (int i = 0; i < kPairs; ++i) {
                JcpAccum a;
                if (fast) {
                  a = evaluator(targets[i], nodes[i], scalar[i], normal[0][i],
                                normal[1][i], normal[2][i], 39.25, 0.1235);
                } else {
                  const std::array<NodeMoments, 3> copied = {
                      normal[0][i], normal[1][i], normal[2][i]};
                  a = taylor_jcp_node_sum(targets[i], nodes[i], scalar[i],
                                          copied, 39.25, 0.1235, coupled,
                                          order, mode);
                }
                sink[fast] += a.phi + a.phin;
              }
            seconds[fast] = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - start)
                                .count();
          }
          best[0] = std::min(best[0], seconds[0]);
          best[1] = std::min(best[1], seconds[1]);
          ratios.push_back(seconds[0] / seconds[1]);
        }
        std::sort(ratios.begin(), ratios.end());
        std::printf("timing order=%d %s: radial %.3f us/pair, fast %.3f "
                    "us/pair, median ratio %.2fx, fast kernel %s "
                    "(checksums %.15g %.15g)\n",
                    order, coupled ? "coupled" : "uncoupled",
                    1e6 * best[0] / (4 * kPairs), 1e6 * best[1] / (4 * kPairs),
                    ratios[ratios.size() / 2],
                    evaluator.fast ? "on" : "off (TITAN_FORCE_SCALAR)",
                    sink[0], sink[1]);
        // One-node kernel against four-node lanes on the same pairs, the
        // target of each group of four shared as in the apply.
        if (evaluator.fast == nullptr || evaluator.fast4 == nullptr)
          continue;
        std::vector<double> lane_ratios;
        double lane_sink[2] = {0, 0}, lane_best[2] = {1e300, 1e300};
        for (int round = 0; round < 9; ++round) {
          double seconds[2] = {0, 0};
          for (const bool four : {round % 2 == 0, round % 2 != 0}) {
            const auto start = std::chrono::steady_clock::now();
            for (int repeat = 0; repeat < 4; ++repeat)
              for (int i = 0; i + 4 <= kPairs; i += 4) {
                const Panel &target = targets[i];
                if (four) {
                  const TreeNode *node_ptrs[4];
                  const double *density_ptrs[16];
                  for (int lane = 0; lane < 4; ++lane) {
                    node_ptrs[lane] = &nodes[i + lane];
                    density_ptrs[4 * lane] = scalar[i + lane].values.data();
                    for (int c = 0; c < 3; ++c)
                      density_ptrs[4 * lane + 1 + c] =
                          normal[c][i + lane].values.data();
                  }
                  JcpAccum got[4];
                  evaluator.fast4(target, node_ptrs, density_ptrs, 39.25,
                                  0.1235, got);
                  for (const JcpAccum &a : got)
                    lane_sink[1] += a.phi + a.phin;
                } else {
                  for (int lane = 0; lane < 4; ++lane) {
                    const JcpAccum a = evaluator.fast(
                        target, nodes[i + lane], scalar[i + lane].values.data(),
                        normal[0][i + lane].values.data(),
                        normal[1][i + lane].values.data(),
                        normal[2][i + lane].values.data(), 39.25, 0.1235);
                    lane_sink[0] += a.phi + a.phin;
                  }
                }
              }
            seconds[four] = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - start)
                                .count();
          }
          lane_best[0] = std::min(lane_best[0], seconds[0]);
          lane_best[1] = std::min(lane_best[1], seconds[1]);
          lane_ratios.push_back(seconds[0] / seconds[1]);
        }
        std::sort(lane_ratios.begin(), lane_ratios.end());
        std::printf("timing order=%d %s: fast %.3f us/pair, four-lane %.3f "
                    "us/pair, median ratio %.2fx (checksums %.15g %.15g)\n",
                    order, coupled ? "coupled" : "uncoupled",
                    1e6 * lane_best[0] / (4 * kPairs),
                    1e6 * lane_best[1] / (4 * kPairs),
                    lane_ratios[lane_ratios.size() / 2], lane_sink[0],
                    lane_sink[1]);
      }
  }

  if (failures != 0) {
    std::printf("FAIL tabi_treecode_radial_oracle (%d)\n", failures);
    return 1;
  }
  std::printf("PASS tabi_treecode_radial_oracle\n");
  return 0;
}
