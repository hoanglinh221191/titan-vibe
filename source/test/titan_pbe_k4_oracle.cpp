#include "tabi/tabi_fmm_m2l_local.h"
#include "tabi/tabi_near.h"
#include "tabi/tabi_preconditioner.h"
#include "tabi/tabi_treecode_internal.h"
#include "titan_executor.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace T = TitanPBE::Tabi;
namespace D = TitanPBE::Tabi::detail;
namespace F = titan::tabi::fmm_m2l;
constexpr double pi = 3.14159265358979323846;

void check(const char *name, double value, double ref, double tol = 2e-7) {
  if (!std::isfinite(value) ||
      std::abs(value - ref) > tol * std::max(1e-8, std::abs(ref))) {
    std::cerr << name << " got=" << value << " expected=" << ref << '\n';
    throw std::runtime_error(name);
  }
}

// Independent four-point derivative of Gk-G0. No production derivative formula.
double fd(const T::Panel &t, const T::Panel &s, double k) {
  auto central = [&](long double h) {
    long double sum = 0.0;
    for (double a : {-1.0, 1.0})
      for (double b : {-1.0, 1.0}) {
        const long double dx = t.x + a * h * t.nx - s.x - b * h * s.nx;
        const long double dy = t.y + a * h * t.ny - s.y - b * h * s.ny;
        const long double dz = t.z + a * h * t.nz - s.z - b * h * s.nz;
        const long double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        sum += a * b * std::expm1(-k * r) / (4 * pi * r);
      }
    return s.area * sum / (4 * h * h);
  };
  return static_cast<double>((4 * central(0.001L) - central(0.002L)) / 3);
}

void pair_oracle() {
  TitanExecutor executor(1);
  T::Panel t, s;
  t.x = 7;
  t.y = -3;
  t.z = 5;
  t.nx = 0.6;
  t.ny = 0;
  t.nz = 0.8;
  s.x = -2;
  s.y = 1;
  s.z = 0.5;
  s.area = 1.7;
  for (const auto normal :
       {std::array<double, 3>{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.6, 0.8, 0}}) {
    s.nx = normal[0];
    s.ny = normal[1];
    s.nz = normal[2];
    for (double k : {0.0, 0.127, 0.3, 0.8}) {
      const double ref = fd(t, s, k);
      const auto near = T::evaluate_jcp_centroid_monopole_kernel(t, s, 80, k);
      check("near K4 vs FD", near.k4, ref);
      double scalar = 0, series = 0;
      for (int c = 0; c < 3; ++c) {
        scalar += normal[c] *
                  (D::mixed_normal_component_green(t, s.x, s.y, s.z, c, k) -
                   D::mixed_normal_component_green(t, s.x, s.y, s.z, c, 0));
        const auto a =
            D::mixed_normal_component_green_series(t, s.x, s.y, s.z, c, k, 4);
        const auto b =
            D::mixed_normal_component_green_series(t, s.x, s.y, s.z, c, 0, 4);
        series += normal[c] * (a.coeff[0] - b.coeff[0]);
      }
      check("scalar component K4 vs FD", s.area * scalar, ref);
      check("Taylor K4 vs FD", s.area * series, ref);
      F::CartesianJcpMoments moments;
      moments.order = 4;
      moments.term_count = F::cartesian_moment_term_count(4);
      for (int c = 0; c < 3; ++c)
        moments.phi_normal[c][0] = normal[c] * s.area;
      F::PrototypePanel target;
      target.x = t.x;
      target.y = t.y;
      target.z = t.z;
      target.nx = t.nx;
      target.ny = t.ny;
      target.nz = t.nz;
      check("Cartesian M2P K4 vs FD",
            F::evaluate_cartesian_jcp_m2p(target, s.x, s.y, s.z, moments, 80,
                                          k)[3],
            ref);
      const auto op = F::build_cartesian_jcp_m2l_operator_p4(s.x, s.y, s.z, t.x,
                                                             t.y, t.z, 80, k);
      const auto local =
          F::apply_cartesian_jcp_m2l_operator_p4_order(moments, op, 4);
      check("Cartesian M2L K4 vs FD",
            F::evaluate_cartesian_jcp_l2p_order(target, t.x, t.y, t.z, local,
                                                4)[3],
            ref);
      const auto compressed =
          F::build_cartesian_jcp_compressed_green_operator_p4(s.x, s.y, s.z,
                                                              t.x, t.y, t.z, k);
      const auto packed =
          F::apply_cartesian_jcp_compressed_green_operator_p4_order(
              moments, compressed, 80, 4);
      check("compressed/BIE M2L K4 vs FD",
            F::evaluate_cartesian_jcp_l2p_order(target, t.x, t.y, t.z, packed,
                                                4)[3],
            ref);
      F::PrototypePanel source;
      source.x = s.x;
      source.y = s.y;
      source.z = s.z;
      source.area = s.area;
      source.nx = s.nx;
      source.ny = s.ny;
      source.nz = s.nz;
      for (int leaf_size : {1, 64}) {
        F::PrototypeConfig cfg;
        cfg.max_leaf_size = leaf_size;
        const auto result = F::evaluate_jcp_offline_operator_analytic(
            cfg, {target, source}, {0.0, 1.0}, {0.0, 0.0}, 80, k, executor);
        check("FMM analytic operator K4 vs FD", result.phin.at(0), ref);
      }
      // Exercise the independent centroid formula used to build LU blocks.
      auto pt = t, ps = s;
      pt.has_triangle_geometry = ps.has_triangle_geometry = true;
      pt.v0_index = ps.v0_index = 0;
      const auto precond = T::build_jcp_block_preconditioner_from_panel_groups(
          {pt, ps}, true, 80, k, 40.5, 0.50625, {40.5, 40.5, 0.50625, 0.50625},
          "centroid", "off", "singular", {{0, 1}}, executor);
      if (precond->valid_blocks != 1)
        throw std::runtime_error("preconditioner block not exercised");
      std::vector<double> recovered;
      precond->apply({-near.k2, 40.5, -ref, 0.0}, recovered);
      for (std::size_t i = 0; i < 4; ++i)
        if (std::abs(recovered[i] - (i == 1 ? 1.0 : 0.0)) > 1e-9)
          throw std::runtime_error("preconditioner K4 vs independent FD RHS");
      std::vector<T::Panel> panels(9, s);
      panels[0] = t;
      std::vector<std::size_t> order(9);
      std::iota(order.begin(), order.end(), 0);
      std::vector<double> unknowns(18, 0);
      for (int j = 1; j < 9; ++j)
        unknowns[j] = 1;
      const auto pc = D::build_tree_panel_cache(panels, order);
      const auto uc =
          D::build_tree_unknown_cache(panels, unknowns, true, order);
      for (const std::string kernel : {"scalar", "avx2"}) {
        std::size_t count = 0;
        const auto leaf = D::direct_jcp_ordered_sum(t, 0, pc, uc, 80, k, true,
                                                    order, 1, 8, kernel, count);
        check("leaf K4 vs FD", leaf.phin, 8 * ref);
      }
      T::TreeNode node;
      node.first_panel = 1;
      node.panel_count = 8;
      for (const std::string kernel : {"tile-scalar", "tile-avx2"}) {
        std::size_t count = 0;
        const auto leaf = D::direct_jcp_ordered_tile_sum(
            t, 0, pc, uc, 80, k, true, node, kernel, count);
        check("leaf tile K4 vs FD", leaf.phin, 8 * ref);
      }
      const std::vector<std::size_t> targets{0, 0, 0, 0};
      for (const std::string kernel :
           {"scalar", "target-avx2", "target-avx2-v2", "target-avx2-4x4",
            "soa-avx2", "fused-soa-avx2"}) {
        std::vector<D::JcpAccum> accum;
        std::size_t count = 0, vectorized = 0, lanes = 0, soa = 0;
        D::direct_jcp_target_batch_tile_sum(panels, targets, pc, uc, 80, k,
                                            true, node, kernel, accum, count,
                                            vectorized, lanes, true, soa);
        for (const auto &v : accum)
          check("batch AVX2 K4 vs FD", v.phin, 8 * ref);
        if (kernel != "scalar" && vectorized == 0 && soa == 0)
          throw std::runtime_error(
              "AVX2 batch oracle did not exercise vector path");
      }
      for (bool vector_exp : {false, true}) {
        std::vector<D::JcpAccum> accum;
        std::size_t count = 0, vectorized = 0, lanes = 0;
        D::direct_jcp_packed_target_source_span_avx2(
            panels, targets, pc, uc, 80, k, true, 1, 8, accum, count,
            vectorized, lanes, vector_exp);
        for (const auto &v : accum)
          check("packed AVX2 K4 vs FD", v.phin, 8 * ref);
        if (vectorized == 0)
          throw std::runtime_error("packed AVX2 not exercised");
      }
    }
  }
}

// Substitute exact Born boundary traces into BOTH coupled Juffer rows, using
// production near kernels and sphere quadrature. A parity test cannot do this.
void born_rows() {
  double worst = 0;
  for (double epsin : {1.0, 4.0})
    for (double salt : {0.0, 0.15, 0.5}) {
      const double a = 6, epsout = 80, eps = epsout / epsin,
                   k = 0.3291 * std::sqrt(salt);
      const double phi = 1 / (4 * pi * epsout * a * (1 + k * a));
      const double phin = -1 / (4 * pi * epsin * a * a);
      T::Panel t;
      t.z = a;
      t.nz = 1;
      t.nx = 0;
      t.ny = 0;
      std::array<double, 4> integral{};
      constexpr int n = 8192;
      for (int i = 0; i < n; ++i) {
        const double th = (i + 0.5) * pi / n;
        T::Panel s;
        s.x = a * std::sin(th);
        s.z = a * std::cos(th);
        s.y = 0;
        s.nx = std::sin(th);
        s.nz = std::cos(th);
        s.ny = 0;
        s.area = 2 * pi * a * a * std::sin(th) * pi / n;
        auto v = T::evaluate_jcp_centroid_monopole_kernel(t, s, eps, k);
        integral[0] += v.k1;
        integral[1] += v.k2;
        integral[2] += v.k3;
        integral[3] += v.k4;
      }
      const double s1 = 1 / (4 * pi * epsin * a), s2 = phin;
      const double r1 = (0.5 * (1 + eps) * phi - integral[0] * phin -
                         integral[1] * phi - s1) /
                        s1;
      const double r2 = (0.5 * (1 + 1 / eps) * phin - integral[2] * phin -
                         integral[3] * phi - s2) /
                        s2;
      worst = std::max({worst, std::abs(r1), std::abs(r2)});
      if (std::abs(r1) > 1e-6 || std::abs(r2) > 1e-6)
        throw std::runtime_error("Born coupled row oracle");
    }
  std::cout << "PASS TABI-01 Born coupled rows worst_relative_residual="
            << worst << '\n';
}

int main() {
  try {
    pair_oracle();
    born_rows();
    std::cout
        << "PASS TABI-01 FD near/leaf/all AVX2 "
           "variants/scalar/Taylor/Cartesian/BIE/FMM/preconditioner kernels\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
  }
}
