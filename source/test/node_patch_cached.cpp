#include "tabi/tabi_near.h"
#include "tabi/tabi_node_patch.h"
#include "tabi/tabi_node_patch_curved.h"
#include "titan_executor.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
using namespace TitanPBE::Tabi;
int main(int argc, char **argv) {
  try {
    TitanExecutor executor(4);
    Options o;
    o.primary_surface = "sphere-patch";
    o.surface_mode = "vdw";
    o.sphere_mesh_mode = "icosphere";
    o.surface_point_count = 96;
    o.discretization = "node-patch";
    o.salt_m = .15;
    o.eps_in = 1;
    o.eps_out = 78.5;
    o.near_integration = "analytic_singular";
    if (argc > 1) {
      o.surface_point_count = std::stoi(argv[1]);
      double offset = argc > 2 ? std::stod(argv[2]) : .8;
      o.salt_m = argc > 3 ? std::stod(argv[3]) : .15;
      o.residual_tol = 1e-10;
      o.max_iter = 300;
      o.preconditioner = "diagonal";
      o.sum_mode = "direct";
      if (argc > 4) {
        o.sum_mode = "treecode";
        o.formulation = "jcp2013-treecode";
        o.tree_opening_theta = std::stod(argv[4]);
        o.tree_opening_theta_explicit = true;
        o.tree_leaf_size = 16;
        o.tree_leaf_size_explicit = true;
        o.tree_builder = "topdown-aabb";
        o.moment_audit = "sampled";
        o.leaf_mac = "on";
        o.taylor_order = 5;
      }
      auto start = std::chrono::steady_clock::now();
      auto r = solve_linear_pb(
          {{0, 0, 0, 0, 1.7}, {offset * 1.7, 0, 0, 1, .017}}, o, executor);
      std::cout << std::setprecision(17)
                << "RESULT {\"energy\":" << r.energy_kcal_mol
                << ",\"residual\":" << r.true_residual_rel
                << ",\"iterations\":" << r.iterations
                << ",\"panels\":" << r.panel_count << ",\"discretization\":\""
                << r.discretization << "\",\"far\":" << r.far_interactions
                << ",\"moment_error\":" << r.moment_translation_rel_error_max
                << ",\"moment_audit_nodes\":" << r.moment_audit_nodes
                << ",\"gmres_ms\":" << r.gmres_ms
                << ",\"matvec_ms_total\":" << r.matvec_ms_total
                << ",\"matvec_ms_avg\":" << r.matvec_ms_avg
                << ",\"mesh_ms\":" << r.surface_mesh_ms
                << ",\"projection_ms\":" << r.energy_projection_ms
                << ",\"memory_bytes\":" << r.memory_bytes << ",\"seconds\":"
                << std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count()
                << "}\n";
      if (!r.converged || r.discretization != "node_patch" ||
          r.near_kernel_policy != "curved_dual")
        throw std::runtime_error("curved production solve failed");
      return 0;
    }
    int checks = 0;
    auto require = [&](bool value, const char *name) {
      ++checks;
      if (!value)
        throw std::runtime_error(name);
    };
    const auto triangles =
        build_surface({{0, 0, 0, 0, 1.7}}, o, nullptr, &executor);
    auto cached = build_node_patch_panels(triangles),
         streamed = build_node_patch_panels(triangles);
    double kappa = .3291 * std::sqrt(.15);
    require(prepare_curved_node_patches(triangles, cached, 78.5, kappa,
                                        64 * 1024 * 1024, executor),
            "sphere detected");
    require(prepare_curved_node_patches(triangles, streamed, 78.5, kappa, 0,
                                        executor),
            "streamed sphere detected");
    require(!cached[0].node_patch_geometry->curved->pair_blocks.empty(),
            "pair cache built");
    require(streamed[0].node_patch_geometry->curved->pair_blocks.empty(),
            "zero budget streams");
    double area = 0, nx = 0, ny = 0, nz = 0, xnx = 0;
    for (const auto &p : cached) {
      const auto &c = *curved_dual_cell(p);
      area += c.area;
      nx += c.normal_integral[0];
      ny += c.normal_integral[1];
      nz += c.normal_integral[2];
      // Index through kMomentSlot: the moment array is packed to the 56 monomials of degree <= 5, not the 6x6x6 lattice.
      xnx += c.moments[1][kMomentSlot[1 * 36]] + c.center[0] * c.moments[1][0];
      require(std::abs(c.moments[0][0] - c.area) < 1e-12, "zeroth moment area");
    }
    const double spherearea = 4 * std::numbers::pi * 1.7 * 1.7;
    require(std::abs(area / spherearea - 1) < 1e-8,
            "curved area partitions sphere");
    require(std::abs(nx) + std::abs(ny) + std::abs(nz) < 1e-11,
            "normal moment conservation");
    require(std::abs(xnx / (spherearea * 1.7 / 3) - 1) < 1e-8,
            "analytic first normal moment");
    for (std::size_t i = 0; i < cached.size(); i += 9)
      for (std::size_t j = 0; j < cached.size(); j += 7) {
        auto a = curved_node_pair(cached[i], i, cached[j], j, 78.5, kappa),
             b = curved_node_pair(streamed[i], i, streamed[j], j, 78.5, kappa);
        for (int k = 0; k < 4; ++k)
          require(a[k] == b[k], "cached and streamed blocks identical");
        auto reordered =
            curved_node_pair(cached[i], 0, cached[j], 1, 78.5, kappa);
        require(a == reordered,
                "cache uses geometry identity, not local indexing");
      }
    bool rejected = false;
    try {
      curved_node_pair(cached[0], 0, cached[1], 1, 2, kappa);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    require(rejected, "cache rejects dielectric mismatch");
    auto lazy = build_node_patch_panels(triangles);
    const std::size_t tiny_budget = 64 * sizeof(CurvedKernelCacheEntry);
    require(prepare_curved_node_patches(triangles, lazy, 78.5, kappa,
                                        tiny_budget, executor, false),
            "lazy cache prepared");
    require(lazy[0].node_patch_geometry->curved->sparse_capacity *
                    sizeof(CurvedKernelCacheEntry) <=
                tiny_budget,
            "sparse cache obeys budget");
    std::vector<unsigned char> matched(512);
    executor.parallel_for(
        matched.size(), 4, 1,
        [&](std::size_t, std::size_t begin, std::size_t end) {
          for (auto task = begin; task < end; ++task) {
            // First contend for one previously empty non-neighbor entry,
            // then overfill the deliberately tiny table with distinct pairs.
            std::size_t i = task < 128 ? 0 : task % 7;
            std::size_t j = task < 128 ? 1 : (task * 37) % cached.size();
            auto a = curved_node_pair(lazy[i], i, lazy[j], j, 78.5, kappa),
                 b = curved_node_pair(cached[i], i, cached[j], j, 78.5, kappa);
            matched[task] = a == b;
          }
        });
    require(std::all_of(matched.begin(), matched.end(),
                        [](auto x) { return x != 0; }),
            "concurrent lazy cache and saturated fallback preserve kernels");
    std::swap(cached[1], cached[2]);
    std::swap(streamed[1], streamed[2]);
    std::vector<std::size_t> order(cached.size());
    std::vector<double> values(2 * cached.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
      order[i] = order.size() - 1 - i;
      values[i] = std::sin(double(i));
      values[order.size() + i] = std::cos(double(i));
    }
    for (int repeat = 0; repeat < 2; ++repeat) {
      auto a = accumulate_jcp_near_ordered(cached[0], 0, cached, values, 78.5,
                                           kappa, true, order, 0, order.size(),
                                           "analytic_singular", "off");
      auto b = accumulate_jcp_near_ordered(
          streamed[0], 0, streamed, values, 78.5, kappa, true, order, 0,
          order.size(), "analytic_singular", "off");
      require(std::abs(a.phi - b.phi) < 1e-12 &&
                  std::abs(a.phin - b.phin) < 1e-12,
              "reordered cached matvec and changed trace match streamed "
              "evaluation");
      for (auto &v : values)
        v = v * .31 + .17;
    }
    auto flat_triangles = triangles;
    for (auto &p : flat_triangles)
      p.sphere_patch_quadrature.reset();
    auto flat = build_node_patch_panels(flat_triangles);
    require(!prepare_curved_node_patches(flat_triangles, flat, 78.5, kappa, 0,
                                         executor),
            "flat legacy surface preserved");
    Options alias = o;
    alias.discretization = "centroid";
    alias.formulation = "tabi2-nodepatch";
    alias.surface_only = true;
    auto alias_parents =
        build_surface({{0, 0, 0, 0, 1.7}}, alias, nullptr, &executor);
    require(
        std::all_of(alias_parents.begin(), alias_parents.end(),
                    [](const Panel &p) {
                      return p.sphere_patch_quadrature &&
                             p.sphere_patch_quadrature->nodes.empty();
                    }),
        "node alias retains parent geometry without triangle self preparation");
    alias.fmm_m2l = "single-matvec-analytic";
    auto surface_only = solve_linear_pb({{0, 0, 0, 1, 1.7}}, alias, executor);
    require(surface_only.converged && surface_only.surface_only,
            "surface-only does not invoke unsupported FMM kernels");
    alias.surface_only = false;
    rejected = false;
    try {
      solve_linear_pb({{0, 0, 0, 1, 1.7}}, alias, executor);
    } catch (const std::runtime_error &error) {
      rejected = std::string(error.what()) ==
                 "curved_node_patch_requires_direct_or_classic_treecode";
    }
    require(rejected, "unsupported curved FMM backend is explicit");
    std::cout << "PASS " << checks << " curved-cache checks\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
