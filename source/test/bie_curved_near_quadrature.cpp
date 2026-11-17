// Near quadrature of curved node-patch cells across a narrow gap.
//
// Two vdW spheres of R = 1.7 A (sphere-patch icosphere, 96 points, exact
// spherical parents), the second placed so that the gap along one node's
// direction is 0.2, 0.05, 0.01 or 0.002 A. The K2 coefficient between the two
// facing nodes is compared with an independent order-256 reference (from the
// September 2026 review, kept with bie-exp). The stored order-4 samples miss
// it by 0.017%, 1.5%, 18% and 71%; the adaptive near quadrature
// (geometric_pair in tabi_node_patch_curved.cpp) must be within 3e-6.
//
// Build: link with the engine objects (test/pbe_p0_oracles.py --case near).
#include "tabi/tabi_node_patch.h"
#include "tabi/tabi_node_patch_curved.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace TitanPBE::Tabi;

int main() {
  try {
    TitanExecutor executor(4);
    Options o;
    o.primary_surface = "sphere-patch";
    o.surface_mode = "vdw";
    o.sphere_mesh_mode = "icosphere";
    o.surface_point_count = 96;
    o.discretization = "node-patch";
    o.surface_only = true;
    constexpr double kEps = 78.5, kKappa = 0.0, kRadius = 1.7;

    const auto single = build_surface({{0, 0, 0, 1, kRadius}}, o, nullptr,
                                      &executor);
    auto single_nodes = build_node_patch_panels(single);
    prepare_curved_node_patches(single, single_nodes, kEps, kKappa, 0,
                                executor, false);
    const Panel source = single_nodes.front();

    const double gaps[] = {0.2, 0.05, 0.01, 0.002};
    const double references[] = {12.563530879298456, 28.245176266499676,
                                 34.4301082962227, 35.72561815795781};
    double worst = 0;
    for (int k = 0; k < 4; ++k) {
      const double factor = (2 * kRadius + gaps[k]) / kRadius;
      const auto triangles = build_surface(
          {{0, 0, 0, 1, kRadius},
           {factor * source.x, factor * source.y, factor * source.z, 1,
            kRadius}},
          o, nullptr, &executor);
      auto nodes = build_node_patch_panels(triangles);
      if (!prepare_curved_node_patches(triangles, nodes, kEps, kKappa, 0,
                                       executor, false))
        throw std::runtime_error("curved cells not built");
      const auto nearest = [&](double f) {
        std::size_t best = 0;
        double best_d = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < nodes.size(); ++i) {
          const double x = nodes[i].x - f * source.x,
                       y = nodes[i].y - f * source.y,
                       z = nodes[i].z - f * source.z;
          const double d = x * x + y * y + z * z;
          if (d < best_d) {
            best_d = d;
            best = i;
          }
        }
        return best;
      };
      // Target on the second sphere facing the gap, source on the first.
      const std::size_t i = nearest(1 + gaps[k] / kRadius), j = nearest(1);
      const auto block = curved_node_pair(nodes[i], i, nodes[j], j, kEps, kKappa);
      const double error = std::abs(block[1] / references[k] - 1);
      worst = std::max(worst, error);
      std::printf("gap %.3f A: K2 %.15g reference %.15g relative error %.3e\n",
                  gaps[k], block[1], references[k], error);
    }
    if (!(worst < 3e-6)) {
      std::printf("FAIL bie_curved_near_quadrature worst %.3e\n", worst);
      return 1;
    }
    std::printf("PASS bie_curved_near_quadrature worst %.3e\n", worst);
    return 0;
  } catch (const std::exception &e) {
    std::printf("FAIL %s\n", e.what());
    return 1;
  }
}
