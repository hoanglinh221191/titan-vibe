#include "tabi/tabi_near.h"
#include "tabi/tabi_node_patch.h"
#include "tabi/tabi_node_patch_curved.h"
#include "titan_executor.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <set>
using namespace TitanPBE::Tabi;
int main(int argc, char **argv) {
  try {
    if (argc != 5)
      return 2;
    int points = std::stoi(argv[1]), direction = std::stoi(argv[4]);
    double offset = std::stod(argv[2]), salt = std::stod(argv[3]), R = 1.7,
           kappa = .3291 * std::sqrt(salt);
    const bool coupled = kappa > 1e-14;
    const double pi = std::numbers::pi;
    std::array<double, 3> axis =
        direction ? std::array<double, 3>{1 / std::sqrt(3.), 1 / std::sqrt(3.),
                                          1 / std::sqrt(3.)}
                  : std::array<double, 3>{1, 0, 0};
    Options o;
    o.primary_surface = "sphere-patch";
    o.surface_mode = "vdw";
    o.surface_only = true;
    o.discretization = "node-patch";
    o.sphere_mesh_mode = "icosphere";
    o.surface_point_count = points;
    TitanExecutor executor(4);
    auto triangles = build_surface({{0, 0, 0, 0, R}}, o, nullptr, &executor);
    auto nodes = build_node_patch_panels(triangles);
    if (!prepare_curved_node_patches(triangles, nodes, 78.5, kappa, 0,
                                     executor))
      return 3;
    std::vector<std::array<double, 2>> trace(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const auto &p = nodes[i];
      double cosine = (p.x * axis[0] + p.y * axis[1] + p.z * axis[2]) / R;
      long double power = 1, P = 1, prev = 1, z = kappa * R,
                  ratio = z > 0 ? 1 + 1 / z : 0;
      for (int l = 0; l < 200; ++l) {
        if (l == 1)
          P = cosine;
        else if (l > 1) {
          auto next = ((2 * l - 1) * cosine * P - (l - 1) * prev) / l;
          prev = P;
          P = next;
        }
        auto gamma = z > 0 ? l - z * ratio : -(l + 1);
        auto a = (2 * l + 1) * power / (4 * pi * R * (l - 78.5 * gamma));
        trace[i][0] += a * P;
        trace[i][1] += a * 78.5 * gamma / R * P;
        power *= offset;
        if (z > 0)
          ratio = (2 * l + 3) / z + 1 / ratio;
        if (l > 2 && std::abs(power) < 1e-20)
          break;
      }
    }
    std::set<std::size_t> targets;
    for (std::size_t i = 0; i < 16; ++i)
      targets.insert(i * nodes.size() / 16);
    std::size_t closest = 0;
    double distance = 1e300;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const auto &p = nodes[i];
      double d = std::pow(p.x - offset * R * axis[0], 2) +
                 std::pow(p.y - offset * R * axis[1], 2) +
                 std::pow(p.z - offset * R * axis[2], 2);
      if (d < distance) {
        distance = d;
        closest = i;
      }
    }
    targets.insert(closest);
    std::cout << std::setprecision(17);
    for (auto i : targets) {
      std::array<std::array<double, 4>, 3> sum{};
      for (std::size_t j = 0; j < nodes.size(); ++j) {
        auto b =
            evaluate_jcp_near_pair_kernel(nodes[i], i, nodes[j], j, 78.5, kappa,
                                          coupled, "analytic_singular", "off");
        if (!b.valid)
          return 4;
        auto &g = sum[i == j ? 0 : b.one_ring ? 1 : 2];
        g[0] += b.k1 * trace[j][1];
        g[1] += b.k2 * trace[j][0];
        g[2] += b.k3 * trace[j][1];
        g[3] += b.k4 * trace[j][0];
      }
      for (int g = 0; g < 3; ++g)
        for (int b = 0; b < 4; ++b)
          std::cout << i << ',' << g << ',' << b << ',' << sum[g][b] << '\n';
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
