#include "tabi/tabi_radial_taylor.h"
#include <cmath>
#include <iostream>
using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::detail;
int main() {
  int checks = 0;
  double maximum = 0;
  for (int order = 2; order <= 5; ++order)
    for (double kappa : {0., 1e-8, .13, 2.})
      for (double eps : {1e-15, 1., 20., 78.5})
        for (int seed = 0; seed < 8; ++seed) {
          Panel t;
          t.x = 3 + .17 * seed;
          t.y = -1.4;
          t.z = 2;
          t.nx = .3;
          t.ny = .4;
          t.nz = std::sqrt(.75);
          TreeNode c;
          c.cx = .3;
          c.cy = .7;
          c.cz = -.2;
          NodeMoments scalar;
          std::array<NodeMoments, 3> normal{};
          const auto &indices = moment_indices(order);
          for (std::size_t m = 0; m < indices.size(); ++m) {
            const auto &a = indices[m];
            double scale = std::pow(.2, a.ax + a.ay + a.az);
            scalar.values[m] = scale * std::cos(m + seed);
            for (int j = 0; j < 3; ++j)
              normal[j].values[m] = scale * std::sin(m * .7 + j + seed);
          }
          for (bool coupled : {false, true}) {
            auto ref = sum_jcp_channels(taylor_jcp_node_channel_sum_series(
                t, c, scalar, normal, eps, kappa, coupled, order, "analytic"));
            auto fast = radial_jcp_node_sum(t, c, scalar, normal, eps, kappa,
                                            coupled, order);
            double error =
                std::max(std::abs(fast.phi - ref.phi),
                         std::abs(fast.phin - ref.phin)) /
                std::max({1., std::abs(ref.phi), std::abs(ref.phin)});
            maximum = std::max(maximum, error);
            ++checks;
            if (!std::isfinite(error) || error > 1e-10) {
              std::cerr << "FAIL " << order << ' ' << kappa << ' ' << eps << ' '
                        << error << ' ' << fast.phi << ' ' << ref.phi << '\n';
              return 1;
            }
          }
        }
  std::cout << "PASS " << checks
            << " radial tensor comparisons; max=" << maximum << '\n';
}
