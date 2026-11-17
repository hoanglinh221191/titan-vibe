// p1 is degree one in all four kernels (manuscript section 2.3): the far
// field of a source cell keeps its monopole and first Cartesian moments, so
// its error falls as (s/r)^2 with the size s of the cell, in K1, K2, K3 and
// K4 alike (a monopole alone falls as s/r).  One cluster of sources, off its
// box centre, with mixed normals and densities, shrunk about that centre,
// against the direct pair sum.  The error is the RMS over 64 target
// directions and normals, against the cancellation-free size sum_j |K x_j|.
// The moments are built both in one leaf and from shifted child moments
// (leaf 3), which must agree to rounding.
//
// Built against the engine objects without tabi_fmm_m2l_prototype.o (this
// file includes that source to reach its far-cell kernel).
#include "../src/tabi/tabi_fmm_m2l_prototype.cpp"
#include <iomanip>
#include <random>
using namespace titan::tabi::fmm_m2l;
int main() {
  std::mt19937 rng(20261018);
  std::uniform_real_distribution<double> u(-1.0, 1.0), w(0.2, 1.0);
  const auto unit = [&] {
    double x, y, z, n;
    do {
      x = u(rng);
      y = u(rng);
      z = u(rng);
      n = std::sqrt(x * x + y * y + z * z);
    } while (n < 0.1 || n > 1.0);
    return std::array<double, 3>{x / n, y / n, z / n};
  };
  std::vector<PrototypePanel> panels(24);
  std::vector<double> phi(panels.size()), phin(panels.size());
  for (std::size_t j = 0; j < panels.size(); ++j) {
    auto &p = panels[j];
    // Skewed so the box centre is not the centroid.
    p.x = 0.6 * u(rng) + 0.25 * std::abs(u(rng));
    p.y = 0.4 * u(rng);
    p.z = 0.5 * u(rng) * std::abs(u(rng));
    const auto n = unit();
    p.nx = n[0];
    p.ny = n[1];
    p.nz = n[2];
    p.area = w(rng);
    phi[j] = u(rng) + 0.3;
    phin[j] = u(rng) - 0.2;
  }
  std::vector<std::array<double, 3>> dirs(64), normals(64);
  for (std::size_t k = 0; k < dirs.size(); ++k) {
    dirs[k] = unit();
    normals[k] = unit();
  }
  // The cluster shrunk by h about its box centre q (which stays q), the
  // targets fixed at r from q: the far-field error of degree one falls as
  // h^2 at any kappa r, so the order does not mix in the kernel's shape.
  const PointSoA base = make_panel_source_points(panels);
  const Tree base_tree = build_tree(base, panels.size());
  const Cell &base_root = base_tree.cells.front();
  const std::array<double, 3> q{base_root.qx, base_root.qy, base_root.qz};
  const double s = base_root.expansion_radius, r = 6.0;
  int failures = 0;
  std::cout << std::setprecision(4) << "P1_DEGREE cluster of " << panels.size()
            << " sources, expansion radius " << s << " A shrunk by h, targets at "
            << r << " A\n";
  for (const double kappa : {0.0, 0.02, 0.125, 0.5}) {
    const double eps = 80.0;
    std::array<double, 4> previous{}, order{};
    double shift = 0.0;
    for (const double h : {1.0, 0.5, 0.25, 0.125}) {
      std::vector<PrototypePanel> shrunk = panels;
      for (auto &p : shrunk) {
        p.x = q[0] + h * (p.x - q[0]);
        p.y = q[1] + h * (p.y - q[1]);
        p.z = q[2] + h * (p.z - q[2]);
      }
      const PointSoA points = make_panel_source_points(shrunk);
      const std::array<Tree, 2> trees{build_tree(points, shrunk.size()),
                                      build_tree(points, 3)};
      std::array<int, 2> roots{-1, -1};
      std::array<std::vector<PanelSourceMoments>, 2> moments;
      for (int t = 0; t < 2; ++t) {
        for (std::size_t c = 0; c < trees[t].cells.size() && roots[t] < 0; ++c)
          if (trees[t].cells[c].begin == 0 &&
              trees[t].cells[c].end == shrunk.size())
            roots[t] = static_cast<int>(c);
        compute_panel_source_moments_into(trees[t], shrunk, phi, phin,
                                          moments[t]);
      }
      const Cell &root = trees[0].cells[static_cast<std::size_t>(roots[0])];
      const Cell &split = trees[1].cells[static_cast<std::size_t>(roots[1])];
      std::array<double, 4> err2{}, size2{};
      for (std::size_t k = 0; k < dirs.size(); ++k) {
        PrototypePanel target;
        target.x = q[0] + r * dirs[k][0];
        target.y = q[1] + r * dirs[k][1];
        target.z = q[2] + r * dirs[k][2];
        target.nx = normals[k][0];
        target.ny = normals[k][1];
        target.nz = normals[k][2];
        JcpChannelValues exact{}, size{};
        for (std::size_t j = 0; j < shrunk.size(); ++j) {
          const auto pair = jcp_direct_pair_channels(target, shrunk[j], phi[j],
                                                     phin[j], eps, kappa);
          for (int c = 0; c < 4; ++c) {
            exact[c] += pair[c];
            size[c] += std::abs(pair[c]);
          }
        }
        const auto approx = jcp_far_cell_channels(
            target, root, moments[0][static_cast<std::size_t>(roots[0])], eps,
            kappa);
        const auto shifted = jcp_far_cell_channels(
            target, split, moments[1][static_cast<std::size_t>(roots[1])], eps,
            kappa);
        for (int c = 0; c < 4; ++c) {
          err2[c] += (approx[c] - exact[c]) * (approx[c] - exact[c]);
          size2[c] += size[c] * size[c];
          if (!(kappa == 0.0 && (c == 0 || c == 3)))
            shift = std::max(shift, std::abs(shifted[c] - approx[c]) / size[c]);
        }
      }
      for (int c = 0; c < 4; ++c) {
        const double e = size2[c] > 0.0 ? std::sqrt(err2[c] / size2[c]) : 0.0;
        if (h < 1.0 && e > 0.0)
          order[c] = std::log2(previous[c] / e);
        previous[c] = e;
      }
    }
    std::cout << "  kappa=" << kappa << ": error at h = 1/8";
    for (int c = 0; c < 4; ++c)
      std::cout << " K" << c + 1 << " " << previous[c];
    std::cout << "; order h 1/4 -> 1/8";
    for (int c = 0; c < 4; ++c)
      std::cout << " K" << c + 1 << " " << order[c];
    std::cout << "; child-shift difference " << shift << "\n";
    // K1 and K4 vanish identically at kappa 0 (G0 - Gk = 0).
    for (int c = 0; c < 4; ++c) {
      if (kappa == 0.0 && (c == 0 || c == 3))
        continue;
      if (!(order[c] > 1.9 && order[c] < 2.1))
        ++failures;
    }
    if (!(shift < 1e-12))
      ++failures;
  }
  std::cout << (failures ? "FAILED " : "PASSED ") << failures << "\n";
  return failures ? 1 : 0;
}
