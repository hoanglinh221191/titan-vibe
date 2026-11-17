// GRID-03 oracle for src/pbe/grid/galerkin_mg.*:
//  1. the stored coarse operators equal P^T A P computed densely, on grids of
//     odd and even size;
//  2. the preconditioned CG solution equals a dense direct solve of the same
//     Dirichlet problem;
//  3. convergence: stationary V-cycle and CG on a uniform Laplacian and on a
//     dielectric sphere (eps 2 inside, 80 outside), 33^3 and 65^3;
//  4. the result does not depend on the thread count (bitwise).
#include "pbe/grid/galerkin_mg.h"
#include "pbe/grid/grid_common.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace TitanPBE;
namespace mg = TitanPBE::galerkin_mg;

namespace {

int failures = 0;

void check(bool ok, const char *what, double value) {
  std::printf("%s %s (%.3e)\n", ok ? "PASS" : "FAIL", what, value);
  if (!ok)
    ++failures;
}

struct Grid {
  std::size_t nx, ny, nz, n;
  std::vector<double> ex, ey, ez, diag;
  mg::FineOperator op() const {
    return {nx, ny, nz, ex.data(), ey.data(), ez.data(), diag.data()};
  }
};

// h^2-scaled 7-point operator from node permittivities (harmonic face means)
// plus an optional screening term.
Grid make_grid(std::size_t nx, std::size_t ny, std::size_t nz,
               const std::vector<double> &eps, const std::vector<double> &screen) {
  Grid g{nx, ny, nz, nx * ny * nz, {}, {}, {}, {}};
  g.ex.assign(g.n, 0.0);
  g.ey.assign(g.n, 0.0);
  g.ez.assign(g.n, 0.0);
  g.diag.assign(g.n, 0.0);
  auto face = [](double a, double b) { return 2.0 * a * b / (a + b); };
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t id = idx3(i, j, k, nx, ny);
        if (i + 1 < nx) g.ex[id] = face(eps[id], eps[id + 1]);
        if (j + 1 < ny) g.ey[id] = face(eps[id], eps[id + nx]);
        if (k + 1 < nz) g.ez[id] = face(eps[id], eps[id + nx * ny]);
      }
  for (std::size_t k = 1; k + 1 < nz; ++k)
    for (std::size_t j = 1; j + 1 < ny; ++j)
      for (std::size_t i = 1; i + 1 < nx; ++i) {
        const std::size_t id = idx3(i, j, k, nx, ny);
        g.diag[id] = g.ex[id] + g.ex[id - 1] + g.ey[id] + g.ey[id - nx] +
                     g.ez[id] + g.ez[id - nx * ny] + screen[id];
      }
  return g;
}

Grid random_grid(std::size_t nx, std::size_t ny, std::size_t nz, std::mt19937 &rng) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<double> eps(nx * ny * nz), screen(nx * ny * nz);
  for (auto &e : eps) e = u(rng) < 0.5 ? 2.0 + u(rng) : 80.0 * (0.5 + u(rng));
  for (auto &s : screen) s = u(rng) < 0.5 ? 0.0 : u(rng);
  return make_grid(nx, ny, nz, eps, screen);
}

Grid sphere_grid(std::size_t n, double eps_in, double eps_out) {
  std::vector<double> eps(n * n * n), screen(n * n * n, 0.0);
  const double c = 0.5 * static_cast<double>(n - 1), r = 0.3 * static_cast<double>(n - 1);
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t i = 0; i < n; ++i) {
        const double d2 = (i - c) * (i - c) + (j - c) * (j - c) + (k - c) * (k - c);
        eps[idx3(i, j, k, n, n)] = d2 < r * r ? eps_in : eps_out;
      }
  return make_grid(n, n, n, eps, screen);
}

// dense interior matrix of a level: the fine 7-point operator or a stored
// 27-point level
Eigen::MatrixXd dense_fine(const Grid &g) {
  const std::size_t mx = g.nx - 2, my = g.ny - 2, mz = g.nz - 2;
  const Eigen::Index m = static_cast<Eigen::Index>(mx * my * mz);
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(m, m);
  auto row = [&](std::size_t i, std::size_t j, std::size_t k) {
    return static_cast<Eigen::Index>(((k - 1) * my + (j - 1)) * mx + (i - 1));
  };
  auto inside = [&](std::size_t i, std::size_t j, std::size_t k) {
    return i >= 1 && i + 1 < g.nx && j >= 1 && j + 1 < g.ny && k >= 1 && k + 1 < g.nz;
  };
  for (std::size_t k = 1; k + 1 < g.nz; ++k)
    for (std::size_t j = 1; j + 1 < g.ny; ++j)
      for (std::size_t i = 1; i + 1 < g.nx; ++i) {
        const std::size_t id = idx3(i, j, k, g.nx, g.ny);
        const Eigen::Index r = row(i, j, k);
        a(r, r) = g.diag[id];
        if (inside(i + 1, j, k)) a(r, row(i + 1, j, k)) = -g.ex[id];
        if (inside(i - 1, j, k)) a(r, row(i - 1, j, k)) = -g.ex[id - 1];
        if (inside(i, j + 1, k)) a(r, row(i, j + 1, k)) = -g.ey[id];
        if (inside(i, j - 1, k)) a(r, row(i, j - 1, k)) = -g.ey[id - g.nx];
        if (inside(i, j, k + 1)) a(r, row(i, j, k + 1)) = -g.ez[id];
        if (inside(i, j, k - 1)) a(r, row(i, j, k - 1)) = -g.ez[id - g.nx * g.ny];
      }
  return a;
}

Eigen::MatrixXd dense_level(const mg::Hierarchy::Level &level) {
  const std::size_t mx = level.nx - 2, my = level.ny - 2, mz = level.nz - 2;
  const Eigen::Index m = static_cast<Eigen::Index>(mx * my * mz);
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(m, m);
  auto row = [&](std::size_t i, std::size_t j, std::size_t k) {
    return static_cast<Eigen::Index>(((k - 1) * my + (j - 1)) * mx + (i - 1));
  };
  for (std::size_t k = 1; k + 1 < level.nz; ++k)
    for (std::size_t j = 1; j + 1 < level.ny; ++j)
      for (std::size_t i = 1; i + 1 < level.nx; ++i) {
        const std::size_t id = idx3(i, j, k, level.nx, level.ny);
        const Eigen::Index r = row(i, j, k);
        a(r, r) = level.coef[id * 14];
        for (int s = 1; s <= 13; ++s) {
          int dx, dy, dz;
          mg::Hierarchy::forward_offset(s, dx, dy, dz);
          const std::size_t gi = i + dx, gj = j + dy, gk = k + dz;
          if (gi < 1 || gi + 1 >= level.nx || gj < 1 || gj + 1 >= level.ny ||
              gk < 1 || gk + 1 >= level.nz)
            continue;
          const double v = level.coef[id * 14 + s];
          a(r, row(gi, gj, gk)) = v;
          a(row(gi, gj, gk), r) = v;
        }
      }
  return a;
}

// trilinear prolongation from coarse interior to fine interior
Eigen::MatrixXd dense_p(std::size_t fnx, std::size_t fny, std::size_t fnz,
                        std::size_t cnx, std::size_t cny, std::size_t cnz) {
  const std::size_t fx = fnx - 2, fy = fny - 2, fz = fnz - 2;
  const std::size_t cx = cnx - 2, cy = cny - 2, cz = cnz - 2;
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(fx * fy * fz),
                                            static_cast<Eigen::Index>(cx * cy * cz));
  auto weights = [](std::size_t g, std::size_t cn, std::vector<std::pair<std::size_t, double>> &out) {
    out.clear();
    if (g % 2 == 0) out.push_back({g / 2, 1.0});
    else { out.push_back({(g - 1) / 2, 0.5}); out.push_back({(g + 1) / 2, 0.5}); }
    out.erase(std::remove_if(out.begin(), out.end(),
                             [&](const auto &e) { return e.first < 1 || e.first + 1 >= cn; }),
              out.end());
  };
  std::vector<std::pair<std::size_t, double>> wi, wj, wk;
  for (std::size_t k = 1; k + 1 < fnz; ++k)
    for (std::size_t j = 1; j + 1 < fny; ++j)
      for (std::size_t i = 1; i + 1 < fnx; ++i) {
        weights(i, cnx, wi); weights(j, cny, wj); weights(k, cnz, wk);
        const Eigen::Index r = static_cast<Eigen::Index>(((k - 1) * fy + (j - 1)) * fx + (i - 1));
        for (auto [K, a] : wk)
          for (auto [J, b] : wj)
            for (auto [I, c] : wi)
              p(r, static_cast<Eigen::Index>(((K - 1) * cy + (J - 1)) * cx + (I - 1))) += a * b * c;
      }
  return p;
}

void galerkin_case(std::size_t nx, std::size_t ny, std::size_t nz, std::mt19937 &rng,
                   TitanExecutor &executor) {
  const Grid g = random_grid(nx, ny, nz, rng);
  mg::Options options;
  options.direct_unknowns = 8; // force coarsening of the small test grid
  mg::Hierarchy h;
  h.build(g.op(), options, executor);
  Eigen::MatrixXd fine = dense_fine(g);
  double worst = 0.0;
  std::size_t fnx = nx, fny = ny, fnz = nz;
  for (int l = 1; l < h.level_count(); ++l) {
    const auto &c = h.level(l);
    const Eigen::MatrixXd p = dense_p(fnx, fny, fnz, c.nx, c.ny, c.nz);
    const Eigen::MatrixXd expected = p.transpose() * fine * p;
    const Eigen::MatrixXd stored = dense_level(c);
    worst = std::max(worst, (expected - stored).cwiseAbs().maxCoeff() /
                                expected.cwiseAbs().maxCoeff());
    fine = expected;
    fnx = c.nx; fny = c.ny; fnz = c.nz;
  }
  char label[96];
  std::snprintf(label, sizeof label, "Galerkin A_H = P^T A P on %zux%zux%zu (%d levels)",
                nx, ny, nz, h.level_count());
  check(h.level_count() >= 2 && worst < 1e-13, label, worst);
}

void direct_case(std::size_t nx, std::size_t ny, std::size_t nz, std::mt19937 &rng,
                 TitanExecutor &executor) {
  const Grid g = random_grid(nx, ny, nz, rng);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::vector<double> b(g.n, 0.0), x(g.n, 0.0);
  for (auto &v : b) v = u(rng);
  // Dirichlet data on the outer layer of x
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i)
        if (i == 0 || j == 0 || k == 0 || i + 1 == nx || j + 1 == ny || k + 1 == nz)
          x[idx3(i, j, k, nx, ny)] = u(rng);
  // dense reference: interior system with the boundary moved to the rhs
  const std::size_t mx = nx - 2, my = ny - 2;
  Eigen::VectorXd rhs(static_cast<Eigen::Index>(mx * my * (nz - 2)));
  for (std::size_t k = 1; k + 1 < nz; ++k)
    for (std::size_t j = 1; j + 1 < ny; ++j)
      for (std::size_t i = 1; i + 1 < nx; ++i) {
        const std::size_t id = idx3(i, j, k, nx, ny);
        double v = b[id];
        if (i == 1) v += g.ex[id - 1] * x[id - 1];
        if (i + 2 == nx) v += g.ex[id] * x[id + 1];
        if (j == 1) v += g.ey[id - nx] * x[id - nx];
        if (j + 2 == ny) v += g.ey[id] * x[id + nx];
        if (k == 1) v += g.ez[id - nx * ny] * x[id - nx * ny];
        if (k + 2 == nz) v += g.ez[id] * x[id + nx * ny];
        rhs[static_cast<Eigen::Index>(((k - 1) * my + (j - 1)) * mx + (i - 1))] = v;
      }
  const Eigen::VectorXd reference = dense_fine(g).llt().solve(rhs);
  mg::Options options;
  options.direct_unknowns = 8;
  mg::Hierarchy h;
  h.build(g.op(), options, executor);
  const auto result = mg::solve(h, b.data(), x.data(), 200, true,
                                [](const mg::IterationReport &r) { return r.residual_rms < 1e-13; },
                                executor);
  double worst = 0.0;
  for (std::size_t k = 1; k + 1 < nz; ++k)
    for (std::size_t j = 1; j + 1 < ny; ++j)
      for (std::size_t i = 1; i + 1 < nx; ++i)
        worst = std::max(worst, std::abs(x[idx3(i, j, k, nx, ny)] -
                                          reference[static_cast<Eigen::Index>(
                                              ((k - 1) * my + (j - 1)) * mx + (i - 1))]));
  char label[96];
  std::snprintf(label, sizeof label, "PCG equals dense solve on %zux%zux%zu, %d iterations",
                nx, ny, nz, result.iterations);
  check(result.converged && worst < 1e-10, label, worst);
}

int iterations_to(const Grid &g, bool use_cg, double target, TitanExecutor &executor,
                  std::vector<double> *solution = nullptr, std::size_t workers = 1) {
  std::vector<double> b(g.n, 0.0), x(g.n, 0.0);
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  for (auto &v : b) v = u(rng);
  mg::Options options;
  options.workers = workers;
  mg::Hierarchy h;
  h.build(g.op(), options, executor);
  const double r0 = mg::residual_rms(h, b.data(), x.data(), executor);
  const auto result = mg::solve(
      h, b.data(), x.data(), 100, use_cg,
      [&](const mg::IterationReport &r) { return r.residual_rms <= target * r0; },
      executor);
  if (solution != nullptr) *solution = x;
  return result.converged ? result.iterations : 1000;
}

} // namespace

int main() {
  std::mt19937 rng(20261023);
  TitanExecutor serial(1);
  TitanExecutor threads(6);

  galerkin_case(9, 9, 9, rng, serial);
  galerkin_case(10, 7, 12, rng, serial);
  galerkin_case(8, 11, 9, rng, serial);
  direct_case(9, 10, 11, rng, serial);
  direct_case(12, 9, 8, rng, serial);

  for (std::size_t n : {33u, 65u}) {
    const Grid uniform = sphere_grid(n, 80.0, 80.0);
    const Grid sphere = sphere_grid(n, 2.0, 80.0);
    const int vu = iterations_to(uniform, false, 1e-8, serial);
    const int cu = iterations_to(uniform, true, 1e-8, serial);
    const int vs = iterations_to(sphere, false, 1e-8, serial);
    const int cs = iterations_to(sphere, true, 1e-8, serial);
    std::printf("n=%zu  1e-8 reduction: uniform V %d, CG %d; sphere 2/80 V %d, CG %d\n",
                n, vu, cu, vs, cs);
    check(vu <= 10, "uniform Laplacian: stationary V-cycle within 10 cycles", vu);
    check(cu <= 8, "uniform Laplacian: CG within 8 iterations", cu);
    check(vs <= 25, "dielectric sphere: stationary V-cycle within 25 cycles", vs);
    check(cs <= 12, "dielectric sphere: CG within 12 iterations", cs);
  }

  const Grid sphere = sphere_grid(65, 2.0, 80.0);
  std::vector<double> one, six;
  iterations_to(sphere, true, 1e-10, serial, &one, 1);
  iterations_to(sphere, true, 1e-10, threads, &six, 6);
  check(one == six, "CG solution identical with 1 and 6 threads", 0.0);

  std::printf("%s\n", failures == 0 ? "ALL PASS" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
