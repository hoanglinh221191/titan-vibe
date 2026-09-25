#include "tabi/analyzer_impl.h"
#include "titan_executor.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

// Native 3-D FD solve with exact Born Dirichlet data. The boundary is supplied
// analytically to isolate the volume operator from finite-box DH
// approximations.
void born_grid() {
  TitanExecutor executor(4);
  const double radius = 3.0, epsout = 80.0;
  double worst_energy = 0.0, worst_salt = 0.0;
  for (double epsin : {1.0, 80.0}) {
    std::array<double, 3> potentials{};
    double vacuum = 0.0;
    int case_index = 0;
    for (double salt : {0.0, 0.15, 0.5, -1.0}) {
      TitanPBE::AnalyzerImpl a;
      a.args_.salt = std::max(0.0, salt);
      a.args_.ion_radius = 0.0;
      a.args_.probe_radius = 0.0;
      a.args_.rms_tol = 1e-11;
      a.args_.dphi_tol = 1e-10;
      a.args_.max_iter = 5000;
      TitanPBE::AnalyzerImpl::GridShape g;
      g.h = 0.125;
      g.nx = g.ny = g.nz = 145;
      g.n = g.nx * g.ny * g.nz;
      g.ox = g.oy = g.oz = -9.0;
      a.ensure_workspace(g);
      TitanPBE::AnalyzerImpl::Atom atom;
      atom.q = 1.0;
      atom.radius = radius;
      const std::vector<TitanPBE::AnalyzerImpl::Atom> atoms{atom};
      a.classify_surface(g, atoms, executor);
      a.classify_ion_accessibility(g, atoms, executor);
      a.assign_charges(g, atoms);
      const double eo = salt < 0 ? epsin : epsout;
      const double k = 0.3291 * std::sqrt(std::max(0.0, salt));
      a.build_coefficients(g, atoms, epsin, eo, k, executor);
      for (std::size_t z = 0; z < g.nz; ++z)
        for (std::size_t y = 0; y < g.ny; ++y)
          for (std::size_t x = 0; x < g.nx; ++x) {
            const double xx = g.ox + x * g.h, yy = g.oy + y * g.h,
                         zz = g.oz + z * g.h;
            const double r = std::sqrt(xx * xx + yy * yy + zz * zz);
            const double reaction =
                (1 / (eo * (1 + k * radius)) - 1 / epsin) / radius;
            const double phi =
                r < radius
                    ? (r > 0 ? 1 / (epsin * r) + reaction : 0.0)
                    : std::exp(-k * (r - radius)) / (eo * (1 + k * radius) * r);
            a.grid_.phi[(z * g.ny + y) * g.nx + x] = phi;
          }
      const auto stats = a.solve_sor(g, executor, 0);
      if (!stats.converged)
        throw std::runtime_error("Born grid solve did not converge");
      const double center =
          a.grid_.phi[((g.nz / 2) * g.ny + g.ny / 2) * g.nx + g.nx / 2];
      if (salt < 0)
        vacuum = center;
      else
        potentials[case_index++] = center;
    }
    for (int i = 0; i < 3; ++i) {
      const double k =
          0.3291 * std::sqrt(std::array<double, 3>{0, 0.15, 0.5}[i]);
      const double ref = (1 / (epsout * (1 + k * radius)) - 1 / epsin) / radius;
      if (std::abs(ref) > 1e-12) {
        const double err = std::abs((potentials[i] - vacuum) / ref - 1);
        std::cout << "Born eps_in=" << epsin << " salt_index=" << i
                  << " phi_rxn=" << potentials[i] - vacuum << " exact=" << ref
                  << " relative_error=" << err << '\n';
        worst_energy = std::max(worst_energy, err);
      }
      if (i > 0) {
        const double salt_ref = -k / (epsout * (1 + k * radius));
        const double err =
            std::abs((potentials[i] - potentials[0]) / salt_ref - 1);
        std::cout << "Born salt_shift_relative_error=" << err << '\n';
        worst_salt = std::max(worst_salt, err);
      }
    }
  }
  if (worst_energy > 0.01 || worst_salt > 0.01)
    throw std::runtime_error("Born grid energy/salt shift exceeds 1 percent");
  std::cout << "PASS GRID-01 3D Born salt=0/0.15/0.5 eps_in=1/80 exact "
               "boundary; worst_energy="
            << worst_energy << " worst_salt_shift=" << worst_salt << '\n';
}

// Compiled with -fno-access-control: exercise the production coefficient
// without exposing private analyzer state in the public API or copying its
// formula.
int run() {
  TitanPBE::AnalyzerImpl analyzer;
  TitanPBE::AnalyzerImpl::GridShape shape;
  analyzer.grid_.ion_accessible = {1, 0};
  analyzer.args_.temp = 298.15;
  double worst = 0.0;
  for (double eps : {2.0, 40.0, 80.0}) {
    for (double salt : {0.0, 0.15, 0.5}) {
      analyzer.args_.salt = salt;
      const double kappa = 0.3291 * std::sqrt(salt);
      double previous = 1.0;
      for (double h : {0.4, 0.2, 0.1}) {
        shape.h = h;
        const double screen =
            analyzer.effective_screen_term(shape, eps, kappa, 0);
        if (analyzer.effective_screen_term(shape, eps, kappa, 1) != 0.0 ||
            analyzer.effective_screen_term(shape, eps, kappa, 2) != 0.0)
          throw std::runtime_error(
              "screen leaks into excluded/out-of-range cells");
        if (salt == 0.0) {
          if (screen != 0.0)
            throw std::runtime_error("nonzero salt-free screen");
          continue;
        }
        // A homogeneous solvent has the exact charge-free solution exp(-k*x).
        // Apply the seven-point production stencil at x=0 (constant in y,z).
        const double diffusion = eps * (2.0 * std::cosh(kappa * h) - 2.0);
        const double defect = std::abs(diffusion - screen) / diffusion;
        if (defect > 0.001 || defect > 0.26 * previous)
          throw std::runtime_error(
              "Yukawa stencil does not converge at physical Debye length");
        previous = defect;
        worst = std::max(worst, defect);
      }
    }
  }
  // Delphi uses its own Debye constant, but the same dielectric convention.
  analyzer.args_.salt = 0.15;
  analyzer.args_.test_operator_convention = "delphi";
  analyzer.args_.test_pqr = "oracle.pqr";
  const double debye = 0.01990076478 * std::sqrt(298.15 * 80.0 / 0.15);
  const double expected = 80.0 * shape.h * shape.h / (debye * debye);
  if (std::abs(analyzer.effective_screen_term(shape, 80.0, 0.0, 0) - expected) >
      1e-14)
    throw std::runtime_error("Delphi convention changed");
  std::cout << "PASS GRID-01 physical Debye stencil, masks, zero salt, Delphi; "
               "worst_defect="
            << worst << '\n';
  return 0;
}

int main() {
  try {
    run();
    born_grid();
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
  }
}
