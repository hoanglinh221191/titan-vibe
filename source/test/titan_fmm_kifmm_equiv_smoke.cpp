#include "tabi/tabi_fmm_kifmm_solve.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmScalarSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::direct_scalar_potential;
using titan::tabi::fmm_kifmm::equivalent_scalar_potential;
using titan::tabi::fmm_kifmm::make_cube_surface_rule;
using titan::tabi::fmm_kifmm::solve_upward_equivalent;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

std::vector<KifmmScalarSource> deterministic_sources(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {{box.center.x - 0.16 * h, box.center.y + 0.09 * h,
        box.center.z - 0.06 * h},
       0.70},
      {{box.center.x + 0.12 * h, box.center.y - 0.15 * h,
        box.center.z + 0.10 * h},
       -0.55},
      {{box.center.x + 0.04 * h, box.center.y + 0.14 * h,
        box.center.z - 0.17 * h},
       0.41},
      {{box.center.x - 0.09 * h, box.center.y - 0.04 * h,
        box.center.z + 0.13 * h},
       -0.23},
      {{box.center.x + 0.01 * h, box.center.y + 0.01 * h,
        box.center.z + 0.02 * h},
       0.37},
  };
}

std::vector<Vec3> external_targets(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {box.center.x + 7.60 * h, box.center.y + 0.31 * h,
       box.center.z - 0.26 * h},
      {box.center.x - 7.75 * h, box.center.y - 0.18 * h,
       box.center.z + 0.44 * h},
      {box.center.x + 0.27 * h, box.center.y + 7.90 * h,
       box.center.z + 0.19 * h},
      {box.center.x - 0.36 * h, box.center.y - 7.65 * h,
       box.center.z - 0.21 * h},
      {box.center.x + 0.16 * h, box.center.y + 0.42 * h,
       box.center.z + 8.10 * h},
      {box.center.x - 0.28 * h, box.center.y - 0.24 * h,
       box.center.z - 7.95 * h},
      {box.center.x + 6.80 * h, box.center.y + 5.90 * h,
       box.center.z - 0.37 * h},
      {box.center.x - 6.50 * h, box.center.y + 0.48 * h,
       box.center.z + 6.70 * h},
  };
}

double relative_rms(const std::vector<double> &approx,
                    const std::vector<double> &direct) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const double diff = approx[i] - direct[i];
    err_sq += diff * diff;
    ref_sq += direct[i] * direct[i];
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

void run_parity_case(int order, KifmmKernel kernel, double kappa,
                     const std::string &solver, double required_rel_rms) {
  const KifmmBox box{{0.15, -0.25, 0.40}, 1.0, 0, 7};
  const auto rule = make_cube_surface_rule(order, 1.05, 2.40);
  const auto sources = deterministic_sources(box);
  const auto targets = external_targets(box);
  const auto result =
      solve_upward_equivalent(rule, box, sources, kernel, kappa, solver);
  require(result.finite, "equivalent solve finite for order " +
                             std::to_string(order) + " solver " + solver);
  require(result.solver == solver, "equivalent solve records solver");
  require(result.charges.size() == rule.unit_equiv.size(),
          "equivalent charge count matches surface");
  require(std::isfinite(result.relative_residual),
          "equivalent residual finite");
  require(std::isfinite(result.max_abs_residual),
          "equivalent max residual finite");

  std::vector<double> approx;
  std::vector<double> direct;
  for (const Vec3 target : targets) {
    approx.push_back(equivalent_scalar_potential(rule, box, result, kernel,
                                                 kappa, target));
    direct.push_back(direct_scalar_potential(sources, kernel, kappa, target));
  }
  const double rms = relative_rms(approx, direct);
  require(std::isfinite(rms), "external relative RMS finite");
  std::cout << std::setprecision(10) << "KIFMM_EQUIV order=" << order
            << " kernel="
            << (kernel == KifmmKernel::Laplace ? "laplace" : "yukawa")
            << " kappa=" << kappa << " solver=" << solver
            << " equiv_residual=" << result.relative_residual
            << " max_abs_residual=" << result.max_abs_residual
            << " external_rel_rms=" << rms << '\n';
  if (required_rel_rms > 0.0) {
    require(rms <= required_rel_rms,
            "external relative RMS gate order " + std::to_string(order) +
                " solver " + solver + " rms=" + std::to_string(rms));
  }
}

} // namespace

int main() {
  run_parity_case(4, KifmmKernel::Laplace, 0.0, "qr", 0.0);
  run_parity_case(6, KifmmKernel::Laplace, 0.0, "qr", 1.0e-6);
  run_parity_case(8, KifmmKernel::Laplace, 0.0, "qr", 0.0);

  for (const double kappa : {0.05, 0.12746, 0.25}) {
    run_parity_case(6, KifmmKernel::Yukawa, kappa, "qr", 1.0e-6);
  }
  run_parity_case(6, KifmmKernel::Yukawa, 0.12746, "svd", 1.0e-6);
  run_parity_case(6, KifmmKernel::Laplace, 0.0, "normal-tikhonov", 0.0);

  std::cout << "PASS titan_fmm_kifmm_equiv_smoke\n";
  return EXIT_SUCCESS;
}
