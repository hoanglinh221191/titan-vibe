#include "tabi/tabi_fmm_kifmm_fields.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmFieldFamily;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::direct_field_grad_target;
using titan::tabi::fmm_kifmm::direct_field_potential;
using titan::tabi::fmm_kifmm::equivalent_field_grad_target;
using titan::tabi::fmm_kifmm::equivalent_field_potential;
using titan::tabi::fmm_kifmm::make_cube_surface_rule;
using titan::tabi::fmm_kifmm::solve_upward_equivalent_field;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

double dot(Vec3 a, Vec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(Vec3 a) {
  return std::sqrt(dot(a, a));
}

double distance(Vec3 a, Vec3 b) {
  return norm({a.x - b.x, a.y - b.y, a.z - b.z});
}

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

std::vector<KifmmPanelSource> deterministic_sources(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {{box.center.x - 0.16 * h, box.center.y + 0.09 * h,
        box.center.z - 0.06 * h},
       unit({0.70, -0.20, 0.42}),
       0.86,
       0.35,
       -0.52},
      {{box.center.x + 0.12 * h, box.center.y - 0.15 * h,
        box.center.z + 0.10 * h},
       unit({-0.30, 0.84, -0.18}),
       1.04,
       -0.41,
       0.67},
      {{box.center.x + 0.04 * h, box.center.y + 0.14 * h,
        box.center.z - 0.17 * h},
       unit({0.25, 0.36, 0.91}),
       0.73,
       0.58,
       -0.21},
      {{box.center.x - 0.09 * h, box.center.y - 0.04 * h,
        box.center.z + 0.13 * h},
       unit({-0.61, -0.31, 0.55}),
       0.92,
       -0.27,
       0.43},
      {{box.center.x + 0.01 * h, box.center.y + 0.01 * h,
        box.center.z + 0.02 * h},
       unit({0.10, -0.74, -0.49}),
       0.64,
       0.22,
       -0.33},
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

double rel_rms(const std::vector<double> &approx,
               const std::vector<double> &direct) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const double d = approx[i] - direct[i];
    err_sq += d * d;
    ref_sq += direct[i] * direct[i];
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

double rel_rms_vec(const std::vector<Vec3> &approx,
                   const std::vector<Vec3> &direct) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const Vec3 d{approx[i].x - direct[i].x, approx[i].y - direct[i].y,
                 approx[i].z - direct[i].z};
    err_sq += dot(d, d);
    ref_sq += dot(direct[i], direct[i]);
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

const char *family_name(KifmmFieldFamily family) {
  switch (family) {
  case KifmmFieldFamily::PhinCharge:
    return "phin_charge";
  case KifmmFieldFamily::PhiNormalX:
    return "phi_normal_x";
  case KifmmFieldFamily::PhiNormalY:
    return "phi_normal_y";
  case KifmmFieldFamily::PhiNormalZ:
    return "phi_normal_z";
  }
  return "unknown";
}

const char *kernel_name(KifmmKernel kernel) {
  return kernel == KifmmKernel::Laplace ? "laplace" : "yukawa";
}

void run_field_case(KifmmFieldFamily family, KifmmKernel kernel,
                    double kappa) {
  const KifmmBox box{{0.15, -0.25, 0.40}, 1.0, 0, 7};
  const auto rule = make_cube_surface_rule(6, 1.05, 2.40);
  const auto sources = deterministic_sources(box);
  const auto targets = external_targets(box);
  double target_min_distance = std::numeric_limits<double>::infinity();
  double target_max_distance = 0.0;
  for (const Vec3 target : targets) {
    const double d = distance(target, box.center);
    target_min_distance = std::min(target_min_distance, d);
    target_max_distance = std::max(target_max_distance, d);
  }

  const auto result =
      solve_upward_equivalent_field(rule, box, sources, family, kernel, kappa,
                                    "qr");
  require(result.finite, std::string("field solve finite for ") +
                             family_name(family) + " " + kernel_name(kernel));
  require(result.status == "pass", "field solve status pass");
  require(result.rank > 0, "field solve rank recorded");
  require(std::isfinite(result.condition_estimate),
          "field solve condition finite");

  std::vector<double> approx_phi;
  std::vector<double> direct_phi;
  std::vector<Vec3> approx_grad;
  std::vector<Vec3> direct_grad;
  for (const Vec3 target : targets) {
    approx_phi.push_back(equivalent_field_potential(rule, box, result, family,
                                                    kernel, kappa, target));
    direct_phi.push_back(
        direct_field_potential(sources, family, kernel, kappa, target));
    approx_grad.push_back(equivalent_field_grad_target(
        rule, box, result, family, kernel, kappa, target));
    direct_grad.push_back(
        direct_field_grad_target(sources, family, kernel, kappa, target));
  }
  const double potential_rms = rel_rms(approx_phi, direct_phi);
  const double gradient_rms = rel_rms_vec(approx_grad, direct_grad);
  std::cout << std::setprecision(10)
            << "KIFMM_FIELD order=" << rule.order
            << " equiv_scale=" << rule.equiv_scale
            << " check_scale=" << rule.check_scale
            << " external_target_min_distance=" << target_min_distance
            << " external_target_max_distance=" << target_max_distance
            << " family=" << family_name(family)
            << " kernel=" << kernel_name(kernel) << " kappa=" << kappa
            << " solver=" << result.solver
            << " condition_estimate=" << result.condition_estimate
            << " rank=" << result.rank << " rhs_norm=" << result.rhs_norm
            << " charge_norm=" << result.charge_norm
            << " equiv_residual=" << result.relative_residual
            << " potential_rel_rms=" << potential_rms
            << " gradient_rel_rms=" << gradient_rms << '\n';
  require(std::isfinite(potential_rms) && potential_rms <= 1.0e-5,
          "field potential RMS gate");
  require(std::isfinite(gradient_rms) && gradient_rms <= 5.0e-5,
          "field gradient RMS gate");
}

} // namespace

int main() {
  for (const auto family : {KifmmFieldFamily::PhinCharge,
                           KifmmFieldFamily::PhiNormalX,
                           KifmmFieldFamily::PhiNormalY,
                           KifmmFieldFamily::PhiNormalZ}) {
    run_field_case(family, KifmmKernel::Laplace, 0.0);
    run_field_case(family, KifmmKernel::Yukawa, 0.12746);
  }

  const KifmmBox box{{0.0, 0.0, 0.0}, 1.0, 0, 1};
  const auto rule = make_cube_surface_rule(6, 1.05, 2.40);
  const auto sources = deterministic_sources(box);
  const auto bad = solve_upward_equivalent_field(
      rule, box, sources, KifmmFieldFamily::PhinCharge, KifmmKernel::Yukawa,
      -0.1, "qr");
  require(!bad.finite && bad.status == "invalid_kappa",
          "negative kappa rejected for field solve");

  std::cout << "PASS titan_fmm_kifmm_field_smoke\n";
  return EXIT_SUCCESS;
}
