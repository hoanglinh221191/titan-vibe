#include "tabi/tabi_fmm_kifmm_solve.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmScalarSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::direct_scalar_grad_target;
using titan::tabi::fmm_kifmm::direct_scalar_potential;
using titan::tabi::fmm_kifmm::equivalent_scalar_grad_target;
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
double norm(Vec3 v) {
  return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double rel_grad_error(const std::vector<Vec3> &approx,
                      const std::vector<Vec3> &direct) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const Vec3 d{approx[i].x - direct[i].x, approx[i].y - direct[i].y,
                 approx[i].z - direct[i].z};
    err_sq += d.x * d.x + d.y * d.y + d.z * d.z;
    ref_sq += direct[i].x * direct[i].x + direct[i].y * direct[i].y +
              direct[i].z * direct[i].z;
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

std::vector<KifmmScalarSource> sources_inside(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {{box.center.x - 0.16 * h, box.center.y + 0.09 * h,
        box.center.z - 0.06 * h},
       0.70},
      {{box.center.x + 0.12 * h, box.center.y - 0.15 * h,
        box.center.z + 0.10 * h},
       -0.55},
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
  };
}

} // namespace

int main() {
  const KifmmBox box{{0.15, -0.25, 0.40}, 1.0, 0, 7};
  const auto rule = make_cube_surface_rule(6, 1.05, 2.40);
  const auto sources = sources_inside(box);

  const auto invalid_kappa = solve_upward_equivalent(
      rule, box, sources, KifmmKernel::Yukawa, -0.10, "qr");
  require(!invalid_kappa.finite, "negative kappa must fail");
  require(invalid_kappa.status == "invalid_kappa",
          "negative kappa status is invalid_kappa");

  const auto unknown_solver = solve_upward_equivalent(
      rule, box, sources, KifmmKernel::Laplace, 0.0, "made-up-solver");
  require(!unknown_solver.finite, "unknown solver must fail");
  require(unknown_solver.status == "unknown_solver",
          "unknown solver status is recorded");

  auto collision_rule = rule;
  collision_rule.unit_check = collision_rule.unit_equiv;
  const auto collision = solve_upward_equivalent(
      collision_rule, box, sources, KifmmKernel::Laplace, 0.0, "qr");
  require(!collision.finite, "check/equiv collision must fail");
  require(collision.status == "surface_collision",
          "collision status is recorded");
  require(collision.diagnostics.min_check_equiv_distance <= 1.0e-12,
          "collision diagnostic distance is small");

  auto outside_sources = sources;
  outside_sources.push_back({{box.center.x + 1.20 * box.half_width,
                              box.center.y, box.center.z},
                             0.10});
  const auto source_outside = solve_upward_equivalent(
      rule, box, outside_sources, KifmmKernel::Laplace, 0.0, "qr");
  require(!source_outside.finite, "source outside source box must fail");
  require(source_outside.status == "source_outside_box",
          "source outside status is recorded");

  const auto pass =
      solve_upward_equivalent(rule, box, sources, KifmmKernel::Yukawa, 0.12746,
                              "qr");
  require(pass.finite, "valid solve must pass");
  require(pass.status == "pass", "valid solve status is pass");
  require(pass.diagnostics.rank > 0, "rank diagnostic recorded");
  require(pass.diagnostics.rhs_norm > 0.0, "rhs norm recorded");
  require(pass.diagnostics.charge_norm > 0.0, "charge norm recorded");
  require(pass.diagnostics.condition_estimate > 1.0,
          "condition estimate recorded");
  require(pass.diagnostics.min_check_equiv_distance > 0.0,
          "min check/equiv distance recorded");
  require(pass.diagnostics.min_source_equiv_distance > 0.0,
          "min source/equiv distance recorded");

  std::vector<Vec3> approx_grad;
  std::vector<Vec3> direct_grad;
  for (const Vec3 target : external_targets(box)) {
    approx_grad.push_back(
        equivalent_scalar_grad_target(rule, box, pass, KifmmKernel::Yukawa,
                                      0.12746, target));
    direct_grad.push_back(
        direct_scalar_grad_target(sources, KifmmKernel::Yukawa, 0.12746,
                                  target));
  }
  const double grad_rms = rel_grad_error(approx_grad, direct_grad);
  require(std::isfinite(grad_rms), "gradient RMS finite");
  require(grad_rms <= 5.0e-6, "gradient RMS gate");

  const Vec3 target = external_targets(box).front();
  const double direct = direct_scalar_potential(sources, KifmmKernel::Yukawa,
                                                0.12746, target);
  const double approx = equivalent_scalar_potential(
      rule, box, pass, KifmmKernel::Yukawa, 0.12746, target);
  require(std::isfinite(direct) && std::isfinite(approx),
          "potential evaluators remain finite");
  require(norm(approx_grad.front()) > 0.0, "equivalent gradient nonzero");

  std::cout << "PASS titan_fmm_kifmm_hardening_smoke\n";
  return EXIT_SUCCESS;
}
