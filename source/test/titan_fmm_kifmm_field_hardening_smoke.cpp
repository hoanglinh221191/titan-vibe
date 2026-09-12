#include "tabi/tabi_fmm_kifmm_fields.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmFieldFamily;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::assemble_kifmm_jcp_channels;
using titan::tabi::fmm_kifmm::grad_source_laplace_G;
using titan::tabi::fmm_kifmm::grad_source_yukawa_G;
using titan::tabi::fmm_kifmm::grad_target_of_source_gradient_component;
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

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

Vec3 add(Vec3 a, Vec3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 sub(Vec3 a, Vec3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

double component(Vec3 v, int c) {
  if (c == 0) {
    return v.x;
  }
  if (c == 1) {
    return v.y;
  }
  return v.z;
}

Vec3 grad_source_component_fd(KifmmKernel kernel, Vec3 target, Vec3 source,
                              double kappa, int source_component) {
  const double h = 2.0e-6 * std::max(1.0, norm(target));
  const Vec3 dx{h, 0.0, 0.0};
  const Vec3 dy{0.0, h, 0.0};
  const Vec3 dz{0.0, 0.0, h};
  auto source_grad_component = [&](Vec3 t) {
    const Vec3 g = kernel == KifmmKernel::Laplace
                       ? grad_source_laplace_G(t, source)
                       : grad_source_yukawa_G(t, source, kappa);
    return component(g, source_component);
  };
  return {(source_grad_component(add(target, dx)) -
           source_grad_component(sub(target, dx))) /
              (2.0 * h),
          (source_grad_component(add(target, dy)) -
           source_grad_component(sub(target, dy))) /
              (2.0 * h),
          (source_grad_component(add(target, dz)) -
           source_grad_component(sub(target, dz))) /
              (2.0 * h)};
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
  };
}

void check_mixed_derivative(KifmmKernel kernel, double kappa) {
  const Vec3 target{2.4, -1.7, 3.2};
  const Vec3 source{-0.2, 0.3, -0.4};
  for (int c = 0; c < 3; ++c) {
    const Vec3 analytic =
        grad_target_of_source_gradient_component(kernel, target, source, kappa,
                                                 c);
    const Vec3 fd = grad_source_component_fd(kernel, target, source, kappa, c);
    const Vec3 diff{sub(analytic, fd)};
    const double rel = norm(diff) / std::max(norm(fd), 1.0e-300);
    std::cout << std::setprecision(10)
              << "KIFMM_MIXED_DERIV kernel="
              << (kernel == KifmmKernel::Laplace ? "laplace" : "yukawa")
              << " component=" << c << " rel=" << rel << '\n';
    require(std::isfinite(rel) && rel <= 2.0e-8,
            "analytic mixed derivative matches finite-difference regression");
  }
}

} // namespace

int main() {
  check_mixed_derivative(KifmmKernel::Laplace, 0.0);
  check_mixed_derivative(KifmmKernel::Yukawa, 0.12746);

  const KifmmBox box{{0.15, -0.25, 0.40}, 1.0, 0, 7};
  const auto rule = make_cube_surface_rule(6, 1.05, 2.40);
  auto sources = deterministic_sources(box);

  auto bad_area = sources;
  bad_area[0].area = 0.0;
  const auto bad_area_result = solve_upward_equivalent_field(
      rule, box, bad_area, KifmmFieldFamily::PhinCharge,
      KifmmKernel::Laplace, 0.0, "qr");
  require(!bad_area_result.finite &&
              bad_area_result.status == "invalid_panel_source",
          "nonpositive panel area rejected");

  auto bad_normal = sources;
  bad_normal[0].normal = {};
  const auto bad_normal_result = solve_upward_equivalent_field(
      rule, box, bad_normal, KifmmFieldFamily::PhiNormalX,
      KifmmKernel::Laplace, 0.0, "qr");
  require(!bad_normal_result.finite &&
              bad_normal_result.status == "invalid_panel_source",
          "zero panel normal rejected");

  const auto channels =
      assemble_kifmm_jcp_channels(rule, box, sources, 0.12746, 80.0,
                                  1.0 / 80.0, "qr");
  require(channels.finite && channels.status == "pass",
          "hardened channel solve passes");
  require(channels.equivalent_system_builds == 2,
          "Laplace and Yukawa equivalent systems are factorized once each");
  require(channels.max_condition_estimate < 1.0e12,
          "condition gate remains below development limit");

  std::cout << "PASS titan_fmm_kifmm_field_hardening_smoke\n";
  return EXIT_SUCCESS;
}
