#include "tabi/tabi_fmm_kifmm_kernels.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::grad_source_laplace_G;
using titan::tabi::fmm_kifmm::grad_source_yukawa_G;
using titan::tabi::fmm_kifmm::grad_target_laplace_G;
using titan::tabi::fmm_kifmm::grad_target_yukawa_G;
using titan::tabi::fmm_kifmm::laplace_G;
using titan::tabi::fmm_kifmm::yukawa_G;

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

Vec3 add(Vec3 a, Vec3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 scale(Vec3 a, double s) {
  return {a.x * s, a.y * s, a.z * s};
}

void require_close(double actual, double expected, double tol,
                   const std::string &message) {
  const double denom = std::max(1.0, std::abs(expected));
  if (!std::isfinite(actual) ||
      std::abs(actual - expected) > tol * denom) {
    std::cerr << "FAIL: " << message << " actual=" << actual
              << " expected=" << expected << " tol=" << tol << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void require_vec_close(Vec3 actual, Vec3 expected, double tol,
                       const std::string &message) {
  require_close(actual.x, expected.x, tol, message + ".x");
  require_close(actual.y, expected.y, tol, message + ".y");
  require_close(actual.z, expected.z, tol, message + ".z");
}

template <class Fn>
double centered_difference(Fn &&fn, Vec3 x, Vec3 direction) {
  constexpr double h = 1.0e-6;
  return (fn(add(x, scale(direction, h))) -
          fn(add(x, scale(direction, -h)))) /
         (2.0 * h);
}

} // namespace

int main() {
  const Vec3 target{1.3, -0.7, 2.1};
  const Vec3 source{-0.4, 0.25, -0.9};
  const Vec3 probe_dir{0.3, -0.2, 0.5};
  const double kappa = 0.12746;

  const double lap = laplace_G(target, source);
  const double yuk = yukawa_G(target, source, kappa);
  require(std::isfinite(lap) && lap > 0.0, "Laplace potential finite");
  require(std::isfinite(yuk) && yuk > 0.0, "Yukawa potential finite");
  require_close(yukawa_G(target, source, 0.0), lap, 1.0e-13,
                "Yukawa kappa=0 matches Laplace");

  const Vec3 grad_t_lap = grad_target_laplace_G(target, source);
  const Vec3 grad_s_lap = grad_source_laplace_G(target, source);
  require_vec_close(grad_s_lap, scale(grad_t_lap, -1.0), 1.0e-12,
                    "Laplace source gradient is negative target gradient");

  const double fd_lap = centered_difference(
      [&](Vec3 x) { return laplace_G(x, source); }, target, probe_dir);
  require_close(dot(grad_t_lap, probe_dir), fd_lap, 5.0e-8,
                "Laplace target gradient finite difference");

  const Vec3 grad_t_yuk = grad_target_yukawa_G(target, source, kappa);
  const Vec3 grad_s_yuk = grad_source_yukawa_G(target, source, kappa);
  require_vec_close(grad_s_yuk, scale(grad_t_yuk, -1.0), 1.0e-12,
                    "Yukawa source gradient is negative target gradient");

  const double fd_yuk = centered_difference(
      [&](Vec3 x) { return yukawa_G(x, source, kappa); }, target, probe_dir);
  require_close(dot(grad_t_yuk, probe_dir), fd_yuk, 5.0e-8,
                "Yukawa target gradient finite difference");

  std::cout << "PASS titan_fmm_kifmm_kernel_smoke\n";
  return EXIT_SUCCESS;
}
