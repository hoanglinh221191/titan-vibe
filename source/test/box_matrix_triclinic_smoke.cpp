#include "BoxMatrix.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void require_close(double actual, double expected, const char *label) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > 1.0e-10) {
    std::cerr << "FAIL box_matrix_triclinic_smoke: " << label
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_triclinic_box() {
  box value;
  value.set_a(20.0);
  value.set_b(20.0);
  value.set_c(20.0);
  value.set_alpha(90.0);
  value.set_beta(90.0);
  value.set_gamma(60.0);
  return value;
}

} // namespace

int main() {
  const BoxMatrix matrix(make_triclinic_box());
  double dx = 11.0;
  double dy = 10.0 * std::sqrt(3.0);
  double dz = 0.0;
  matrix.apply_mic(dx, dy, dz);
  require_close(dx, 1.0, "scalar x");
  require_close(dy, 0.0, "scalar y");
  require_close(dz, 0.0, "scalar z");

#if TITAN_COMPILED_WITH_AVX
  alignas(32) double xs[4]{11.0, -11.0, 1.0, 0.0};
  alignas(32) double ys[4]{10.0 * std::sqrt(3.0),
                           -10.0 * std::sqrt(3.0), 0.0, 0.0};
  alignas(32) double zs[4]{};
  __m256d vx = _mm256_load_pd(xs);
  __m256d vy = _mm256_load_pd(ys);
  __m256d vz = _mm256_load_pd(zs);
  matrix.apply_mic_avx(vx, vy, vz);
  _mm256_store_pd(xs, vx);
  _mm256_store_pd(ys, vy);
  _mm256_store_pd(zs, vz);
  require_close(xs[0], 1.0, "AVX lane 0 x");
  require_close(ys[0], 0.0, "AVX lane 0 y");
  require_close(xs[1], -1.0, "AVX lane 1 x");
  require_close(ys[1], 0.0, "AVX lane 1 y");
  require_close(xs[2], 1.0, "AVX lane 2 x");
#endif

  std::cout << "PASS box_matrix_triclinic_smoke\n";
  return 0;
}
