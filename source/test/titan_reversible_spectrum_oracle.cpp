// KIN-09 oracle: titan_markov::reversible_spectrum against Eigen's general
// eigensolver on random reversible transition matrices, and its refusal of
// matrices that are not in detailed balance.
#include "titan_markov_estimator.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char *what, int n) {
  std::printf("%s n=%d %s\n", ok ? "PASS" : "FAIL", n, what);
  if (!ok)
    ++failures;
}

void reversible_case(int n, std::mt19937 &rng) {
  // symmetric positive flux x_ij (some pairs absent), T = x / row sum,
  // pi = row sum / total: detailed balance by construction
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j) {
      const double v = (j == i + 1 || j == i || unit(rng) < 0.3)
                           ? std::exp(6.0 * unit(rng) - 3.0)
                           : 0.0;
      x(i, j) = x(j, i) = v;
    }
  const Eigen::VectorXd rows = x.rowwise().sum();
  const Eigen::MatrixXd t = x.array().colwise() / rows.array();
  const Eigen::VectorXd pi = rows / rows.sum();

  Eigen::VectorXd values;
  Eigen::MatrixXd vectors;
  const bool solved =
      titan_markov::reversible_spectrum(t, pi, values, &vectors);
  check(solved, "reversible matrix accepted", n);
  if (!solved)
    return;

  Eigen::EigenSolver<Eigen::MatrixXd> general(t, false);
  std::vector<double> reference;
  double max_imag = 0.0;
  for (int k = 0; k < n; ++k) {
    reference.push_back(general.eigenvalues()[k].real());
    max_imag = std::max(max_imag, std::abs(general.eigenvalues()[k].imag()));
  }
  std::vector<double> mine(values.data(), values.data() + n);
  std::sort(reference.begin(), reference.end());
  std::sort(mine.begin(), mine.end());
  double value_error = 0.0;
  for (int k = 0; k < n; ++k)
    value_error = std::max(value_error, std::abs(reference[k] - mine[k]));
  check(value_error < 1e-11 && max_imag < 1e-9,
        "eigenvalues equal the general solver's", n);

  double residual = 0.0, orthogonality = 0.0, norm_error = 0.0;
  for (int k = 0; k < n; ++k) {
    residual = std::max(
        residual, (t * vectors.col(k) - values[k] * vectors.col(k)).norm());
    norm_error = std::max(norm_error, std::abs(vectors.col(k).norm() - 1.0));
    for (int l = 0; l < k; ++l) {
      const double inner =
          vectors.col(k).dot(pi.asDiagonal() * vectors.col(l));
      const double scale =
          std::sqrt(vectors.col(k).dot(pi.asDiagonal() * vectors.col(k)) *
                    vectors.col(l).dot(pi.asDiagonal() * vectors.col(l)));
      orthogonality = std::max(orthogonality, std::abs(inner) / scale);
    }
  }
  check(residual < 1e-10, "right eigenvectors satisfy T psi = lambda psi", n);
  check(orthogonality < 1e-10, "right eigenvectors are pi-orthogonal", n);
  check(norm_error < 1e-12, "right eigenvectors have unit 2-norm", n);

  Eigen::VectorXd only_values;
  check(titan_markov::reversible_spectrum(t, pi, only_values, nullptr) &&
            (only_values - values).cwiseAbs().maxCoeff() < 1e-12,
        "eigenvalue-only call agrees", n);
}

void refusal_cases(std::mt19937 &rng) {
  std::uniform_real_distribution<double> unit(0.1, 1.0);
  const int n = 6;
  Eigen::MatrixXd t(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      t(i, j) = unit(rng);
  t = (t.array().colwise() / t.rowwise().sum().array()).matrix();
  Eigen::VectorXd pi;
  std::string error;
  titan_markov::stationary(t, pi, error);
  Eigen::VectorXd values;
  check(!titan_markov::reversible_spectrum(t, pi, values, nullptr),
        "non-reversible matrix refused (caller uses the general solver)", n);

  Eigen::VectorXd zero_pi = Eigen::VectorXd::Constant(n, 1.0 / (n - 1));
  zero_pi[2] = 0.0;
  check(!titan_markov::reversible_spectrum(t, zero_pi, values, nullptr),
        "non-positive stationary vector refused", n);
}

} // namespace

int main() {
  std::mt19937 rng(20261023);
  for (int n : {2, 3, 10, 60, 200})
    reversible_case(n, rng);
  refusal_cases(rng);
  std::printf("%s\n", failures == 0 ? "ALL PASS" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
