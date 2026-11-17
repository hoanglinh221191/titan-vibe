// KIN-13 oracle for TitanTram (src/titan_tram.*), the -dham estimator.
//
// Model: 9 microstates on a ring with energies E_m, lumped into 3 Markov
// states (m / 3). At temperature T the microstate dynamics is Metropolis on
// the ring, in detailed balance with exp(-E_m / kT). Within a Markov state
// the microstate energies differ by several kT, so the bias of one
// temperature relative to another is not a function of the Markov state.
//
//  1. Exact data: expected lag-1 counts C^k_ij = N_k sum_{m in i, n in j}
//     pi^k_m P^k_mn and one sample per microstate with multiplicity
//     N_k pi^k_m. TRAM must return the exact target distribution at a
//     simulated and at an unsimulated temperature, the exact lumped
//     transition matrix of each ensemble, and sample weights that add up to
//     the exact microstate populations (1e-8). The state-mean small-bias
//     estimate that -dham used before (Rosta-Hummer with b = (beta_k -
//     beta_t) <E>_k(i)) is evaluated on the same data and must miss.
//  2. Sampled data at three temperatures: TRAM within 0.01 of the exact
//     target populations.
//  3. Invariance to an energy offset; identical results with 1 and 4
//     threads; one ensemble at the target gives the reversible MSM.
//  4. Refusals: no transitions, a transition into a state without samples,
//     states not connected by transitions.
#include "titan_tram.h"
#include "titan_executor.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const std::string &what, double value) {
  std::printf("%s %s (%.3e)\n", ok ? "PASS" : "FAIL", what.c_str(), value);
  if (!ok)
    ++failures;
}

constexpr int kMicro = 9;
constexpr int kStates = 3;
const double kEnergy[kMicro] = {0.0, 6.0, 2.0, 3.0, -1.0, 8.0, 5.0, 1.5, 9.0};

int lump(int m) { return m / 3; }

double beta(double t) { return 1.0 / (TitanTram::kBoltzmannKJ * t); }

Eigen::VectorXd micro_pi(double t) {
  Eigen::VectorXd p(kMicro);
  for (int m = 0; m < kMicro; ++m)
    p[m] = std::exp(-beta(t) * kEnergy[m]);
  return p / p.sum();
}

// Metropolis on a ring: move to m +- 1 with probability 1/2 min(1, ...)
Eigen::MatrixXd micro_p(double t) {
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(kMicro, kMicro);
  for (int m = 0; m < kMicro; ++m) {
    for (int d : {-1, 1}) {
      const int n = (m + d + kMicro) % kMicro;
      p(m, n) = 0.5 * std::min(1.0, std::exp(-beta(t) * (kEnergy[n] - kEnergy[m])));
    }
    p(m, m) = 1.0 - p.row(m).sum();
  }
  return p;
}

Eigen::VectorXd exact_state_pi(double t) {
  const Eigen::VectorXd p = micro_pi(t);
  Eigen::VectorXd s = Eigen::VectorXd::Zero(kStates);
  for (int m = 0; m < kMicro; ++m)
    s[lump(m)] += p[m];
  return s;
}

Eigen::MatrixXd exact_counts(double t, double n) {
  const Eigen::VectorXd p = micro_pi(t);
  const Eigen::MatrixXd tm = micro_p(t);
  Eigen::MatrixXd c = Eigen::MatrixXd::Zero(kStates, kStates);
  for (int m = 0; m < kMicro; ++m)
    for (int q = 0; q < kMicro; ++q)
      c(lump(m), lump(q)) += n * p[m] * tm(m, q);
  return c;
}

void load_exact(TitanTram &tram, const std::vector<double> &temps,
                double offset = 0.0) {
  for (std::size_t k = 0; k < temps.size(); ++k) {
    const double n = 1.0e6 * (1.0 + 0.3 * static_cast<double>(k));
    const Eigen::VectorXd p = micro_pi(temps[k]);
    for (int m = 0; m < kMicro; ++m)
      tram.add_sample(static_cast<int>(k), lump(m), kEnergy[m] + offset, n * p[m]);
    const Eigen::MatrixXd c = exact_counts(temps[k], n);
    for (int i = 0; i < kStates; ++i)
      for (int j = 0; j < kStates; ++j)
        tram.add_transition(static_cast<int>(k), i, j, c(i, j));
  }
}

// The small-bias DHAM estimate with state-mean energies that -dham used
// before KIN-13, on the same exact data.
Eigen::VectorXd state_mean_dham(const std::vector<double> &temps, double target) {
  Eigen::MatrixXd numerator = Eigen::MatrixXd::Zero(kStates, kStates);
  Eigen::MatrixXd denominator = Eigen::MatrixXd::Zero(kStates, kStates);
  for (double t : temps) {
    const double n = 1.0e6;
    const Eigen::VectorXd p = micro_pi(t);
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(kStates), pop = Eigen::VectorXd::Zero(kStates);
    for (int m = 0; m < kMicro; ++m) {
      mean[lump(m)] += p[m] * kEnergy[m];
      pop[lump(m)] += p[m];
    }
    mean = mean.cwiseQuotient(pop);
    const Eigen::MatrixXd c = exact_counts(t, n);
    for (int i = 0; i < kStates; ++i)
      for (int j = 0; j < kStates; ++j) {
        numerator(i, j) += c(i, j);
        const double bi = (beta(t) - beta(target)) * mean[i];
        const double bj = (beta(t) - beta(target)) * mean[j];
        denominator(i, j) += c.row(i).sum() * std::exp(-0.5 * (bj - bi));
      }
  }
  Eigen::MatrixXd m = numerator.cwiseQuotient(denominator);
  for (int i = 0; i < kStates; ++i)
    m.row(i) /= m.row(i).sum();
  Eigen::VectorXd pi = Eigen::VectorXd::Constant(kStates, 1.0 / kStates);
  for (int it = 0; it < 100000; ++it)
    pi = (m.transpose() * pi).eval();
  return pi / pi.sum();
}

void exact_case(double target) {
  const std::vector<double> temps{290.0, 310.0, 340.0};
  TitanTram tram(kStates, temps, target);
  load_exact(tram, temps);
  const auto report = tram.solve(200000, 1e-13);
  const std::string at = " (target " + std::to_string(static_cast<int>(target)) + " K)";
  check(report.converged, "exact data: TRAM converges in " +
                              std::to_string(report.iterations) + " iterations" + at,
        report.max_change);
  if (!report.converged)
    return;
  const double pi_error =
      (tram.target_stationary() - exact_state_pi(target)).cwiseAbs().maxCoeff();
  check(pi_error < 1e-8, "exact data: target populations exact" + at, pi_error);

  double tpm_error = 0.0;
  for (int k = 0; k < 3; ++k) {
    Eigen::MatrixXd expected = exact_counts(temps[static_cast<std::size_t>(k)], 1.0);
    for (int i = 0; i < kStates; ++i)
      expected.row(i) /= expected.row(i).sum();
    tpm_error = std::max(tpm_error,
                         (tram.transition_matrix(k) - expected).cwiseAbs().maxCoeff());
  }
  check(tpm_error < 1e-8, "exact data: lumped transition matrix of every ensemble" + at,
        tpm_error);

  // sample weights: one sample per (ensemble, microstate), insertion order
  const std::vector<double> w = tram.sample_weights();
  Eigen::VectorXd micro = Eigen::VectorXd::Zero(kMicro);
  for (std::size_t x = 0; x < w.size(); ++x)
    micro[static_cast<int>(x % kMicro)] += w[x];
  const double weight_error = (micro - micro_pi(target)).cwiseAbs().maxCoeff();
  check(weight_error < 1e-8, "exact data: sample weights give microstate populations" + at,
        weight_error);

  const double old_error =
      (state_mean_dham(temps, target) - exact_state_pi(target)).cwiseAbs().maxCoeff();
  check(old_error > 1e-5, "exact data: the old state-mean DHAM estimate misses by" + at,
        old_error);
}

void sampled_case() {
  const std::vector<double> temps{290.0, 310.0, 340.0};
  const double target = 300.0;
  TitanTram tram(kStates, temps, target);
  std::mt19937_64 rng(20261023);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (std::size_t k = 0; k < temps.size(); ++k) {
    const Eigen::MatrixXd p = micro_p(temps[k]);
    int m = 0;
    int previous = -1;
    for (int step = 0; step < 400000; ++step) {
      tram.add_sample(static_cast<int>(k), lump(m), kEnergy[m]);
      if (previous >= 0)
        tram.add_transition(static_cast<int>(k), lump(previous), lump(m));
      previous = m;
      const double r = u(rng);
      double acc = 0.0;
      int next = m;
      for (int q = 0; q < kMicro; ++q) {
        acc += p(m, q);
        if (r < acc) {
          next = q;
          break;
        }
      }
      m = next;
    }
  }
  const auto report = tram.solve(100000, 1e-10);
  const double error = report.converged
                           ? (tram.target_stationary() - exact_state_pi(target))
                                 .cwiseAbs()
                                 .maxCoeff()
                           : 1.0;
  check(report.converged && error < 0.01,
        "sampled data (3 x 4e5 steps): target populations within 0.01", error);
}

void invariance_cases() {
  const std::vector<double> temps{290.0, 310.0, 340.0};
  TitanTram a(kStates, temps, 300.0), b(kStates, temps, 300.0);
  load_exact(a, temps);
  load_exact(b, temps, 1.0e5);
  a.solve(200000, 1e-12);
  b.solve(200000, 1e-12);
  const double offset_error =
      (a.target_stationary() - b.target_stationary()).cwiseAbs().maxCoeff();
  check(offset_error < 1e-10, "energy offset of 1e5 kJ/mol changes nothing", offset_error);

  TitanTram serial(kStates, temps, 300.0), threaded(kStates, temps, 300.0);
  load_exact(serial, temps);
  load_exact(threaded, temps);
  TitanExecutor one(1), four(4);
  serial.solve(5000, 1e-12, &one);
  threaded.solve(5000, 1e-12, &four);
  check(serial.target_stationary() == threaded.target_stationary() &&
            serial.sample_weights() == threaded.sample_weights(),
        "identical results with 1 and 4 threads", 0.0);

  // one ensemble at the target: the reversible MSM of the counts
  TitanTram single(2, {300.0}, 300.0);
  const double c[2][2] = {{90, 10}, {10, 40}};
  for (int i = 0; i < 2; ++i) {
    single.add_sample(0, i, 3.0 * i, i == 0 ? 100.0 : 50.0);
    for (int j = 0; j < 2; ++j)
      single.add_transition(0, i, j, c[i][j]);
  }
  const auto r = single.solve(10000, 1e-13);
  const Eigen::MatrixXd t = single.transition_matrix(0);
  const double msm_error = std::max({std::abs(t(0, 1) - 0.1), std::abs(t(1, 0) - 0.2),
                                     std::abs(single.target_stationary()[0] - 2.0 / 3.0)});
  check(r.converged && msm_error < 1e-10,
        "one ensemble at the target gives the reversible MSM", msm_error);
}

void refusal_cases() {
  {
    TitanTram tram(2, {300.0, 320.0}, 300.0);
    tram.add_sample(0, 0, 1.0);
    tram.add_sample(0, 1, 2.0);
    const auto r = tram.solve(100, 1e-8);
    check(!r.converged && r.message.find("transition") != std::string::npos,
          "no transitions refused: " + r.message, 0.0);
  }
  {
    TitanTram tram(2, {300.0}, 300.0);
    tram.add_sample(0, 0, 1.0);
    tram.add_transition(0, 0, 1);
    const auto r = tram.solve(100, 1e-8);
    check(!r.converged && r.message.find("samples") != std::string::npos,
          "transition into a state without samples refused: " + r.message, 0.0);
  }
  {
    TitanTram tram(3, {300.0}, 300.0);
    for (int i = 0; i < 3; ++i)
      tram.add_sample(0, i, 1.0, 2.0);
    tram.add_transition(0, 0, 1);
    tram.add_transition(0, 1, 0);
    tram.add_transition(0, 2, 2);
    const auto r = tram.solve(100, 1e-8);
    check(!r.converged && r.message.find("connected") != std::string::npos,
          "states not connected by transitions refused: " + r.message, 0.0);
  }
}

} // namespace

int main() {
  exact_case(310.0); // a simulated temperature
  exact_case(300.0); // between two simulated temperatures
  sampled_case();
  invariance_cases();
  refusal_cases();
  std::printf("%s\n", failures == 0 ? "ALL PASS" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
