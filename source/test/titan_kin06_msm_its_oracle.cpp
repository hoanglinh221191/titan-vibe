// KIN-06 oracle: implied timescales (ITS) when |lambda| = 1 or lambda is
// complex, in every ITS path, rank by rank.
//
// Convention: eigenvalues sorted by |lambda| (descending) after the one
// nearest to +1 (the stationary mode, rank 0); rank k >= 1 gives
// t_k = -lag dt / ln|lambda_k|, or NaN when |lambda_k - 1| <= 1e-8, |lambda_k|
// >= 1, |Im lambda_k| > 1e-6 or |lambda_k| < 1e-15. A NaN keeps its rank.
//
//   periodic  production MSM (TitanMSM, reversible MLE) of a bipartite 4-state
//             chain: every step goes from {0,1} to {2,3} and back, so the
//             eigenvalues are 1, -1, mu, -mu. Expected ITS: NaN, t(mu), t(mu).
//             compute_implied_timescales, msm_timescales.dat and the ITS scan
//             (msm_its.dat, lag 1) must all give that, rank by rank.
//   blockdiag two disconnected 2-state blocks (eigenvalues 1, 1, 0.9, 0.7)
//             through the validation path compute_implied_timescales_from_tpm:
//             NaN, t(0.9), t(0.7); two eigenvalues within 1e-8 of 1.
//   complex   (new helper only) a 3-cycle with complex eigenvalues: NaN, NaN.
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "info.h"
#include "titan_executor.h"
#include "titan_msm.h"
#include "titan_validation_core.h"
#if __has_include("titan_its.h")
#include "titan_its.h"
#define KIN06_HAVE_ITS_HELPER 1
#endif

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

bool same(const std::vector<double> &got, const std::vector<double> &want) {
  if (got.size() != want.size())
    return false;
  for (std::size_t i = 0; i < got.size(); ++i) {
    if (std::isnan(want[i]) != std::isnan(got[i]))
      return false;
    if (!std::isnan(want[i]) &&
        std::abs(got[i] - want[i]) > 1e-6 * std::max(1.0, std::abs(want[i])))
      return false;
  }
  return true;
}

std::string show(const std::vector<double> &v) {
  std::ostringstream o;
  o << std::setprecision(6) << '[';
  for (std::size_t i = 0; i < v.size(); ++i)
    o << (i ? ", " : "") << v[i];
  o << ']';
  return o.str();
}

// second column of every non-comment row
std::vector<double> read_column(const std::string &path, int column) {
  std::ifstream in(path);
  std::string line;
  std::vector<double> out;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    std::string token;
    for (int c = 0; c <= column && (row >> token); ++c)
      if (c == column)
        out.push_back(token == "nan" || token == "-nan" ? NAN : std::stod(token));
  }
  return out;
}

std::vector<double> its_row(const std::string &path, int lag) {
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    int l;
    double lag_ns;
    if (!(row >> l >> lag_ns) || l != lag)
      continue;
    std::vector<double> out;
    std::string token;
    while (row >> token)
      out.push_back(token == "nan" || token == "-nan" ? NAN : std::stod(token));
    return out;
  }
  return {};
}

} // namespace

int main() {
  // ---------------------------------------------------------------- periodic
  constexpr double dt_ns = 0.1;
  constexpr int n_frames = 4000;
  const double centre[4] = {-15.0, -5.0, 5.0, 15.0};
  auto msm = std::make_unique<TitanMSM>(4, 1, dt_ns);
  std::vector<int> truth(n_frames);
  unsigned x = 12345u;
  int state = 0;
  for (int i = 0; i < n_frames; ++i) {
    truth[static_cast<std::size_t>(i)] = state;
    Eigen::VectorXd p(1);
    p(0) = centre[state] + 0.01 * std::sin(0.7 * i);
    msm->add_projected_frame(p, 0, i);
    msm->record_frame_time(100.0 * i, 0);
    x = 1103515245u * x + 12345u;
    const double u = ((x >> 8) & 0xffffu) / 65536.0;
    // A = {0,1} -> B = {2,3} -> A; first member preferred from first member
    const bool first = (state % 2 == 0) ? (u < 0.8) : (u < 0.3);
    state = (state < 2) ? (first ? 2 : 3) : (first ? 0 : 1);
  }
  msm->add_trajectory_break();
  TitanExecutor executor(1);
  msm->train_kmeans(executor, 100, 42u);
  msm->build_markov_model(1);
  if (!msm->is_built()) {
    std::cout << "FAIL KIN-06 synthetic periodic MSM was not built\n";
    return 1;
  }
  // k-means must recover the four states
  std::vector<int> label_of(4, -1);
  bool clusters_ok = true;
  for (int i = 0; i < n_frames; ++i) {
    const int t = truth[static_cast<std::size_t>(i)];
    const int m = msm->discrete_traj()[static_cast<std::size_t>(i)];
    if (label_of[t] < 0)
      label_of[t] = m;
    clusters_ok = clusters_ok && label_of[t] == m;
  }
  check(clusters_ok, "periodic: k-means recovers the 4 states");

  Eigen::EigenSolver<Eigen::MatrixXd> es(msm->tpm(), false);
  std::vector<std::complex<double>> ev(es.eigenvalues().data(),
                                       es.eigenvalues().data() + 4);
  std::sort(ev.begin(), ev.end(), [](auto a, auto b) {
    return std::abs(a) > std::abs(b);
  });
  double mu = 0.0;
  for (auto l : ev)
    if (std::abs(std::abs(l) - 1.0) > 1e-6)
      mu = std::max(mu, std::abs(l));
  const double t_mu = -dt_ns / std::log(mu);
  const std::vector<double> want_periodic = {NAN, t_mu, t_mu};
  std::cout << std::setprecision(6) << "case periodic: eigenvalues";
  for (auto l : ev)
    std::cout << ' ' << l.real();
  std::cout << "; expected ITS " << show(want_periodic) << " ns\n";

  const auto its_direct = msm->compute_implied_timescales(3);
  std::cout << "  compute_implied_timescales: " << show(its_direct) << '\n';
  check(same(its_direct, want_periodic),
        "periodic: compute_implied_timescales keeps the rank of lambda = -1 (NaN)");

  arguments args;
  info info_output(args);
  std::fstream log("log.txt", std::ios::out);
  msm->export_results(executor, info_output, log, 3, true, false, -1, {1},
                      TitanMSMItsSolverMode::Full, 0, 256, true);
  const auto table = read_column(info_output.get_msm_timescales_filename(), 1);
  const auto scan = its_row(info_output.get_msm_its_filename(), 1);
  std::cout << "  msm_timescales.dat: " << show(table) << '\n'
            << "  ITS scan lag 1:     " << show(scan) << '\n';
  check(same(table, want_periodic), "periodic: msm_timescales.dat rank by rank");
  check(same(scan, want_periodic), "periodic: ITS scan rank by rank");
#ifdef KIN06_HAVE_ITS_HELPER
  check(msm->unit_eigenvalue_count() == 1,
        "periodic: one eigenvalue within 1e-8 of 1 (-1 is not counted)");
#endif

  // --------------------------------------------------------------- blockdiag
  Eigen::MatrixXd block = Eigen::MatrixXd::Zero(4, 4);
  block << 0.9, 0.1, 0.0, 0.0, //
      0.2, 0.8, 0.0, 0.0,      //
      0.0, 0.0, 0.95, 0.05,    //
      0.0, 0.0, 0.05, 0.95;
  const std::vector<double> want_block = {NAN, -dt_ns / std::log(0.9),
                                          -dt_ns / std::log(0.7)};
  const auto its_block =
      TitanValidation::compute_implied_timescales_from_tpm(block, dt_ns, 3);
  std::cout << "case blockdiag: expected " << show(want_block)
            << ", validation path " << show(its_block) << '\n';
  check(same(its_block, want_block),
        "blockdiag: second unit eigenvalue keeps its rank (NaN)");
#ifdef KIN06_HAVE_ITS_HELPER
  {
    Eigen::EigenSolver<Eigen::MatrixXd> eb(block, false);
    check(titan_its::count_unit_eigenvalues(eb.eigenvalues()) == 2,
          "blockdiag: two eigenvalues within 1e-8 of 1 reported");
    const auto via_helper =
        titan_its::implied_timescales(eb.eigenvalues(), dt_ns, 3);
    check(same(via_helper, want_block), "blockdiag: shared helper agrees");
  }
  // ----------------------------------------------------------------- complex
  {
    Eigen::MatrixXd cyc(3, 3);
    cyc << 0.1, 0.9, 0.0, //
        0.0, 0.1, 0.9,    //
        0.9, 0.0, 0.1;
    Eigen::EigenSolver<Eigen::MatrixXd> ec(cyc, false);
    const auto its_c = titan_its::implied_timescales(ec.eigenvalues(), dt_ns, 2);
    std::cout << "case complex: " << show(its_c) << '\n';
    check(its_c.size() == 2 && std::isnan(its_c[0]) && std::isnan(its_c[1]),
          "complex: NaN for the complex pair, ranks kept");
    check(titan_its::count_unit_eigenvalues(ec.eigenvalues()) == 1,
          "complex: one eigenvalue within 1e-8 of 1");
  }
#else
  check(false, "no shared ITS helper reporting the multiplicity of eigenvalue 1");
#endif

  if (failures) {
    std::cout << "FAIL KIN-06 MSM ITS: " << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-06 MSM ITS: NaN with the rank kept for |lambda| = 1 "
               "and complex modes in every path; multiplicity of 1 reported\n";
  return 0;
}
