// KIN-20 oracle: MBAR bootstrap errors of autocorrelated energies.
//
// One replica at 300 K, N = 20000 frames, energy E_t = E0 + s x_t with x an
// AR(1) process x_t = phi x_{t-1} + e_t (unit Gaussian innovations,
// stationary start). At the sampled temperature MBAR gives the plain mean, so
// the standard error of <E> is that of the mean of an AR(1) series:
//   SE^2 = s^2 Var(x) / N [1 + 2 sum_{k=1}^{N-1} (1 - k/N) phi^k],
//   Var(x) = 1/(1 - phi^2),
// i.e. sqrt(g) times the i.i.d. error, g = (1 + phi)/(1 - phi) (19 for
// phi = 0.9, 3 for phi = 0.5). The bootstrap StdE (400 replicates) must be
// within [0.80, 1.20] of SE. Justification of the window: the moving-block
// bootstrap with block length L underestimates the variance by about
// tau/L = g/(2L) (block length 4 g: -12 % in variance, -6 % in SE), and its
// sampling scatter for one series is about 5 % (n/L ~ 200 blocks, 400
// replicates); a 20000-frame simulation of 40 series gave 0.944 +- 0.056.
// Frame-wise i.i.d. resampling gives 1/sqrt(g) of SE (0.23 and 0.58).
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "arguments.h"
#include "titan_executor.h"
#include "titan_mbar.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kFrames = 20000;
constexpr double kE0 = -5000.0;
constexpr double kScale = 30.0;

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

double analytic_se(double phi) {
  long double sum = 0.0L, pk = 1.0L;
  for (int k = 1; k < kFrames; ++k) {
    pk *= phi;
    sum += (1.0L - static_cast<long double>(k) / kFrames) * pk;
  }
  const long double var = kScale * kScale / (1.0L - phi * phi);
  return static_cast<double>(std::sqrt(var / kFrames * (1.0L + 2.0L * sum)));
}

void run_case(const std::string &tag, double phi, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  {
    std::ofstream out("e_" + tag + ".xvg");
    out << "# synthetic AR(1) energies\n@    title \"Potential\"\n"
        << std::setprecision(15);
    double x = normal(rng) / std::sqrt(1.0 - phi * phi);
    for (int i = 0; i < kFrames; ++i) {
      if (i > 0)
        x = phi * x + normal(rng);
      out << static_cast<double>(i) << ' ' << kE0 + kScale * x << '\n';
    }
  }
  std::ofstream("ener_" + tag + ".txt") << "e_" + tag + ".xvg\n";
  std::ofstream("temps.txt") << "300\n";
  arguments args;
  args.mbar.flag_mbar = true;
  args.mbar.ener_list_file = "ener_" + tag + ".txt";
  args.mbar.temps_file = "temps.txt";
  args.mbar.temp_mode = "fixed";
  args.mbar.target_temp = 300.0;
  args.mbar.window_match = "auto";
  args.mbar.bootstrap_samples = 400;
  TitanMBAR::FrameMetadata meta; // energy-only path
  TitanExecutor executor(2);
  std::fstream log("log_" + tag + ".txt", std::ios::out);
  const auto r =
      TitanMBAR::run_temperature_ladder_reweighting(args, meta, executor, log);
  TitanMBAR::write_summary_file("summary_" + tag + ".txt", r, args, log);
  const double se = analytic_se(phi);
  const double g = (1.0 + phi) / (1.0 - phi);
  const double ratio = r.target_std_energy / se;
  std::cout << std::setprecision(5) << std::fixed << "case phi=" << phi
            << " (g=" << g << "): bootstrap StdE = " << r.target_std_energy
            << " kJ/mol, analytic SE = " << se << ", ratio = " << ratio
            << " (i.i.d. ratio " << 1.0 / std::sqrt(g) << ")\n";
  check(r.solved && r.converged, tag + ": MBAR solved and converged");
  check(ratio > 0.80 && ratio < 1.20,
        tag + ": bootstrap StdE within [0.80, 1.20] of the AR(1) SE");
  std::ifstream in("summary_" + tag + ".txt");
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string summary = ss.str();
  check(summary.find("Bootstrap Dropped") != std::string::npos,
        tag + ": summary reports dropped (non-converged) replicates");
  check(summary.find("Bootstrap Block Length") != std::string::npos,
        tag + ": summary reports the block length per replica");
}

} // namespace

int main() {
  run_case("phi090", 0.9, 20261022u);
  run_case("phi050", 0.5, 777u);
  if (failures) {
    std::cout << "FAIL KIN-20 MBAR bootstrap: " << failures
              << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-20 MBAR bootstrap: block bootstrap error of the "
               "mean matches the AR(1) standard error\n";
  return 0;
}
