// KIN-22 oracle: -entropy (Duan-Zhang interaction entropy) on Gaussian
// interaction energies, through impl_entropy_duan_zhang (the -entropy output
// path: entropy_result.dat and entropy_result_time.dat).
//
// For Gaussian fluctuations of the interaction energy with standard deviation
// sigma the interaction entropy is exact:
//   -T dS = kT ln <exp(beta (E - <E>))> = sigma^2 / (2 kT).
// Sample: N = 20000 deterministic Gaussian quantiles z_i = Phi^-1((i-1/2)/N),
// E_elec = -60 + 5 z_i kJ/mol in a fixed shuffled order, E_vdW = -50 kJ/mol;
// T = 298.15 K (default -temperature), energies in kJ/mol.
//   analytic sigma^2/(2kT) = 5.0424 kJ/mol; the exponential average of this
//   finite quantile sample is 5.0239 (-0.37 %, the missing tail beyond
//   z = 4.06). Tolerance 1 % of the analytic value.
//   The 2-sigma trimmed estimate of the old default is 3.4056 (-32 %).
// The time series must hold the same estimator on every prefix [0, i]:
// checked at 12 prefixes against a direct long-double evaluation (5e-5
// relative: the file prints 6 significant digits), and its last row must
// equal the total. entropy_result.dat must state the trimming used.
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "arguments.h"
#include "folder.h"
#include "functions.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kN = 20000;
constexpr double kSigma = 5.0;
constexpr double kT = 298.15;
constexpr double kBoltz = 1.38064e-3 * 6.02214; // kJ/mol/K as in the engine

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

// Acklam's inverse normal CDF refined by one Halley step (|err| < 1e-13).
double inv_norm(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                             -2.759285104469687e+02, 1.383577518672690e+02,
                             -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                             -1.556989798598866e+02, 6.680131188771972e+01,
                             -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                             -2.400758277161838e+00, -2.549732539343734e+00,
                             4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                             2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p > 1 - 0.02425) {
    const double q = std::sqrt(-2 * std::log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double u = e * 2.5066282746310002 * std::exp(x * x / 2); // sqrt(2 pi)
  return x - u / (1 + x * u / 2);
}

// direct estimate on the prefix [0, n): kT ln mean exp(beta (E - mean))
double direct(const std::vector<double> &e, std::size_t n) {
  long double mean = 0.0L;
  for (std::size_t i = 0; i < n; ++i)
    mean += e[i];
  mean /= n;
  const long double rt = kBoltz * kT;
  long double mx = -1e300L;
  for (std::size_t i = 0; i < n; ++i)
    mx = std::max(mx, (e[i] - mean) / rt);
  long double s = 0.0L;
  for (std::size_t i = 0; i < n; ++i)
    s += std::exp((e[i] - mean) / rt - mx);
  return static_cast<double>(rt * (mx + std::log(s / n)));
}

std::string read_file(const std::string &path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

} // namespace

int main() {
  std::vector<double> z(kN);
  for (int i = 0; i < kN; ++i)
    z[static_cast<std::size_t>(i)] = inv_norm((i + 0.5) / kN);
  // fixed shuffle (LCG Fisher-Yates)
  unsigned x = 20261022u;
  for (int i = kN - 1; i > 0; --i) {
    x = 1664525u * x + 1013904223u;
    std::swap(z[static_cast<std::size_t>(i)], z[x % static_cast<unsigned>(i + 1)]);
  }
  std::vector<double> e_ele(kN), e_vdw(kN, -50.0), total(kN), time(kN);
  std::vector<std::vector<double>> e_ele_v(kN), e_vdw_v(kN);
  for (int i = 0; i < kN; ++i) {
    e_ele[static_cast<std::size_t>(i)] = -60.0 + kSigma * z[static_cast<std::size_t>(i)];
    total[static_cast<std::size_t>(i)] = e_ele[static_cast<std::size_t>(i)] - 50.0;
    time[static_cast<std::size_t>(i)] = 10.0 * i;
    e_ele_v[static_cast<std::size_t>(i)] = {e_ele[static_cast<std::size_t>(i)]};
    e_vdw_v[static_cast<std::size_t>(i)] = {-50.0};
  }

  arguments args;
  args.entropy = true;
  args.gbsa_temperature = kT;
  args.eu_par = 1.0;
  args.eu = "kJ/mol";
  args.tu = "ps";
  folders paths;
  paths.emplace_back(std::filesystem::path("."));
  const auto t0 = std::chrono::steady_clock::now();
  impl_entropy_duan_zhang(paths, e_ele_v, e_vdw_v, time, args, e_ele, e_vdw);
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  const double analytic = kSigma * kSigma / (2.0 * kBoltz * kT);
  const double sample = direct(total, kN);
  const std::string result = read_file("entropy_result.dat");
  double reported = NAN;
  {
    const std::string key = "Interaction entropy (-T*DeltaS)= ";
    const auto pos = result.find(key);
    if (pos != std::string::npos)
      reported = std::stod(result.substr(pos + key.size()));
  }
  std::cout << std::setprecision(6) << std::fixed
            << "case gaussian sigma=5 kJ/mol, T=298.15 K, N=20000: -TdS = "
            << reported << " kJ/mol; analytic " << analytic
            << ", quantile sample " << sample << "; impl_entropy_duan_zhang "
            << seconds << " s\n";
  check(std::abs(reported / analytic - 1.0) < 0.01,
        "-TdS = sigma^2/(2kT) within 1 % at default settings");
  check(result.find("filter_sigma") != std::string::npos,
        "entropy_result.dat states the sample trimming (filter_sigma)");

  // time series: prefix estimates
  std::ifstream ts("entropy_result_time.dat");
  std::string line;
  std::getline(ts, line); // header
  std::vector<double> series;
  while (std::getline(ts, line)) {
    std::istringstream row(line);
    double t, v;
    if (row >> t >> v)
      series.push_back(v);
  }
  check(series.size() == static_cast<std::size_t>(kN),
        "entropy_result_time.dat has one row per frame");
  if (series.size() == static_cast<std::size_t>(kN)) {
    double worst = 0.0;
    for (std::size_t n : {2u, 3u, 5u, 10u, 31u, 100u, 317u, 1000u, 3162u,
                          10000u, 15000u, static_cast<unsigned>(kN)}) {
      const double ref = direct(total, n);
      const double rel = std::abs(series[n - 1] - ref) / std::max(1e-3, std::abs(ref));
      worst = std::max(worst, rel);
      if (rel > 5e-5)
        std::cout << "  prefix " << n << ": file " << series[n - 1]
                  << ", direct " << ref << '\n';
    }
    std::cout << "  time series: worst relative deviation at 12 prefixes "
              << std::scientific << worst << std::fixed << '\n';
    check(worst <= 5e-5, "every checked prefix holds the untrimmed estimate");
    check(std::abs(series.back() - reported) <= 5e-5 * std::abs(reported) + 1e-5,
          "last row of the time series equals the total");
  }

  if (failures) {
    std::cout << "FAIL KIN-22 interaction entropy: " << failures
              << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-22 interaction entropy: untrimmed estimate equals "
               "sigma^2/(2kT), trimming reported, prefix series exact\n";
  return 0;
}
