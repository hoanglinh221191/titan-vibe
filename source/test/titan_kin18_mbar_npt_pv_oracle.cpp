// KIN-18 oracle: MBAR in the NPT ensemble adds p*V to the reduced potential.
//
// Two replicas at 300 K and 320 K, both at p = 1 bar, 20000 frames each.
// Configurational energy U and volume V are independent with densities of
// states U^(a-1) and V^(b-1), so at temperature T and pressure p
//   U ~ Gamma(a, kT),   p V c ~ Gamma(b, kT)   (c = 0.0602214 kJ/mol/(bar nm^3))
// and the enthalpy H = U + p V c ~ Gamma(a + b, kT):
//   <H>(T) = U0 + (a + b) k T,   Cp = (a + b) k.
// The .xvg column holds U0 + U (U0 = -5000 kJ/mol), the trajectory metadata
// the box volume of each frame. MBAR at 310 K must return <H> and Cp of the
// enthalpy (a = 60, b = 40: 257.75 kJ/mol above U0 and 0.8314 kJ/mol/K).
// Without the p V term it returns <U> and Cv of U (154.65, 0.4989): 40 % low.
// Tolerance: <H> - U0 within 1.0 kJ/mol (the standard error of the
// reweighted mean is about 0.2 kJ/mol: SD(H) = 10 kT = 25.8 kJ/mol over
// 4e4 frames), Cp within 5 % (relative SE of a variance ~ sqrt(2/N) = 1 %).
//
// Controls: the same data with -mbar-pressure-bar 0 and with ensemble nvt
// must NOT add p V (result = U statistics) and the summary must say so.
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
#include <tuple>
#include <vector>

namespace {

constexpr double kB = 0.008314462618;     // kJ/mol/K
constexpr double kPV = 0.06022140857;     // kJ/mol per bar nm^3
constexpr double kU0 = -5000.0;
constexpr double kA = 60.0, kBdof = 40.0; // shape of U and of pV
constexpr int kFrames = 20000;
const double kTemps[2] = {300.0, 320.0};
constexpr double kTarget = 310.0;

struct Replica {
  std::vector<double> u, v; // U (kJ/mol, without U0), V (nm^3)
};

Replica sample(int r, double pressure_bar) {
  std::mt19937_64 rng(1000u + static_cast<unsigned>(r));
  const double kT = kB * kTemps[r];
  std::gamma_distribution<double> gu(kA, kT), gv(kBdof, kT);
  Replica out;
  for (int i = 0; i < kFrames; ++i) {
    out.u.push_back(gu(rng));
    out.v.push_back(gv(rng) / (pressure_bar * kPV));
  }
  return out;
}

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

TitanMBAR::CoreResult run(const std::string &tag, const Replica *reps,
                          const std::string &ensemble, double pressure) {
  std::ofstream list("ener_" + tag + ".txt");
  for (int r = 0; r < 2; ++r) {
    const std::string path = "e" + std::to_string(r) + "_" + tag + ".xvg";
    list << path << '\n';
    std::ofstream out(path);
    out << "# synthetic NPT energies\n@    title \"Potential\"\n"
        << std::setprecision(15);
    for (int i = 0; i < kFrames; ++i)
      out << static_cast<double>(i) << ' ' << kU0 + reps[r].u[i] << '\n';
  }
  list.close();
  std::ofstream("temps.txt") << "300\n320\n";
  arguments args;
  args.mbar.flag_mbar = true;
  args.mbar.ener_list_file = "ener_" + tag + ".txt";
  args.mbar.temps_file = "temps.txt";
  args.mbar.temp_mode = "fixed";
  args.mbar.ensemble = ensemble;
  args.mbar.target_pressure_bar = pressure;
  args.mbar.target_temp = kTarget;
  args.mbar.window_match = "auto";
  TitanMBAR::FrameMetadata meta;
  long g = 0;
  for (int r = 0; r < 2; ++r)
    for (int i = 0; i < kFrames; ++i) {
      meta.replica_ids.push_back(r);
      meta.timestamps_ps.push_back(static_cast<double>(i));
      meta.global_indices.push_back(g++);
      meta.box_volumes_nm3.push_back(reps[r].v[i]);
    }
  TitanExecutor executor(2);
  std::fstream log("log_" + tag + ".txt", std::ios::out);
  auto result =
      TitanMBAR::run_temperature_ladder_reweighting(args, meta, executor, log);
  TitanMBAR::write_summary_file("summary_" + tag + ".txt", result, args, log);
  return result;
}

std::string read_file(const std::string &path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

} // namespace

int main() {
  const double pressure = 1.0;
  const Replica reps[2] = {sample(0, pressure), sample(1, pressure)};
  const double h_expected = (kA + kBdof) * kB * kTarget;
  const double cp_expected = (kA + kBdof) * kB;
  const double u_expected = kA * kB * kTarget;
  const double cv_u_expected = kA * kB;
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << "case npt p=1 bar: expected <H>-U0 = " << h_expected
            << " kJ/mol, Cp = " << cp_expected << " kJ/mol/K (U only: "
            << u_expected << ", " << cv_u_expected << ")\n";

  const auto npt = run("npt", reps, "npt", pressure);
  check(npt.solved && npt.converged, "npt: MBAR solved and converged");
  const double h = npt.target_mean_energy - kU0;
  std::cout << "  npt: <H>-U0 = " << h << ", Cp = " << npt.target_heat_capacity
            << '\n';
  check(std::abs(h - h_expected) < 1.0, "npt: <H> = U0 + (a+b) k T within 1 kJ/mol");
  check(std::abs(npt.target_heat_capacity / cp_expected - 1.0) < 0.05,
        "npt: Cp = (a+b) k within 5 %");
  const std::string s_npt = read_file("summary_npt.txt");
  check(s_npt.find("Box Volume Metadata") == std::string::npos,
        "npt summary: no 'Box Volume Metadata: available' line");
  check(s_npt.find("pV Term") != std::string::npos &&
            s_npt.find("included (p = 1") != std::string::npos,
        "npt summary: states that p V is included at 1 bar");

  // Controls: p V must not be added without a pressure or in NVT.
  for (const auto &[tag, ensemble, p] :
       {std::tuple<std::string, std::string, double>{"p0", "npt", 0.0},
        std::tuple<std::string, std::string, double>{"nvt", "nvt", 1.0}}) {
    const auto r = run(tag, reps, ensemble, p);
    const double e = r.target_mean_energy - kU0;
    std::cout << "  " << tag << ": <E>-U0 = " << e << ", C = "
              << r.target_heat_capacity << '\n';
    check(r.solved && r.converged, tag + ": MBAR solved and converged");
    check(std::abs(e - u_expected) < 1.0,
          tag + ": no p V term, <U> = a k T within 1 kJ/mol");
    check(std::abs(r.target_heat_capacity / cv_u_expected - 1.0) < 0.05,
          tag + ": C = a k within 5 %");
    const std::string s = read_file("summary_" + tag + ".txt");
    check(s.find("not included") != std::string::npos,
          tag + " summary: states that p V is not included");
  }

  if (failures) {
    std::cout << "FAIL KIN-18 MBAR NPT: " << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-18 MBAR NPT: reduced potential beta (U + p V), "
               "<H> and Cp of the analytic Gamma ensemble; no p V without "
               "pressure or in NVT\n";
  return 0;
}
