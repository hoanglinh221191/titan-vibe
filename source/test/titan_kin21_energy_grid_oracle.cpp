// KIN-21 oracle: energy grid coarser than the trajectory grid.
//
// Two replicas (300 K, 310 K), trajectory frames (TICA metadata) every 1 ps
// over 0..1999 ps, energy .xvg rows every 2 ps over 0..1998 ps. Each energy
// row is one sample: MBAR must use 1000 samples per replica (the frame nearest
// to each energy row), not 2000 (every odd frame used to copy the energy of
// its neighbour row at the 0.5 dt_energy tolerance), and the window audit must
// flag the replica as energy_coarser_than_traj.
// Control: energy every 0.5 ps (finer than the trajectory): one sample per
// frame (2000 per replica), no flag, unchanged.
// The matcher (TitanDHAM::window_match_energy_to_dtraj, shared by MBAR and
// DHAM) is also checked directly: 1000 matched frames, 1000 left without an
// energy value.
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "arguments.h"
#include "titan_dham.h"
#include "titan_executor.h"
#include "titan_mbar.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr int kFrames = 2000;

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

double energy(int r, double t) {
  return -5000.0 + 15.0 * r + 6.0 * std::sin(0.37 * t + 1.3 * r) +
         4.0 * std::cos(1.1 * t + 0.2 * r);
}

TitanMBAR::CoreResult run(const std::string &tag, double dt_energy) {
  std::ofstream list("ener_" + tag + ".txt");
  for (int r = 0; r < 2; ++r) {
    const std::string path = "e" + std::to_string(r) + "_" + tag + ".xvg";
    list << path << '\n';
    std::ofstream out(path);
    out << "# synthetic energies\n@    title \"Potential\"\n"
        << std::setprecision(12);
    for (double t = 0.0; t <= kFrames - 1 + 1e-9; t += dt_energy)
      out << t << ' ' << energy(r, t) << '\n';
  }
  list.close();
  std::ofstream("temps.txt") << "300\n310\n";
  arguments args;
  args.mbar.flag_mbar = true;
  args.mbar.ener_list_file = "ener_" + tag + ".txt";
  args.mbar.temps_file = "temps.txt";
  args.mbar.temp_mode = "fixed";
  args.mbar.target_temp = 300.0;
  args.mbar.window_match = "auto";
  TitanMBAR::FrameMetadata meta;
  long g = 0;
  for (int r = 0; r < 2; ++r)
    for (int i = 0; i < kFrames; ++i) {
      meta.replica_ids.push_back(r);
      meta.timestamps_ps.push_back(static_cast<double>(i));
      meta.global_indices.push_back(g++);
    }
  TitanExecutor executor(2);
  std::fstream log("log_" + tag + ".txt", std::ios::out);
  auto result =
      TitanMBAR::run_temperature_ladder_reweighting(args, meta, executor, log);
  TitanMBAR::write_window_audit_file("audit_" + tag + ".dat", result.window_audit);
  return result;
}

} // namespace

int main() {
  // the shared matcher on its own
  std::vector<double> traj_t;
  for (int i = 0; i < kFrames; ++i)
    traj_t.push_back(static_cast<double>(i));
  std::vector<std::pair<double, double>> coarse;
  for (int i = 0; i < kFrames; i += 2)
    coarse.emplace_back(static_cast<double>(i), energy(0, i));
  const auto m = TitanDHAM::window_match_energy_to_dtraj(traj_t, coarse, 0.0);
  std::size_t finite = 0;
  for (double v : m.aligned_values)
    finite += std::isfinite(v) ? 1u : 0u;
  std::cout << "case matcher: " << m.kept_frames << " kept, " << finite
            << " with an energy value\n";
  check(finite == kFrames / 2, "matcher: one frame per energy row (1000)");

  const auto coarse_run = run("coarse", 2.0);
  std::cout << "case coarse (energy 2 ps, frames 1 ps): samples "
            << coarse_run.total_usable_samples << ", N_k";
  for (int n : coarse_run.sampled_counts)
    std::cout << ' ' << n;
  std::cout << '\n';
  check(coarse_run.solved, "coarse: MBAR solved");
  check(coarse_run.total_usable_samples == kFrames,
        "coarse: 2000 usable samples (one per energy row)");
  check(coarse_run.sampled_counts.size() == 2 &&
            coarse_run.sampled_counts[0] == kFrames / 2 &&
            coarse_run.sampled_counts[1] == kFrames / 2,
        "coarse: N_k = 1000 per state");
#ifdef TITAN_MBAR_AUDIT_ENERGY_COARSER_FLAG
  bool flagged = coarse_run.window_audit.size() == 2;
  for (const auto &row : coarse_run.window_audit)
    flagged = flagged && row.energy_coarser_than_traj;
  check(flagged, "coarse: window audit flags energy_coarser_than_traj");
#else
  check(false, "coarse: window audit has no energy_coarser_than_traj flag");
#endif

  const auto fine_run = run("fine", 0.5);
  std::cout << "case fine (energy 0.5 ps, frames 1 ps): samples "
            << fine_run.total_usable_samples << '\n';
  check(fine_run.solved, "fine: MBAR solved");
  check(fine_run.total_usable_samples == 2 * kFrames,
        "fine: 4000 usable samples (one per frame)");
#ifdef TITAN_MBAR_AUDIT_ENERGY_COARSER_FLAG
  bool unflagged = true;
  for (const auto &row : fine_run.window_audit)
    unflagged = unflagged && !row.energy_coarser_than_traj;
  check(unflagged, "fine: no energy_coarser_than_traj flag");
#endif

  if (failures) {
    std::cout << "FAIL KIN-21 energy grid: " << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-21 energy grid: one frame per energy row, N_k halved "
               "and flagged when the energy grid is coarser\n";
  return 0;
}
