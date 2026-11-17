// KIN-17 oracle: MBAR demux-map states resolved by time.
//
// Two T-REMD replicas (300 K, 330 K). Exchanges every 20 ps; the true state
// of replica 0 at time t is s0(ceil(t / 20)) (block pattern from an LCG),
// replica 1 holds the other state. Trajectory frames (TICA metadata) are at
// t = 100, 120, ..., 580 ps for replica 0 and 100..400 ps for replica 1, so
// N_k depends on the labels. The energy .xvg rows hold time, energy and the
// temperature of the state sampled at that time.
//
//   A  aligned: energy rows and a map without time column at the frame times
//   B  energy every 10 ps from t = 0, map without time column (row = energy row)
//   C  energy as B, map with a time column, rows every 20 ps (replica_index.xvg)
//   D  as B without frame metadata, -b 100 -e 580 -dt 20
//   E  as C, ladder inferred from the map and the temperature column
//   F  as C, map ending at t = 300 ps: must refuse
//   G  as B, map without time column but one row per FRAME: must refuse
//
// Every sample must carry the state of its own time and N_k the counts of the
// true states; A, B, C and E (same 41 samples) must give the same N_k, frame
// weights and <E>(300 K). Linked against the engine objects.
#include "arguments.h"
#include "titan_executor.h"
#include "titan_mbar.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double kTemps[2] = {300.0, 330.0};

int block_state0(long q) {
  unsigned x = 2463534242u ^ static_cast<unsigned>(q * 2654435761u);
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  return static_cast<int>((x >> 7) & 1u);
}

// true state of replica r at time t (t a multiple of 10 ps)
int true_state(int r, double t) {
  const long q = static_cast<long>(std::ceil(t / 20.0 - 1e-9));
  const int s0 = block_state0(q);
  return r == 0 ? s0 : 1 - s0;
}

double energy(int r, double t) {
  const int s = true_state(r, t);
  return -5000.0 + 12.0 * s + 6.0 * std::sin(0.37 * t + 1.3 * r) +
         4.0 * std::cos(1.1 * t + 0.2 * r);
}

void write_energy(const std::string &path, int r, double t0, double t1,
                  double dt) {
  std::ofstream out(path);
  out << "# synthetic T-REMD energies\n@    title \"Energies\"\n"
      << std::setprecision(12);
  for (double t = t0; t <= t1 + 1e-9; t += dt)
    out << t << ' ' << energy(r, t) << ' ' << kTemps[true_state(r, t)] << '\n';
}

// row with the replica at state 0 and at state 1
void write_map_row(std::ofstream &out, double t, bool time_column) {
  const int replica_at_state0 = true_state(0, t) == 0 ? 0 : 1;
  if (time_column)
    out << std::fixed << std::setprecision(2) << t << "  ";
  out << replica_at_state0 << ' ' << 1 - replica_at_state0 << '\n';
}

void write_map(const std::string &path, double t0, double t1, double dt,
               bool time_column) {
  std::ofstream out(path);
  out << "# replica index per state\n";
  for (double t = t0; t <= t1 + 1e-9; t += dt)
    write_map_row(out, t, time_column);
}

struct Case {
  std::string name;
  TitanMBAR::CoreResult result;
};

TitanMBAR::CoreResult run(const std::string &tag, const std::string &map,
                          bool metadata, bool temps_file, double b, double e,
                          int dt) {
  arguments args;
  args.mbar.flag_mbar = true;
  args.mbar.ener_list_file = "ener_" + tag + ".txt";
  args.mbar.temps_file = temps_file ? "temps.txt" : "";
  args.mbar.temp_mode = "per_frame";
  args.mbar.demux_map_file = map;
  args.mbar.target_temp = 300.0;
  args.mbar.window_match = "auto";
  if (b >= 0.0) {
    args.flag_begin_set = true;
    args.begin_time = b;
  }
  if (e >= 0.0) {
    args.flag_end_set = true;
    args.end_time = e;
  }
  args.dt = dt;
  TitanMBAR::FrameMetadata meta;
  if (metadata) {
    // replica 1's trajectory is shorter (100..400 ps), so N_k depends on the
    // labels and a wrong label changes N_k and f_k
    long g = 0;
    for (int r = 0; r < 2; ++r)
      for (int i = 0; i < (r == 0 ? 25 : 16); ++i) {
        meta.replica_ids.push_back(r);
        meta.timestamps_ps.push_back(100.0 + 20.0 * i);
        meta.global_indices.push_back(g++);
      }
  }
  TitanExecutor executor(2);
  std::fstream log("log_" + tag + ".txt", std::ios::out);
  return TitanMBAR::run_temperature_ladder_reweighting(args, meta, executor,
                                                       log);
}

int failures = 0;
void check(bool ok, const std::string &what) {
  if (!ok) {
    ++failures;
    std::cout << "  FAIL " << what << '\n';
  }
}

// ref == nullptr: no reference with the same sample set (D)
void compare(const std::string &name, const TitanMBAR::CoreResult &r,
             const TitanMBAR::CoreResult *ref, std::size_t expected_samples) {
  check(r.solved, name + ": MBAR solved (abort reason: " +
                      (r.abort_reason.empty() ? "none" : r.abort_reason) + ")");
  if (!r.solved)
    return;
  int wrong_state = 0;
  std::vector<int> true_counts(2, 0);
  std::map<std::pair<int, long>, double> weights;
  for (const auto &row : r.frame_weights) {
    const int truth = true_state(row.replica_id, row.time_ps);
    ++true_counts[static_cast<std::size_t>(truth)];
    if (row.sampled_state_index != truth)
      ++wrong_state;
    weights[{row.replica_id, std::lround(row.time_ps * 1000.0)}] = row.weight;
  }
  check(r.frame_weights.size() == expected_samples,
        name + ": " + std::to_string(r.frame_weights.size()) + " samples, " +
            std::to_string(expected_samples) + " expected");
  check(wrong_state == 0, name + ": " + std::to_string(wrong_state) + "/" +
                              std::to_string(r.frame_weights.size()) +
                              " samples carry the state of another time");
  check(r.sampled_counts == true_counts,
        name + ": N_k equal to the counts of the true states (" +
            std::to_string(true_counts[0]) + "/" +
            std::to_string(true_counts[1]) + ")");
  std::cout << std::setprecision(12) << "  " << name << ": N_k = "
            << (r.sampled_counts.size() > 1 ? r.sampled_counts[0] : -1) << '/'
            << (r.sampled_counts.size() > 1 ? r.sampled_counts[1] : -1)
            << ", samples " << r.frame_weights.size() << ", <E>(300 K) = "
            << r.target_mean_energy << ", wrong states " << wrong_state
            << '\n';
  if (!ref)
    return;
  check(r.sampled_counts == ref->sampled_counts,
        name + ": N_k equal to the aligned case");
  double max_dw = 0.0;
  for (const auto &row : ref->frame_weights) {
    const auto it =
        weights.find({row.replica_id, std::lround(row.time_ps * 1000.0)});
    if (it == weights.end()) {
      max_dw = INFINITY;
      break;
    }
    max_dw = std::max(max_dw, std::abs(it->second - row.weight) /
                                  std::max(1e-300, std::abs(row.weight)));
  }
  check(weights.size() == ref->frame_weights.size() && max_dw <= 1e-10,
        name + ": frame weights equal to the aligned case (max rel diff " +
            std::to_string(max_dw) + ")");
  const double de = std::abs(r.target_mean_energy - ref->target_mean_energy);
  check(de <= 1e-10 * std::abs(ref->target_mean_energy),
        name + ": <E>(300 K) equal to the aligned case");
}

} // namespace

int main() {
  try {
    std::ofstream("temps.txt") << "300\n330\n";
    for (int r = 0; r < 2; ++r) {
      write_energy("a" + std::to_string(r) + ".xvg", r, 100.0, 580.0, 20.0);
      write_energy("b" + std::to_string(r) + ".xvg", r, 0.0, 600.0, 10.0);
    }
    std::ofstream("ener_a.txt") << "a0.xvg\na1.xvg\n";
    for (const char *tag : {"b", "c", "d", "e", "f", "g"})
      std::ofstream(std::string("ener_") + tag + ".txt") << "b0.xvg\nb1.xvg\n";
    write_map("map_a.txt", 100.0, 580.0, 20.0, false);
    write_map("map_b.txt", 0.0, 600.0, 10.0, false);
    write_map("map_c.xvg", 0.0, 600.0, 20.0, true);
    write_map("map_f.xvg", 0.0, 300.0, 20.0, true);
    write_map("map_g.txt", 100.0, 580.0, 20.0, false);

    // frames: replica 0 at 100..580 ps, replica 1 at 100..400 ps (41)
    const auto a = run("a", "map_a.txt", true, true, -1, -1, 0);
    // independent check of the aligned case: states of their own time
    compare("A aligned", a, nullptr, 41);
    if (!a.solved)
      throw std::runtime_error("aligned case did not solve");
    compare("B energy every 10 ps from t=0, map = energy rows",
            run("b", "map_b.txt", true, true, -1, -1, 0), &a, 41);
    compare("C map with time column every 20 ps",
            run("c", "map_c.xvg", true, true, -1, -1, 0), &a, 41);
    // no metadata: every energy row of both replicas in 100..580 ps, step 20
    compare("D no frame metadata, -b 100 -e 580 -dt 20",
            run("d", "map_b.txt", false, true, 100.0, 580.0, 20), nullptr, 50);
    compare("E ladder inferred from the map",
            run("e", "map_c.xvg", true, false, -1, -1, 0), &a, 41);
    const auto f = run("f", "map_f.xvg", true, true, -1, -1, 0);
    check(!f.solved && f.abort_reason.find("demux") != std::string::npos,
          "F map ending at 300 ps refused (abort reason: " +
              (f.abort_reason.empty() ? std::string("none") : f.abort_reason) +
              ")");
    const auto g = run("g", "map_g.txt", true, true, -1, -1, 0);
    check(!g.solved && g.abort_reason.find("demux") != std::string::npos,
          "G map without time column and one row per frame refused (abort "
          "reason: " +
              (g.abort_reason.empty() ? std::string("none") : g.abort_reason) +
              ")");
    if (failures) {
      std::cout << "FAIL KIN-17 " << failures << " check(s) failed\n";
      return 1;
    }
    std::cout << "PASS KIN-17 MBAR demux map: states by time for nstenergy != "
                 "xtc stride, t0 != 0, -b/-e/-dt, time-column map, inferred "
                 "ladder; uncovered and ambiguous maps refused\n";
    return 0;
  } catch (const std::exception &error) {
    std::cout << "FAIL KIN-17 " << error.what() << '\n';
    return 1;
  }
}
