// KIN-19 oracle: DHAM trajectory-break detection with a repeated first frame.
//
// Two replicas (300 K, 315 K) of 40 projected frames every 10 ps in two well
// separated 1-D clusters; the energy depends on the microstate only, so the
// state-mean energies are exact whatever frames are matched.
//   A  plain run
//   B  the same data, but replica 0 starts with its t = 0 frame twice (two
//      concatenated .xtc/.xvg parts that both contain t = 0), in the
//      trajectory and in the energy file
//   C  the repeated t = 0 frame in the trajectory only
// The repeated frame adds one transition of time gap 0, which is a break; all
// other transitions are regular (10 ps). Under TRAM (KIN-13) a transition
// needs a sample with an energy at both ends; one of the two t = 0 frames
// gets no energy row of its own (KIN-21), so B and C count the transitions
// of A less the one that starts at that frame (DHAM TRAM Transitions in the
// log), and the reweighted model may move by about one frame's weight
// (checked < 0.05). With dt taken from the first two timestamps
// (0 ps) every regular transition of replica 0 was declared a break (B, C);
// the repeated energy row also stopped the nearest-row scan of the energy
// matcher at t = 0, leaving 39 of 41 frames without energy (B).
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "info.h"
#include "postproc/engine_postproc.h"
#include "titan_executor.h"
#include "titan_msm.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kFrames = 40;
constexpr double kDt = 10.0;

std::vector<int> state_sequence(unsigned seed) {
  std::vector<int> s(kFrames);
  unsigned x = seed;
  int state = static_cast<int>(seed & 1u);
  for (int i = 0; i < kFrames; ++i) {
    x = 1103515245u * x + 12345u;
    if (((x >> 16) & 0x7fffu) % 10 < 3)
      state = 1 - state;
    s[static_cast<std::size_t>(i)] = state;
  }
  return s;
}

struct Run {
  Eigen::VectorXd pi;
  Eigen::MatrixXd tpm;
  bool converged{false};
  int label_plus{-1}; // microstate label of the +5 cluster
  long transitions{-1}; // "DHAM TRAM Transitions" in the log
};

long logged_transitions(const std::string &path) {
  std::ifstream in(path);
  std::string line;
  long value = -1;
  while (std::getline(in, line)) {
    if (line.find("DHAM TRAM Transitions") == std::string::npos)
      continue;
    const std::size_t colon = line.rfind(':');
    if (colon != std::string::npos)
      value = std::stol(line.substr(colon + 1));
  }
  return value;
}

Run run(const std::string &tag, bool duplicate_first,
        bool duplicate_energy_row = true) {
  const std::vector<int> s[2] = {state_sequence(7u), state_sequence(1234u)};
  std::ofstream list("ener_" + tag + ".txt");
  for (int r = 0; r < 2; ++r) {
    const std::string path = "e" + std::to_string(r) + "_" + tag + ".xvg";
    list << path << '\n';
    std::ofstream out(path);
    out << "# synthetic energy\n@    title \"energy\"\n" << std::setprecision(10);
    for (int i = 0; i < kFrames; ++i) {
      const double e = -100.0 + 4.0 * s[r][static_cast<std::size_t>(i)];
      out << kDt * i << ' ' << e << '\n';
      if (i == 0 && r == 0 && duplicate_first && duplicate_energy_row)
        out << kDt * i << ' ' << e << '\n';
    }
  }
  list.close();
  std::ofstream("temps.txt") << "300\n315\n";

  auto msm = std::make_unique<TitanMSM>(2, 1, 0.01);
  int64_t g = 0;
  for (int r = 0; r < 2; ++r) {
    for (int i = 0; i < kFrames; ++i) {
      const int copies = (i == 0 && r == 0 && duplicate_first) ? 2 : 1;
      for (int c = 0; c < copies; ++c) {
        Eigen::VectorXd x(1);
        x(0) = (s[r][static_cast<std::size_t>(i)] ? 5.0 : -5.0) +
               0.01 * std::sin(1.0 + i + 3.0 * r);
        msm->add_projected_frame(x, r, g++);
        msm->record_frame_time(kDt * i, r);
      }
    }
    msm->add_trajectory_break();
  }
  TitanExecutor executor(1);
  msm->train_kmeans(executor, 100, 42u);
  msm->build_markov_model(1);
  if (!msm->is_built())
    throw std::runtime_error(tag + ": synthetic MSM was not built");

  arguments args;
  args.dham.flag_dham = true;
  args.dham.temps_file = "temps.txt";
  args.dham.ener_list_file = "ener_" + tag + ".txt";
  args.dham.target_temp = 300.0;
  args.dham.window_match = "auto";
  info info_output(args);
  std::fstream log("log_" + tag + ".txt", std::ios::out);
  std::vector<double> weights;
  Run out;
  run_dham_postproc(args, msm, executor, info_output, log, weights, &out.tpm,
                    &out.pi, &out.converged);
  log.close();
  out.transitions = logged_transitions("log_" + tag + ".txt");
  // label of the +5 cluster: replica 1 frame of state 1
  const auto &micro = msm->discrete_traj();
  const std::size_t offset = duplicate_first ? kFrames + 1 : kFrames;
  for (int i = 0; i < kFrames; ++i)
    if (s[1][static_cast<std::size_t>(i)] == 1) {
      out.label_plus = micro[offset + static_cast<std::size_t>(i)];
      break;
    }
  return out;
}

} // namespace

void show(const char *name, const Run &r) {
  std::cout << "  " << name << ": converged " << r.converged << " pi = "
            << (r.pi.size() == 2 ? r.pi(0) : NAN) << ' '
            << (r.pi.size() == 2 ? r.pi(1) : NAN) << '\n';
}

// max |a - b| over pi and T, states matched through the +5 cluster label
double difference(const Run &a, const Run &b) {
  const int pa = a.label_plus, pb = b.label_plus;
  double worst = 0.0;
  for (int i = 0; i < 2; ++i) {
    const int ia = i == 0 ? pa : 1 - pa, ib = i == 0 ? pb : 1 - pb;
    worst = std::max(worst, std::abs(a.pi(ia) - b.pi(ib)));
    for (int j = 0; j < 2; ++j) {
      const int ja = j == 0 ? pa : 1 - pa, jb = j == 0 ? pb : 1 - pb;
      worst = std::max(worst, std::abs(a.tpm(ia, ja) - b.tpm(ib, jb)));
    }
  }
  return worst;
}

int main() {
  try {
    const Run a = run("A", false);
    // B: repeated t = 0 frame in the trajectory and in the energy file
    const Run b = run("B", true, true);
    // C: repeated t = 0 frame in the trajectory only
    const Run c = run("C", true, false);
    std::cout << std::setprecision(15);
    show("A", a);
    show("B", b);
    show("C", c);
    if (!a.converged || a.pi.size() != 2 || a.tpm.rows() != 2)
      throw std::runtime_error("reference run A did not converge");
    int bad = 0;
    for (const auto &[name, r] : {std::pair<const char *, const Run *>{"B", &b},
                                  std::pair<const char *, const Run *>{"C", &c}}) {
      if (!r->converged || r->pi.size() != 2 || r->tpm.rows() != 2) {
        std::cout << "  FAIL run " << name
                  << " did not converge: the count matrix lost its "
                     "transitions\n";
        ++bad;
        continue;
      }
      const double worst = difference(a, *r);
      std::cout << "  transitions A " << a.transitions << ", " << name << ' '
                << r->transitions << "; max |A - " << name
                << "| over pi and T = " << worst << '\n';
      if (r->transitions < a.transitions - 1 || r->transitions > a.transitions ||
          a.transitions <= 0) {
        std::cout << "  FAIL run " << name
                  << ": the repeated t = 0 frame changed the transition count\n";
        ++bad;
      }
      if (worst > 0.05) {
        std::cout << "  FAIL run " << name
                  << ": the reweighted model moved by more than one frame's "
                     "weight\n";
        ++bad;
      }
    }
    if (bad)
      throw std::runtime_error(std::to_string(bad) + " run(s) differ from A");
    std::cout << "PASS KIN-19 DHAM breaks: a repeated t = 0 frame is one "
                 "break; the transitions are those of the plain run less the "
                 "one that starts at the frame without an energy row\n";
    return 0;
  } catch (const std::exception &error) {
    std::cout << "FAIL KIN-19 " << error.what() << '\n';
    return 1;
  }
}
