// KIN-16 oracle: DHAM frame weights after window matching.
//
// Two replicas (300 K, 315 K) of 40 projected frames each, every 10 ps from
// t = 0, in two well separated 1-D clusters. The energy .xvg of replica 0
// starts at t = 100 ps (10 frames trimmed at the start); that of replica 1
// covers 50..330 ps (5 frames trimmed at the start, 6 at the end). After the
// real run_dham_postproc, every frame inside its energy window must carry
// its OWN TRAM sample weight and every frame outside must carry weight 0;
// the weights sum to 1. The expected weights come from a TitanTram built
// here from the same frames (energies read back from the .xvg files, lag-1
// transitions between consecutive kept frames): a weight shifted by the
// trimmed start of a window lands on another frame with another energy.
// (Before KIN-13 the expected weight was p_m / N_m of the microstate.)
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "info.h"
#include "postproc/engine_postproc.h"
#include "titan_executor.h"
#include "titan_msm.h"
#include "titan_tram.h"

#include <sstream>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kFrames = 40;
constexpr double kDt = 10.0;

// Two-state sequence with persistence (deterministic LCG).
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

void write_energy(const std::string &path, const std::vector<int> &states,
                  int replica, double t_begin, double t_end) {
  std::ofstream out(path);
  out << "# synthetic energy\n@    title \"energy\"\n";
  out << std::setprecision(10);
  for (int i = 0; i < kFrames; ++i) {
    const double t = kDt * i;
    if (t < t_begin - 1e-9 || t > t_end + 1e-9)
      continue;
    const double e = -100.0 + 4.0 * states[static_cast<std::size_t>(i)] +
                     0.3 * std::sin(0.1 * t + replica);
    out << t << ' ' << e << '\n';
  }
}

} // namespace

int main() {
  try {
    const std::vector<int> s0 = state_sequence(7u), s1 = state_sequence(1234u);
    const double window[2][2] = {{100.0, 390.0}, {50.0, 330.0}};
    write_energy("e0.xvg", s0, 0, window[0][0], window[0][1]);
    write_energy("e1.xvg", s1, 1, window[1][0], window[1][1]);
    std::ofstream("temps.txt") << "300\n315\n";
    std::ofstream("ener_list.txt") << "e0.xvg\ne1.xvg\n";

    auto msm = std::make_unique<TitanMSM>(2, 1, 0.01);
    for (int r = 0; r < 2; ++r) {
      const auto &s = r == 0 ? s0 : s1;
      for (int i = 0; i < kFrames; ++i) {
        Eigen::VectorXd x(1);
        x(0) = (s[static_cast<std::size_t>(i)] ? 5.0 : -5.0) +
               0.01 * std::sin(1.0 + i + 3.0 * r);
        msm->add_projected_frame(x, r, static_cast<int64_t>(r * kFrames + i));
        msm->record_frame_time(kDt * i, r);
      }
      msm->add_trajectory_break();
    }
    TitanExecutor executor(1);
    msm->train_kmeans(executor, 100, 42u);
    msm->build_markov_model(1);
    if (!msm->is_built())
      throw std::runtime_error("synthetic MSM was not built");

    arguments args;
    args.dham.flag_dham = true;
    args.dham.temps_file = "temps.txt";
    args.dham.ener_list_file = "ener_list.txt";
    args.dham.target_temp = 300.0;
    args.dham.window_match = "auto";
    info info_output(args);
    std::fstream log("log.txt", std::ios::out);
    std::vector<double> weights;
    Eigen::VectorXd pi;
    bool converged = false;
    run_dham_postproc(args, msm, executor, info_output, log, weights, nullptr,
                      &pi, &converged);
    if (!converged || pi.size() != 2)
      throw std::runtime_error("DHAM did not converge on the synthetic case");

    const auto &micro = msm->discrete_traj();
    auto kept = [&](int r, int i) {
      const double t = kDt * i;
      return t >= window[r][0] - 1e-9 && t <= window[r][1] + 1e-9;
    };
    // energies as written (10 significant digits), read back by time
    auto read_energy = [](const std::string &path) {
      std::vector<std::pair<double, double>> rows;
      std::ifstream in(path);
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '@')
          continue;
        std::istringstream fields(line);
        double t = 0.0, e = 0.0;
        if (fields >> t >> e)
          rows.push_back({t, e});
      }
      return rows;
    };
    const std::vector<std::pair<double, double>> energy[2] = {
        read_energy("e0.xvg"), read_energy("e1.xvg")};
    auto energy_at = [&](int r, int i) {
      for (const auto &[t, e] : energy[r])
        if (std::abs(t - kDt * i) < 1e-9)
          return e;
      throw std::runtime_error("no energy row for a kept frame");
    };
    TitanTram reference(2, {300.0, 315.0}, 300.0);
    std::vector<std::size_t> sample_frame;
    for (int r = 0; r < 2; ++r) {
      for (int i = 0; i < kFrames; ++i) {
        if (!kept(r, i))
          continue;
        const std::size_t g = static_cast<std::size_t>(r * kFrames + i);
        reference.add_sample(r, micro[g], energy_at(r, i));
        sample_frame.push_back(g);
        if (i + 1 < kFrames && kept(r, i + 1))
          reference.add_transition(r, micro[g], micro[g + 1]);
      }
    }
    if (!reference.solve(args.dham.max_iter, args.dham.tolerance).converged)
      throw std::runtime_error("reference TRAM did not converge");
    std::vector<double> expected_weight(static_cast<std::size_t>(2 * kFrames), 0.0);
    const std::vector<double> w = reference.sample_weights();
    for (std::size_t x = 0; x < w.size(); ++x)
      expected_weight[sample_frame[x]] = w[x];

    if (weights.size() != static_cast<std::size_t>(2 * kFrames))
      throw std::runtime_error("weight vector has " +
                               std::to_string(weights.size()) + " entries");
    int wrong = 0, wrong_outside = 0;
    double sum = 0.0;
    for (int r = 0; r < 2; ++r) {
      for (int i = 0; i < kFrames; ++i) {
        const std::size_t g = static_cast<std::size_t>(r * kFrames + i);
        const int m = micro[g];
        const double expected = expected_weight[g];
        sum += weights[g];
        if (std::abs(weights[g] - expected) > 1e-10 * std::max(1e-3, expected)) {
          ++wrong;
          if (!kept(r, i))
            ++wrong_outside;
          if (wrong <= 6)
            std::cout << "  frame replica " << r << " t=" << kDt * i
                      << " ps state " << m << ": weight " << weights[g]
                      << ", expected " << expected << '\n';
        }
      }
    }
    std::cout << std::setprecision(15) << "pi = " << pi(0) << ' ' << pi(1)
              << ", reference pi = " << reference.target_stationary()(0) << ' '
              << reference.target_stationary()(1) << ", weight sum = " << sum
              << '\n';
    if (wrong > 0)
      throw std::runtime_error(
          std::to_string(wrong) + "/" + std::to_string(2 * kFrames) +
          " frames carry a weight other than their own TRAM sample weight (" +
          std::to_string(wrong_outside) +
          " of them outside the energy window)");
    if (std::abs(sum - 1.0) > 1e-12)
      throw std::runtime_error("weights do not sum to 1");
    std::cout << "PASS KIN-16 DHAM weights: trimmed start/end windows, each "
                 "frame its own TRAM sample weight, 0 outside the window\n";
    return 0;
  } catch (const std::exception &error) {
    std::cout << "FAIL KIN-16 " << error.what() << '\n';
    return 1;
  }
}
