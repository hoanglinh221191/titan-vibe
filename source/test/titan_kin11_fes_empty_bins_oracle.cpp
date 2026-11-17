// KIN-11 oracle: empty bins of the 2-D FES (pca/tica/mbar/umap *_fes_2d.dat).
//
// TitanPCA::export_projection_fes_analysis on 42 synthetic projected frames,
// a 10 x 10 grid over [0, 10]^2 (bin width 1):
//   (0, 0) and (10, 10)   one frame each (corner bins, they fix the range)
//   (2.5, 2.5)            30 frames: basin A, bin (2, 2), F = 0
//   (6.5, 6.5)            10 frames: basin B, bin (6, 6), F = kT ln 3,
//                         every neighbour bin empty
// At T = 3000 K (kT = 24.94 kJ/mol), F_B = 27.40 kJ/mol > 20 kJ/mol. Checks:
//   - the 96 empty bins are written as "nan" (old: 20.0000) and hold NaN;
//   - the 4 sampled bins hold -kT ln(P / P_max), scaled by 1 / eu_par;
//   - B is a basin: an empty neighbour is not lower (old: the 20 kJ/mol fill
//     beat F_B, so only basin A was found);
//   - the same with eu_par = 4.184 (kcal/mol output): empty bins stay "nan".
// Linked against the engine objects (test/engine_linked_oracles.py).
#include "arguments.h"
#include "titan_executor.h"
#include "titan_pca.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

struct Row {
  double x, y;
  std::string f;
};

std::vector<Row> read_fes(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    Row r{};
    if (iss >> r.x >> r.y >> r.f)
      rows.push_back(r);
  }
  return rows;
}

void run_case(double eu_par, const std::string &tag) {
  arguments args;
  args.fes_grid = 10;
  args.eu_par = eu_par;
  const double temp = 3000.0;
  const double kT = 0.008314463 * temp;

  std::vector<std::pair<double, double>> proj;
  proj.push_back({0.0, 0.0});
  proj.push_back({10.0, 10.0});
  for (int i = 0; i < 30; ++i)
    proj.push_back({2.5, 2.5});
  for (int i = 0; i < 10; ++i)
    proj.push_back({6.5, 6.5});
  AlignedDoubleVector weights(proj.size(), 1.0);

  TitanExecutor executor(2);
  std::fstream log("kin11_" + tag + ".log", std::ios::out | std::ios::trunc);
  const std::string fes_path = "kin11_" + tag + "_fes_2d.dat";
  const std::string rep_path = "kin11_" + tag + "_representatives.dat";
  TitanPCA::FES2DResult fes;
  std::vector<TitanPCA::ClusterNode> clusters;
  const bool ok = TitanPCA::export_projection_fes_analysis(
      "KIN11", "PCA", args, proj, weights, false, temp, "uniform", executor,
      log, fes_path, rep_path, {}, {}, &fes, &clusters);
  check(ok, tag + ": export_projection_fes_analysis returned true");
  if (!ok)
    return;

  const double scale = 1.0 / eu_par;
  const double f_b = kT * std::log(3.0) * scale;
  const double f_corner = kT * std::log(30.0) * scale;
  int empty_nan = 0, empty_other = 0, sampled_ok = 0;
  for (int ix = 0; ix < fes.grid_size; ++ix)
    for (int iy = 0; iy < fes.grid_size; ++iy) {
      const auto &bin = fes.bin(ix, iy);
      if (bin.count_w <= 0.0) {
        (std::isnan(bin.free_energy) ? empty_nan : empty_other)++;
        continue;
      }
      double want = f_corner;
      if (ix == 2 && iy == 2)
        want = 0.0;
      else if (ix == 6 && iy == 6)
        want = f_b;
      if (std::abs(bin.free_energy - want) <= 1e-9 * (1.0 + want))
        ++sampled_ok;
      else
        std::cout << "    bin (" << ix << "," << iy << ") F " << bin.free_energy
                  << " expected " << want << '\n';
    }
  check(empty_nan == 96 && empty_other == 0,
        tag + ": 96 empty bins hold NaN in memory (" + std::to_string(empty_nan) +
            " NaN, " + std::to_string(empty_other) + " other)");
  check(sampled_ok == 4, tag + ": 4 sampled bins hold -kT ln(P/Pmax)/eu_par");

  const auto rows = read_fes(fes_path);
  int nan_rows = 0, finite_rows = 0, twenty_rows = 0;
  for (const auto &r : rows) {
    if (r.f == "nan")
      ++nan_rows;
    else {
      ++finite_rows;
      if (r.f == "20.0000")
        ++twenty_rows;
    }
  }
  check(rows.size() == 100 && nan_rows == 96 && finite_rows == 4,
        tag + ": file has 96 \"nan\" rows and 4 numbers (" +
            std::to_string(nan_rows) + " nan, " + std::to_string(finite_rows) +
            " numbers, " + std::to_string(twenty_rows) + " of them 20.0000)");

  bool found_b = false, found_a = false;
  for (const auto &c : clusters) {
    if (std::abs(c.min_pc1 - 6.5) < 1e-9 && std::abs(c.min_pc2 - 6.5) < 1e-9 &&
        std::abs(c.min_energy - f_b) <= 1e-9 * (1.0 + f_b))
      found_b = true;
    if (std::abs(c.min_pc1 - 2.5) < 1e-9 && std::abs(c.min_pc2 - 2.5) < 1e-9)
      found_a = true;
  }
  std::cout << "    " << clusters.size() << " basin(s); F_B = " << f_b << '\n';
  check(found_a, tag + ": basin A at (2.5, 2.5)");
  check(found_b, tag + ": basin B at (6.5, 6.5), all neighbours empty (F_B = " +
                     std::to_string(f_b) + " in the output unit)");
}

} // namespace

int main() {
  std::cout << "case kJ (eu_par 1)\n";
  run_case(1.0, "kj");
  std::cout << "case kcal (eu_par 4.184)\n";
  run_case(4.184, "kcal");
  if (failures) {
    std::cout << "FAIL KIN-11 FES empty bins: " << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-11 FES empty bins: nan in file and memory, basins next "
               "to empty bins found, any energy unit\n";
  return 0;
}
