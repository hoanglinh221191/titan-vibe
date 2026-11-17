// Per-atom Shrake-Rupley SASA of one structure at several sphere-point
// counts, with the time per call (MDB-15). Driven by
// test/sasa_point_convergence.py.
//
//   titan_sasa_point_convergence atoms.txt out.csv 96 384 960 ...
//
// atoms.txt: one atom per line, "x y z radius" in Angstrom. out.csv: one row
// per atom, one SASA column per point count. The timing (best of three calls,
// one thread) goes to stdout as "points,ms".
#include "titan_sasa.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: " << argv[0] << " atoms.txt out.csv points...\n";
    return EXIT_FAILURE;
  }
  std::ifstream in(argv[1]);
  std::vector<TitanSASA::Atom> atoms;
  TitanSASA::Atom atom;
  while (in >> atom.x >> atom.y >> atom.z >> atom.radius)
    atoms.push_back(atom);
  if (atoms.empty()) {
    std::cerr << "no atoms read\n";
    return EXIT_FAILURE;
  }
  std::vector<int> counts;
  for (int i = 3; i < argc; ++i)
    counts.push_back(std::atoi(argv[i]));

  std::vector<std::vector<double>> areas;
  std::cout << "points,ms\n";
  for (const int count : counts) {
    TitanSASA::Options options;
    options.point_count = count;
    std::vector<double> area;
    double best_ms = 1.0e300;
    for (int repeat = 0; repeat < 3; ++repeat) {
      const auto begin = std::chrono::steady_clock::now();
      area = TitanSASA::compute_per_atom_sasa(
          std::span<const TitanSASA::Atom>(atoms.data(), atoms.size()), options);
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - begin)
                            .count();
      best_ms = std::min(best_ms, ms);
    }
    std::cout << count << ',' << std::setprecision(6) << best_ms << '\n';
    areas.push_back(std::move(area));
  }

  std::ofstream out(argv[2]);
  out << std::setprecision(12);
  for (std::size_t c = 0; c < counts.size(); ++c)
    out << (c ? "," : "") << "points_" << counts[c];
  out << '\n';
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    for (std::size_t c = 0; c < counts.size(); ++c)
      out << (c ? "," : "") << areas[c][a];
    out << '\n';
  }
  return EXIT_SUCCESS;
}
