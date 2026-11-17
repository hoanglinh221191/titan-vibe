// Where the fused GBSA kernel spends its time (ENE-15, ENE-14, ENE-04):
// one thread, in process, two synthetic 627-atom molecules in contact, a SASA
// cache so that no surface is computed. Times the cross pairs with and
// without the per-residue decomposition and the one-body (self, nonpolar,
// intramolecular Delta-GB) calls of each molecule, 15 repetitions, min and
// median. Build it against two trees and alternate the two executables
// pinned to one P-core (on this machine separate-process timings drift):
//   python test/engine_linked_oracles.py --cases gbsa_timing --work-root <dir>
//   <dir>/gbsa_timing.exe [CHARMM36|AMBER99SB] [repetitions]
#include "gbsa_intrin.h"
#include "MOL.h"
#include "MOL_info.h"
#include "MOL_para.h"
#include "cache_sigma_epsilon.h"
#include "info.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {
constexpr std::size_t kAtoms = 627;
struct Molecule {
  MOL mol;
  std::vector<mol_info> info;
  Multi_MOL_parameters parameters;
  AlignedDoubleVector free_radii, complex_radii;
  std::vector<double> lcpo;
};
double hash01(unsigned a, unsigned b) {
  unsigned x = a * 2654435761u ^ (b + 0x9e3779b9u + (a << 6) + (a >> 2));
  x ^= x >> 15; x *= 2246822519u; x ^= x >> 13; x *= 3266489917u; x ^= x >> 16;
  return (x & 0xffffffu) / double(0x1000000);
}
Molecule make(unsigned seed, double shift_x, const char *chain) {
  Molecule m;
  m.mol.pos.resize(kAtoms);
  m.parameters.resize(kAtoms);
  m.free_radii.resize(kAtoms);
  m.complex_radii.resize(kAtoms);
  m.lcpo.assign(kAtoms, 0.0);
  for (std::size_t a = 0; a < kAtoms; ++a) {
    const std::size_t ix = a % 9, iy = (a / 9) % 9, iz = a / 81;
    m.mol.pos.x[a] = shift_x + 2.2 * ix + 0.5 * (hash01(seed, 3 * a) - 0.5);
    m.mol.pos.y[a] = 2.2 * iy + 0.5 * (hash01(seed, 3 * a + 1) - 0.5);
    m.mol.pos.z[a] = 2.2 * iz + 0.5 * (hash01(seed, 3 * a + 2) - 0.5);
    m.parameters.charge[a] = 1.6 * (hash01(seed + 7, a) - 0.5);
    m.parameters.radius[a] = 1.2 + 0.8 * hash01(seed + 11, a);
    m.parameters.sigma[a] = 0.25 + 0.1 * hash01(seed + 13, a);
    m.parameters.epsilon[a] = 0.2 + 0.3 * hash01(seed + 17, a);
    m.free_radii[a] = 0.13 + 0.22 * hash01(seed + 19, a);
    m.complex_radii[a] = m.free_radii[a] * (1.0 + 0.15 * hash01(seed + 23, a));
    m.info.emplace_back(1 + a / 10, "C", "C" + std::to_string(a + 1), chain, "T");
  }
  return m;
}
} // namespace

int main(int argc, char **argv) {
  const std::string ff = argc > 1 ? argv[1] : "CHARMM36";
  const int reps = argc > 2 ? std::atoi(argv[2]) : 15;
  Molecule rec = make(1u, 0.0, "A"), lig = make(3u, 2.2 * 9, "B");
  arguments args;
  args.igb = 8;
  args.saltcon = 0.0;
  args.force_f = ff;
  info inf(args);
  GBSASasaCache sasa;
  sasa.sasa_rec_free.assign(kAtoms, 0.0);
  sasa.sasa_lig_free.assign(kAtoms, 0.0);
  sasa.sasa_complex.assign(2 * kAtoms, 0.0);
  cache_sigma_epsilon cache;
  auto run = [&](const char *what, std::size_t rb, std::size_t re,
                 std::size_t lb, std::size_t le, bool ro, bool lo,
                 bool decomp) {
    std::vector<double> t;
    double g = 0;
    for (int r = 0; r < reps; ++r) {
      GBSA_Aggregator agg(args, inf);
      GBSA_energy_snap e;
      const auto t0 = std::chrono::steady_clock::now();
      cal_energy_gbsa_avx2(rec.mol, lig.mol, rec.info, lig.info, cache, rb,
                           re, lb, le, e, 0, args, rec.parameters,
                           lig.parameters, rec.lcpo, rec.lcpo, lig.lcpo,
                           lig.lcpo, rec.free_radii, lig.free_radii,
                           rec.complex_radii, lig.complex_radii,
                           decomp ? &agg : nullptr, &sasa, true, ro, lo);
      t.push_back(std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t0)
                      .count());
      g = e.G_polar + e.E_vdw + e.E_elec;
    }
    std::sort(t.begin(), t.end());
    std::printf("  %s %-36s min %7.3f  median %7.3f ms  (sum %.10g)\n",
                ff.c_str(), what, t.front(), t[t.size() / 2], g);
  };
  run("cross pairs, decomposition", 0, kAtoms, 0, kAtoms, false, false, true);
  run("cross pairs, no decomposition", 0, kAtoms, 0, kAtoms, false, false,
      false);
  run("receptor one-body, decomposition", 0, kAtoms, 0, 0, true, false, true);
  run("ligand one-body, decomposition", 0, 0, 0, kAtoms, false, true, true);
  return 0;
}
