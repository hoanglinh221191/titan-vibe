// ENE-05: the fused free/complex Born-radius passes must return exactly what
// the separate calls returned.
//
// For every GB model (igb 1/2/5/8) and rgbmax (off, 25, 12 A) the four
// binding endpoints from cal_binding_born_radii_parallel_avx2() are compared
// bit for bit with four cal_effective_born_radii_parallel_avx2() calls, and
// cal_free_and_complex_born_radii_avx2() over sub-ranges with the two serial
// calls it replaces (the CAS path).  The inputs are barnase/barstar from the
// JPCB binding fixture, with CAS deletion sentinels on one residue and a
// coincident-atom pair in a second, small case to cover the NaN path.  A
// timing block then alternates old and new on a ~13k-atom receptor (igb 5,
// or GB_DRIVER_TIMING_IGB).
//
// The pair loops are vectorised (AVX2): HCT for igb 1/2/5, HCT, far-field
// series and neck term for igb 8 (GBn2). With
// GB_DRIVER_DUMP=<file> every endpoint array is written out; running once
// with TITAN_FORCE_SCALAR=1 (the scalar loop) and once without, then
// "gb_binding_radii_fused_driver --compare scalar.bin simd.bin" checks that
// the vectorised radii agree with the scalar ones to 1e-12 relative and that
// NaN radii sit at the same atoms.
//
// Build: link with the engine objects of the same profile (all but Main.o).
// Usage: gb_binding_radii_fused_driver <barnase.pqr> <barstar.pqr>
//            [timing_threads=1] [timing_rounds=10]
//        gb_binding_radii_fused_driver --compare <a.bin> <b.bin>
// e.g. cmd /c start /affinity 1 /wait /b gb_binding_radii_fused_driver ...
// pins the run to logical CPU 0 (a performance core on the dev machine).
#include "MOL.h"
#include "MOL_info.h"
#include "gbsa_intrin.h"
#include "titan_executor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::ofstream *g_dump = nullptr;

void dump_array(const AlignedDoubleVector &values) {
  const std::uint64_t n = values.size();
  g_dump->write(reinterpret_cast<const char *>(&n), sizeof(n));
  g_dump->write(reinterpret_cast<const char *>(values.data()),
                static_cast<std::streamsize>(n * sizeof(double)));
}

// Reads two dumps record by record; returns the number of failures.
int compare_dumps(const char *a_path, const char *b_path) {
  std::ifstream a(a_path, std::ios::binary), b(b_path, std::ios::binary);
  if (!a || !b) {
    std::cerr << "cannot open dumps\n";
    return 1;
  }
  double worst = 0;
  std::size_t arrays = 0, values = 0, nan_mismatch = 0;
  std::uint64_t na = 0, nb = 0;
  while (a.read(reinterpret_cast<char *>(&na), sizeof(na)) &&
         b.read(reinterpret_cast<char *>(&nb), sizeof(nb))) {
    if (na != nb) {
      std::cerr << "array size mismatch\n";
      return 1;
    }
    std::vector<double> x(na), y(nb);
    a.read(reinterpret_cast<char *>(x.data()),
           static_cast<std::streamsize>(na * sizeof(double)));
    b.read(reinterpret_cast<char *>(y.data()),
           static_cast<std::streamsize>(nb * sizeof(double)));
    for (std::size_t i = 0; i < na; ++i) {
      if (std::isnan(x[i]) || std::isnan(y[i])) {
        nan_mismatch += std::isnan(x[i]) != std::isnan(y[i]);
        continue;
      }
      worst = std::max(worst, std::fabs(x[i] - y[i]) /
                                  std::max(std::fabs(x[i]), 1e-300));
      ++values;
    }
    ++arrays;
  }
  std::cout << "compare: arrays=" << arrays << " radii=" << values
            << " max relative difference=" << worst
            << " nan position mismatches=" << nan_mismatch << "\n";
  const bool ok = arrays > 0 && worst <= 1e-12 && nan_mismatch == 0;
  std::cout << (ok ? "PASS" : "FAIL") << " gb radii scalar vs vectorised\n";
  return ok ? 0 : 1;
}

struct Molecule {
  MOL mol;
  std::vector<mol_info> info;
  AlignedDoubleVector vdw;
};

Molecule read_pqr(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  Molecule out;
  std::vector<double> x, y, z;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0)
      continue;
    std::istringstream fields(line);
    std::string record, serial, atom, residue, chain;
    long resnr = 0;
    double px = 0, py = 0, pz = 0, q = 0, radius = 0;
    fields >> record >> serial >> atom >> residue >> chain >> resnr >> px >>
        py >> pz >> q >> radius;
    x.push_back(px);
    y.push_back(py);
    z.push_back(pz);
    out.info.emplace_back(static_cast<size_t>(resnr), atom, atom, chain,
                          residue);
    out.vdw.push_back(radius);
  }
  out.mol.pos.resize(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    out.mol.pos.x[i] = x[i];
    out.mol.pos.y[i] = y[i];
    out.mol.pos.z[i] = z[i];
  }
  return out;
}

Molecule replicate(const Molecule &base, int copies, double shift) {
  Molecule out;
  const size_t n = base.mol.size();
  out.mol.pos.resize(n * copies);
  for (int c = 0; c < copies; ++c)
    for (size_t i = 0; i < n; ++i) {
      const size_t k = c * n + i;
      out.mol.pos.x[k] = base.mol.pos.x[i] + shift * c;
      out.mol.pos.y[k] = base.mol.pos.y[i];
      out.mol.pos.z[k] = base.mol.pos.z[i];
      out.info.push_back(base.info[i]);
      out.vdw.push_back(base.vdw[i]);
    }
  return out;
}

bool same_bits(const AlignedDoubleVector &a, const AlignedDoubleVector &b) {
  return a.size() == b.size() &&
         (a.empty() ||
          std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

struct Endpoints {
  AlignedDoubleVector rec_free, lg_free, rec_complex, lg_complex;
};

Endpoints old_binding(const Molecule &rec, const Molecule &lg, int igb,
                      double rgbmax, TitanExecutor &executor) {
  Endpoints e;
  e.rec_free.resize(rec.mol.size());
  e.lg_free.resize(lg.mol.size());
  e.rec_complex.resize(rec.mol.size());
  e.lg_complex.resize(lg.mol.size());
  cal_effective_born_radii_parallel_avx2(rec.mol, rec.vdw, rec.info,
                                         e.rec_free, igb, executor, nullptr,
                                         nullptr, nullptr, rgbmax);
  cal_effective_born_radii_parallel_avx2(lg.mol, lg.vdw, lg.info, e.lg_free,
                                         igb, executor, nullptr, nullptr,
                                         nullptr, rgbmax);
  cal_effective_born_radii_parallel_avx2(rec.mol, rec.vdw, rec.info,
                                         e.rec_complex, igb, executor,
                                         &lg.mol, &lg.vdw, &lg.info, rgbmax);
  cal_effective_born_radii_parallel_avx2(lg.mol, lg.vdw, lg.info,
                                         e.lg_complex, igb, executor,
                                         &rec.mol, &rec.vdw, &rec.info, rgbmax);
  return e;
}

Endpoints new_binding(const Molecule &rec, const Molecule &lg, int igb,
                      double rgbmax, TitanExecutor &executor) {
  Endpoints e;
  cal_binding_born_radii_parallel_avx2(rec.mol, rec.vdw, rec.info, lg.mol,
                                       lg.vdw, lg.info, igb, executor, rgbmax,
                                       e.rec_free, e.lg_free, e.rec_complex,
                                       e.lg_complex);
  return e;
}

bool same_endpoints(const Endpoints &a, const Endpoints &b) {
  return same_bits(a.rec_free, b.rec_free) && same_bits(a.lg_free, b.lg_free) &&
         same_bits(a.rec_complex, b.rec_complex) &&
         same_bits(a.lg_complex, b.lg_complex);
}

size_t count_nan(const Endpoints &e) {
  size_t n = 0;
  for (const auto *v : {&e.rec_free, &e.lg_free, &e.rec_complex,
                        &e.lg_complex})
    n += std::count_if(v->begin(), v->end(),
                       [](double r) { return std::isnan(r); });
  return n;
}

int compare_case(const char *label, const Molecule &rec, const Molecule &lg,
                 TitanExecutor &executor) {
  int failures = 0;
  const double inf = std::numeric_limits<double>::infinity();
  for (const int igb : {1, 2, 5, 8})
    for (const double rgbmax : {inf, 25.0, 12.0}) {
      const Endpoints ref = old_binding(rec, lg, igb, rgbmax, executor);
      const Endpoints got = new_binding(rec, lg, igb, rgbmax, executor);
      if (g_dump != nullptr)
        for (const auto *v : {&got.rec_free, &got.lg_free, &got.rec_complex,
                              &got.lg_complex})
          dump_array(*v);
      bool ok = same_endpoints(ref, got);

      // CAS path: serial, over sub-ranges, receptor and ligand as "mol".
      for (const bool receptor_side : {true, false}) {
        const Molecule &m = receptor_side ? rec : lg;
        const Molecule &env = receptor_side ? lg : rec;
        const size_t n = m.mol.size();
        const size_t split = n / 3;
        AlignedDoubleVector ref_free(n, -1.0), ref_complex(n, -1.0);
        AlignedDoubleVector got_free(n, -1.0), got_complex(n, -1.0);
        for (const auto [b, e] : {std::pair<size_t, size_t>{0, split},
                                  std::pair<size_t, size_t>{split, n}}) {
          cal_effective_born_radii_avx2(m.mol, m.vdw, m.info, ref_free, igb,
                                        b, e, nullptr, nullptr, nullptr,
                                        rgbmax);
          cal_effective_born_radii_avx2(m.mol, m.vdw, m.info, ref_complex, igb,
                                        b, e, &env.mol, &env.vdw, &env.info,
                                        rgbmax);
          cal_free_and_complex_born_radii_avx2(m.mol, m.vdw, m.info, got_free,
                                               got_complex, igb, b, e, env.mol,
                                               env.vdw, env.info, rgbmax);
        }
        ok = ok && same_bits(ref_free, got_free) &&
             same_bits(ref_complex, got_complex);
      }
      std::cout << label << " igb=" << igb << " rgbmax=" << rgbmax
                << " nan=" << count_nan(ref) << (ok ? " identical" : " MISMATCH")
                << "\n";
      failures += ok ? 0 : 1;
    }
  return failures;
}

} // namespace

int main(int argc, char **argv) {
  if (argc == 4 && std::string(argv[1]) == "--compare")
    return compare_dumps(argv[2], argv[3]);
  std::ofstream dump_stream;
  if (const char *dump_path = std::getenv("GB_DRIVER_DUMP")) {
    dump_stream.open(dump_path, std::ios::binary);
    g_dump = &dump_stream;
  }
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: gb_binding_radii_fused_driver <rec.pqr> <lig.pqr> "
                 "[timing_threads=1] [timing_rounds=10]\n";
    return 2;
  }
  const size_t timing_threads = argc > 3 ? std::stoul(argv[3]) : 1;
  const int timing_rounds = argc > 4 ? std::stoi(argv[4]) : 10;
  TitanExecutor executor(8);
  Molecule rec = read_pqr(argv[1]);
  const Molecule lg = read_pqr(argv[2]);

  // CAS deletion sentinels: one receptor residue's side chain stops
  // descreening but keeps a finite target radius.
  const size_t sentinel_res = rec.info[rec.info.size() / 2].get_resnr();
  size_t sentinels = 0;
  for (size_t i = 0; i < rec.info.size(); ++i) {
    const std::string name = rec.info[i].get_atom();
    if (rec.info[i].get_resnr() == sentinel_res && name != "N" &&
        name != "CA" && name != "C" && name != "O" && name != "H" &&
        name != "HA" && name != "CB") {
      rec.vdw[i] = 0.001;
      ++sentinels;
    }
  }
  std::cout << "receptor=" << rec.mol.size() << " ligand=" << lg.mol.size()
            << " cas_sentinels=" << sentinels << "\n";

  int failures = compare_case("barnase_barstar", rec, lg, executor);

  // Coincident atoms make the integral invalid: every radius of that atom
  // must be NaN on both paths, and only the complex one for a partner clash.
  Molecule small_rec = rec, small_lg = lg;
  small_rec.mol.pos.x[5] = small_rec.mol.pos.x[4];
  small_rec.mol.pos.y[5] = small_rec.mol.pos.y[4];
  small_rec.mol.pos.z[5] = small_rec.mol.pos.z[4];
  small_lg.mol.pos.x[7] = small_rec.mol.pos.x[9];
  small_lg.mol.pos.y[7] = small_rec.mol.pos.y[9];
  small_lg.mol.pos.z[7] = small_rec.mol.pos.z[9];
  failures += compare_case("coincident_atoms", small_rec, small_lg, executor);

  // Timing on a ~12k-atom receptor (seven barnase copies 45 A apart) and
  // barstar; igb 5 unless GB_DRIVER_TIMING_IGB says otherwise, no cutoff.  Old and new alternate within each round, in
  // flipped order every round; the median of the per-round ratios is robust
  // to the clock drift of this machine.  Run single-threaded and pinned to a
  // performance core for a stable number.
  const char *timing_igb_env = std::getenv("GB_DRIVER_TIMING_IGB");
  const int timing_igb = timing_igb_env != nullptr ? std::atoi(timing_igb_env) : 5;
  TitanExecutor timing_executor(timing_threads);
  const Molecule big = replicate(rec, 7, 45.0);
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> ratios;
  double best_old = 1e300, best_new = 1e300;
  for (int round = 0; round < timing_rounds; ++round) {
    double seconds[2] = {0.0, 0.0};
    for (const bool use_new : {round % 2 == 0, round % 2 != 0}) {
      const auto start = std::chrono::steady_clock::now();
      const Endpoints e =
          use_new ? new_binding(big, lg, timing_igb, inf, timing_executor)
                  : old_binding(big, lg, timing_igb, inf, timing_executor);
      seconds[use_new] = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();
      if (e.rec_free.empty())
        return 3;
    }
    best_old = std::min(best_old, seconds[0]);
    best_new = std::min(best_new, seconds[1]);
    ratios.push_back(seconds[0] / seconds[1]);
  }
  std::sort(ratios.begin(), ratios.end());
  std::cout << "timing igb=" << timing_igb << " receptor=" << big.mol.size()
            << " ligand="
            << lg.mol.size() << " threads=" << timing_threads
            << " rounds=" << timing_rounds << " best_old_s=" << best_old
            << " best_new_s=" << best_new
            << " min_ratio=" << best_old / best_new
            << " median_round_ratio=" << ratios[ratios.size() / 2]
            << " range=" << ratios.front() << ".." << ratios.back() << "\n";

  if (failures != 0) {
    std::cout << "FAIL gb_binding_radii_fused_driver (" << failures << ")\n";
    return 1;
  }
  std::cout << "PASS gb_binding_radii_fused_driver\n";
  return 0;
}
