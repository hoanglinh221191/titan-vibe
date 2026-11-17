// GBSA kernel oracle (ENE-14, ENE-15): the polar terms of the fused kernel
// against the scalar titan_gb::polar_pair_energy_kj_mol.
//
// Synthetic molecules of 627 atoms (jittered lattice, charges in
// [-0.8, 0.8], free Born radii 0.13..0.35 nm, complex radii larger). For
// Still GB without salt, with 0.15 M salt and with ALPB, through the real
// cal_energy_gbsa_avx2:
//   1. ENE-15: receptor-only and ligand-only calls (the other range empty,
//      so no cross pair) give self + intramolecular Delta-GB; they must
//      equal the long-double scalar sum to 1e-12;
//   2. ENE-15: the whole kernel in one call and split in the engine's tiles
//      (64-atom receptor tiles with every ligand atom, plus ligand one-body
//      tiles) must agree to 1e-12 relative (molecules 22 A apart);
//   3. ENE-14: the cross pairs of two molecules in contact (one-body terms
//      off) must equal the long-double scalar sum to 1e-12; the 3.2e-6
//      exponential of the old cross loop misses this by orders of magnitude.
// It prints every sum to 17 digits, so the values of two builds can be put
// side by side. Linked against the engine objects.
#include "gbsa_intrin.h"
#include "MOL.h"
#include "MOL_info.h"
#include "MOL_para.h"
#include "cache_sigma_epsilon.h"
#include "titan_electrostatics.h"
#include "titan_gb_model.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
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
  x ^= x >> 15;
  x *= 2246822519u;
  x ^= x >> 13;
  x *= 3266489917u;
  x ^= x >> 16;
  return (x & 0xffffffu) / double(0x1000000);
}

Molecule make(unsigned seed, double shift_x) {
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
    m.parameters.radius[a] = 1.2 + 0.8 * hash01(seed + 11, a); // Angstrom
    m.parameters.sigma[a] = 0.25 + 0.1 * hash01(seed + 13, a);
    m.parameters.epsilon[a] = 0.2 + 0.3 * hash01(seed + 17, a);
    m.free_radii[a] = 0.13 + 0.22 * hash01(seed + 19, a);
    m.complex_radii[a] = m.free_radii[a] * (1.0 + 0.15 * hash01(seed + 23, a));
    m.info.emplace_back(1 + a / 10, "C", "C" + std::to_string(a + 1), "A",
                        "T");
  }
  return m;
}

long double reference_onebody(const Molecule &m, const arguments &args,
                              double kappa, double size_complex,
                              double size_free) {
  const bool alpb = args.flag_alpb;
  long double sum = 0.0L;
  for (std::size_t i = 0; i < kAtoms; ++i) {
    const double q = m.parameters.charge[i];
    sum += 0.5L * (static_cast<long double>(titan_gb::polar_pair_energy_kj_mol(
                       q, q, 0.0, m.complex_radii[i], m.complex_radii[i],
                       kappa, args.int_diel, args.ext_diel, alpb,
                       size_complex)) -
                   titan_gb::polar_pair_energy_kj_mol(
                       q, q, 0.0, m.free_radii[i], m.free_radii[i], kappa,
                       args.int_diel, args.ext_diel, alpb, size_free));
    for (std::size_t k = i + 1; k < kAtoms; ++k) {
      const double dx = m.mol.pos.x[i] - m.mol.pos.x[k];
      const double dy = m.mol.pos.y[i] - m.mol.pos.y[k];
      const double dz = m.mol.pos.z[i] - m.mol.pos.z[k];
      const double r2 = 0.01 * (dx * dx + dy * dy + dz * dz);
      sum += static_cast<long double>(titan_gb::polar_pair_energy_kj_mol(
                 q, m.parameters.charge[k], r2, m.complex_radii[i],
                 m.complex_radii[k], kappa, args.int_diel, args.ext_diel,
                 alpb, size_complex)) -
             titan_gb::polar_pair_energy_kj_mol(
                 q, m.parameters.charge[k], r2, m.free_radii[i],
                 m.free_radii[k], kappa, args.int_diel, args.ext_diel, alpb,
                 size_free);
    }
  }
  return sum;
}

long double reference_cross(const Molecule &rec, const Molecule &lig,
                            const arguments &args, double kappa,
                            double size_complex) {
  long double sum = 0.0L;
  for (std::size_t j = 0; j < kAtoms; ++j)
    for (std::size_t i = 0; i < kAtoms; ++i) {
      const double dx = rec.mol.pos.x[i] - lig.mol.pos.x[j];
      const double dy = rec.mol.pos.y[i] - lig.mol.pos.y[j];
      const double dz = rec.mol.pos.z[i] - lig.mol.pos.z[j];
      const double r2 = 0.01 * (dx * dx + dy * dy + dz * dz);
      sum += titan_gb::polar_pair_energy_kj_mol(
          rec.parameters.charge[i], lig.parameters.charge[j], r2,
          rec.complex_radii[i], lig.complex_radii[j], kappa, args.int_diel,
          args.ext_diel, args.flag_alpb, size_complex);
    }
  return sum;
}

// A_det of each molecule and of the complex, as the kernel computes them
// without a SASA cache
void adet_sizes(const Molecule &rec, const Molecule &lig, double &size_rec,
                double &size_lig, double &size_complex) {
  std::vector<double> xs, ys, zs, rs;
  auto add = [&](const Molecule &m) {
    for (std::size_t a = 0; a < kAtoms; ++a) {
      xs.push_back(m.mol.pos.x[a]);
      ys.push_back(m.mol.pos.y[a]);
      zs.push_back(m.mol.pos.z[a]);
      rs.push_back(m.parameters.radius[a]);
    }
  };
  add(rec);
  add(lig);
  size_rec = titan_gb::electrostatic_size_adet_nm(
      {xs.data(), kAtoms}, {ys.data(), kAtoms}, {zs.data(), kAtoms},
      {rs.data(), kAtoms});
  size_lig = titan_gb::electrostatic_size_adet_nm(
      {xs.data() + kAtoms, kAtoms}, {ys.data() + kAtoms, kAtoms},
      {zs.data() + kAtoms, kAtoms}, {rs.data() + kAtoms, kAtoms});
  size_complex = titan_gb::electrostatic_size_adet_nm(
      {xs.data(), 2 * kAtoms}, {ys.data(), 2 * kAtoms},
      {zs.data(), 2 * kAtoms}, {rs.data(), 2 * kAtoms});
}

GBSA_energy_snap call(Molecule &rec, Molecule &lig, const arguments &args,
                      std::size_t rb, std::size_t re, std::size_t lb,
                      std::size_t le, bool rec_onebody, bool lig_onebody) {
  cache_sigma_epsilon cache;
  GBSA_energy_snap e;
  cal_energy_gbsa_avx2(rec.mol, lig.mol, rec.info, lig.info, cache, rb, re,
                       lb, le, e, 0, args, rec.parameters, lig.parameters,
                       rec.lcpo, rec.lcpo, lig.lcpo, lig.lcpo, rec.free_radii,
                       lig.free_radii, rec.complex_radii, lig.complex_radii,
                       nullptr, nullptr, false, rec_onebody, lig_onebody);
  return e;
}

int failures = 0;
void check(double got, long double ref, double tol, const std::string &what) {
  const double rel = static_cast<double>(std::abs((got - ref) / ref));
  std::printf("    %-34s %.17g  (ref %.17Lg, rel %.2e)\n", what.c_str(), got,
              ref, rel);
  if (!(rel <= tol)) {
    ++failures;
    std::printf("  FAIL %s: relative difference %.3e > %.0e\n", what.c_str(),
                rel, tol);
  }
}

} // namespace

int main() {
  Molecule rec = make(1u, 0.0), lig = make(2u, 20.0 + 2.2 * 9);
  // in contact: the lattices are 2.2 A apart along x
  Molecule lig_near = make(3u, 2.2 * 9);
  struct Case {
    const char *name;
    double salt;
    bool alpb;
  } cases[] = {{"still_nosalt", 0.0, false},
               {"still_salt0.15", 0.15, false},
               {"alpb_salt0.15", 0.15, true}};
  for (const auto &c : cases) {
    arguments args;
    args.igb = 5;
    args.int_diel = 1.0;
    args.ext_diel = 78.5;
    args.saltcon = c.salt;
    args.gbsa_temperature = 298.15;
    args.flag_alpb = c.alpb;
    args.force_f = "AMBER99SB";
    args.gbsa_sasa_points = 12;
    const double kappa = titan_electrostatics::amber_gb_kappa_nm_inverse(
        args.saltcon, args.ext_diel, args.gbsa_temperature);
    std::printf("case %s\n", c.name);

    double size_rec = 0, size_lig = 0, size_complex = 0;
    if (c.alpb)
      adet_sizes(rec, lig, size_rec, size_lig, size_complex);

    const auto rec_only = call(rec, lig, args, 0, kAtoms, 0, 0, true, false);
    const auto lig_only = call(rec, lig, args, 0, 0, 0, kAtoms, false, true);
    check(rec_only.G_polar,
          reference_onebody(rec, args, kappa, size_complex, size_rec), 1e-12,
          "receptor self + intramolecular");
    check(lig_only.G_polar,
          reference_onebody(lig, args, kappa, size_complex, size_lig), 1e-12,
          "ligand self + intramolecular");

    const auto whole = call(rec, lig, args, 0, kAtoms, 0, kAtoms, true, true);
    GBSA_energy_snap tiled;
    for (std::size_t b = 0; b < kAtoms; b += 64) {
      const auto t = call(rec, lig, args, b, std::min(kAtoms, b + 64), 0,
                          kAtoms, true, false);
      tiled.G_polar += t.G_polar;
      tiled.E_vdw += t.E_vdw;
      tiled.E_elec += t.E_elec;
    }
    for (std::size_t b = 0; b < kAtoms; b += 64) {
      const auto t = call(rec, lig, args, 0, 0, b, std::min(kAtoms, b + 64),
                          false, true);
      tiled.G_polar += t.G_polar;
      tiled.E_vdw += t.E_vdw;
      tiled.E_elec += t.E_elec;
    }
    check(tiled.G_polar, whole.G_polar, 1e-12, "G_polar tiled vs one call");
    check(tiled.E_vdw, whole.E_vdw, 1e-12, "E_vdw tiled vs one call");
    // E_elec here is a ~50 kJ/mol net of ~1e4 kJ/mol of pair terms of both
    // signs: its summation-order noise is 1e-13..1e-12 relative
    check(tiled.E_elec, whole.E_elec, 1e-11, "E_elec tiled vs one call");
    std::printf("    G_polar whole %.17g, E_vdw %.17g, E_elec %.17g\n",
                whole.G_polar, whole.E_vdw, whole.E_elec);

    double near_rec = 0, near_lig = 0, near_complex = 0;
    if (c.alpb)
      adet_sizes(rec, lig_near, near_rec, near_lig, near_complex);
    const auto cross =
        call(rec, lig_near, args, 0, kAtoms, 0, kAtoms, false, false);
    check(cross.G_polar,
          reference_cross(rec, lig_near, args, kappa, near_complex), 1e-12,
          "ENE-14 cross G_polar in contact");
  }
  if (failures) {
    std::printf("FAIL GBSA kernel %d check(s)\n", failures);
    return 1;
  }
  std::printf("PASS GBSA kernel: intramolecular Delta-GB and cross G_polar "
              "equal to the scalar long-double sums to 1e-12 (ENE-15, "
              "ENE-14); tiled kernel equal to one call\n");
  return 0;
}
