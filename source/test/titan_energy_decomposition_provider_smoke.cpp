#include "titan_energy_decomposition.h"
#include "titan_electrostatics.h"

#include "arguments.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_energy_decomposition_provider_smoke: " << message
              << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-12) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_energy_decomposition_provider_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box(double length_a) {
  box out;
  out.set_a(length_a);
  out.set_b(length_a);
  out.set_c(length_a);
  out.set_alpha(90.0);
  out.set_beta(90.0);
  out.set_gamma(90.0);
  return out;
}

Multi_MOL make_molecule(const std::vector<double> &xs, double charge_scale,
                        double sigma_scale, double epsilon_scale,
                        double box_length_a = 80.0) {
  Multi_MOL out;
  out.info.resize(xs.size());
  out.parameters.resize(xs.size());
  for (std::size_t i = 0; i < xs.size(); ++i) {
    out.parameters[i].set_data(charge_scale * static_cast<double>(i + 1), 0.1,
                               sigma_scale * static_cast<double>(i + 1),
                               epsilon_scale * static_cast<double>(i + 1),
                               1.0);
  }
  MOL frame(make_box(box_length_a), 0.0, "AMBER99SB");
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(xs.size(), 0.0);
  frame.pos.z.assign(xs.size(), 0.0);
  out.m_mol.push_back(std::move(frame));
  return out;
}

EnergyDecompositionProvider make_provider(const Multi_MOL &rec,
                                          const Multi_MOL &lig,
                                          cache_sigma_epsilon &cache,
                                          arguments &args, double cutoff_nm,
                                          bool cutoff_active,
                                          const std::string &force_field =
                                              "AMBER99SB") {
  cache.force_field_name = force_field;
  args.force_f = force_field;
  args.epsilon_in = 1.0;
  EnergyDecompositionProvider provider;
  provider.configure(rec, lig, cache, args, cutoff_nm, cutoff_active);
  return provider;
}

void require_batch_matches_scalar(const EnergyDecompositionProvider &provider,
                                  const MOL &rec_frame, const MOL &lig_frame,
                                  std::size_t rec_begin, std::size_t rec_end,
                                  std::size_t lig_begin, std::size_t lig_end,
                                  const std::string &label) {
  const std::size_t rec_count = rec_end - rec_begin;
  const std::size_t lig_count = lig_end - lig_begin;
  std::vector<double> elec(rec_count * lig_count, -1.0);
  std::vector<double> vdw(rec_count * lig_count, -1.0);
  require(provider.inter_pair_energies_batch(rec_frame, lig_frame, rec_begin,
                                             rec_end, lig_begin, lig_end, elec,
                                             vdw),
          label + " batch accepted");

  for (std::size_t lig_pos = lig_begin; lig_pos < lig_end; ++lig_pos) {
    for (std::size_t rec_pos = rec_begin; rec_pos < rec_end; ++rec_pos) {
      const auto scalar =
          provider.pair_energy(rec_frame, lig_frame, rec_pos, lig_pos);
      const std::size_t row = lig_pos - lig_begin;
      const std::size_t col = rec_pos - rec_begin;
      const std::size_t index =
          titan_energy_cache_index(row, rec_count, col);
      require_close(elec[index], scalar.elec, label + " elec parity");
      require_close(vdw[index], scalar.vdw, label + " vdw parity");
    }
  }
}

void run_cache_index_and_minimal_amber_contract() {
  require(titan_energy_cache_index(2, 5, 3) == 13,
          "cache index is mol2 row-major over mol1 columns");
  Multi_MOL rec =
      make_molecule({0.0, 3.0, 7.0}, 0.1, 0.12, 0.20);
  Multi_MOL lig = make_molecule({11.0, 17.0}, -0.2, 0.15, 0.25);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  require_batch_matches_scalar(provider, rec.m_mol[0], lig.m_mol[0], 1, 3, 0,
                               2, "minimal AMBER subset");
}

void run_amber_lennard_jones_scale_contract() {
  Multi_MOL rec = make_molecule({0.0}, 0.1, 0.30, 0.20);
  Multi_MOL lig = make_molecule({5.0}, -0.2, 0.40, 0.45);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  const auto pair = provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double sigma = 0.5 * (0.30 + 0.40);
  const double epsilon = std::sqrt(0.20 * 0.45);
  const double sr6 = std::pow(sigma / 0.5, 6.0);
  const double expected_vdw = 4.0 * epsilon * (sr6 * sr6 - sr6);
  require_close(pair.vdw, expected_vdw,
                "AMBER Lennard-Jones uses 4*epsilon convention");
}

void run_force_field_formula_contracts() {
  Multi_MOL rec = make_molecule({0.0}, 0.1, 0.30, 0.20);
  Multi_MOL lig = make_molecule({5.0}, -0.2, 0.40, 0.45);
  arguments args;
  cache_sigma_epsilon cache;

  auto opls = make_provider(rec, lig, cache, args, -1.0, false, "OPLSAA");
  const auto opls_pair =
      opls.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double opls_sigma = std::sqrt(0.30 * 0.40);
  const double opls_epsilon = std::sqrt(0.20 * 0.45);
  const double opls_sr6 = std::pow(opls_sigma / 0.5, 6.0);
  require_close(opls_pair.vdw,
                4.0 * opls_epsilon *
                    (opls_sr6 * opls_sr6 - opls_sr6),
                "OPLS geometric mixing contract");

  rec.parameters[0].set_epsilon_sigma_radius(4.0e-6, 2.0e-3, 0.1);
  lig.parameters[0].set_epsilon_sigma_radius(9.0e-6, 8.0e-3, 0.1);
  auto gromos =
      make_provider(rec, lig, cache, args, -1.0, false, "GROMOS54A7");
  const auto gromos_pair =
      gromos.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double c6_ij = std::sqrt(2.0e-3 * 8.0e-3);
  const double c12_ij = std::sqrt(4.0e-6 * 9.0e-6);
  const double inv_r6 = std::pow(1.0 / 0.5, 6.0);
  require_close(gromos_pair.vdw,
                c12_ij * inv_r6 * inv_r6 - c6_ij * inv_r6,
                "GROMOS C12/C6 coefficient contract");

  rec.parameters[0].set_epsilon_sigma_radius(0.20, 0.30, 0.1);
  lig.parameters[0].set_epsilon_sigma_radius(0.45, 0.40, 0.1);
  rec.parameters[0].set_type_id(0);
  lig.parameters[0].set_type_id(1);
  cache.charmm_type_count = 2;
  cache.charmm_pair_sigma.assign(4, -1.0);
  cache.charmm_pair_epsilon.assign(4, -1.0);
  cache.charmm_pair_sigma[1] = cache.charmm_pair_sigma[2] = 0.25;
  cache.charmm_pair_epsilon[1] = cache.charmm_pair_epsilon[2] = 0.80;
  auto charmm =
      make_provider(rec, lig, cache, args, -1.0, false, "CHARMM36");
  const auto charmm_pair =
      charmm.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double charmm_sr6 = std::pow(0.25 / 0.5, 6.0);
  require_close(charmm_pair.vdw,
                4.0 * 0.80 * (charmm_sr6 * charmm_sr6 - charmm_sr6),
                "CHARMM pair override contract");
}

void run_screened_electrostatics_contract() {
  Multi_MOL rec = make_molecule({0.0}, 0.5, 0.30, 0.20);
  Multi_MOL lig = make_molecule({10.0}, -0.25, 0.40, 0.45);
  arguments args;
  cache_sigma_epsilon cache;
  args.flag_elec_screening = true;
  args.salt_concentration = 0.15;
  args.epsilon_in = 2.0;
  auto provider =
      make_provider(rec, lig, cache, args, -1.0, false, "AMBER99SB");
  args.flag_elec_screening = true;
  args.salt_concentration = 0.15;
  args.epsilon_in = 2.0;
  provider.configure(rec, lig, cache, args, -1.0, false);
  const auto pair =
      provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double lambda_debye = titan_electrostatics::debye_length_nm(0.15);
  const double expected = 138.935457839 * 0.5 * -0.25 *
                          std::exp(-1.0 / lambda_debye) / (2.0 * 1.0);
  require_close(pair.elec, expected, "AMBER Debye screening contract");
}

void run_cutoff_and_raw_coordinate_contract() {
  Multi_MOL rec = make_molecule({1.0}, 0.1, 0.12, 0.20, 20.0);
  Multi_MOL lig = make_molecule({3.0, 19.0}, -0.2, 0.15, 0.25, 20.0);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, 0.3, true);
  std::vector<double> elec(2, -1.0);
  std::vector<double> vdw(2, -1.0);
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 1,
                                             0, 2, elec, vdw),
          "raw-coordinate/cutoff batch accepted");
  require_close(elec[0], provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0)
                             .elec,
                "raw-coordinate pair elec matches scalar");
  require(std::fabs(elec[0]) > 0.0, "raw pair inside cutoff is retained");
  require_close(elec[1], 0.0, "pair beyond cutoff has zero elec");
  require_close(vdw[1], 0.0, "pair beyond cutoff has zero vdw");
}

void run_cross_boundary_raw_coordinate_contract() {
  Multi_MOL rec = make_molecule({1.0}, 0.5, 0.30, 0.20, 20.0);
  Multi_MOL lig = make_molecule({19.0}, -0.25, 0.40, 0.45, 20.0);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  const auto pair =
      provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double expected = 138.935457839 * 0.5 * -0.25 / 1.8;
  require_close(pair.elec, expected,
                "orthorhombic cross-boundary coordinates remain raw", 1e-9);
}

void run_triclinic_raw_coordinate_contract() {
  Multi_MOL rec = make_molecule({0.0}, 0.5, 0.30, 0.20, 20.0);
  Multi_MOL lig = make_molecule({11.0}, -0.25, 0.40, 0.45, 20.0);
  box triclinic = make_box(20.0);
  triclinic.set_gamma(60.0);
  rec.m_mol[0].Box = triclinic;
  lig.m_mol[0].Box = triclinic;
  lig.m_mol[0].pos.y[0] = 10.0 * std::sqrt(3.0);

  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  const auto pair =
      provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0);
  const double raw_distance_nm =
      0.1 * std::sqrt(11.0 * 11.0 + 300.0);
  const double expected =
      138.935457839 * 0.5 * -0.25 / raw_distance_nm;
  require_close(pair.elec, expected, "triclinic coordinates remain raw", 1e-9);
}

void run_tail_and_determinism_contract() {
  Multi_MOL rec =
      make_molecule({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 0.1, 0.12, 0.20);
  Multi_MOL lig = make_molecule({10.0, 14.0}, -0.2, 0.15, 0.25);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  for (std::size_t rec_count : {5u, 6u, 7u}) {
    require_batch_matches_scalar(provider, rec.m_mol[0], lig.m_mol[0], 0,
                                 rec_count, 0, 2,
                                 "AMBER vector block plus tail");
  }

  std::vector<double> first_elec(14, -1.0);
  std::vector<double> first_vdw(14, -1.0);
  std::vector<double> second_elec(14, -1.0);
  std::vector<double> second_vdw(14, -1.0);
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 7,
                                             0, 2, first_elec, first_vdw),
          "first deterministic batch accepted");
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 7,
                                             0, 2, second_elec, second_vdw),
          "second deterministic batch accepted");
  require(first_elec == second_elec && first_vdw == second_vdw,
          "repeated batch output is deterministic");
}

} // namespace

int main() {
  run_cache_index_and_minimal_amber_contract();
  run_amber_lennard_jones_scale_contract();
  run_force_field_formula_contracts();
  run_screened_electrostatics_contract();
  run_cutoff_and_raw_coordinate_contract();
  run_cross_boundary_raw_coordinate_contract();
  run_triclinic_raw_coordinate_contract();
  run_tail_and_determinism_contract();
  std::cout << "PASS titan_energy_decomposition_provider_smoke\n";
  return 0;
}
