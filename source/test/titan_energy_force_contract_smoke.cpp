#include "functions.h"
#include "titan_electrostatics.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

void print_force_map(size_t, std::vector<std::vector<data_out>> &,
                     std::fstream &, bool, const arguments &,
                     std::vector<energy_res> &, size_t, info &) {}

namespace {

void require_close(double actual, double expected, const std::string &label,
                   double tolerance = 1.0e-11) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_energy_force_contract_smoke: " << label
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box() { return box(80.0, 80.0, 80.0, 90.0, 90.0, 90.0); }

MOL make_frame(double x, const std::string &force_field) {
  MOL frame(make_box(), 0.0, force_field);
  frame.pos.x = {x};
  frame.pos.y = {0.0};
  frame.pos.z = {0.0};
  return frame;
}

Multi_MOL_parameters make_parameters(double charge, double sigma,
                                     double epsilon) {
  Multi_MOL_parameters parameters(1);
  parameters[0].set_data(charge, 0.1, sigma, epsilon, 1.0);
  return parameters;
}

void run_amber_dielectric_and_screening_contract() {
  const MOL rec = make_frame(0.0, "AMBER99SB");
  const MOL lig = make_frame(10.0, "AMBER99SB");
  const auto rec_parameters = make_parameters(1.0, 0.0, 0.0);
  const auto lig_parameters = make_parameters(-0.5, 0.0, 0.0);
  const std::vector<mol_info> rec_info(1), lig_info(1);
  cache_sigma_epsilon cache;
  cache_energy energy;
  cache_force force_value;
  arguments args;
  args.force_f = "AMBER99SB";
  args.flag_cal_force = true;
  args.epsilon_in = 4.0;

  cal_energy_force_amber(rec, lig, rec_info, lig_info, cache, 0, 1, energy,
                         force_value, 0, args, rec_parameters, lig_parameters);
  const double prefactor = 138.935457839 * 1.0 * -0.5 / 4.0;
  require_close(energy.elec_ener, prefactor,
                "AMBER force path honors epsilon_in");
  require_close(force_value.forces[0].Fx, -prefactor,
                "AMBER unscreened force derivative");

  args.flag_elec_screening = true;
  args.salt_concentration = 0.15;
  cal_energy_force_amber(rec, lig, rec_info, lig_info, cache, 0, 1, energy,
                         force_value, 0, args, rec_parameters, lig_parameters);
  const double lambda = titan_electrostatics::debye_length_nm(0.15);
  require_close(lambda, 0.78456058871820455,
                "Debye length oracle at 0.15 M", 1.0e-13);
  const double screened_energy = prefactor * std::exp(-1.0 / lambda);
  const double screened_force =
      -prefactor * std::exp(-1.0 / lambda) * (1.0 + 1.0 / lambda);
  require_close(energy.elec_ener, screened_energy,
                "AMBER force path honors Debye screening");
  require_close(force_value.forces[0].Fx, screened_force,
                "AMBER screened force derivative");
}

void run_opls_and_gromos_dielectric_contract() {
  const MOL rec = make_frame(0.0, "OPLSAA");
  const MOL lig = make_frame(10.0, "OPLSAA");
  const auto rec_parameters = make_parameters(1.0, 0.0, 0.0);
  const auto lig_parameters = make_parameters(-0.5, 0.0, 0.0);
  const std::vector<mol_info> rec_info(1), lig_info(1);
  cache_sigma_epsilon cache;
  cache_energy energy;
  cache_force force_value;
  arguments args;
  args.flag_cal_force = true;
  args.epsilon_in = 5.0;
  const double expected = 138.935457839 * 1.0 * -0.5 / 5.0;

  args.force_f = "OPLSAA";
  cal_energy_force_opls(rec, lig, rec_info, lig_info, cache, 0, 1, energy,
                        force_value, 0, args, rec_parameters, lig_parameters);
  require_close(energy.elec_ener, expected,
                "OPLS force path honors epsilon_in");

  args.force_f = "GROMOS54A7";
  cal_energy_force_gromos(rec, lig, rec_info, lig_info, cache, 0, 1, energy,
                          force_value, 0, args, rec_parameters, lig_parameters);
  require_close(energy.elec_ener, expected,
                "GROMOS force path honors epsilon_in");
}

void run_filtered_cache_zero_contract() {
  const MOL rec = make_frame(0.0, "AMBER99SB");
  const MOL lig = make_frame(1.0, "AMBER99SB");
  const auto rec_parameters = make_parameters(0.5, 0.3, 0.2);
  const auto lig_parameters = make_parameters(-0.5, 0.3, 0.2);
  const std::vector<mol_info> rec_info(1), lig_info(1);
  arguments args;
  args.force_f = "AMBER99SB";
  args.flag_cal_force = true;
  args.self_int = true;
  cache_sigma_epsilon cache(1, args);
  cache.set_elec_vdw(0, 123.0, 456.0);
  cache.set_force(0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  cache_energy energy;
  cache_force force_value;
  cal_energy_force_amber(rec, lig, rec_info, lig_info, cache, 0, 1, energy,
                         force_value, 0, args, rec_parameters, lig_parameters);
  require_close(cache.elec_ener[0], 0.0,
                "filtered pair clears stale electrostatic cache");
  require_close(cache.vdw_ener[0], 0.0,
                "filtered pair clears stale Lennard-Jones cache");
  const auto cached_force = cache.get_force(0);
  if (!cached_force.has_value()) {
    std::cerr << "FAIL titan_energy_force_contract_smoke: force cache missing\n";
    std::exit(1);
  }
  for (const auto &component : cached_force->forces) {
    require_close(component.Fx, 0.0,
                  "filtered pair clears stale cached force x");
    require_close(component.Fy, 0.0,
                  "filtered pair clears stale cached force y");
    require_close(component.Fz, 0.0,
                  "filtered pair clears stale cached force z");
  }
}

} // namespace

int main() {
  run_amber_dielectric_and_screening_contract();
  run_opls_and_gromos_dielectric_contract();
  run_filtered_cache_zero_contract();
  std::cout << "PASS titan_energy_force_contract_smoke\n";
  return 0;
}
