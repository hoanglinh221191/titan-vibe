#include "functions.h"
#include "titan_energy_decomposition.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require_close(double actual, double expected, const std::string &label,
                   double tolerance = 2.0e-10) {
  const double scale = std::max(1.0, std::fabs(expected));
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance * scale) {
    std::cerr << "FAIL titan_nonbonded_kernel_parity_smoke: " << label
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box() { return box(80.0, 80.0, 80.0, 90.0, 90.0, 90.0); }

Multi_MOL make_molecule(const std::vector<double> &xs,
                        const std::string &force_field, double charge_scale,
                        double sigma_scale, double epsilon_scale) {
  Multi_MOL out;
  out.info.resize(xs.size());
  out.parameters.resize(xs.size());
  for (std::size_t i = 0; i < xs.size(); ++i) {
    out.parameters[i].set_data(charge_scale * static_cast<double>(i + 1), 0.1,
                               sigma_scale * (1.0 + 0.05 * i),
                               epsilon_scale * (1.0 + 0.10 * i), 1.0);
  }
  MOL frame(make_box(), 0.0, force_field);
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(xs.size(), 0.0);
  frame.pos.z.assign(xs.size(), 0.0);
  out.m_mol.push_back(std::move(frame));
  return out;
}

void run_force_field(const std::string &force_field) {
  const bool gromos = force_field.find("GROMOS") != std::string::npos;
  Multi_MOL rec = make_molecule({0.0, 2.0, 4.0, 6.0, 8.0}, force_field, 0.08,
                                gromos ? 2.0e-3 : 0.30,
                                gromos ? 4.0e-6 : 0.20);
  Multi_MOL lig = make_molecule({15.0, 19.0}, force_field, -0.12,
                                gromos ? 8.0e-3 : 0.38,
                                gromos ? 9.0e-6 : 0.45);
  arguments args;
  args.force_f = force_field;
  args.epsilon_in = 2.5;
  cache_sigma_epsilon scalar_cache;
  scalar_cache.force_field_name = force_field;
  cache_sigma_epsilon avx_cache;
  avx_cache.force_field_name = force_field;
  cache_energy scalar, avx;

  if (force_field.find("AMBER") != std::string::npos ||
      force_field.find("CHARMM") != std::string::npos) {
    cal_energy_amber(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                     scalar_cache, 0, lig.m_mol[0].size(), scalar, 0, args,
                     rec.parameters, lig.parameters);
    cal_energy_amber_intrin(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                            avx_cache, 0, lig.m_mol[0].size(), avx, 0, args,
                            rec.parameters, lig.parameters);
  } else if (force_field.find("OPLS") != std::string::npos) {
    cal_energy_opls(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                    scalar_cache, 0, lig.m_mol[0].size(), scalar, 0, args,
                    rec.parameters, lig.parameters);
    cal_energy_opls_intrin(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                           avx_cache, 0, lig.m_mol[0].size(), avx, 0, args,
                           rec.parameters, lig.parameters);
  } else if (force_field.find("GROMOS") != std::string::npos) {
    cal_energy_gromos(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                      scalar_cache, 0, lig.m_mol[0].size(), scalar, 0, args,
                      rec.parameters, lig.parameters);
    cal_energy_gromos_intrin(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                             avx_cache, 0, lig.m_mol[0].size(), avx, 0, args,
                             rec.parameters, lig.parameters);
  } else {
    args.epsilon_in = 15.0;
    cal_energy_martini(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                       scalar_cache, 0, lig.m_mol[0].size(), scalar, 0, args,
                       rec.parameters, lig.parameters);
    cal_energy_martini_intrin(rec.m_mol[0], lig.m_mol[0], rec.info, lig.info,
                              avx_cache, 0, lig.m_mol[0].size(), avx, 0, args,
                              rec.parameters, lig.parameters);
  }

  require_close(avx.elec_ener, scalar.elec_ener,
                force_field + " scalar/AVX electrostatics");
  require_close(avx.vdw_ener, scalar.vdw_ener,
                force_field + " scalar/AVX Lennard-Jones");

  cache_sigma_epsilon provider_cache;
  provider_cache.force_field_name = force_field;
  EnergyDecompositionProvider provider;
  provider.configure(rec, lig, provider_cache, args, -1.0, false);
  double provider_elec = 0.0;
  double provider_vdw = 0.0;
  for (std::size_t j = 0; j < lig.m_mol[0].size(); ++j) {
    for (std::size_t i = 0; i < rec.m_mol[0].size(); ++i) {
      const auto pair =
          provider.pair_energy(rec.m_mol[0], lig.m_mol[0], i, j);
      provider_elec += pair.elec;
      provider_vdw += pair.vdw;
    }
  }
  require_close(provider_elec, scalar.elec_ener,
                force_field + " provider electrostatics");
  require_close(provider_vdw, scalar.vdw_ener,
                force_field + " provider Lennard-Jones");
}

} // namespace

int main() {
  para duplicate_types;
  duplicate_types.base.emplace_back(0.30, 0.20, 12.0, 1.0,
                                    std::string("CT"));
  duplicate_types.register_atom_type("CT");
  duplicate_types.base.emplace_back(0.31, 0.21, 12.0, 1.0,
                                    std::string("CT"));
  duplicate_types.register_atom_type("CT");
  if (duplicate_types.atom_type_to_id.at("CT") != 1) {
    std::cerr << "FAIL titan_nonbonded_kernel_parity_smoke: duplicate atom "
                 "types must use the final parameter record\n";
    return 1;
  }
  run_force_field("AMBER19SB");
  run_force_field("CHARMM36");
  run_force_field("OPLSAA");
  run_force_field("GROMOS54A7");
  run_force_field("MARTINI");
  std::cout << "PASS titan_nonbonded_kernel_parity_smoke\n";
  return 0;
}
