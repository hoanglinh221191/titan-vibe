// ENE-12 unit oracle for src/topology/gmx_topology_match.cpp:
//  1. match_by_names equals the former O(N_pdb x N_topo) search on random
//     structures (subsets, shuffles, N/C-prefixed residue names, repeats);
//  2. atom_name_key / residue_names_compatible on fixed cases;
//  3. compare_positional_names counts a rotated atom order as mismatches.
#include "topology/gmx_topology_match.h"

#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

int failures = 0;

void expect(bool condition, const std::string &what) {
  if (!condition) {
    ++failures;
    std::printf("FAIL %s\n", what.c_str());
  }
}

// The search post_process_atom_from_topo used before ENE-12.
std::vector<long long> reference_match(const std::vector<mol_info> &info,
                                       const std::vector<topo_info_temp> &topo) {
  std::vector<long long> out(info.size(), -1);
  std::vector<bool> consumed(topo.size(), false);
  for (std::size_t i = 0; i < info.size(); ++i) {
    for (std::size_t j = 0; j < topo.size(); ++j) {
      if (consumed[j])
        continue;
      bool names_match = info[i].get_res_name() == topo[j].res_name;
      if (!names_match) {
        const std::string pdb = info[i].get_res_name();
        const std::string top = topo[j].res_name;
        if (pdb.length() == top.length() + 1 && (pdb[0] == 'N' || pdb[0] == 'C') &&
            pdb.substr(1) == top)
          names_match = true;
      }
      if (names_match && info[i].get_atom() == topo[j].atom_name) {
        out[i] = static_cast<long long>(j);
        consumed[j] = true;
        break;
      }
    }
  }
  return out;
}

} // namespace

int main() {
  std::mt19937 rng(20261022);
  const std::vector<std::string> residues{"ALA", "GLY", "ASP", "LYS", "HIS"};
  const std::vector<std::string> atoms{"N", "CA", "C", "O", "CB", "H", "HA"};
  std::size_t compared = 0;
  for (int trial = 0; trial < 400; ++trial) {
    std::vector<topo_info_temp> topo;
    const int nres = 1 + static_cast<int>(rng() % 12);
    for (int r = 0; r < nres; ++r) {
      std::string res = residues[rng() % residues.size()];
      if (rng() % 7 == 0)
        res = (rng() % 2 ? "N" : "C") + res;
      const int natoms = 1 + static_cast<int>(rng() % atoms.size());
      for (int a = 0; a < natoms; ++a)
        topo.emplace_back(atoms[rng() % atoms.size()], "T", res, 0.1 * a, 12.0, r + 1);
    }
    std::vector<mol_info> info;
    for (const auto &t : topo) {
      if (rng() % 5 == 0)
        continue; // structure lacks this atom
      std::string res = t.res_name;
      if (rng() % 6 == 0)
        res = (rng() % 2 ? "N" : "C") + res; // terminal-prefixed in the PDB
      std::string atom = t.atom_name;
      if (rng() % 15 == 0)
        atom = "XX"; // atom absent from the topology
      info.emplace_back(static_cast<std::size_t>(t.res_nr), "?", atom, "A", res);
    }
    std::shuffle(info.begin(), info.end(), rng);
    const auto got = TitanTopo::match_by_names(info, topo);
    const auto want = reference_match(info, topo);
    expect(got == want, "match_by_names differs from the former search in trial " +
                            std::to_string(trial));
    compared += info.size();
  }

  expect(TitanTopo::atom_name_key(" 1hd1") == "HD11", "atom_name_key 1HD1");
  expect(TitanTopo::atom_name_key("O5*") == "O5'", "atom_name_key O5*");
  expect(TitanTopo::residue_names_compatible("HIS", "HSD"), "HIS ~ HSD");
  expect(TitanTopo::residue_names_compatible("HIE", "HISB") , "HIE ~ HISB");
  expect(TitanTopo::residue_names_compatible("NALA", "ALA"), "NALA ~ ALA");
  expect(TitanTopo::residue_names_compatible("ASP", "ASH"), "ASP ~ ASH");
  expect(TitanTopo::residue_names_compatible("CHIE", "HSE"), "CHIE ~ HSE");
  expect(!TitanTopo::residue_names_compatible("ALA", "GLY"), "ALA !~ GLY");
  expect(!TitanTopo::residue_names_compatible("ASP", "GLU"), "ASP !~ GLU");

  std::vector<topo_info_temp> topo;
  std::vector<mol_info> same, rotated;
  for (int r = 1; r <= 4; ++r)
    for (const char *atom : {"N", "CA", "C", "O"})
      topo.emplace_back(atom, "T", "ALA", 0.0, 12.0, r);
  for (const auto &t : topo)
    same.emplace_back(static_cast<std::size_t>(t.res_nr), "?", t.atom_name, "A", "ALA");
  for (int r = 1; r <= 4; ++r)
    for (const char *atom : {"CA", "C", "O", "N"})
      rotated.emplace_back(static_cast<std::size_t>(r), "?", atom, "A", "ALA");
  const auto ok = TitanTopo::compare_positional_names(same, topo);
  const auto bad = TitanTopo::compare_positional_names(rotated, topo);
  expect(ok.atom_mismatches == 0 && ok.residue_mismatches == 0, "same order has no mismatch");
  expect(bad.atom_mismatches == 16 && bad.examples.size() == 10,
         "rotated order: every atom mismatched, ten examples");

  if (failures == 0)
    std::printf("PASS topology_name_match_oracle: %zu atoms matched like the former search, "
                "name rules and positional check\n", compared);
  return failures == 0 ? 0 : 1;
}
