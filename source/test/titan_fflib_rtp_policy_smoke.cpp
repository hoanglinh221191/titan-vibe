#include "topology/fflib_charge_policy.h"
#include "topology/fflib_rtp_parser.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

const info_aa *find_atom(const amino_acid &amine, const std::string &residue,
                         const std::string &atom) {
  const auto residue_it = amine.amino_acids.find(residue);
  if (residue_it == amine.amino_acids.end())
    return nullptr;
  for (const auto &entry : residue_it->second)
    if (entry.atom == atom)
      return &entry;
  return nullptr;
}

} // namespace

int main(int argc, char **argv) {
  std::istringstream rtp(R"RTP(
[ bondedtypes ]
1 1 1 2 1 3 1 0

[ MET ]
  [ atoms ]
    N   NH1  -0.470  0
   HN     H   0.310  1
   CA   CT1   0.070  2
    C     C   0.510  3
    O     O  -0.510  4
  [ bonds ]
    N    HN
    N    CA
    C    CA
    C    +N
  [ impropers ]
    N    -C    CA    HN
  [ cmap ]
   -C     N    CA     C    +N

[ BENZ ]
  [ atoms ]
   CG CG2R61 -0.115 0
   HG  HGR61  0.115 1
  [ bonds ]
   CG HG
)RTP");
  amino_acid amine;
  const auto parsed = TitanFfLib::parse_gromacs_rtp_atoms(rtp, amine);
  if (!parsed || parsed->residue_count != 2 || parsed->atom_count != 7)
    return 1;
  if (amine.amino_acids.at("MET").size() != 5 ||
      amine.amino_acids.at("BENZ").size() != 2)
    return 2;
  const info_aa *n = find_atom(amine, "MET", "N");
  const info_aa *ca = find_atom(amine, "MET", "CA");
  const info_aa *benz = find_atom(amine, "BENZ", "CG");
  if (n == nullptr || n->atom_type != "NH1" ||
      std::fabs(n->charge + 0.470) > 1.0e-12)
    return 3;
  if (ca == nullptr || ca->atom_type != "CT1" ||
      std::fabs(ca->charge - 0.070) > 1.0e-12)
    return 4;
  if (benz == nullptr || benz->atom_type != "CG2R61")
    return 5;
  if (find_atom(amine, "MET", "-C") != nullptr ||
      find_atom(amine, "MET", "+N") != nullptr)
    return 6;

  std::istringstream malformed(R"RTP(
[ ALA ]
  [ atoms ]
  N NH1
)RTP");
  amino_acid rejected;
  if (TitanFfLib::parse_gromacs_rtp_atoms(malformed, rejected))
    return 7;

  using TitanFfLib::resolve_charge_source_policy;
  const auto none = resolve_charge_source_policy(false, false, false);
  const auto ligand_only = resolve_charge_source_policy(false, true, false);
  const auto receptor_only = resolve_charge_source_policy(false, false, true);
  const auto both = resolve_charge_source_policy(false, true, true);
  const auto common = resolve_charge_source_policy(true, false, false);
  if (none.molecule1_uses_topology || none.molecule2_uses_topology)
    return 8;
  if (ligand_only.molecule1_uses_topology || !ligand_only.molecule2_uses_topology)
    return 9;
  if (!receptor_only.molecule1_uses_topology || receptor_only.molecule2_uses_topology)
    return 10;
  if (!both.molecule1_uses_topology || !both.molecule2_uses_topology ||
      !common.molecule1_uses_topology || !common.molecule2_uses_topology)
    return 11;

  if (argc == 2) {
    std::ifstream actual(argv[1]);
    amino_acid actual_amine;
    const auto actual_result =
        TitanFfLib::parse_gromacs_rtp_atoms(actual, actual_amine);
    if (!actual_result)
      return 12;
    const info_aa *actual_n = find_atom(actual_amine, "MET", "N");
    const info_aa *actual_ca = find_atom(actual_amine, "MET", "CA");
    const info_aa *actual_benz = find_atom(actual_amine, "BENZ", "CG");
    if (actual_n == nullptr || actual_n->atom_type != "NH1" ||
        actual_ca == nullptr || actual_ca->atom_type != "CT1" ||
        actual_benz == nullptr || actual_benz->atom_type != "CG2R61")
      return 13;
    std::cout << " ACTUAL_RTP residues=" << actual_result->residue_count
              << " atoms=" << actual_result->atom_count;
  }

  std::cout << " PASS titan_fflib_rtp_policy_smoke residues=2 atoms=7\n";
  return 0;
}
