param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\mol_topology_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'mol_topology_smoke.cpp'
$exe = Join-Path $WorkRoot 'mol_topology_smoke.exe'
$top = Join-Path $WorkRoot 'mini.top'
$itp = Join-Path $WorkRoot 'mini.itp'
$standaloneTop = Join-Path $WorkRoot 'standalone_ligand.top'

@'
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
CT 6 12.011 0.0 A 0.400000 0.250000
HZ 1  1.008 0.0 A 0.000000 0.000000

[ moleculetype ]
MOLA 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0
3 HZ 1 RES A3 1 0.0 1.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
MOLB 3

[ atoms ]
1 CT 1 RES B1 1 0.0 12.0
2 CT 1 RES B2 1 0.0 12.0

[ bonds ]
1 2 1
'@ | Set-Content -LiteralPath $itp -Encoding ASCII

@"
#include `"mini.itp`"

[ system ]
mini

[ molecules ]
MOLA 2
MOLB 1
"@ | Set-Content -LiteralPath $top -Encoding ASCII

@'
[ moleculetype ]
LIG 3

[ atoms ]
1 CT 1 LIG L1 1 0.0 12.0
2 CT 1 LIG L2 1 0.0 12.0

[ bonds ]
1 2 1

[ molecules ]
LIG 1
'@ | Set-Content -LiteralPath $standaloneTop -Encoding ASCII

@'
#include "mol_topology.h"
#include "MOL_info.h"

#include <iostream>
#include <cmath>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL mol_topology_smoke: " << message << "\n";
    std::exit(1);
  }
}

void expect_eq_size(size_t actual, size_t expected, const char *message) {
  if (actual != expected) {
    std::cerr << "FAIL mol_topology_smoke: " << message << ": got "
              << actual << " expected " << expected << "\n";
    std::exit(1);
  }
}

std::vector<mol_info> make_selection(std::initializer_list<size_t> source_ids) {
  std::vector<mol_info> info(source_ids.size());
  size_t i = 0;
  for (const size_t source_id : source_ids) {
    info[i].set_original_atom_id(source_id);
    ++i;
  }
  return info;
}

} // namespace

int main(int argc, char **argv) {
  expect(argc == 3, "expected topology path arguments");
  const std::string top_path = argv[1];

  auto full_info = make_selection({0, 1, 2, 3, 4, 5, 6, 7});
  auto full_result = build_gromacs_mol_topology(top_path, full_info);
  expect(full_result.has_value(), "full topology parse/project failed");
  const MolTopology &full = *full_result;
  expect(full.available(), "full topology not marked available");
  expect_eq_size(full.atom_count(), 8, "full atom_count mismatch");
  expect_eq_size(full.bonds().size(), 5, "full expanded bond count mismatch");
  expect_eq_size(full.neighbors(1).size(), 2,
                 "CSR neighbor count for atom 1 mismatch");
  expect_eq_size(full.neighbors(3).size(), 1,
                 "CSR neighbor count for atom 3 mismatch");
  expect_eq_size(full.molecule_ranges().size(), 3,
                 "molecule range expansion mismatch");
  expect_eq_size(full.vdw_radius_count(), 8,
                 "topology vdW radius expansion mismatch");
  expect(full.vdw_radius(0).has_value(), "topology vdW radius unavailable");
  const double expected_radius = 0.5 * std::pow(2.0, 1.0 / 6.0) * 4.0;
  expect(std::abs(*full.vdw_radius(0) - expected_radius) < 1.0e-10,
         "topology vdW radius value mismatch");
  expect(full.vdw_radius_source(0) ==
             MolTopologyVdwRadiusSource::GromacsTopology,
         "topology vdW radius source mismatch");
  expect(full.vdw_radius(2).has_value(),
         "explicit zero-LJ topology radius treated as unavailable");
  expect(*full.vdw_radius(2) == 0.0,
         "explicit zero-LJ topology radius was not preserved");
  expect(full.vdw_radius_source(2) ==
             MolTopologyVdwRadiusSource::GromacsTopologyExplicitZeroLJ,
         "zero-LJ topology radius lost force-field provenance");

  auto subset_info = make_selection({3, 4});
  auto subset_result = build_gromacs_mol_topology(top_path, subset_info);
  expect(subset_result.has_value(), "subset topology project failed");
  expect_eq_size(subset_result->atom_count(), 2, "subset atom_count mismatch");
  expect_eq_size(subset_result->bonds().size(), 1,
                 "subset projected bond count mismatch");
  expect_eq_size(subset_result->neighbors(0).size(), 1,
                 "subset CSR neighbor mismatch");
  expect_eq_size(subset_result->vdw_radius_count(), 2,
                 "subset projected vdW radius count mismatch");

  MolTopology parameter_cache;
  const std::vector<double> parameter_radii{1.7, 1.55};
  parameter_cache.cache_missing_vdw_radii(parameter_radii);
  parameter_cache.cache_missing_vdw_radii(parameter_radii);
  expect(!parameter_cache.available(),
         "radius-only cache unexpectedly marked graph available");
  expect_eq_size(parameter_cache.vdw_radius_count(), 2,
                 "parameter vdW radius cache count mismatch");
  expect(parameter_cache.vdw_radius_source(1) ==
             MolTopologyVdwRadiusSource::MultiMolParameters,
         "parameter vdW radius source mismatch");
  MolTopology explicit_zero_parameter_cache;
  const std::vector<double> zero_parameter_radii{0.0};
  const std::vector<double> zero_parameter_sigma{0.0};
  const std::vector<double> zero_parameter_epsilon{0.0};
  const std::vector<std::int32_t> assigned_type_ids{7};
  explicit_zero_parameter_cache.cache_missing_vdw_radii(
      zero_parameter_radii, zero_parameter_sigma, zero_parameter_epsilon,
      assigned_type_ids);
  expect(explicit_zero_parameter_cache.vdw_radius(0).has_value() &&
             *explicit_zero_parameter_cache.vdw_radius(0) == 0.0 &&
             explicit_zero_parameter_cache.vdw_radius_source(0) ==
                 MolTopologyVdwRadiusSource::
                     MultiMolParametersExplicitZeroLJ,
         "assigned parameter zero-LJ site was treated as missing");
  MolTopology missing_zero_parameter_cache;
  missing_zero_parameter_cache.cache_missing_vdw_radii(
      zero_parameter_radii, zero_parameter_sigma, zero_parameter_epsilon);
  expect(!missing_zero_parameter_cache.vdw_radius(0).has_value(),
         "unassigned zero radius was treated as explicit zero-LJ");
  missing_zero_parameter_cache.cache_missing_vdw_radii(
      zero_parameter_radii, zero_parameter_sigma, zero_parameter_epsilon,
      assigned_type_ids);
  expect(missing_zero_parameter_cache.vdw_radius(0).has_value() &&
             *missing_zero_parameter_cache.vdw_radius(0) == 0.0 &&
             missing_zero_parameter_cache.vdw_radius_source(0) ==
                 MolTopologyVdwRadiusSource::
                     MultiMolParametersExplicitZeroLJ,
         "late parameter assignment did not refresh zero-LJ provenance");

  auto source_zero_info = make_selection({0});
  auto source_zero_result =
      build_gromacs_mol_topology(top_path, source_zero_info);
  expect(source_zero_result.has_value(),
         "source atom ID zero was treated as unavailable");
  expect_eq_size(source_zero_result->source_atom_ids().front(), 0,
                 "source atom ID zero was not preserved");

  std::vector<mol_info> partial_info(2);
  partial_info[0].set_original_atom_id(3);
  auto partial_result = build_gromacs_mol_topology(top_path, partial_info);
  expect(!partial_result.has_value(),
         "partially populated source IDs unexpectedly projected");

  std::vector<mol_info> ambiguous_subset_info(2);
  auto ambiguous_subset_result =
      build_gromacs_mol_topology(top_path, ambiguous_subset_info);
  expect(!ambiguous_subset_result.has_value(),
         "subset without source provenance unexpectedly projected");

  std::vector<mol_info> sequential_info(8);
  auto sequential_result = build_gromacs_mol_topology(top_path, sequential_info);
  expect(sequential_result.has_value(), "sequential fallback failed");
  expect_eq_size(sequential_result->bonds().size(), 5,
                 "sequential fallback bond count mismatch");

  auto standalone_ligand_info = make_selection({6, 7});
  auto standalone_result =
      build_gromacs_mol_topology(argv[2], standalone_ligand_info);
  expect(standalone_result.has_value(),
         "standalone topology local fallback failed");
  expect_eq_size(standalone_result->bonds().size(), 1,
                 "standalone topology local fallback bond count mismatch");

  auto missing_result = build_gromacs_mol_topology("missing_topology.top",
                                                  full_info);
  expect(!missing_result.has_value(), "missing topology unexpectedly parsed");

  std::cout << "PASS mol_topology_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\mol_topology.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile mol_topology smoke'
}

& $exe $top $standaloneTop
if ($LASTEXITCODE -ne 0) {
    throw "mol_topology smoke failed with exit code $LASTEXITCODE"
}
