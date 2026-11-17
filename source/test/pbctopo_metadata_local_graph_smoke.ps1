param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_metadata_local_graph_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$source = Join-Path $WorkRoot 'metadata_local_graph_smoke.cpp'
$exe = Join-Path $WorkRoot 'metadata_local_graph_smoke.exe'

@'
#include "pbctopo/local_graph.h"
#include "pbctopo/reference_strain.h"
#include "Multi_MOL.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

mol_info make_atom(std::size_t raw_resnr, const char *atom_name,
                   const char *chain, const char *res_name,
                   std::size_t chainnum, std::size_t internal_resnr) {
  mol_info info(raw_resnr, "C", atom_name, chain, res_name);
  info.set_chainnum(chainnum);
  info.set_internal_resnr(internal_resnr);
  info.set_original_atom_id(raw_resnr * 100 + internal_resnr);
  return info;
}

} // namespace

int main() {
  Multi_MOL system;
  system.info.push_back(make_atom(1, "CA", "A", "ALA", 7, 0));
  system.info.push_back(make_atom(1, "CB", "A", "ALA", 7, 0));
  system.info.push_back(make_atom(1, "CA", "A", "GLY", 7, 1));
  system.info.push_back(make_atom(1, "O", "A", "GLY", 7, 1));
  system.info.push_back(make_atom(1, "CA", "A", "SER", 7, 10));
  system.info.push_back(make_atom(1, "CA", "A", "THR", 7, 11));
  system.info.push_back(make_atom(1, "CA", "B", "ALA", 8, 2));

  MOL mol;
  mol.pos.resize(system.info.size());
  mol.pos.x[0] = 1.0;
  mol.pos.y[0] = 1.0;
  mol.pos.z[0] = 1.0;
  mol.pos.x[1] = 19.0;
  mol.pos.y[1] = 1.0;
  mol.pos.z[1] = 1.0;
  mol.pos.x[2] = 3.0;
  mol.pos.y[2] = 3.0;
  mol.pos.z[2] = 3.0;
  mol.pos.x[3] = 4.0;
  mol.pos.y[3] = 3.0;
  mol.pos.z[3] = 3.0;
  mol.pos.x[4] = 3.2;
  mol.pos.y[4] = 3.0;
  mol.pos.z[4] = 3.0;
  mol.pos.x[5] = 12.0;
  mol.pos.y[5] = 12.0;
  mol.pos.z[5] = 12.0;
  mol.pos.x[6] = 15.0;
  mol.pos.y[6] = 15.0;
  mol.pos.z[6] = 15.0;
  mol.Box = box(20.0, 20.0, 20.0, 90.0, 90.0, 90.0);

  titan_pbctopo::PbctopoSystemStaticData data =
      titan_pbctopo::build_pbctopo_system_static_data(system, mol);

  expect(data.metadata_graph_available, "metadata graph not available");
  expect(data.metadata_residue_groups == 5, "wrong metadata residue group count");
  expect(data.metadata_intra_edges == 2, "wrong metadata intra-edge count");
  expect(data.metadata_inter_edges == 3, "wrong metadata inter-edge count");
  expect(data.metadata_inter_hard_edges == 1 &&
             data.metadata_inter_soft_edges == 2,
         "metadata gap or implausible distance was not downgraded");
  expect(data.projected_metadata_edges.size() == 5,
         "wrong projected metadata edge count");
  expect(data.continuity_pairs_local.size() == 5,
         "metadata edges were not promoted to continuity pairs");
  std::size_t hard_continuity = 0;
  std::size_t soft_continuity = 0;
  for (const auto &pair : data.continuity_pairs_local) {
    if (titan_pbctopo::pbctopo_continuity_is_hard(pair.kind))
      ++hard_continuity;
    else
      ++soft_continuity;
  }
  expect(hard_continuity == 3 && soft_continuity == 2,
         "metadata continuity hard/soft kind was not preserved");
  titan_pbctopo::PbctopoContinuityPair soft_probe;
  soft_probe.kind = titan_pbctopo::PbctopoContinuityKind::MetadataSoft;
  expect(!titan_pbctopo::pbctopo_continuity_hard_violation(
             soft_probe, 10000.0, 25.0),
         "soft metadata continuity was treated as a hard rejection");
  titan_pbctopo::PbctopoContinuityPair hard_probe;
  hard_probe.kind = titan_pbctopo::PbctopoContinuityKind::MetadataHard;
  expect(titan_pbctopo::pbctopo_continuity_hard_violation(
             hard_probe, 36.0, 25.0),
         "hard metadata continuity did not enforce its cutoff");
  expect(data.continuity_anchor_local.size() == 5,
         "metadata residue anchors were not recorded");

  std::vector<titan_pbctopo::PbctopoAtom> atoms;
  titan_pbctopo::append_atoms_for_system(mol, &system.info, &data, 0, atoms);
  expect(atoms[0].chain_key == 7, "chain key does not use internal chainnum");
  expect(atoms[0].residue_id == 0, "residue id does not use internal resnr");
  expect(atoms[2].residue_id == 1,
         "raw residue labels were used instead of internal resnr");

  std::vector<titan_pbctopo::PbctopoLocalEdge> edges;
  std::unordered_set<std::uint64_t> edge_set;
  titan_pbctopo::PbctopoTimingStats timing;
  titan_pbctopo::append_metadata_graph_edges_for_system(data, 0, atoms, edges,
                                                        edge_set, timing);
  expect(edges.size() == 3, "soft metadata gap entered hard edge graph");
  expect(timing.pbctopo_metadata_residue_groups == 5,
         "metadata residue groups were not counted");
  expect(timing.pbctopo_metadata_intra_edges == 2,
         "metadata intra edges were not counted");
  expect(timing.pbctopo_metadata_inter_edges == 3,
         "metadata inter edges were not counted");
  expect(timing.pbctopo_metadata_inter_hard_edges == 1 &&
             timing.pbctopo_metadata_inter_soft_edges == 2,
         "metadata hard/soft inter-edge counters mismatch");
  for (const auto &edge : edges) {
    expect(edge.kind == titan_pbctopo::PbctopoLocalEdgeKind::MetadataSkeleton,
           "metadata edge kind was not preserved");
    expect(edge.ref_d2 == 0.0, "metadata edge ref distance is not neutral");
    expect(edge.inv_ref_d2 == 0.0, "metadata edge inverse ref is not neutral");
    expect(std::isinf(edge.break_cutoff_d2),
           "metadata edge break cutoff is not neutral");
    const double strain = titan_pbctopo::pbctopo_symmetric_reference_strain(
        4.0, edge.ref_d2);
    expect(std::isfinite(strain) && strain == 0.0,
           "metadata edge contributed non-finite strain");
  }

  edges.push_back({5, 6, 400.0, 1.0 / 400.0, 900.0,
                   titan_pbctopo::PbctopoLocalEdgeKind::GeometryHeuristic});
  const auto hard_components =
      titan_pbctopo::build_hard_connected_components_from_local_edges(
          atoms.size(), edges);
  expect(hard_components.size() == 4 && hard_components[0].size() == 4 &&
             hard_components[1].size() == 1 &&
             hard_components[2].size() == 1 &&
             hard_components[3].size() == 1,
         "isolated atom was rigid-merged through a geometry fallback edge");

  edges.back().kind = titan_pbctopo::PbctopoLocalEdgeKind::TopologyBond;
  const auto topology_components =
      titan_pbctopo::build_hard_connected_components_from_local_edges(
          atoms.size(), edges);
  expect(topology_components.size() == 3,
         "topology edge did not participate in hard connectivity");

  std::vector<titan_pbctopo::PbctopoAtom> geometry_atoms(3);
  for (auto &atom : geometry_atoms) {
    atom.is_heavy = true;
    atom.owner = 0;
  }
  geometry_atoms[0].sx = 0.10;
  geometry_atoms[1].sx = 0.12;
  geometry_atoms[2].sx = 0.11;
  std::vector<titan_pbctopo::PbctopoLocalEdge> geometry_edges;
  std::unordered_set<std::uint64_t> geometry_edge_set;
  const std::vector<char> eligible{1, 1, 0};
  titan_pbctopo::build_local_neighborhood_graph(
      geometry_atoms, mol.Box, geometry_edges, geometry_edge_set, 1.0, false,
      eligible);
  expect(geometry_edges.size() == 1 && geometry_edges[0].a != 2 &&
             geometry_edges[0].b != 2,
         "geometry fallback included an atom outside its eligibility mask");

  geometry_atoms.resize(2);
  geometry_atoms[0].sx = 0.10;
  geometry_atoms[1].sx = 0.40;
  geometry_edges.clear();
  geometry_edge_set.clear();
  titan_pbctopo::build_local_neighborhood_graph(
      geometry_atoms, mol.Box, geometry_edges, geometry_edge_set, 1.0);
  expect(geometry_edges.empty(),
         "unbounded nearest-neighbor fallback created a long edge");

  std::cout << "PASS pbctopo_metadata_local_graph_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $source -Encoding ASCII

$gxx = (Get-Command g++ -ErrorAction Stop).Source
$compileArgs = @(
    '-std=c++23',
    '-O2',
    '-Isrc',
    $source,
    'src/pbctopo/local_graph.cpp',
    'src/pbctopo/lattice_math.cpp',
    'src/topology/mol_topology.cpp',
    '-o',
    $exe
)

Push-Location $repo
try {
    & $gxx @compileArgs
    if ($LASTEXITCODE -ne 0) {
        throw "g++ exited with code $LASTEXITCODE"
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "metadata local graph smoke exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
