param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_virtual_discrete_search_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_virtual_discrete_search_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_virtual_discrete_search_smoke.exe'

@'
#include "pbctopo/virtual_discrete_search.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void require_true(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL pbctopo_virtual_discrete_search_smoke: " << message
              << '\n';
    std::exit(1);
  }
}

box make_box() {
  box value;
  value.set_a(10.0);
  value.set_b(10.0);
  value.set_c(10.0);
  return value;
}

titan_pbctopo::PbctopoAtom make_atom(double sx, double sy, double sz) {
  titan_pbctopo::PbctopoAtom atom;
  atom.sx = sx;
  atom.sy = sy;
  atom.sz = sz;
  return atom;
}

int virtual_shift(double coord, double center) {
  const double low = center - 0.5;
  const double high = center + 0.5;
  for (int shift = 0; shift <= 2; ++shift) {
    const double image = coord + static_cast<double>(shift);
    if (image + 1.0e-12 >= low && image < high - 1.0e-12)
      return shift;
  }
  return -99;
}

bool same_state(const titan_pbctopo::PbctopoVirtualDiscreteState &lhs,
                const titan_pbctopo::PbctopoVirtualDiscreteState &rhs) {
  return std::fabs(lhs.cx - rhs.cx) <= 1.0e-12 &&
         std::fabs(lhs.cy - rhs.cy) <= 1.0e-12 &&
         std::fabs(lhs.cz - rhs.cz) <= 1.0e-12 &&
         std::fabs(lhs.seam_sum - rhs.seam_sum) <= 1.0e-12 &&
         lhs.source == rhs.source;
}

} // namespace

int main() {
  using namespace titan_pbctopo;
  const box Box = make_box();

  std::vector<PbctopoAtom> narrow_atoms{
      make_atom(0.1000, 0.20, 0.30), make_atom(0.1005, 0.20, 0.30)};
  const auto narrow_plan = build_pbctopo_virtual_discrete_search_plan(
      narrow_atoms, {}, Box, {});
  const bool kept_narrow_interval = std::any_of(
      narrow_plan.x_states.begin(), narrow_plan.x_states.end(),
      [](const PbctopoVirtualAxisState &state) {
        const double width = state.interval_end - state.interval_begin;
        return width > 1.0e-5 && width < 8.3e-3;
      });
  require_true(kept_narrow_interval,
               "exact sweep must preserve a sub-grid-width interval");

  std::vector<PbctopoAtom> seam_atoms{
      make_atom(0.90, 0.25, 0.35), make_atom(0.10, 0.25, 0.35)};
  PbctopoLocalEdge topology_edge;
  topology_edge.a = 0;
  topology_edge.b = 1;
  topology_edge.kind = PbctopoLocalEdgeKind::TopologyBond;
  const auto seam_plan = build_pbctopo_virtual_discrete_search_plan(
      seam_atoms, {topology_edge}, Box, {});
  require_true(!seam_plan.x_states.empty(), "x axis states missing");
  require_true(seam_plan.x_states.front().topology_broken == 0,
               "topology seam minimum should preserve the wrapped bond");
  for (const auto &state : seam_plan.x_states) {
    const int delta = virtual_shift(seam_atoms[1].sx, state.center) -
                      virtual_shift(seam_atoms[0].sx, state.center);
    const std::size_t brute_broken = delta == 1 ? 0u : 1u;
    require_true(state.topology_broken == brute_broken,
                 "sweep seam count must match brute force");
  }

  require_true(narrow_plan.search_exhaustive,
               "small state space should be exhaustive");
  require_true(narrow_plan.states.size() == narrow_plan.total_state_space,
               "exhaustive state count should equal Cartesian product");

  std::vector<PbctopoAtom> large_atoms;
  for (int idx = 0; idx < 20; ++idx) {
    const double value = 0.02 + 0.045 * static_cast<double>(idx);
    large_atoms.push_back(make_atom(value, value + 0.001, value + 0.002));
  }
  PbctopoVirtualDiscreteSearchOptions shortlist_options;
  shortlist_options.top_k_per_axis = 4;
  shortlist_options.max_exhaustive_states = 4096;
  const auto shortlist_plan = build_pbctopo_virtual_discrete_search_plan(
      large_atoms, {}, Box, {}, shortlist_options);
  require_true(!shortlist_plan.search_exhaustive,
               "large state space should use the axis shortlist");
  require_true(shortlist_plan.states.size() == 64,
               "4-by-4-by-4 shortlist should evaluate 64 states");
  require_true(std::isfinite(shortlist_plan.omitted_seam_lower_bound),
               "shortlist should report an omitted-state lower bound");
  require_true(!shortlist_plan.branch_and_bound_used &&
                   !std::isfinite(shortlist_plan.lower_bound_gap),
               "P1 shortlist must not claim branch-and-bound certification");

  std::vector<PbctopoVirtualDiscreteSeedHint> hints{
      {2.35, 2.35, 2.35, 10.0, PbctopoVirtualSeedSource::Previous}};
  const auto hinted_plan = build_pbctopo_virtual_discrete_search_plan(
      narrow_atoms, {}, Box, hints);
  require_true(hinted_plan.x_states.size() == narrow_plan.x_states.size(),
               "seed hints must not alter the exact axis state set");
  require_true(hinted_plan.states.front().source ==
                   PbctopoVirtualSeedSource::Previous,
               "hint source should be carried for diagnostics");

  const auto repeat_plan = build_pbctopo_virtual_discrete_search_plan(
      narrow_atoms, {}, Box, hints);
  require_true(repeat_plan.states.size() == hinted_plan.states.size(),
               "repeat plan size changed");
  for (std::size_t idx = 0; idx < hinted_plan.states.size(); ++idx) {
    require_true(same_state(hinted_plan.states[idx], repeat_plan.states[idx]),
                 "discrete state ordering must be deterministic");
  }

  PbctopoVirtualDiscreteSearchWorkspace workspace;
  PbctopoVirtualDiscreteSearchPlan workspace_plan_first;
  PbctopoVirtualDiscreteSearchPlan workspace_plan_second;
  build_pbctopo_virtual_discrete_search_plan(
      large_atoms, {}, Box, hints, shortlist_options, workspace,
      workspace_plan_first);
  build_pbctopo_virtual_discrete_search_plan(
      large_atoms, {}, Box, hints, shortlist_options, workspace,
      workspace_plan_second);
  require_true(!workspace_plan_first.workspace_reused,
               "first workspace build must allocate its sweep capacity");
  require_true(workspace_plan_second.workspace_reused,
               "second workspace build must reuse sweep capacity");
  require_true(workspace_plan_second.workspace_bytes_estimate > 0,
               "workspace memory estimate must be reported");
  require_true(workspace_plan_second.event_count_x > 0 &&
                   workspace_plan_second.event_count_y > 0 &&
                   workspace_plan_second.event_count_z > 0,
               "per-axis event counts must be reported");
  require_true(workspace_plan_first.states.size() ==
                   workspace_plan_second.states.size(),
               "workspace reuse changed plan size");
  for (std::size_t idx = 0; idx < workspace_plan_first.states.size(); ++idx) {
    require_true(same_state(workspace_plan_first.states[idx],
                            workspace_plan_second.states[idx]),
                 "workspace reuse changed deterministic state ordering");
  }

  std::cout << "PASS pbctopo_virtual_discrete_search_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\virtual_discrete_search.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile pbctopo virtual discrete search smoke'
    }
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "pbctopo virtual discrete search smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
