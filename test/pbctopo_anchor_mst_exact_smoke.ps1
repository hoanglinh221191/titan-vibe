param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_anchor_mst_exact_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'pbctopo_anchor_mst_exact_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_anchor_mst_exact_smoke.exe'

@'
#include "pbctopo/anchor_graph.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_anchor_mst_exact_smoke: " << message << "\n";
    std::exit(1);
  }
}

} // namespace

int main() {
  using namespace titan_pbctopo;
  constexpr std::size_t count = 600;
  std::vector<Vec3> points;
  points.reserve(count);
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double t = static_cast<double>(idx);
    const double layer = static_cast<double>(idx % 7);
    points.push_back({0.13 * t + 0.017 * layer,
                      3.0 * std::sin(0.071 * t) + 0.19 * layer,
                      2.0 * std::cos(0.047 * t) + 0.11 * (idx % 11)});
  }
  std::vector<std::size_t> anchors(count);
  std::iota(anchors.begin(), anchors.end(), 0);
  std::vector<PbctopoContinuityPair> continuity;

  double dense_mst2 = 0.0;
  double dense_max_d2 = 0.0;
  PbctopoAnchorGraphScore dense_graph;
  PbctopoTimingStats dense_timing;
  require(anchor_graph_metrics_enhanced(
              points, anchors, continuity, dense_mst2, dense_max_d2,
              dense_graph, nullptr, PbctopoRefineScorePhase::Exact,
              AnchorBackend::Prim, &dense_timing, nullptr),
          "dense Prim failed");

  double spatial_mst2 = 0.0;
  double spatial_max_d2 = 0.0;
  PbctopoAnchorGraphScore spatial_graph;
  PbctopoTimingStats spatial_timing;
  require(anchor_graph_metrics_enhanced(
              points, anchors, continuity, spatial_mst2, spatial_max_d2,
              spatial_graph, nullptr, PbctopoRefineScorePhase::Exact,
              AnchorBackend::SpatialAll, &spatial_timing, nullptr),
          "certified spatial path failed");

  const double mst_tol = 1.0e-9 * std::max(1.0, std::abs(dense_mst2));
  const double max_tol = 1.0e-9 * std::max(1.0, std::abs(dense_max_d2));
  require(std::abs(spatial_mst2 - dense_mst2) <= mst_tol,
          "spatial MST weight differs from exact Prim");
  require(std::abs(spatial_max_d2 - dense_max_d2) <= max_tol,
          "spatial MST bottleneck differs from exact Prim");
  require(spatial_timing.component_graph_refine_anchor_spatial_calls == 1,
          "spatial backend was not exercised");
  require(
      spatial_timing.component_graph_refine_anchor_spatial_exact_certified +
              spatial_timing.component_graph_refine_anchor_spatial_fallbacks ==
          1,
      "spatial result was returned without an exact certificate or fallback");

  std::cout << "PASS pbctopo_anchor_mst_exact_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src `
    src\pbctopo\anchor_graph.cpp `
    src\pbctopo\candidate_order.cpp `
    src\pbctopo\lattice_math.cpp `
    -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile PBCTOPO exact anchor MST smoke'
}

& $exe
if ($LASTEXITCODE -ne 0) {
    throw "PBCTOPO exact anchor MST smoke failed with exit code $LASTEXITCODE"
}
