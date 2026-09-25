param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_evidence_connectivity_margin_oracle')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'oracle.cpp'
$exe = Join-Path $WorkRoot 'oracle.exe'

@'
#include "pbctopo/evidence_connectivity_margin.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using namespace titan_pbctopo;

void fail(const char *message, std::size_t node_count,
          std::uint64_t assignment, double threshold) {
  std::cerr << "FAIL pbctopo_evidence_connectivity_margin_oracle: "
            << message << " node_count=" << node_count
            << " assignment=" << assignment
            << " threshold=" << threshold << "\n";
  std::exit(1);
}

std::size_t component_count(
    std::size_t node_count,
    const std::vector<PbctopoEvidenceSupportEdge> &edges,
    bool apply_threshold, double threshold) {
  std::vector<std::vector<std::size_t>> adjacency(node_count);
  for (const auto &edge : edges) {
    if (apply_threshold && edge.support < threshold)
      continue;
    adjacency[edge.a].push_back(edge.b);
    adjacency[edge.b].push_back(edge.a);
  }

  std::vector<bool> visited(node_count, false);
  std::size_t count = 0;
  for (std::size_t root = 0; root < node_count; ++root) {
    if (visited[root])
      continue;
    ++count;
    std::vector<std::size_t> stack{root};
    visited[root] = true;
    while (!stack.empty()) {
      const std::size_t node = stack.back();
      stack.pop_back();
      for (const std::size_t neighbor : adjacency[node]) {
        if (visited[neighbor])
          continue;
        visited[neighbor] = true;
        stack.push_back(neighbor);
      }
    }
  }
  return count;
}

struct OracleResult {
  std::size_t thresholded_edge_count = 0;
  std::size_t thresholded_component_count = 0;
  bool thresholded_connected = false;
  bool supplied_relation_spans = false;
  bool margin_defined = false;
  double maximum_bottleneck_support =
      std::numeric_limits<double>::quiet_NaN();
};

OracleResult brute_force_oracle(
    std::size_t node_count,
    const std::vector<PbctopoEvidenceSupportEdge> &edges,
    double threshold) {
  OracleResult result;
  result.thresholded_edge_count = static_cast<std::size_t>(std::count_if(
      edges.begin(), edges.end(),
      [threshold](const auto &edge) { return edge.support >= threshold; }));
  result.thresholded_component_count =
      component_count(node_count, edges, true, threshold);
  result.thresholded_connected = result.thresholded_component_count == 1;
  result.supplied_relation_spans =
      node_count == 1 || component_count(node_count, edges, false, 0.0) == 1;

  if (node_count == 1 || !result.supplied_relation_spans)
    return result;

  bool found_tree = false;
  double best_bottleneck = -std::numeric_limits<double>::infinity();
  const std::uint64_t subset_count = std::uint64_t{1} << edges.size();
  for (std::uint64_t mask = 0; mask < subset_count; ++mask) {
    if (static_cast<std::size_t>(std::popcount(mask)) != node_count - 1)
      continue;

    std::vector<PbctopoEvidenceSupportEdge> selected;
    selected.reserve(node_count - 1);
    double bottleneck = std::numeric_limits<double>::infinity();
    for (std::size_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
      if ((mask & (std::uint64_t{1} << edge_index)) == 0)
        continue;
      selected.push_back(edges[edge_index]);
      bottleneck = std::min(bottleneck, edges[edge_index].support);
    }
    if (component_count(node_count, selected, false, 0.0) != 1)
      continue;
    found_tree = true;
    best_bottleneck = std::max(best_bottleneck, bottleneck);
  }

  if (!found_tree)
    return result;
  result.margin_defined = true;
  result.maximum_bottleneck_support = best_bottleneck;
  return result;
}

bool edge_is_present(const PbctopoEvidenceSupportEdge &query,
                     const std::vector<PbctopoEvidenceSupportEdge> &edges) {
  return std::any_of(edges.begin(), edges.end(), [&](const auto &edge) {
    return query.a == edge.a && query.b == edge.b &&
           query.support == edge.support;
  });
}

int main() {
  constexpr double thresholds[] = {-0.5, 0.0, 0.5};
  std::uint64_t checked_cases = 0;

  for (std::size_t node_count = 1; node_count <= 4; ++node_count) {
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t a = 0; a < node_count; ++a) {
      for (std::size_t b = a + 1; b < node_count; ++b)
        pairs.emplace_back(a, b);
    }

    std::uint64_t assignment_count = 1;
    for (std::size_t i = 0; i < pairs.size(); ++i)
      assignment_count *= 4;

    for (std::uint64_t assignment = 0; assignment < assignment_count;
         ++assignment) {
      std::uint64_t encoded = assignment;
      std::vector<PbctopoEvidenceSupportEdge> edges;
      for (const auto &[a, b] : pairs) {
        const int state = static_cast<int>(encoded % 4);
        encoded /= 4;
        if (state == 0)
          continue;
        const double support = state == 1 ? -1.0 : state == 2 ? 0.0 : 1.0;
        edges.push_back({a, b, support});
      }

      for (const double threshold : thresholds) {
        ++checked_cases;
        const auto expected =
            brute_force_oracle(node_count, edges, threshold);
        const auto actual = compute_pbctopo_evidence_connectivity_margin(
            node_count, edges, threshold);

        if (actual.canonical_edge_count != edges.size())
          fail("canonical edge count mismatch", node_count, assignment,
               threshold);
        if (actual.thresholded_edge_count != expected.thresholded_edge_count)
          fail("thresholded edge count mismatch", node_count, assignment,
               threshold);
        if (actual.thresholded_component_count !=
            expected.thresholded_component_count)
          fail("thresholded component count mismatch", node_count, assignment,
               threshold);
        if (actual.thresholded_connected != expected.thresholded_connected)
          fail("thresholded connectedness mismatch", node_count, assignment,
               threshold);
        if (actual.supplied_relation_spans !=
            expected.supplied_relation_spans)
          fail("supplied-domain spanning mismatch", node_count, assignment,
               threshold);
        if (actual.margin_defined != expected.margin_defined)
          fail("margin definedness mismatch", node_count, assignment,
               threshold);

        if (!expected.margin_defined) {
          if (!std::isnan(actual.maximum_bottleneck_support) ||
              !std::isnan(actual.signed_margin) ||
              !actual.spanning_tree_witness.empty())
            fail("undefined margin exposed numeric/tree output", node_count,
                 assignment, threshold);
          continue;
        }

        const double expected_margin =
            expected.maximum_bottleneck_support - threshold;
        if (actual.maximum_bottleneck_support !=
                expected.maximum_bottleneck_support ||
            actual.signed_margin != expected_margin)
          fail("brute-force bottleneck or signed margin mismatch", node_count,
               assignment, threshold);
        if ((actual.signed_margin >= 0.0) != actual.thresholded_connected)
          fail("margin sign disagrees with connectedness", node_count,
               assignment, threshold);
        if (actual.spanning_tree_witness.size() != node_count - 1 ||
            component_count(node_count, actual.spanning_tree_witness, false,
                            0.0) != 1)
          fail("invalid spanning-tree witness", node_count, assignment,
               threshold);

        double witness_bottleneck = std::numeric_limits<double>::infinity();
        for (const auto &edge : actual.spanning_tree_witness) {
          if (edge.a >= edge.b || !edge_is_present(edge, edges))
            fail("tree witness contains a noncanonical or absent edge",
                 node_count, assignment, threshold);
          witness_bottleneck = std::min(witness_bottleneck, edge.support);
        }
        if (witness_bottleneck != expected.maximum_bottleneck_support ||
            !edge_is_present(actual.bottleneck_witness,
                             actual.spanning_tree_witness) ||
            actual.bottleneck_witness.support !=
                expected.maximum_bottleneck_support)
          fail("bottleneck witness mismatch", node_count, assignment,
               threshold);
      }
    }
  }

  if (checked_cases != 12495)
    fail("unexpected exhaustive case count", 0, checked_cases, 0.0);
  std::cout << "PASS pbctopo_evidence_connectivity_margin_oracle cases="
            << checked_cases << "\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ascii

& $Compiler -std=c++23 -O2 -Wall -Wextra -pedantic -Isrc $src `
    src/pbctopo/evidence_connectivity_margin.cpp `
    -o $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile evidence connectivity-margin oracle"
}
& $exe
if ($LASTEXITCODE -ne 0) {
    throw "Evidence connectivity-margin oracle failed"
}
