param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_evidence_connectivity_margin_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
#include "pbctopo/evidence_connectivity_margin.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_evidence_connectivity_margin_smoke: "
              << message << "\n";
    std::exit(1);
  }
}

void expect_near(double actual, double expected, const char *message) {
  expect(std::abs(actual - expected) <= 1.0e-12, message);
}

template <class Function>
void expect_invalid(Function &&function, const char *message) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return;
  }
  expect(false, message);
}

template <class Function>
void expect_overflow(Function &&function, const char *message) {
  try {
    function();
  } catch (const std::overflow_error &) {
    return;
  }
  expect(false, message);
}

bool same_edge(const PbctopoEvidenceSupportEdge &lhs,
               const PbctopoEvidenceSupportEdge &rhs) {
  return lhs.a == rhs.a && lhs.b == rhs.b && lhs.support == rhs.support;
}

bool same_tree(const std::vector<PbctopoEvidenceSupportEdge> &lhs,
               const std::vector<PbctopoEvidenceSupportEdge> &rhs) {
  return lhs.size() == rhs.size() &&
         std::equal(lhs.begin(), lhs.end(), rhs.begin(), same_edge);
}

int main() {
  const std::vector<PbctopoEvidenceSupportEdge> positive_edges{
      {0, 1, 0.9}, {1, 2, 0.8}, {2, 3, 0.7}, {0, 3, 0.2}};
  const auto positive = compute_pbctopo_evidence_connectivity_margin(
      4, positive_edges, 0.6);
  expect(positive.canonical_edge_count == 4 &&
             positive.thresholded_edge_count == 3 &&
             positive.thresholded_component_count == 1 &&
             positive.thresholded_connected &&
             positive.supplied_relation_spans && positive.margin_defined &&
             positive.spanning_tree_witness.size() == 3,
         "connected positive-margin graph was summarized incorrectly");
  expect_near(positive.maximum_bottleneck_support, 0.7,
              "positive graph bottleneck support was incorrect");
  expect_near(positive.signed_margin, 0.1,
              "positive graph signed margin was incorrect");
  expect(same_edge(positive.bottleneck_witness, {2, 3, 0.7}),
         "positive graph bottleneck witness was incorrect");

  const std::vector<PbctopoEvidenceSupportEdge> negative_edges{
      {0, 1, 0.9}, {1, 2, 0.4}};
  const auto negative = compute_pbctopo_evidence_connectivity_margin(
      3, negative_edges, 0.5);
  expect(!negative.thresholded_connected &&
             negative.thresholded_component_count == 2 &&
             negative.supplied_relation_spans && negative.margin_defined,
         "threshold-disconnected spanning domain was summarized incorrectly");
  expect_near(negative.maximum_bottleneck_support, 0.4,
              "negative graph bottleneck support was incorrect");
  expect_near(negative.signed_margin, -0.1,
              "negative graph signed margin was incorrect");

  const std::vector<PbctopoEvidenceSupportEdge> triangle_edges{
      {0, 2, 0.8}, {1, 2, 0.7}, {0, 1, 0.9}};
  const auto triangle = compute_pbctopo_evidence_connectivity_margin(
      3, triangle_edges, 0.65);
  expect(triangle.thresholded_connected && triangle.margin_defined &&
             triangle.spanning_tree_witness.size() == 2 &&
             same_edge(triangle.spanning_tree_witness[0], {0, 1, 0.9}) &&
             same_edge(triangle.spanning_tree_witness[1], {0, 2, 0.8}),
         "cycle edge was not rejected by deterministic Kruskal selection");
  expect_near(triangle.maximum_bottleneck_support, 0.8,
              "triangle bottleneck support was incorrect");

  const std::vector<PbctopoEvidenceSupportEdge> duplicate_edges{
      {1, 0, 0.3}, {0, 1, 0.9}, {1, 0, 0.6}, {1, 2, 0.8}};
  const auto duplicate = compute_pbctopo_evidence_connectivity_margin(
      3, duplicate_edges, 0.7);
  expect(duplicate.canonical_edge_count == 2 &&
             duplicate.thresholded_edge_count == 2 &&
             duplicate.thresholded_connected && duplicate.margin_defined,
         "duplicate endpoint pairs were not reduced to maximum support");
  expect_near(duplicate.maximum_bottleneck_support, 0.8,
              "duplicate maximum support did not reach the spanning tree");

  const auto singleton = compute_pbctopo_evidence_connectivity_margin(
      1, std::span<const PbctopoEvidenceSupportEdge>{}, 0.5);
  expect(singleton.thresholded_connected &&
             singleton.thresholded_component_count == 1 &&
             singleton.supplied_relation_spans && !singleton.margin_defined &&
             std::isnan(singleton.maximum_bottleneck_support) &&
             std::isnan(singleton.signed_margin) &&
             singleton.spanning_tree_witness.empty(),
         "singleton connectedness or undefined margin was incorrect");

  const std::vector<PbctopoEvidenceSupportEdge> unspanned_edges{
      {0, 1, 0.9}, {2, 3, 0.8}};
  const auto unspanned = compute_pbctopo_evidence_connectivity_margin(
      4, unspanned_edges, 0.5);
  expect(!unspanned.thresholded_connected &&
             unspanned.thresholded_component_count == 2 &&
             !unspanned.supplied_relation_spans &&
             !unspanned.margin_defined &&
             std::isnan(unspanned.maximum_bottleneck_support) &&
             std::isnan(unspanned.signed_margin) &&
             unspanned.spanning_tree_witness.empty(),
         "unspanned relation domain acquired a margin or tree witness");

  const std::vector<PbctopoEvidenceSupportEdge> boundary_edges{
      {0, 1, 0.8}, {1, 2, 0.5}};
  const auto boundary = compute_pbctopo_evidence_connectivity_margin(
      3, boundary_edges, 0.5);
  expect(boundary.thresholded_connected && boundary.margin_defined,
         "edge at the inclusive threshold was excluded");
  expect_near(boundary.signed_margin, 0.0,
              "threshold-boundary margin was not zero");

  const std::vector<PbctopoEvidenceSupportEdge> deterministic_edges{
      {0, 1, 0.9}, {0, 2, 0.8}, {1, 2, 0.8},
      {1, 3, 0.7}, {2, 3, 0.7}};
  std::vector<PbctopoEvidenceSupportEdge> reordered_edges{
      {3, 2, 0.7}, {2, 1, 0.8}, {3, 1, 0.7},
      {2, 0, 0.8}, {1, 0, 0.9}};
  const auto deterministic = compute_pbctopo_evidence_connectivity_margin(
      4, deterministic_edges, 0.65);
  const auto reordered = compute_pbctopo_evidence_connectivity_margin(
      4, reordered_edges, 0.65);
  expect(deterministic.thresholded_component_count ==
                 reordered.thresholded_component_count &&
             deterministic.maximum_bottleneck_support ==
                 reordered.maximum_bottleneck_support &&
             deterministic.signed_margin == reordered.signed_margin &&
             same_tree(deterministic.spanning_tree_witness,
                       reordered.spanning_tree_witness) &&
             same_edge(deterministic.bottleneck_witness,
                       reordered.bottleneck_witness),
         "input order or endpoint orientation changed deterministic output");

  const std::vector<PbctopoEvidenceSupportEdge> relabeled_edges{
      {3, 1, 0.9}, {3, 0, 0.8}, {1, 0, 0.8},
      {1, 2, 0.7}, {0, 2, 0.7}};
  const auto relabeled = compute_pbctopo_evidence_connectivity_margin(
      4, relabeled_edges, 0.65);
  expect(deterministic.thresholded_component_count ==
                 relabeled.thresholded_component_count &&
             deterministic.thresholded_connected ==
                 relabeled.thresholded_connected &&
             deterministic.maximum_bottleneck_support ==
                 relabeled.maximum_bottleneck_support &&
             deterministic.signed_margin == relabeled.signed_margin,
         "node relabeling changed scalar connectivity-margin results");

  expect_invalid(
      [] {
        (void)compute_pbctopo_evidence_connectivity_margin(
            0, std::span<const PbctopoEvidenceSupportEdge>{}, 0.5);
      },
      "zero-node graph was accepted");
  expect_invalid(
      [] {
        const std::vector<PbctopoEvidenceSupportEdge> edges{{0, 0, 1.0}};
        (void)compute_pbctopo_evidence_connectivity_margin(1, edges, 0.5);
      },
      "self edge was accepted");
  expect_invalid(
      [] {
        const std::vector<PbctopoEvidenceSupportEdge> edges{{0, 2, 1.0}};
        (void)compute_pbctopo_evidence_connectivity_margin(2, edges, 0.5);
      },
      "out-of-range endpoint was accepted");
  expect_invalid(
      [] {
        const std::vector<PbctopoEvidenceSupportEdge> edges{
            {0, 1, std::numeric_limits<double>::quiet_NaN()}};
        (void)compute_pbctopo_evidence_connectivity_margin(2, edges, 0.5);
      },
      "NaN support was accepted");
  expect_invalid(
      [] {
        const std::vector<PbctopoEvidenceSupportEdge> edges{
            {0, 1, std::numeric_limits<double>::infinity()}};
        (void)compute_pbctopo_evidence_connectivity_margin(2, edges, 0.5);
      },
      "infinite support was accepted");
  expect_invalid(
      [] {
        (void)compute_pbctopo_evidence_connectivity_margin(
            1, std::span<const PbctopoEvidenceSupportEdge>{},
            std::numeric_limits<double>::quiet_NaN());
      },
      "NaN threshold was accepted");
  expect_invalid(
      [] {
        (void)compute_pbctopo_evidence_connectivity_margin(
            1, std::span<const PbctopoEvidenceSupportEdge>{},
            std::numeric_limits<double>::infinity());
      },
      "infinite threshold was accepted");
  expect_overflow(
      [] {
        const std::vector<PbctopoEvidenceSupportEdge> edges{
            {0, 1, std::numeric_limits<double>::max()}};
        (void)compute_pbctopo_evidence_connectivity_margin(
            2, edges, -std::numeric_limits<double>::max());
      },
      "overflowing signed margin was accepted");

  std::cout << "PASS pbctopo_evidence_connectivity_margin_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ascii

& $Compiler -std=c++23 -O2 -Wall -Wextra -pedantic -Isrc $src `
    src/pbctopo/evidence_connectivity_margin.cpp `
    -o $exe
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile evidence connectivity-margin smoke"
}
& $exe
if ($LASTEXITCODE -ne 0) {
    throw "Evidence connectivity-margin smoke failed"
}
