param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_sparse_pair_shift_table_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_sparse_pair_shift_table_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_sparse_pair_shift_table_smoke.exe'

@'
#include "pbctopo/component_search_pair_table.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_sparse_pair_shift_table_smoke: " << message
              << "\n";
    std::exit(1);
  }
}

void expect_eq_size(std::size_t actual, std::size_t expected,
                    const char *message) {
  if (actual != expected) {
    std::cerr << "FAIL pbctopo_sparse_pair_shift_table_smoke: " << message
              << ": got " << actual << " expected " << expected << "\n";
    std::exit(1);
  }
}

void expect_masks_equal(const titan_pbctopo::PbctopoShiftMask &actual,
                        const titan_pbctopo::PbctopoShiftMask &expected,
                        const char *message) {
  if (actual.shift_count() != expected.shift_count()) {
    std::cerr << "FAIL pbctopo_sparse_pair_shift_table_smoke: " << message
              << ": shift counts differ\n";
    std::exit(1);
  }
  for (std::size_t shift_idx = 0; shift_idx < actual.shift_count();
       ++shift_idx) {
    if (actual.test(shift_idx) != expected.test(shift_idx)) {
      std::cerr << "FAIL pbctopo_sparse_pair_shift_table_smoke: " << message
                << ": bit " << shift_idx << " differs\n";
      std::exit(1);
    }
  }
}

} // namespace

int main() {
  using titan_pbctopo::PbctopoComponentPairShiftEntry;
  using titan_pbctopo::PbctopoComponentPairShiftTable;
  using titan_pbctopo::PbctopoShift;
  using titan_pbctopo::PbctopoShiftMask;
  using titan_pbctopo::ClashCellKey;
  using titan_pbctopo::ClashCellKeyHash;
  using titan_pbctopo::pbctopo_enforce_pair_arc_consistency;

  PbctopoShiftMask dynamic_mask(40);
  dynamic_mask.set_all();
  expect_eq_size(dynamic_mask.count(), 40, "dynamic full mask count mismatch");
  expect(dynamic_mask.test(0), "dynamic full mask missing bit 0");
  expect(dynamic_mask.test(39), "dynamic full mask missing bit 39");

  PbctopoShiftMask only_tail(40);
  only_tail.set_only(39);
  const bool changed = dynamic_mask.intersect_with(only_tail);
  expect(changed, "intersect_with should report changed dynamic mask");
  expect_eq_size(dynamic_mask.count(), 1, "dynamic intersect count mismatch");
  expect(!dynamic_mask.test(0), "dynamic intersect kept cleared bit 0");
  expect(dynamic_mask.test(39), "dynamic intersect dropped bit 39");

  std::vector<std::size_t> visited;
  dynamic_mask.for_each_set_bit([&](std::size_t bit) { visited.push_back(bit); });
  expect_eq_size(visited.size(), 1, "dynamic iteration count mismatch");
  expect_eq_size(visited.front(), 39, "dynamic iteration order mismatch");

  PbctopoShiftMask small_mask(27);
  small_mask.set_all();
  small_mask.clear(13);
  expect_eq_size(small_mask.count(), 26, "small clear count mismatch");
  expect(!small_mask.test(13), "small clear did not clear requested bit");

  PbctopoShiftMask inline_125(125);
  inline_125.set_all();
  expect_eq_size(inline_125.count(), 125,
                 "125-shift inline mask count mismatch");
  expect_eq_size(inline_125.storage_bytes_estimate(),
                 sizeof(PbctopoShiftMask),
                 "125-shift mask unexpectedly allocated heap storage");
  PbctopoShiftMask heap_129(129);
  heap_129.set_all();
  expect(heap_129.storage_bytes_estimate() > sizeof(PbctopoShiftMask),
         "129-shift mask did not use dynamic fallback storage");

  std::unordered_map<ClashCellKey, int, ClashCellKeyHash> cell_map;
  cell_map[{0, 0, 0}] = 1;
  cell_map[{1ll << 30, 0, 0}] = 2;
  expect_eq_size(cell_map.size(), 2,
                 "wide spatial cell keys aliased distinct cells");
  ClashCellKey checked_cell;
  expect(!titan_pbctopo::pbctopo_make_cell_key(
             std::numeric_limits<double>::infinity(), 0.0, 0.0, 1.0,
             checked_cell),
         "non-finite spatial cell coordinate was accepted");
  const ClashCellKey max_cell{std::numeric_limits<std::int64_t>::max(), 0, 0};
  expect(!titan_pbctopo::pbctopo_offset_cell_key(
             max_cell, 1, 0, 0, checked_cell),
         "spatial neighbor overflow was not rejected");

  PbctopoComponentPairShiftTable table;
  table.reset(2, 40);
  expect_eq_size(table.dense_entries_equivalent(), 6400,
                 "dense equivalent mismatch");
  expect_eq_size(table.materialized_entries(), 0,
                 "new sparse table should not materialize entries");
  expect_eq_size(table.compat_mask(0, 3, 1).count(), 40,
                 "default compat mask should allow all shifts");

  const auto &implicit = table.entry_or_default(0, 3, 1, 37);
  expect(!implicit.hard_clash, "implicit sparse entry should be compatible");
  expect(std::isinf(implicit.min_pair_d2),
         "implicit sparse entry should have infinite min pair distance");
  expect_eq_size(implicit.clash_atom_a,
                 std::numeric_limits<std::size_t>::max(),
                 "implicit sparse entry clash atom mismatch");

  PbctopoComponentPairShiftEntry clash;
  clash.hard_clash = true;
  clash.min_pair_d2 = 0.25;
  clash.bbox_d2 = 0.5;
  clash.clash_atom_a = 11;
  clash.clash_atom_b = 29;
  table.store_directed_pair(0, 3, 1, 37, clash);
  table.store_directed_pair(1, 37, 0, 3, clash);

  expect(!table.compat_mask(0, 3, 1).test(37),
         "hard clash did not clear forward compat mask");
  expect(!table.compat_mask(1, 37, 0).test(3),
         "hard clash did not clear reverse compat mask");
  expect_eq_size(table.materialized_entries(), 2,
                 "materialized sparse entry count mismatch");
  expect_eq_size(table.sparse_rows(), 2, "non-empty sparse row count mismatch");
  expect(table.peak_bytes_estimate() > 0, "peak byte estimate should be positive");
  expect(table.dense_entries_equivalent() >= table.materialized_entries(),
         "dense equivalent should dominate materialized entries");

  const auto &stored = table.entry_or_default(0, 3, 1, 37);
  expect(stored.hard_clash, "stored sparse entry lost hard clash flag");
  expect_eq_size(stored.clash_atom_a, 11, "stored clash atom a mismatch");
  expect_eq_size(stored.clash_atom_b, 29, "stored clash atom b mismatch");

  PbctopoComponentPairShiftEntry compatible;
  compatible.min_pair_d2 = 4.0;
  table.store_directed_pair(0, 3, 1, 37, compatible);
  expect(table.compat_mask(0, 3, 1).test(37),
         "compatible overwrite did not restore compat mask bit");
  expect_eq_size(table.materialized_entries(), 2,
                 "overwrite should not add another sparse entry");

  std::vector<PbctopoShift> shifts(3);
  for (std::size_t idx = 0; idx < shifts.size(); ++idx) {
    shifts[idx].shift_index = idx;
    shifts[idx].ix = static_cast<int>(idx);
    shifts[idx].tx = static_cast<double>(idx);
  }
  PbctopoComponentPairShiftTable relative_table;
  relative_table.reset(2, shifts);
  relative_table.store_relative_pair(0, 1, {1, 0, 0}, clash);
  expect(relative_table.at(0, 0, 1, 1).hard_clash,
         "relative entry missing for absolute pair 0->1");
  expect(relative_table.at(0, 1, 1, 2).hard_clash,
         "equivalent absolute pair did not reuse relative entry");
  expect(!relative_table.at(0, 0, 1, 2).hard_clash,
         "different relative shift resolved to stored entry");
  const PbctopoShiftMask relative_forward_first =
      relative_table.compat_mask(0, 0, 1);
  const PbctopoShiftMask relative_forward_second =
      relative_table.compat_mask(0, 0, 1);
  expect_masks_equal(relative_forward_first, relative_forward_second,
                     "repeated relative compat calls changed result");
  expect(!relative_forward_first.test(1),
         "relative clash did not clear first equivalent pair");
  expect(relative_forward_first.test(0) && relative_forward_first.test(2),
         "relative missing entries lost default compatibility");
  expect(!relative_table.compat_mask(0, 1, 1).test(2),
         "relative clash did not clear second equivalent pair");
  expect(!relative_table.compat_mask(1, 1, 0).test(0),
         "canonical relative clash did not clear reverse orientation");
  expect_eq_size(relative_table.materialized_entries(), 1,
                 "relative table should store one canonical pair entry");
  expect_eq_size(relative_table.sparse_rows(), 1,
                 "relative table should use one triangular row");

  PbctopoComponentPairShiftTable cache_table;
  cache_table.reset(2, shifts);
  cache_table.store_relative_pair(0, 1, {1, 0, 0}, clash);
  const PbctopoShiftMask cached_before_insert =
      cache_table.compat_mask(0, 0, 1);
  expect(!cached_before_insert.test(1) && cached_before_insert.test(2),
         "initial relative cache contents mismatch");

  cache_table.store_relative_pair(0, 1, {2, 0, 0}, clash);
  const PbctopoShiftMask cached_after_insert =
      cache_table.compat_mask(0, 0, 1);
  expect(!cached_after_insert.test(1) && !cached_after_insert.test(2),
         "relative cache stayed stale after inserting an entry");
  const PbctopoShiftMask reverse_after_insert =
      cache_table.compat_mask(1, 2, 0);
  expect(!reverse_after_insert.test(0) && !reverse_after_insert.test(1) &&
             reverse_after_insert.test(2),
         "inserted clashes did not clear reverse cached orientation");

  cache_table.store_relative_pair(0, 1, {1, 0, 0}, compatible);
  const PbctopoShiftMask cached_after_overwrite =
      cache_table.compat_mask(0, 0, 1);
  expect(cached_after_overwrite.test(1) && !cached_after_overwrite.test(2),
         "relative cache stayed stale after overwriting an entry");
  const PbctopoShiftMask reverse_after_overwrite =
      cache_table.compat_mask(1, 2, 0);
  expect(!reverse_after_overwrite.test(0) && reverse_after_overwrite.test(1),
         "reverse relative cache stayed stale after overwrite");
  expect_masks_equal(cached_after_overwrite,
                     cache_table.compat_mask(0, 0, 1),
                     "rebuilt relative cache was not stable");

  PbctopoComponentPairShiftTable arc_table;
  arc_table.reset(3, shifts);
  for (int rhs_shift : {0, 2}) {
    arc_table.store_relative_pair(0, 1, {rhs_shift, 0, 0}, clash);
    arc_table.store_relative_pair(1, 0, {-rhs_shift, 0, 0}, clash);
  }
  for (int rhs_shift : {0, 1}) {
    const int delta = rhs_shift - 1;
    arc_table.store_relative_pair(1, 2, {delta, 0, 0}, clash);
    arc_table.store_relative_pair(2, 1, {-delta, 0, 0}, clash);
  }
  std::vector<PbctopoShiftMask> domains(3, PbctopoShiftMask::full(3));
  domains[0] = PbctopoShiftMask::only(3, 0);
  std::size_t arc_passes = 0;
  std::size_t arc_revisions = 0;
  std::size_t arc_removals = 0;
  std::size_t arc_queue_pushes = 0;
  expect(pbctopo_enforce_pair_arc_consistency(
             arc_table, domains, &arc_passes, &arc_revisions, &arc_removals,
             &arc_queue_pushes),
         "arc-consistent fixture was rejected");
  expect_eq_size(domains[0].count(), 1, "arc root domain changed");
  expect(domains[1].test(1) && domains[1].count() == 1,
         "arc consistency did not propagate root constraint");
  expect(domains[2].test(2) && domains[2].count() == 1,
         "arc consistency did not propagate chained constraint");
  expect(arc_passes == 1, "AC-3 should use one queue execution");
  expect_eq_size(arc_queue_pushes, 4,
                 "AC-3 should enqueue only four constrained directed arcs");
  expect_eq_size(arc_revisions, 4,
                 "AC-3 should visit each constrained directed arc once");
  expect_eq_size(arc_removals, 4,
                 "arc consistency removal count mismatch");

  std::vector<PbctopoShift> large_shifts;
  large_shifts.reserve(125);
  for (int ix = -2; ix <= 2; ++ix) {
    for (int iy = -2; iy <= 2; ++iy) {
      for (int iz = -2; iz <= 2; ++iz) {
        PbctopoShift shift;
        shift.shift_index = large_shifts.size();
        shift.ix = ix;
        shift.iy = iy;
        shift.iz = iz;
        large_shifts.push_back(shift);
      }
    }
  }
  PbctopoComponentPairShiftTable large_table;
  large_table.reset(64, large_shifts);
  expect_eq_size(large_table.materialized_entries(), 0,
                 "large sparse table should start empty");
  expect_eq_size(large_table.sparse_rows(), 0,
                 "large sparse table should not preallocate pair rows");
  expect(large_table.peak_bytes_estimate() < 128 * 1024,
         "empty large sparse table allocated pair-quadratic storage");
  expect(large_table.dense_entries_equivalent() == 64000000,
         "large dense-equivalent telemetry mismatch");

  std::vector<PbctopoShiftMask> unconstrained_domains(
      64, PbctopoShiftMask::full(large_shifts.size()));
  std::size_t unconstrained_passes = 0;
  std::size_t unconstrained_revisions = 0;
  std::size_t unconstrained_removals = 0;
  std::size_t unconstrained_pushes = 0;
  expect(pbctopo_enforce_pair_arc_consistency(
             large_table, unconstrained_domains, &unconstrained_passes,
             &unconstrained_revisions, &unconstrained_removals,
             &unconstrained_pushes),
         "unconstrained sparse table was rejected");
  expect_eq_size(unconstrained_passes, 1,
                 "unconstrained AC-3 pass telemetry mismatch");
  expect_eq_size(unconstrained_revisions, 0,
                 "unconstrained AC-3 scanned component pairs");
  expect_eq_size(unconstrained_removals, 0,
                 "unconstrained AC-3 changed a domain");
  expect_eq_size(unconstrained_pushes, 0,
                 "unconstrained AC-3 populated its queue");

  large_table.store_relative_pair(7, 41, {1, 0, 0}, clash);
  expect_eq_size(large_table.sparse_rows(), 1,
                 "one constrained pair should materialize one row");
  expect_eq_size(large_table.constrained_neighbors(7).size(), 1,
                 "lhs constrained-neighbor count mismatch");
  expect_eq_size(large_table.constrained_neighbors(41).size(), 1,
                 "rhs constrained-neighbor count mismatch");
  std::vector<PbctopoShiftMask> one_pair_domains(
      64, PbctopoShiftMask::full(large_shifts.size()));
  std::size_t one_pair_revisions = 0;
  std::size_t one_pair_pushes = 0;
  expect(pbctopo_enforce_pair_arc_consistency(
             large_table, one_pair_domains, nullptr, &one_pair_revisions,
             nullptr, &one_pair_pushes),
         "single-pair sparse table was rejected");
  expect_eq_size(one_pair_pushes, 2,
                 "single constrained pair should enqueue two arcs");
  expect_eq_size(one_pair_revisions, 2,
                 "single constrained pair should revise two arcs");

  std::cout << "PASS pbctopo_sparse_pair_shift_table_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\pbctopo\component_search_pair_table.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile pbctopo sparse pair-shift table smoke'
}

& $exe
if ($LASTEXITCODE -ne 0) {
    throw "pbctopo sparse pair-shift table smoke failed with exit code $LASTEXITCODE"
}
