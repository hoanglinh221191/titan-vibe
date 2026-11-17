// Check of the PBCTOPO steric scan (evaluate_pbctopo_sterics_impl in
// src/pbctopo/certification.cpp) against the scan it replaced.
//
// The reference below is the previous implementation, copied verbatim: cells
// sorted with std::sort, nine binary searches per atom, every pair tested
// through the shared bookkeeping.  The new scan counting-sorts the cells into
// a dense grid and filters each run with AVX2 integer masks before the same
// per-pair work, in the same order.  On random systems (dimers, dense and
// sparse clouds, metadata-internal pairs, exclusions, a NaN point, points far
// enough to force the sorted fallback, points outside the int64 cell range,
// both scan modes) every field of PbctopoStericAudit and the four vdW pair
// counters must be identical, and the fixed-certificate audit built on the
// metadata-internal mode must report the same internal scan.  A timing block
// then alternates old and new on a 1,300-atom dimer.  Run once more with
// TITAN_FORCE_SCALAR=1 for the scalar run filter.
//
// Build: link with the engine objects of the same profile (all but Main.o).
// Usage: pbctopo_steric_scan_driver [cases] [timing_rounds]
#include "pbctopo/certification.h"
#include "pbctopo/lattice_math.h"
#include "pbctopo/steric_clash.h"
#include "pbctopo/types.h"
#include "titan_simd_dispatch.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <string>
#include <unordered_set>
#include <vector>

using titan_geom::Vec3;
using namespace titan_pbctopo;

namespace reference {

// Copied from src/pbctopo/certification.cpp before the dense-grid scan.
struct PbctopoStericCellKey {
  std::int64_t x = 0;
  std::int64_t y = 0;
  std::int64_t z = 0;

  bool operator==(const PbctopoStericCellKey &) const = default;
};

std::uint64_t steric_pair_key(std::size_t a, std::size_t b) {
  const std::size_t lo = std::min(a, b);
  const std::size_t hi = std::max(a, b);
  return (static_cast<std::uint64_t>(lo) << 32) |
         static_cast<std::uint64_t>(hi);
}

bool point_is_finite(const Vec3 &point) {
  return std::isfinite(point.x) && std::isfinite(point.y) &&
         std::isfinite(point.z);
}

PbctopoStericAudit evaluate_pbctopo_sterics_impl(
    std::span<const PbctopoAtom> atoms,
    std::span<const Vec3> placed_points,
    std::span<const std::size_t> atom_component_index,
    const std::unordered_set<std::uint64_t> &excluded_pairs,
    double fallback_cutoff, PbctopoTimingStats *timing,
    const PbctopoLatticeMetric *wrapped_lattice,
    bool metadata_internal_only) {
  PbctopoStericAudit result;
  if (atoms.empty() || placed_points.size() != atoms.size() ||
      atom_component_index.size() != atoms.size() ||
      !std::isfinite(fallback_cutoff) || fallback_cutoff <= 0.0) {
    return result;
  }
  for (std::size_t component_idx : atom_component_index) {
    if (component_idx >= atoms.size())
      return result;
  }

  double max_radius = 0.0;
  for (const auto &atom : atoms) {
    switch (atom.steric_radius_source) {
    case PbctopoStericRadiusSource::ForceFieldTopology:
    case PbctopoStericRadiusSource::ForceFieldParameter:
      ++result.force_field_radius_atoms;
      break;
    case PbctopoStericRadiusSource::ExplicitZeroLennardJones:
      ++result.explicit_zero_lj_radius_atoms;
      break;
    case PbctopoStericRadiusSource::InferredElement:
      ++result.inferred_radius_atoms;
      break;
    case PbctopoStericRadiusSource::GenericFallback:
      ++result.generic_fallback_radius_atoms;
      break;
    case PbctopoStericRadiusSource::Unavailable:
    default:
      ++result.unavailable_radius_atoms;
      break;
    }
    if (atom.has_vdw_radius && std::isfinite(atom.steric_radius)) {
      max_radius = std::max(max_radius, atom.steric_radius);
    }
  }
  const double grid_cutoff = std::max(fallback_cutoff, 2.0 * max_radius);
  if (!std::isfinite(grid_cutoff) || grid_cutoff <= 0.0)
    return result;
  result.input_valid = true;
  result.termination_reason =
      PbctopoStericTerminationReason::NumericalFailure;
  const double inv_grid = 1.0 / grid_cutoff;
  const double cell_min =
      static_cast<double>(std::numeric_limits<std::int64_t>::min() + 1);
  const double cell_max =
      static_cast<double>(std::numeric_limits<std::int64_t>::max() - 1);
  auto atom_cell = [&](std::size_t atom_idx, PbctopoStericCellKey &cell) {
    const Vec3 &point = placed_points[atom_idx];
    if (!point_is_finite(point))
      return false;
    const double cell_x = std::floor(point.x * inv_grid);
    const double cell_y = std::floor(point.y * inv_grid);
    const double cell_z = std::floor(point.z * inv_grid);
    if (cell_x < cell_min || cell_x > cell_max || cell_y < cell_min ||
        cell_y > cell_max || cell_z < cell_min || cell_z > cell_max) {
      return false;
    }
    cell = {static_cast<std::int64_t>(cell_x),
            static_cast<std::int64_t>(cell_y),
            static_cast<std::int64_t>(cell_z)};
    return true;
  };

  // The scan visits, for each atom in index order, the earlier atoms of the
  // 27 cells around it: cells in (dx, dy, dz) order, atoms of a cell in
  // index order.  That order decides ties of the running minima and which
  // clash ends the scan, so it is kept exactly.  Instead of a hash grid
  // filled as the scan goes (27 lookups and a push per atom), the atoms are
  // sorted once by (cell, index): the three cells along z around a given
  // (x, y) are then one contiguous run, found by one binary search, already
  // in dz order and index order.  An atom whose point or cell is invalid
  // stops the scan when it is reached, so it and every later atom are left
  // out of the index; nothing before it could meet them.
  struct CellEntry {
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;
    std::size_t atom;
  };
  const auto entry_less = [](const CellEntry &lhs, const CellEntry &rhs) {
    if (lhs.x != rhs.x)
      return lhs.x < rhs.x;
    if (lhs.y != rhs.y)
      return lhs.y < rhs.y;
    if (lhs.z != rhs.z)
      return lhs.z < rhs.z;
    return lhs.atom < rhs.atom;
  };
  thread_local std::vector<CellEntry> sorted_cells;
  sorted_cells.clear();
  sorted_cells.reserve(atoms.size());
  for (std::size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    PbctopoStericCellKey cell;
    if (!atom_cell(atom_idx, cell))
      break;
    sorted_cells.push_back({cell.x, cell.y, cell.z, atom_idx});
  }
  std::sort(sorted_cells.begin(), sorted_cells.end(), entry_less);

  for (std::size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const Vec3 &point = placed_points[atom_idx];
    PbctopoStericCellKey cell;
    if (!atom_cell(atom_idx, cell))
      return result;
    const std::int64_t cx = cell.x;
    const std::int64_t cy = cell.y;
    const std::int64_t cz = cell.z;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        const CellEntry first{cx + dx, cy + dy, cz - 1, 0};
        for (auto it = std::lower_bound(sorted_cells.begin(),
                                        sorted_cells.end(), first, entry_less);
             it != sorted_cells.end() && it->x == cx + dx &&
             it->y == cy + dy && it->z <= cz + 1;
             ++it) {
          const std::size_t other_idx = it->atom;
          if (other_idx >= atom_idx)
            continue;
          {
            const bool same_component =
                atom_component_index[atom_idx] ==
                atom_component_index[other_idx];
            const bool metadata_internal_pair =
                same_component && wrapped_lattice != nullptr &&
                wrapped_lattice->valid() &&
                !atoms[atom_idx].has_topology_graph &&
                !atoms[other_idx].has_topology_graph &&
                atoms[atom_idx].has_topology_key &&
                atoms[other_idx].has_topology_key;
            if (metadata_internal_only && !metadata_internal_pair)
              continue;
            if (same_component && !metadata_internal_pair)
              continue;
            // The exclusion lookup comes after the cheap tests above: a pair
            // they skip is skipped whether or not it is excluded, and most
            // neighbouring pairs lie within one component.
            if (excluded_pairs.contains(
                    steric_pair_key(atom_idx, other_idx))) {
              continue;
            }
            const Vec3 &other = placed_points[other_idx];
            const double ddx = point.x - other.x;
            const double ddy = point.y - other.y;
            const double ddz = point.z - other.z;
            const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (!std::isfinite(d2))
              return result;
            ++result.pair_checks;
            if (metadata_internal_pair)
              ++result.metadata_internal_pairs_checked;
            const auto threshold = pbctopo_pair_clash_threshold(
                atoms[atom_idx], atoms[other_idx], fallback_cutoff);
            const double margin = std::sqrt(std::max(0.0, d2)) -
                                  threshold.cutoff;
            // Two different pairs: the one closest to clashing (smallest
            // margin, which decides the certificate) and the one closest in
            // distance (min_pair_d2, which candidate ranking reads).  With
            // unequal radii they need not coincide; reporting the margin
            // pair's distance as the minimum distance was a defect.
            if (margin < result.min_margin_A) {
              result.min_margin_A = margin;
              result.min_margin_pair_d2 = d2;
              result.clash_atom_a = atom_idx;
              result.clash_atom_b = other_idx;
            }
            if (!same_component && d2 < result.min_pair_d2) {
              result.min_pair_d2 = d2;
              result.min_pair_atom_a = atom_idx;
              result.min_pair_atom_b = other_idx;
            }
            const bool placed_clash = pbctopo_atoms_clash(
                d2, atoms[atom_idx], atoms[other_idx], fallback_cutoff,
                timing);
            if (metadata_internal_pair && placed_clash) {
              const auto wrapped_nearest =
                  wrapped_lattice->nearest_image_cartesian(
                      atoms[atom_idx].x - atoms[other_idx].x,
                      atoms[atom_idx].y - atoms[other_idx].y,
                      atoms[atom_idx].z - atoms[other_idx].z);
              if (!wrapped_nearest.valid) {
                result.termination_reason =
                    PbctopoStericTerminationReason::NumericalFailure;
                return result;
              }
              const bool clash_already_present =
                  wrapped_nearest.d2 + 1.0e-12 < threshold.cutoff2;
              if (!clash_already_present) {
                result.metadata_internal_clashes = 1;
                result.hard_clashes = 1;
                result.input_valid = true;
                result.decision_complete = true;
                result.termination_reason =
                    PbctopoStericTerminationReason::FirstHardClash;
                return result;
              }
              continue;
            }
            if (!metadata_internal_pair && placed_clash) {
              result.hard_clashes = 1;
              result.input_valid = true;
              result.decision_complete = true;
              result.termination_reason =
                  PbctopoStericTerminationReason::FirstHardClash;
              return result;
            }
          }
        }
      }
    }
  }
  result.decision_complete = true;
  result.pair_enumeration_complete = true;
  result.scan_complete = true;
  result.termination_reason = PbctopoStericTerminationReason::Completed;
  result.input_valid = true;
  return result;
}

} // namespace reference

namespace {

struct Case {
  std::vector<PbctopoAtom> atoms;
  std::vector<Vec3> points;
  std::vector<std::size_t> component;
  std::unordered_set<std::uint64_t> excluded;
  double fallback_cutoff = 2.2;
  bool with_lattice = false;
  double box = 60.0;
  std::string label;
};

std::uint64_t pair_key(std::size_t a, std::size_t b) {
  const std::size_t lo = std::min(a, b);
  const std::size_t hi = std::max(a, b);
  return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
}

PbctopoLatticeMetric cubic_lattice(double length) {
  const std::array<double, 9> matrix{length, 0.0, 0.0, 0.0, length,
                                     0.0,    0.0, 0.0, length};
  return PbctopoLatticeMetric(std::span<const double, 9>(matrix));
}

// Atoms of `count` components: compact blobs (spacing about 1.5 A) placed so
// that neighbouring blobs touch, or a uniform cloud of the given density.
Case make_case(std::mt19937_64 &rng, std::size_t atoms, std::size_t components,
               const std::string &layout, double radius_scale) {
  Case c;
  c.label = layout + " n=" + std::to_string(atoms) + " k=" +
            std::to_string(components);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::normal_distribution<double> gauss(0.0, 1.0);
  c.atoms.resize(atoms);
  c.points.resize(atoms);
  c.component.resize(atoms);
  const double blob_radius = std::cbrt(static_cast<double>(atoms) /
                                       std::max<std::size_t>(1, components)) *
                             0.9;
  for (std::size_t i = 0; i < atoms; ++i) {
    const std::size_t comp =
        layout == "interleaved" ? i % components : (i * components) / atoms;
    c.component[i] = comp;
    Vec3 p{};
    if (layout == "cloud" || layout == "interleaved") {
      const double side = std::cbrt(static_cast<double>(atoms) * 11.0);
      p = {side * unit(rng), side * unit(rng), side * unit(rng)};
    } else if (layout == "lattice_dimer") {
      // protein-like density without clashes: a jittered 1.55 A lattice,
      // the two components meeting at the plane x = side / 2
      const std::size_t side = static_cast<std::size_t>(
          std::ceil(std::cbrt(static_cast<double>(atoms))));
      const std::size_t ix = i % side, iy = (i / side) % side,
                        iz = i / (side * side);
      p = {1.55 * static_cast<double>(ix) + 0.05 * gauss(rng),
           1.55 * static_cast<double>(iy) + 0.05 * gauss(rng),
           1.55 * static_cast<double>(iz) + 0.05 * gauss(rng)};
      c.component[i] = 2 * ix < side ? 0 : 1;
    } else if (layout == "sparse") {
      p = {6.0 * static_cast<double>(i % 17), 6.0 * static_cast<double>((i / 17) % 17),
           6.0 * static_cast<double>(i / 289)};
    } else {
      // blobs along x, touching
      const double cx = 2.0 * blob_radius * static_cast<double>(comp) * 0.92;
      Vec3 d{gauss(rng), gauss(rng), gauss(rng)};
      const double norm = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z) + 1e-12;
      const double r = blob_radius * std::cbrt(unit(rng));
      p = {cx + d.x / norm * r, d.y / norm * r, d.z / norm * r};
    }
    c.points[i] = p;
    PbctopoAtom &a = c.atoms[i];
    const double pick = unit(rng);
    if (pick < 0.80)
      a.steric_radius_source = PbctopoStericRadiusSource::ForceFieldTopology;
    else if (pick < 0.86)
      a.steric_radius_source = PbctopoStericRadiusSource::ForceFieldParameter;
    else if (pick < 0.90)
      a.steric_radius_source = PbctopoStericRadiusSource::ExplicitZeroLennardJones;
    else if (pick < 0.94)
      a.steric_radius_source = PbctopoStericRadiusSource::InferredElement;
    else if (pick < 0.97)
      a.steric_radius_source = PbctopoStericRadiusSource::GenericFallback;
    else
      a.steric_radius_source = PbctopoStericRadiusSource::Unavailable;
    a.has_vdw_radius =
        a.steric_radius_source != PbctopoStericRadiusSource::Unavailable;
    a.vdw_radius =
        a.steric_radius_source == PbctopoStericRadiusSource::ExplicitZeroLennardJones
            ? 0.0
            : (a.has_vdw_radius ? 1.1 + 0.9 * unit(rng) : 0.0);
    a.steric_radius = radius_scale * a.vdw_radius;
    a.has_topology_graph = unit(rng) < 0.7;
    a.has_topology_key = unit(rng) < 0.8;
    a.x = std::fmod(p.x + 1000.0, c.box);
    a.y = std::fmod(p.y + 1000.0, c.box);
    a.z = std::fmod(p.z + 1000.0, c.box);
  }
  // exclude about a tenth of the close intercomponent pairs
  for (std::size_t i = 0; i < atoms && atoms <= 1500; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      const double dx = c.points[i].x - c.points[j].x;
      const double dy = c.points[i].y - c.points[j].y;
      const double dz = c.points[i].z - c.points[j].z;
      if (dx * dx + dy * dy + dz * dz < 9.0 && unit(rng) < 0.1)
        c.excluded.insert(pair_key(i, j));
    }
  }
  return c;
}

bool same_double(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool same_audit(const PbctopoStericAudit &a, const PbctopoStericAudit &b) {
  return a.pair_checks == b.pair_checks && a.hard_clashes == b.hard_clashes &&
         a.metadata_internal_pairs_checked == b.metadata_internal_pairs_checked &&
         a.metadata_internal_clashes == b.metadata_internal_clashes &&
         a.force_field_radius_atoms == b.force_field_radius_atoms &&
         a.explicit_zero_lj_radius_atoms == b.explicit_zero_lj_radius_atoms &&
         a.inferred_radius_atoms == b.inferred_radius_atoms &&
         a.generic_fallback_radius_atoms == b.generic_fallback_radius_atoms &&
         a.unavailable_radius_atoms == b.unavailable_radius_atoms &&
         same_double(a.min_margin_A, b.min_margin_A) &&
         same_double(a.min_margin_pair_d2, b.min_margin_pair_d2) &&
         a.clash_atom_a == b.clash_atom_a && a.clash_atom_b == b.clash_atom_b &&
         same_double(a.min_pair_d2, b.min_pair_d2) &&
         a.min_pair_atom_a == b.min_pair_atom_a &&
         a.min_pair_atom_b == b.min_pair_atom_b &&
         a.decision_complete == b.decision_complete &&
         a.pair_enumeration_complete == b.pair_enumeration_complete &&
         a.termination_reason == b.termination_reason &&
         a.scan_complete == b.scan_complete && a.input_valid == b.input_valid;
}

bool same_counters(const PbctopoTimingStats &a, const PbctopoTimingStats &b) {
  return a.pbctopo_vdw_pair_checks == b.pbctopo_vdw_pair_checks &&
         a.pbctopo_vdw_radius_pairs == b.pbctopo_vdw_radius_pairs &&
         a.pbctopo_vdw_radius_fallback_pairs ==
             b.pbctopo_vdw_radius_fallback_pairs &&
         a.pbctopo_vdw_explicit_zero_lj_exempt_pairs ==
             b.pbctopo_vdw_explicit_zero_lj_exempt_pairs;
}

struct Tally {
  std::size_t cases = 0, failures = 0, clashes = 0, completed = 0,
              metadata_pairs = 0, invalid = 0;
};

void check(const Case &c, Tally &tally) {
  const PbctopoLatticeMetric lattice = cubic_lattice(c.box);
  const PbctopoLatticeMetric *wrapped = c.with_lattice ? &lattice : nullptr;
  PbctopoTimingStats ta, tb;
  const PbctopoStericAudit ref = reference::evaluate_pbctopo_sterics_impl(
      c.atoms, c.points, c.component, c.excluded, c.fallback_cutoff, &ta,
      wrapped, false);
  const PbctopoStericAudit got = evaluate_pbctopo_final_sterics(
      c.atoms, c.points, c.component, c.excluded, c.fallback_cutoff, &tb,
      wrapped);
  ++tally.cases;
  tally.clashes += ref.hard_clashes;
  tally.completed += ref.scan_complete ? 1 : 0;
  tally.invalid += ref.input_valid ? 0 : 1;
  tally.metadata_pairs += ref.metadata_internal_pairs_checked;
  if (!same_audit(ref, got) || !same_counters(ta, tb)) {
    ++tally.failures;
    std::cout << "FAIL final " << c.label << " pairs " << ref.pair_checks
              << "/" << got.pair_checks << " margin " << ref.min_margin_A
              << "/" << got.min_margin_A << " reason "
              << static_cast<int>(ref.termination_reason) << "/"
              << static_cast<int>(got.termination_reason) << "\n";
  }
  if (c.with_lattice) {
    PbctopoTimingStats tc, td;
    const PbctopoStericAudit ref_internal =
        reference::evaluate_pbctopo_sterics_impl(
            c.atoms, c.points, c.component, c.excluded, c.fallback_cutoff,
            &tc, &lattice, true);
    const PbctopoFixedCertificateAudit fixed =
        evaluate_pbctopo_fixed_certificate_invariants(
            c.atoms, c.points, c.component, {}, {}, c.excluded,
            c.fallback_cutoff, lattice, &td);
    ++tally.cases;
    const bool ok =
        !fixed.evaluated ||
        (fixed.metadata_internal_pairs_checked ==
             ref_internal.metadata_internal_pairs_checked &&
         fixed.metadata_internal_clashes ==
             ref_internal.metadata_internal_clashes &&
         fixed.internal_steric_scan_complete ==
             (ref_internal.input_valid && ref_internal.scan_complete) &&
         same_counters(tc, td));
    if (!ok) {
      ++tally.failures;
      std::cout << "FAIL internal " << c.label << " pairs "
                << ref_internal.metadata_internal_pairs_checked << "/"
                << fixed.metadata_internal_pairs_checked << "\n";
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  const std::size_t cases = argc > 1 ? std::stoul(argv[1]) : 400;
  const int rounds = argc > 2 ? std::stoi(argv[2]) : 8;
  std::cout << "avx2 " << (titan_simd_dispatch::runtime_has_avx2() ? "on" : "off")
            << "\n";
  std::mt19937_64 rng(20261011);
  Tally tally;
  const char *layouts[] = {"dimer", "cloud", "sparse", "interleaved",
                           "lattice_dimer"};
  const std::size_t sizes[] = {0, 1, 2, 7, 60, 300, 1300};
  for (std::size_t k = 0; k < cases; ++k) {
    const std::string layout = layouts[k % 5];
    const std::size_t n = sizes[(k / 5) % 7];
    const std::size_t components = 1 + (k / 35) % 4;
    const double scale = (k % 3 == 0) ? 0.55 : ((k % 3 == 1) ? 0.75 : 0.95);
    Case c = make_case(rng, n, components, layout, scale);
    c.fallback_cutoff = (k % 5 == 0) ? 1.5 : ((k % 5 == 1) ? 3.0 : 2.2);
    c.with_lattice = (k % 2) == 1;
    if (n > 3 && k % 11 == 0) {
      c.points[n / 2].y = std::numeric_limits<double>::quiet_NaN();
      c.label += " nan";
    }
    if (n > 3 && k % 13 == 0) {
      c.points[n / 3].x = 1.0e7;   // spreads the grid past its cell limit
      c.label += " far";
    }
    if (n > 3 && k % 17 == 0) {
      c.points[n - 1].z = 1.0e300; // outside the int64 cell range
      c.label += " huge";
    }
    check(c, tally);
  }
  std::cout << "cases=" << tally.cases << " failures=" << tally.failures
            << " clash_cases=" << tally.clashes
            << " completed_scans=" << tally.completed
            << " invalid_inputs=" << tally.invalid
            << " metadata_pairs_checked=" << tally.metadata_pairs << "\n";

  // Timing: a 1,300-atom dimer at protein density whose scan runs to the end
  // (radius scale 0.4 keeps every pair outside its clash distance).
  std::mt19937_64 timing_rng(7);
  Case bench = make_case(timing_rng, 1300, 2, "lattice_dimer", 0.4);
  for (PbctopoAtom &atom : bench.atoms) {
    atom.steric_radius_source = PbctopoStericRadiusSource::ForceFieldTopology;
    atom.has_vdw_radius = true;
    atom.vdw_radius = 1.5;
    atom.steric_radius = 0.5;
  }
  PbctopoTimingStats sink;
  const PbctopoStericAudit probe = evaluate_pbctopo_final_sterics(
      bench.atoms, bench.points, bench.component, bench.excluded,
      bench.fallback_cutoff, &sink, nullptr);
  std::vector<double> ratios;
  double best_old = 1e300, best_new = 1e300;
  constexpr int kCalls = 200;
  for (int round = 0; round < rounds; ++round) {
    double seconds[2] = {0.0, 0.0};
    for (const bool use_new : {round % 2 == 0, round % 2 != 0}) {
      const auto start = std::chrono::steady_clock::now();
      for (int call = 0; call < kCalls; ++call) {
        const PbctopoStericAudit audit =
            use_new ? evaluate_pbctopo_final_sterics(
                          bench.atoms, bench.points, bench.component,
                          bench.excluded, bench.fallback_cutoff, &sink, nullptr)
                    : reference::evaluate_pbctopo_sterics_impl(
                          bench.atoms, bench.points, bench.component,
                          bench.excluded, bench.fallback_cutoff, &sink,
                          nullptr, false);
        if (audit.pair_checks == 12345678)
          std::cout << "";
      }
      seconds[use_new] = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count() /
                         kCalls;
    }
    best_old = std::min(best_old, seconds[0]);
    best_new = std::min(best_new, seconds[1]);
    ratios.push_back(seconds[0] / seconds[1]);
  }
  std::sort(ratios.begin(), ratios.end());
  std::cout << "timing dimer n=1300 pair_checks=" << probe.pair_checks
            << " complete=" << probe.scan_complete
            << " best_old_ms=" << best_old * 1e3
            << " best_new_ms=" << best_new * 1e3
            << " median_round_ratio=" << ratios[ratios.size() / 2] << "\n";
  std::cout << (tally.failures == 0 ? "PASS" : "FAIL")
            << " pbctopo_steric_scan_driver\n";
  return tally.failures == 0 ? 0 : 1;
}
