#include "tabi/tabi_ses_patch_ray.h"
#include "tabi/tabi_ses_patch_row_worklist.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

PatchMappingGridSpec enclosing_grid(const PatchComplex &complex,
                                    double spacing) {
  std::array<double, 3> lower = {std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::infinity()};
  std::array<double, 3> upper = {-std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()};
  for (const Atom &atom : complex.atoms) {
    const std::array<double, 3> center = {atom.x, atom.y, atom.z};
    const double radius = atom.radius + complex.probe_radius;
    for (std::size_t dimension = 0U; dimension < 3U; ++dimension) {
      lower[dimension] = std::min(lower[dimension], center[dimension] - radius);
      upper[dimension] = std::max(upper[dimension], center[dimension] + radius);
    }
  }
  std::array<std::size_t, 3> count{};
  for (std::size_t dimension = 0U; dimension < 3U; ++dimension) {
    lower[dimension] -= 2.137 * spacing;
    upper[dimension] += 2.0 * spacing;
    count[dimension] = static_cast<std::size_t>(std::ceil(
                           (upper[dimension] - lower[dimension]) / spacing)) +
                       1U;
  }
  return {
      {lower[0], lower[1], lower[2]}, spacing, count[0], count[1], count[2]};
}

PatchComplex production_three_atom_complex() {
  const double root3 = std::sqrt(3.0);
  const std::vector<Atom> atoms = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {0.0, root3, 0.0, 0.0, 1.0}};
  PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  require(complex.status == "pass" && !complex.contacts.empty() &&
              !complex.tori.empty() && !complex.reentrant.empty(),
          "production D1 patch complex is incomplete");
  return complex;
}

std::vector<PatchRefId> canonical_patch_refs(const PatchComplex &complex) {
  std::vector<PatchRefId> refs;
  refs.reserve(complex.contacts.size() + complex.tori.size() +
               complex.reentrant.size());
  const auto append = [&](PatchKind kind, std::size_t count) {
    for (std::size_t index = 0U; index < count; ++index) {
      PatchRefId ref;
      require(PatchRefId::try_make(kind, index, ref),
              "D1 oracle PatchRefId overflow");
      refs.push_back(ref);
    }
  };
  append(PatchKind::contact, complex.contacts.size());
  append(PatchKind::toroidal, complex.tori.size());
  append(PatchKind::reentrant, complex.reentrant.size());
  return refs;
}

bool projected_row_inside(const PatchBounds &bounds,
                          const DecodedPatchRayRow &row) {
  switch (row.axis) {
  case RayAxis::x:
    return row.u >= bounds.lower.y && row.u <= bounds.upper.y &&
           row.v >= bounds.lower.z && row.v <= bounds.upper.z;
  case RayAxis::y:
    return row.u >= bounds.lower.x && row.u <= bounds.upper.x &&
           row.v >= bounds.lower.z && row.v <= bounds.upper.z;
  case RayAxis::z:
    return row.u >= bounds.lower.x && row.u <= bounds.upper.x &&
           row.v >= bounds.lower.y && row.v <= bounds.upper.y;
  }
  return false;
}

PatchBounds audit_worklist_bounds(const PatchComplex &complex, PatchRefId ref,
                                  double spacing) {
  PatchBounds bounds =
      outward_padded_patch_bounds(complex, ref.kind(), ref.index());
  const auto pad = [&](double &lower, double &upper) {
    const double extent = std::abs(upper - lower);
    const double scale =
        std::max({1.0, std::abs(lower), std::abs(upper), extent, spacing});
    const double padding =
        128.0 * std::numeric_limits<double>::epsilon() * scale;
    lower = std::nextafter(lower - padding,
                           -std::numeric_limits<double>::infinity());
    upper = std::nextafter(upper + padding,
                           std::numeric_limits<double>::infinity());
  };
  pad(bounds.lower.x, bounds.upper.x);
  pad(bounds.lower.y, bounds.upper.y);
  pad(bounds.lower.z, bounds.upper.z);
  return bounds;
}

void require_exact_brute_membership(
    const PatchComplex &complex, const PatchMappingGridSpec &spec,
    const PreparedPatchRayRowWorklist &prepared) {
  const std::vector<PatchRefId> refs = canonical_patch_refs(complex);
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    for (std::size_t row_index = 0U;
         row_index < prepared.axes[axis_index].row_count(); ++row_index) {
      DecodedPatchRayRow row;
      require(decode_patch_ray_row(spec, axis, row_index, row),
              "D1 oracle row decode failed");
      std::vector<PatchRefId> expected;
      for (const PatchRefId ref : refs) {
        const PatchBounds bounds =
            audit_worklist_bounds(complex, ref, spec.grid_spacing);
        if (projected_row_inside(bounds, row))
          expected.push_back(ref);
      }
      const std::span<const PatchRefId> actual =
          prepared.axes[axis_index].row(row_index);
      require(std::equal(expected.begin(), expected.end(), actual.begin(),
                         actual.end()),
              "D1 CSR differs from independent AABB row oracle");
    }
  }
}

std::pair<std::size_t, std::size_t>
legacy_grid_interval(double lower_bound, double upper_bound, double lower,
                     double spacing, std::size_t count) {
  if (count == 0U || !std::isfinite(lower_bound) ||
      !std::isfinite(upper_bound) || !std::isfinite(lower) ||
      !std::isfinite(spacing) || !(spacing > 0.0) || lower_bound > upper_bound)
    return {1U, 0U};
  const double grid_upper = lower + spacing * static_cast<double>(count - 1U);
  if (!std::isfinite(grid_upper) || upper_bound < lower ||
      lower_bound > grid_upper)
    return {1U, 0U};
  const double scale =
      std::max({1.0, std::abs(lower_bound), std::abs(upper_bound), spacing});
  const double padding = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  const double first = std::ceil((lower_bound - padding - lower) / spacing);
  const double last = std::floor((upper_bound + padding - lower) / spacing);
  if (!std::isfinite(first) || !std::isfinite(last) || last < 0.0 ||
      first > static_cast<double>(count - 1U))
    return {1U, 0U};
  const double clamped_first = std::max(0.0, first);
  const double clamped_last = std::min(static_cast<double>(count - 1U), last);
  if (clamped_first > clamped_last)
    return {1U, 0U};
  return {static_cast<std::size_t>(clamped_first),
          static_cast<std::size_t>(clamped_last)};
}

std::size_t
require_legacy_range_superset(const PatchComplex &complex,
                              const PatchMappingGridSpec &spec,
                              const PreparedPatchRayRowWorklist &prepared) {
  const std::vector<PatchRefId> refs = canonical_patch_refs(complex);
  std::size_t false_negatives = 0U;
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    for (const PatchRefId ref : refs) {
      const PatchBounds bounds = patch_bounds(complex, ref.kind(), ref.index());
      double lower_u = 0.0;
      double upper_u = 0.0;
      double lower_v = 0.0;
      double upper_v = 0.0;
      double grid_lower_u = 0.0;
      double grid_lower_v = 0.0;
      std::size_t u_count = 0U;
      std::size_t v_count = 0U;
      if (axis == RayAxis::x) {
        lower_u = bounds.lower.y;
        upper_u = bounds.upper.y;
        lower_v = bounds.lower.z;
        upper_v = bounds.upper.z;
        grid_lower_u = spec.lower.y;
        grid_lower_v = spec.lower.z;
        u_count = spec.ny;
        v_count = spec.nz;
      } else if (axis == RayAxis::y) {
        lower_u = bounds.lower.x;
        upper_u = bounds.upper.x;
        lower_v = bounds.lower.z;
        upper_v = bounds.upper.z;
        grid_lower_u = spec.lower.x;
        grid_lower_v = spec.lower.z;
        u_count = spec.nx;
        v_count = spec.nz;
      } else {
        lower_u = bounds.lower.x;
        upper_u = bounds.upper.x;
        lower_v = bounds.lower.y;
        upper_v = bounds.upper.y;
        grid_lower_u = spec.lower.x;
        grid_lower_v = spec.lower.y;
        u_count = spec.nx;
        v_count = spec.ny;
      }
      const auto u = legacy_grid_interval(lower_u, upper_u, grid_lower_u,
                                          spec.grid_spacing, u_count);
      const auto v = legacy_grid_interval(lower_v, upper_v, grid_lower_v,
                                          spec.grid_spacing, v_count);
      if (u.first > u.second || v.first > v.second)
        continue;
      for (std::size_t iv = v.first; iv <= v.second; ++iv) {
        for (std::size_t iu = u.first; iu <= u.second; ++iu) {
          const std::size_t row = iu + u_count * iv;
          const auto candidates = prepared.axes[axis_index].row(row);
          false_negatives += static_cast<std::size_t>(
              !std::binary_search(candidates.begin(), candidates.end(), ref));
        }
      }
    }
  }
  require(false_negatives == 0U,
          "D1 worklist is narrower than the production broad phase");
  return false_negatives;
}

struct IntersectionCoverage {
  std::size_t required_pairs{0U};
  std::size_t false_negatives{0U};
  std::size_t conservative_false_positives{0U};
  std::array<std::size_t, 3> candidate_pairs_by_kind{};
  std::array<std::size_t, 3> required_pairs_by_kind{};
  std::array<std::size_t, 3> false_negatives_by_kind{};
  std::array<std::size_t, 3> conservative_false_positives_by_kind{};
  // Diagnostic only: per-call timing includes timer overhead and is not a
  // production benchmark claim.
  std::array<double, 3> intersection_oracle_ms_by_kind{};
};

IntersectionCoverage
require_intersection_coverage(const PatchComplex &complex,
                              const PatchMappingGridSpec &spec,
                              const PreparedPatchRayRowWorklist &prepared) {
  const std::vector<PatchRefId> refs = canonical_patch_refs(complex);
  const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
  IntersectionCoverage coverage;
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    for (std::size_t row_index = 0U;
         row_index < prepared.axes[axis_index].row_count(); ++row_index) {
      DecodedPatchRayRow row;
      require(decode_patch_ray_row(spec, axis, row_index, row),
              "D1 intersection oracle row decode failed");
      const std::span<const PatchRefId> candidates =
          prepared.axes[axis_index].row(row_index);
      for (const PatchRefId ref : refs) {
        const std::size_t kind = static_cast<std::size_t>(ref.kind());
        require(kind < coverage.candidate_pairs_by_kind.size(),
                "D1 intersection oracle PatchKind is out of range");
        const auto intersection_started = std::chrono::steady_clock::now();
        const PatchAxisRoots roots = intersect_patch_axis(
            complex, probes, ref.kind(), ref.index(), axis, row.u, row.v);
        coverage.intersection_oracle_ms_by_kind[kind] +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - intersection_started)
                .count();
        const bool required = roots.count != 0U || !roots.certified;
        const bool present =
            std::binary_search(candidates.begin(), candidates.end(), ref);
        coverage.required_pairs += static_cast<std::size_t>(required);
        coverage.false_negatives +=
            static_cast<std::size_t>(required && !present);
        coverage.conservative_false_positives +=
            static_cast<std::size_t>(!required && present);
        coverage.candidate_pairs_by_kind[kind] +=
            static_cast<std::size_t>(present);
        coverage.required_pairs_by_kind[kind] +=
            static_cast<std::size_t>(required);
        coverage.false_negatives_by_kind[kind] +=
            static_cast<std::size_t>(required && !present);
        coverage.conservative_false_positives_by_kind[kind] +=
            static_cast<std::size_t>(!required && present);
      }
    }
  }
  require(coverage.false_negatives == 0U,
          "D1 projected AABB missed a certified/uncertified intersection");
  return coverage;
}

std::uint64_t worklist_checksum(const PreparedPatchRayRowWorklist &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  mix(prepared.complex_fingerprint.valid);
  for (const std::uint64_t bits : prepared.lattice_identity.scalar_bits)
    mix(bits);
  for (const std::size_t count : prepared.lattice_identity.counts)
    mix(count);
  mix(prepared.lattice_identity.valid);
  for (const PatchRayAxisRowWorklist &axis : prepared.axes) {
    mix(static_cast<std::uint64_t>(axis.axis));
    mix(axis.u_count);
    mix(axis.v_count);
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const PatchRefId ref : axis.candidates)
      mix(ref.packed());
  }
  mix(prepared.total_rows);
  mix(prepared.candidate_pairs);
  mix(prepared.nonempty_rows);
  mix(prepared.maximum_row_occupancy);
  mix(prepared.projected_aabb_worklist_complete);
  for (const char value : prepared.status)
    mix(static_cast<unsigned char>(value));
  return hash;
}

std::size_t exact_retained_bytes(const PreparedPatchRayRowWorklist &prepared) {
  std::size_t bytes = 0U;
  for (const PatchRayAxisRowWorklist &axis : prepared.axes) {
    bytes += axis.offsets.capacity() * sizeof(std::size_t);
    bytes += axis.candidates.capacity() * sizeof(PatchRefId);
  }
  return bytes;
}

std::size_t
candidate_payload_bytes(const PreparedPatchRayRowWorklist &prepared) {
  std::size_t bytes = 0U;
  for (const PatchRayAxisRowWorklist &axis : prepared.axes)
    bytes += axis.candidates.capacity() * sizeof(PatchRefId);
  return bytes;
}

void check_sphere_boundary_and_rounding() {
  const PatchComplex sphere = build_patch_complex(
      std::vector<Atom>{{0.0, 0.0, 0.0, 0.0, 1.0}}, 1.0, nullptr, false);
  require(sphere.status == "pass" && sphere.contacts.size() == 1U,
          "D1 sphere complex failed");
  const PatchBounds raw = patch_bounds(sphere, PatchKind::contact, 0U);
  const std::array<double, 4> tangent_rows = {
      raw.lower.y,
      std::nextafter(raw.lower.y, -std::numeric_limits<double>::infinity()),
      raw.upper.y,
      std::nextafter(raw.upper.y, std::numeric_limits<double>::infinity())};
  PatchRefId contact;
  require(PatchRefId::try_make(PatchKind::contact, 0U, contact),
          "D1 sphere ref failed");
  for (const double y : tangent_rows) {
    const PatchMappingGridSpec spec{{-2.0, y, 0.0}, 0.5, 9U, 1U, 1U};
    const auto worklist = build_patch_ray_row_worklist(sphere, spec);
    require(validate_patch_ray_row_worklist(sphere, spec, worklist),
            "D1 boundary worklist failed");
    require_exact_brute_membership(sphere, spec, worklist);
    (void)require_legacy_range_superset(sphere, spec, worklist);
    const auto x_candidates = worklist.axes[0].row(0U);
    require(
        std::binary_search(x_candidates.begin(), x_candidates.end(), contact),
        "D1 boundary/nextafter row lost tangent candidate");
    (void)require_intersection_coverage(sphere, spec, worklist);
  }

  const PatchComplex translated = build_patch_complex(
      std::vector<Atom>{{1.0e6, -1.0e6, 2.0e6, 0.0, 1.0}}, 1.0, nullptr, false);
  const PatchMappingGridSpec translated_spec = enclosing_grid(translated, 0.25);
  const auto translated_worklist =
      build_patch_ray_row_worklist(translated, translated_spec);
  require(validate_patch_ray_row_worklist(translated, translated_spec,
                                          translated_worklist),
          "D1 large-coordinate worklist failed");
  require_exact_brute_membership(translated, translated_spec,
                                 translated_worklist);
  (void)require_legacy_range_superset(translated, translated_spec,
                                      translated_worklist);
  (void)require_intersection_coverage(translated, translated_spec,
                                      translated_worklist);
}

void require_invalid(const PatchComplex &complex,
                     const PatchMappingGridSpec &spec,
                     const PreparedPatchRayRowWorklist &prepared,
                     const char *message) {
  require(!validate_patch_ray_row_worklist(complex, spec, prepared), message);
}

void check_mutations_and_overflow() {
  const PatchComplex sphere = build_patch_complex(
      std::vector<Atom>{{0.0, 0.0, 0.0, 0.0, 1.0}}, 1.0, nullptr, false);
  const PatchMappingGridSpec spec{{-2.0, -2.0, -2.0}, 0.5, 9U, 9U, 9U};
  const auto valid = build_patch_ray_row_worklist(sphere, spec);
  require(validate_patch_ray_row_worklist(sphere, spec, valid),
          "D1 mutation baseline failed");
  auto mutated = valid;
  ++mutated.complex_fingerprint.first;
  require_invalid(sphere, spec, mutated, "D1 stale fingerprint accepted");
  mutated = valid;
  ++mutated.lattice_identity.scalar_bits[0];
  require_invalid(sphere, spec, mutated, "D1 stale lattice accepted");
  mutated = valid;
  mutated.axes[0].offsets.back() = 0U;
  require_invalid(sphere, spec, mutated, "D1 bad CSR sentinel accepted");
  mutated = valid;
  require(mutated.axes[0].offsets.size() > 2U,
          "D1 mutation fixture offsets too small");
  mutated.axes[0].offsets[1] = mutated.axes[0].offsets[2] + 1U;
  require_invalid(sphere, spec, mutated, "D1 nonmonotone CSR accepted");
  mutated = valid;
  require(!mutated.axes[0].candidates.empty(),
          "D1 mutation fixture has no candidates");
  mutated.axes[0].candidates.front() = {};
  require_invalid(sphere, spec, mutated, "D1 invalid PatchRef accepted");
  mutated = valid;
  ++mutated.candidate_pairs;
  require_invalid(sphere, spec, mutated, "D1 stale counter accepted");
  mutated = valid;
  mutated.status = "not_built";
  require_invalid(sphere, spec, mutated, "D1 stale status accepted");
  mutated = valid;
  mutated.projected_aabb_worklist_complete = false;
  require_invalid(sphere, spec, mutated, "D1 stale completeness flag accepted");

  PatchRayRowWorklistBuildOptions limit;
  limit.candidate_pair_limit = 0U;
  const auto limited =
      build_patch_ray_row_worklist(sphere, spec, nullptr, limit);
  require(limited.status == "candidate_limit" &&
              !limited.projected_aabb_worklist_complete &&
              limited.axes[0].offsets.empty() &&
              limited.axes[0].candidates.empty(),
          "D1 candidate-pair overflow guard did not fail closed");
  const PatchMappingGridSpec overflow{
      {0.0, 0.0, 0.0}, 1.0, std::numeric_limits<std::size_t>::max(), 2U, 2U};
  const auto invalid_grid = build_patch_ray_row_worklist(sphere, overflow);
  require(invalid_grid.status == "invalid_input" &&
              !invalid_grid.projected_aabb_worklist_complete,
          "D1 row-count overflow did not fail closed");
  PatchComplex nonfinite = sphere;
  nonfinite.contacts.front().center.x =
      std::numeric_limits<double>::quiet_NaN();
  const auto invalid_geometry = build_patch_ray_row_worklist(nonfinite, spec);
  require(invalid_geometry.status == "invalid_patch_geometry" &&
              !invalid_geometry.projected_aabb_worklist_complete,
          "D1 nonfinite patch geometry did not fail closed");
}

void check_determinism_memory_and_ray_comparison() {
  PatchComplex triple = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(triple, 0.5);
  TitanExecutor executor_one(1U);
  TitanExecutor executor_eight(8U);
  const auto one = build_patch_ray_row_worklist(triple, spec, &executor_one);
  const auto eight =
      build_patch_ray_row_worklist(triple, spec, &executor_eight);
  PatchRayRowWorklistBuildOptions reverse_options;
  reverse_options.reverse_count_schedule_for_audit = true;
  const auto reversed = build_patch_ray_row_worklist(
      triple, spec, &executor_eight, reverse_options);
  require(validate_patch_ray_row_worklist(triple, spec, one) &&
              validate_patch_ray_row_worklist(triple, spec, eight) &&
              validate_patch_ray_row_worklist(triple, spec, reversed),
          "D1 production worklist validation failed");
  require(one.requested_workers == 1U && one.bounds_workers_used == 1U &&
              eight.requested_workers == 8U && eight.bounds_workers_used > 1U &&
              eight.bounds_workers_used <= eight.requested_workers &&
              reversed.bounds_workers_used > 1U,
          "D1 executor path did not use deterministic parallel bounds work");
  require_exact_brute_membership(triple, spec, one);
  const std::size_t legacy_false_negatives =
      require_legacy_range_superset(triple, spec, one);
  const IntersectionCoverage coverage =
      require_intersection_coverage(triple, spec, one);
  const std::size_t classified_candidates =
      coverage.candidate_pairs_by_kind[0] +
      coverage.candidate_pairs_by_kind[1] + coverage.candidate_pairs_by_kind[2];
  require(classified_candidates == one.candidate_pairs,
          "D1 per-kind candidate telemetry does not cover the worklist");
  const std::uint64_t checksum = worklist_checksum(one);
  require(checksum == worklist_checksum(eight) &&
              checksum == worklist_checksum(reversed),
          "D1 nt=1/8 or reversed count schedule changed CSR");

  std::array<bool, 3> saw_kind{};
  for (const PatchRayAxisRowWorklist &axis : one.axes) {
    for (const PatchRefId ref : axis.candidates)
      saw_kind[static_cast<std::size_t>(ref.kind())] = true;
  }
  require(std::all_of(saw_kind.begin(), saw_kind.end(),
                      [](bool value) { return value; }),
          "D1 production worklist missed a patch kind");

  PatchRayBuildOptions provenance_options;
  provenance_options.retain_provenance = true;
  const PatchRayGrid ray_grid =
      build_patch_ray_grid(triple, spec.lower, spec.grid_spacing, spec.nx,
                           spec.ny, spec.nz, &executor_one, provenance_options);
  require(ray_grid.status == "pass" && ray_grid.provenance.certified,
          "D1 production ray comparison grid failed");

  const PatchMappingGridSpec memory_spec{
      {-8.0, -3.0, -3.0}, 0.5, 33U, 13U, 13U};
  const PatchComplex single = build_patch_complex(
      std::vector<Atom>{{-4.0, 0.0, 0.0, 0.0, 1.0}}, 1.0, nullptr, false);
  const PatchComplex doubled = build_patch_complex(
      std::vector<Atom>{{-4.0, 0.0, 0.0, 0.0, 1.0}, {4.0, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  const auto one_patch = build_patch_ray_row_worklist(single, memory_spec);
  const auto two_patch = build_patch_ray_row_worklist(doubled, memory_spec);
  require(
      validate_patch_ray_row_worklist(single, memory_spec, one_patch) &&
          validate_patch_ray_row_worklist(doubled, memory_spec, two_patch) &&
          two_patch.candidate_pairs == 2U * one_patch.candidate_pairs &&
          candidate_payload_bytes(two_patch) ==
              2U * candidate_payload_bytes(one_patch) &&
          one_patch.retained_bytes == exact_retained_bytes(one_patch) &&
          two_patch.retained_bytes == exact_retained_bytes(two_patch) &&
          one_patch.accounted_peak_temporary_bytes >=
              one_patch.retained_bytes &&
          two_patch.accounted_peak_temporary_bytes >=
              two_patch.retained_bytes &&
          two_patch.accounted_peak_temporary_bytes <=
              2U * one_patch.accounted_peak_temporary_bytes,
      "D1 candidate or accounted memory is not linear");

  for (int warmup = 0; warmup < 2; ++warmup) {
    const auto ignored = build_patch_ray_row_worklist(triple, spec);
    require(worklist_checksum(ignored) == checksum, "D1 warmup drifted");
  }
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto timed = build_patch_ray_row_worklist(triple, spec);
    require(worklist_checksum(timed) == checksum, "D1 timed build drifted");
    sample = timed.build_ms;
  }
  std::array<double, 7> sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  std::cout
      << "d1_checksum=" << checksum << " d1_rows=" << one.total_rows
      << " d1_candidate_pairs=" << one.candidate_pairs
      << " d1_nonempty_rows=" << one.nonempty_rows
      << " d1_max_occupancy=" << one.maximum_row_occupancy
      << " d1_bounds_workers_nt8=" << eight.bounds_workers_used
      << " d1_false_negatives=" << coverage.false_negatives
      << " d1_legacy_false_negatives=" << legacy_false_negatives
      << " d1_false_positives=" << coverage.conservative_false_positives
      << " d1_required_intersections=" << coverage.required_pairs
      << " d1_contact_candidates=" << coverage.candidate_pairs_by_kind[0]
      << " d1_toroidal_candidates=" << coverage.candidate_pairs_by_kind[1]
      << " d1_reentrant_candidates=" << coverage.candidate_pairs_by_kind[2]
      << " d1_contact_required=" << coverage.required_pairs_by_kind[0]
      << " d1_toroidal_required=" << coverage.required_pairs_by_kind[1]
      << " d1_reentrant_required=" << coverage.required_pairs_by_kind[2]
      << " d1_contact_false_positives="
      << coverage.conservative_false_positives_by_kind[0]
      << " d1_toroidal_false_positives="
      << coverage.conservative_false_positives_by_kind[1]
      << " d1_reentrant_false_positives="
      << coverage.conservative_false_positives_by_kind[2]
      << " d1_contact_oracle_ms=" << coverage.intersection_oracle_ms_by_kind[0]
      << " d1_toroidal_oracle_ms=" << coverage.intersection_oracle_ms_by_kind[1]
      << " d1_reentrant_oracle_ms="
      << coverage.intersection_oracle_ms_by_kind[2]
      << " d1_worklist_bytes=" << one.retained_bytes
      << " d1_peak_bytes=" << one.accounted_peak_temporary_bytes
      << " d1_ray_raw_events=" << ray_grid.raw_events
      << " d1_ray_retained_bytes=" << ray_grid.retained_bytes
      << " d1_provenance_bytes=" << ray_grid.provenance.retained_bytes
      << " d1_memory_pairs_one=" << one_patch.candidate_pairs
      << " d1_memory_pairs_two=" << two_patch.candidate_pairs
      << " d1_memory_candidate_bytes_one=" << candidate_payload_bytes(one_patch)
      << " d1_memory_candidate_bytes_two=" << candidate_payload_bytes(two_patch)
      << " d1_samples_ms=" << samples[0] << ',' << samples[1] << ','
      << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
      << samples[5] << ',' << samples[6] << " d1_median_ms=" << sorted[3]
      << '\n';
}

} // namespace

int main() {
  try {
    check_sphere_boundary_and_rounding();
    check_mutations_and_overflow();
    check_determinism_memory_and_ray_comparison();
    std::cout << "tabi_ses_patch_row_worklist_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_row_worklist_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
