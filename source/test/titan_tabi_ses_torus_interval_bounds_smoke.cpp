#include "tabi/tabi_ses_patch_ray.h"
#include "tabi/tabi_ses_torus_interval_bounds.h"
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
#include <numbers>
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
    for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
      lower[coordinate] =
          std::min(lower[coordinate], center[coordinate] - radius);
      upper[coordinate] =
          std::max(upper[coordinate], center[coordinate] + radius);
    }
  }
  std::array<std::size_t, 3> counts{};
  for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
    lower[coordinate] -= 2.137 * spacing;
    upper[coordinate] += 2.0 * spacing;
    counts[coordinate] =
        static_cast<std::size_t>(
            std::ceil((upper[coordinate] - lower[coordinate]) / spacing)) +
        1U;
  }
  return {
      {lower[0], lower[1], lower[2]}, spacing, counts[0], counts[1], counts[2]};
}

PatchComplex production_three_atom_complex() {
  const double root3 = std::sqrt(3.0);
  const std::vector<Atom> atoms = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {0.0, root3, 0.0, 0.0, 1.0}};
  PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  require(complex.status == "pass" && !complex.contacts.empty() &&
              !complex.tori.empty() && !complex.reentrant.empty(),
          "D2a production complex is incomplete");
  return complex;
}

std::vector<Atom> rotate_translate_atoms(std::span<const Atom> atoms,
                                         double angle,
                                         const Vec3 &translation) {
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const double cx = std::cos(0.43 * angle);
  const double sx = std::sin(0.43 * angle);
  std::vector<Atom> result(atoms.begin(), atoms.end());
  for (Atom &atom : result) {
    const double x = atom.x;
    const double y = atom.y;
    const double z = atom.z;
    const double rotated_x = c * x - s * y;
    const double rotated_y = s * x + c * y;
    atom.x = translation.x + rotated_x;
    atom.y = translation.y + cx * rotated_y - sx * z;
    atom.z = translation.z + sx * rotated_y + cx * z;
  }
  return result;
}

std::vector<Atom> regular_tetrahedron_atoms(double side) {
  const double root3 = std::sqrt(3.0);
  const double root6 = std::sqrt(6.0);
  return {{0.0, 0.0, 0.0, 0.0, 1.0},
          {side, 0.0, 0.0, 0.0, 1.0},
          {0.5 * side, 0.5 * root3 * side, 0.0, 0.0, 1.0},
          {0.5 * side, side / (2.0 * root3), root6 * side / 3.0, 0.0, 1.0}};
}

void require_dense_interval_enclosure(const PatchComplex &complex,
                                      std::size_t torus, std::size_t interval,
                                      double spacing) {
  PatchBounds bounds;
  require(torus_exposed_interval_bound_for_audit(complex, torus, interval,
                                                 spacing, bounds),
          "D2a interval bound construction failed");
  const ToroidalPatch &patch = complex.tori[torus];
  const AngularInterval &range = patch.exposed_intervals[interval];
  constexpr std::size_t theta_samples = 513U;
  constexpr std::size_t beta_samples = 65U;
  const std::array<long double, 3> center = {
      patch.circle_center.x, patch.circle_center.y, patch.circle_center.z};
  const std::array<long double, 3> u = {patch.basis_u.x, patch.basis_u.y,
                                        patch.basis_u.z};
  const std::array<long double, 3> v = {patch.basis_v.x, patch.basis_v.y,
                                        patch.basis_v.z};
  const std::array<long double, 3> axis = {patch.axis.x, patch.axis.y,
                                           patch.axis.z};
  const std::array<double, 3> lower = {bounds.lower.x, bounds.lower.y,
                                       bounds.lower.z};
  const std::array<double, 3> upper = {bounds.upper.x, bounds.upper.y,
                                       bounds.upper.z};
  const long double period = 2.0L * std::numbers::pi_v<long double>;
  const long double expanded_begin =
      static_cast<long double>(range.begin) - kAnalyticSesAngleTolerance;
  const long double expanded_end =
      static_cast<long double>(range.end) + kAnalyticSesAngleTolerance;
  std::array<std::array<long double, 2>, 3> segments{};
  std::size_t segment_count = 0U;
  segments[segment_count++] = {std::max(0.0L, expanded_begin),
                               std::min(period, expanded_end)};
  if (expanded_begin < 0.0L)
    segments[segment_count++] = {period + expanded_begin, period};
  if (expanded_end > period)
    segments[segment_count++] = {0.0L, expanded_end - period};
  const auto check_theta = [&](long double theta) {
    const long double ct = std::cos(theta);
    const long double st = std::sin(theta);
    for (std::size_t ib = 0U; ib < beta_samples; ++ib) {
      const long double beta =
          period * static_cast<long double>(ib) / (beta_samples - 1U);
      const long double cb = std::cos(beta);
      const long double sb = std::sin(beta);
      for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
        const long double radial = u[coordinate] * ct + v[coordinate] * st;
        const long double value =
            center[coordinate] +
            static_cast<long double>(patch.circle_radius) * radial +
            static_cast<long double>(complex.probe_radius) *
                (cb * radial + sb * axis[coordinate]);
        require(value >= static_cast<long double>(lower[coordinate]) &&
                    value <= static_cast<long double>(upper[coordinate]),
                "D2a angular-band long-double sample escaped interval bound");
      }
    }
  };
  for (std::size_t segment = 0U; segment < segment_count; ++segment) {
    for (std::size_t it = 0U; it < theta_samples; ++it) {
      const long double fraction =
          static_cast<long double>(it) / (theta_samples - 1U);
      const long double theta =
          segments[segment][0] +
          fraction * (segments[segment][1] - segments[segment][0]);
      check_theta(theta);
    }
  }
  const std::array<double, 8> boundary_diagnostics = {
      range.begin,
      std::nextafter(range.begin, -std::numeric_limits<double>::infinity()),
      std::nextafter(range.begin, std::numeric_limits<double>::infinity()),
      range.begin - kAnalyticSesAngleTolerance,
      range.end,
      std::nextafter(range.end, -std::numeric_limits<double>::infinity()),
      std::nextafter(range.end, std::numeric_limits<double>::infinity()),
      range.end + kAnalyticSesAngleTolerance};
  for (double theta : boundary_diagnostics) {
    theta = std::fmod(theta, static_cast<double>(period));
    if (theta < 0.0)
      theta += static_cast<double>(period);
    check_theta(theta);
  }
}

void check_extrema_and_boundary_intervals() {
  PatchComplex complex = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  require(complex.status == "pass" && complex.tori.size() == 1U &&
              complex.tori[0].circle_radius < complex.probe_radius,
          "D2a spindle fixture failed");
  const double period = 2.0 * std::numbers::pi;
  complex.tori[0].exposed_intervals = {
      {0.0, period},
      {0.0, 0.0},
      {0.0, std::nextafter(0.0, 1.0)},
      {0.0, 0.5 * std::numbers::pi},
      {0.5 * std::numbers::pi,
       std::nextafter(0.5 * std::numbers::pi,
                      std::numeric_limits<double>::infinity())},
      {5.5, period},
      {0.0, 0.6}};
  for (std::size_t interval = 0U;
       interval < complex.tori[0].exposed_intervals.size(); ++interval)
    require_dense_interval_enclosure(complex, 0U, interval, 0.125);
}

void require_all_interval_enclosures(const PatchComplex &complex,
                                     double spacing) {
  for (std::size_t torus = 0U; torus < complex.tori.size(); ++torus) {
    for (std::size_t interval = 0U;
         interval < complex.tori[torus].exposed_intervals.size(); ++interval)
      require_dense_interval_enclosure(complex, torus, interval, spacing);
  }
}

struct TorusCoverage {
  std::size_t full_candidates{0U};
  std::size_t interval_candidates{0U};
  std::size_t required{0U};
  std::size_t full_false_positives{0U};
  std::size_t interval_false_positives{0U};
  std::size_t false_negatives{0U};
  double oracle_ms{0.0};
};

void accumulate_coverage(TorusCoverage &aggregate,
                         const TorusCoverage &fixture) {
  aggregate.full_candidates += fixture.full_candidates;
  aggregate.interval_candidates += fixture.interval_candidates;
  aggregate.required += fixture.required;
  aggregate.full_false_positives += fixture.full_false_positives;
  aggregate.interval_false_positives += fixture.interval_false_positives;
  aggregate.false_negatives += fixture.false_negatives;
  aggregate.oracle_ms += fixture.oracle_ms;
}

TorusCoverage
require_torus_coverage(const PatchComplex &complex,
                       const PatchMappingGridSpec &spec,
                       const PreparedPatchRayRowWorklist &full,
                       const PreparedTorusIntervalRowWorklist &interval) {
  require(validate_patch_ray_row_worklist(complex, spec, full) &&
              validate_torus_interval_row_worklist(complex, spec, interval),
          "D2a coverage input validation failed");
  const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
  TorusCoverage result;
  const auto started = std::chrono::steady_clock::now();
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    require(full.axes[axis_index].row_count() ==
                interval.axes[axis_index].row_count(),
            "D2a full/interval row counts differ");
    for (std::size_t row = 0U; row < full.axes[axis_index].row_count(); ++row) {
      const auto full_candidates = full.axes[axis_index].row(row);
      const auto interval_candidates = interval.axes[axis_index].row(row);
      require(std::includes(full_candidates.begin(), full_candidates.end(),
                            interval_candidates.begin(),
                            interval_candidates.end()),
              "D2a interval worklist is not a D1 subset");
      for (const PatchRefId ref : full_candidates) {
        if (ref.kind() != PatchKind::toroidal) {
          require(std::binary_search(interval_candidates.begin(),
                                     interval_candidates.end(), ref),
                  "D2a changed contact/reentrant D1 membership");
        }
      }
      DecodedPatchRayRow decoded;
      require(decode_patch_ray_row(spec, axis, row, decoded),
              "D2a row decode failed");
      for (std::size_t torus = 0U; torus < complex.tori.size(); ++torus) {
        PatchRefId ref;
        require(PatchRefId::try_make(PatchKind::toroidal, torus, ref),
                "D2a torus ref overflow");
        const bool in_full = std::binary_search(full_candidates.begin(),
                                                full_candidates.end(), ref);
        const bool in_interval = std::binary_search(
            interval_candidates.begin(), interval_candidates.end(), ref);
        const PatchAxisRoots roots =
            intersect_patch_axis(complex, probes, PatchKind::toroidal, torus,
                                 axis, decoded.u, decoded.v);
        const bool required = roots.count != 0U || !roots.certified;
        result.full_candidates += static_cast<std::size_t>(in_full);
        result.interval_candidates += static_cast<std::size_t>(in_interval);
        result.required += static_cast<std::size_t>(required);
        result.full_false_positives +=
            static_cast<std::size_t>(in_full && !required);
        result.interval_false_positives +=
            static_cast<std::size_t>(in_interval && !required);
        result.false_negatives +=
            static_cast<std::size_t>(required && !in_interval);
      }
    }
  }
  result.oracle_ms = std::chrono::duration<double, std::milli>(
                         std::chrono::steady_clock::now() - started)
                         .count();
  require(result.false_negatives == 0U,
          "D2a interval worklist missed torus root/uncertainty pair");
  return result;
}

std::uint64_t checksum(const PreparedTorusIntervalRowWorklist &prepared) {
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
  for (const PatchRayAxisRowWorklist &axis : prepared.axes) {
    mix(static_cast<std::uint64_t>(axis.axis));
    mix(axis.u_count);
    mix(axis.v_count);
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const PatchRefId ref : axis.candidates)
      mix(ref.packed());
  }
  for (const std::size_t offset : prepared.torus_interval_offsets)
    mix(offset);
  for (const TorusExposedIntervalBound &bound :
       prepared.torus_interval_bounds) {
    mix(bound.ref.packed());
    mix(bound.interval_index);
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.x));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.y));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.z));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.x));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.y));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.z));
    mix(bound.fallback_full_torus);
  }
  for (const TorusIntervalFallbackReason reason :
       prepared.torus_fallback_reasons)
    mix(static_cast<std::uint64_t>(reason));
  mix(prepared.total_rows);
  mix(prepared.raw_rectangle_candidate_pairs);
  mix(prepared.candidate_pairs);
  mix(prepared.nonempty_rows);
  mix(prepared.maximum_row_occupancy);
  mix(prepared.exposed_interval_rectangles);
  mix(prepared.fallback_tori);
  for (const std::size_t count : prepared.fallback_reason_counts)
    mix(count);
  mix(prepared.interval_worklist_complete);
  for (const char value : prepared.status)
    mix(static_cast<unsigned char>(value));
  return hash;
}

std::size_t
exact_retained_bytes(const PreparedTorusIntervalRowWorklist &prepared) {
  std::size_t bytes = 0U;
  for (const PatchRayAxisRowWorklist &axis : prepared.axes) {
    bytes += axis.offsets.capacity() * sizeof(std::size_t);
    bytes += axis.candidates.capacity() * sizeof(PatchRefId);
  }
  bytes += prepared.torus_interval_offsets.capacity() * sizeof(std::size_t);
  bytes += prepared.torus_interval_bounds.capacity() *
           sizeof(TorusExposedIntervalBound);
  bytes += prepared.torus_fallback_reasons.capacity() *
           sizeof(TorusIntervalFallbackReason);
  return bytes;
}

void require_same_csr(const PreparedTorusIntervalRowWorklist &left,
                      const PreparedTorusIntervalRowWorklist &right,
                      const char *message) {
  require(checksum(left) == checksum(right), message);
}

void require_same_patch_membership(
    const PreparedPatchRayRowWorklist &full,
    const PreparedTorusIntervalRowWorklist &interval, PatchRefId ref,
    const char *message) {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(full.axes[axis].row_count() == interval.axes[axis].row_count(),
            message);
    for (std::size_t row = 0U; row < full.axes[axis].row_count(); ++row) {
      const auto full_row = full.axes[axis].row(row);
      const auto interval_row = interval.axes[axis].row(row);
      require(
          std::binary_search(full_row.begin(), full_row.end(), ref) ==
              std::binary_search(interval_row.begin(), interval_row.end(), ref),
          message);
    }
  }
}

void check_fixtures_determinism_and_coverage() {
  PatchComplex production = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(production, 0.5);
  TitanExecutor executor_one(1U);
  TitanExecutor executor_eight(8U);
  const auto full =
      build_patch_ray_row_worklist(production, spec, &executor_one);
  const auto interval =
      build_torus_interval_row_worklist(production, spec, &executor_one);
  const auto threaded =
      build_torus_interval_row_worklist(production, spec, &executor_eight);
  TorusIntervalRowWorklistBuildOptions reverse_options;
  reverse_options.reverse_rectangle_schedule_for_audit = true;
  const auto reversed = build_torus_interval_row_worklist(
      production, spec, &executor_eight, reverse_options);
  require(
      validate_torus_interval_row_worklist(production, spec, interval) &&
          validate_torus_interval_row_worklist(production, spec, threaded) &&
          validate_torus_interval_row_worklist(production, spec, reversed),
      "D2a production interval worklist failed validation");
  require(threaded.bounds_workers_used > 1U &&
              threaded.bounds_workers_used <= threaded.requested_workers,
          "D2a threaded bounds path did not use multiple workers");
  require_same_csr(interval, threaded, "D2a nt1/8 checksum mismatch");
  require_same_csr(interval, reversed, "D2a reverse checksum mismatch");
  const TorusCoverage coverage =
      require_torus_coverage(production, spec, full, interval);
  TorusCoverage aggregate = coverage;
  require_all_interval_enclosures(production, spec.grid_spacing);

  const std::vector<Atom> base_atoms = production.atoms;
  for (std::size_t fixture = 0U; fixture < 4U; ++fixture) {
    const Vec3 translation =
        fixture == 3U ? Vec3{1.0e6, -2.0e6, 0.5e6}
                      : Vec3{4.0 * fixture, -3.0 * fixture, 1.5 * fixture};
    const std::vector<Atom> atoms = rotate_translate_atoms(
        base_atoms, 0.271 * static_cast<double>(fixture + 1U), translation);
    const PatchComplex rotated =
        build_patch_complex(atoms, 1.0, nullptr, false);
    require(rotated.status == "pass" && !rotated.tori.empty() &&
                !rotated.reentrant.empty(),
            "D2a rotated fixture is incomplete");
    const PatchMappingGridSpec rotated_spec = enclosing_grid(rotated, 0.5);
    const auto rotated_full =
        build_patch_ray_row_worklist(rotated, rotated_spec);
    const auto rotated_interval =
        build_torus_interval_row_worklist(rotated, rotated_spec);
    require_all_interval_enclosures(rotated, rotated_spec.grid_spacing);
    accumulate_coverage(aggregate,
                        require_torus_coverage(rotated, rotated_spec,
                                               rotated_full, rotated_interval));
  }

  const PatchComplex tetra =
      build_patch_complex(regular_tetrahedron_atoms(3.2), 1.0, nullptr, false);
  require(tetra.status == "pass" && tetra.tori.size() >= 4U &&
              !tetra.reentrant.empty(),
          "D2a multi-triplet fixture is incomplete");
  const PatchMappingGridSpec tetra_spec = enclosing_grid(tetra, 0.5);
  const auto tetra_full = build_patch_ray_row_worklist(tetra, tetra_spec);
  const auto tetra_interval =
      build_torus_interval_row_worklist(tetra, tetra_spec);
  require_all_interval_enclosures(tetra, tetra_spec.grid_spacing);
  accumulate_coverage(
      aggregate,
      require_torus_coverage(tetra, tetra_spec, tetra_full, tetra_interval));

  const PatchComplex spindle = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  require(spindle.status == "pass" && spindle.tori.size() == 1U &&
              spindle.tori[0].circle_radius < spindle.probe_radius,
          "D2a spindle production fixture failed");
  const PatchMappingGridSpec spindle_spec = enclosing_grid(spindle, 0.25);
  const auto spindle_full = build_patch_ray_row_worklist(spindle, spindle_spec);
  const auto spindle_interval =
      build_torus_interval_row_worklist(spindle, spindle_spec);
  require_all_interval_enclosures(spindle, spindle_spec.grid_spacing);
  accumulate_coverage(aggregate,
                      require_torus_coverage(spindle, spindle_spec,
                                             spindle_full, spindle_interval));
  require(aggregate.false_negatives == 0U,
          "D2a aggregate fixture coverage has a false negative");

  for (int warmup = 0; warmup < 2; ++warmup) {
    const auto ignored = build_torus_interval_row_worklist(production, spec);
    require(checksum(ignored) == checksum(interval), "D2a warmup drifted");
  }
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto timed = build_torus_interval_row_worklist(production, spec);
    require(checksum(timed) == checksum(interval), "D2a timed build drifted");
    sample = timed.build_ms;
  }
  std::array<double, 7> sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  const double candidate_reduction =
      coverage.full_candidates == 0U
          ? 0.0
          : 100.0 *
                static_cast<double>(coverage.full_candidates -
                                    coverage.interval_candidates) /
                static_cast<double>(coverage.full_candidates);
  const double false_positive_reduction =
      coverage.full_false_positives == 0U
          ? 0.0
          : 100.0 *
                static_cast<double>(coverage.full_false_positives -
                                    coverage.interval_false_positives) /
                static_cast<double>(coverage.full_false_positives);
  const double aggregate_candidate_reduction =
      aggregate.full_candidates == 0U
          ? 0.0
          : 100.0 *
                static_cast<double>(aggregate.full_candidates -
                                    aggregate.interval_candidates) /
                static_cast<double>(aggregate.full_candidates);
  const double aggregate_false_positive_reduction =
      aggregate.full_false_positives == 0U
          ? 0.0
          : 100.0 *
                static_cast<double>(aggregate.full_false_positives -
                                    aggregate.interval_false_positives) /
                static_cast<double>(aggregate.full_false_positives);
  std::cout
      << "d2a_checksum=" << checksum(interval)
      << " d2a_full_torus_candidates=" << coverage.full_candidates
      << " d2a_interval_torus_candidates=" << coverage.interval_candidates
      << " d2a_required=" << coverage.required
      << " d2a_full_false_positives=" << coverage.full_false_positives
      << " d2a_interval_false_positives=" << coverage.interval_false_positives
      << " d2a_false_negatives=" << coverage.false_negatives
      << " d2a_candidate_reduction_percent=" << candidate_reduction
      << " d2a_false_positive_reduction_percent=" << false_positive_reduction
      << " d2a_interval_rectangles=" << interval.exposed_interval_rectangles
      << " d2a_interval_rectangle_bytes="
      << interval.torus_interval_bounds.capacity() *
             sizeof(TorusExposedIntervalBound)
      << " d2a_fallback_tori=" << interval.fallback_tori
      << " d2a_fallback_empty="
      << interval.fallback_reason_counts[static_cast<std::size_t>(
             TorusIntervalFallbackReason::empty_intervals)]
      << " d2a_fallback_invalid_interval="
      << interval.fallback_reason_counts[static_cast<std::size_t>(
             TorusIntervalFallbackReason::invalid_interval)]
      << " d2a_fallback_invalid_frame="
      << interval.fallback_reason_counts[static_cast<std::size_t>(
             TorusIntervalFallbackReason::invalid_frame)]
      << " d2a_fallback_sub_bound="
      << interval.fallback_reason_counts[static_cast<std::size_t>(
             TorusIntervalFallbackReason::sub_bound_failure)]
      << " d2a_raw_rectangle_pairs=" << interval.raw_rectangle_candidate_pairs
      << " d2a_bounds_workers_nt8=" << threaded.bounds_workers_used
      << " d2a_retained_bytes=" << interval.retained_bytes
      << " d2a_peak_bytes=" << interval.accounted_peak_temporary_bytes
      << " d2a_oracle_ms=" << coverage.oracle_ms
      << " d2a_aggregate_fixture_count=7"
      << " d2a_aggregate_full_torus_candidates=" << aggregate.full_candidates
      << " d2a_aggregate_interval_torus_candidates="
      << aggregate.interval_candidates
      << " d2a_aggregate_required=" << aggregate.required
      << " d2a_aggregate_full_false_positives="
      << aggregate.full_false_positives
      << " d2a_aggregate_interval_false_positives="
      << aggregate.interval_false_positives
      << " d2a_aggregate_false_negatives=" << aggregate.false_negatives
      << " d2a_aggregate_candidate_reduction_percent="
      << aggregate_candidate_reduction
      << " d2a_aggregate_false_positive_reduction_percent="
      << aggregate_false_positive_reduction
      << " d2a_aggregate_oracle_ms=" << aggregate.oracle_ms
      << " d2a_samples_ms=" << samples[0] << ',' << samples[1] << ','
      << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
      << samples[5] << ',' << samples[6] << " d2a_median_ms=" << sorted[3]
      << '\n';
}

void check_fallback_mutations_and_memory() {
  const PatchComplex production = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(production, 0.5);
  const auto valid = build_torus_interval_row_worklist(production, spec);
  require(validate_torus_interval_row_worklist(production, spec, valid) &&
              valid.retained_bytes == exact_retained_bytes(valid),
          "D2a mutation baseline failed");
  auto mutated = valid;
  ++mutated.complex_fingerprint.first;
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a stale fingerprint accepted");
  mutated = valid;
  ++mutated.lattice_identity.scalar_bits[0];
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a stale lattice accepted");
  mutated = valid;
  mutated.axes[0].offsets.back() = 0U;
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a malformed CSR accepted");
  mutated = valid;
  require(!mutated.axes[0].candidates.empty(),
          "D2a mutation candidate fixture is empty");
  mutated.axes[0].candidates.front() = {};
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a invalid ref accepted");
  mutated = valid;
  ++mutated.candidate_pairs;
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a stale counter accepted");
  mutated = valid;
  mutated.interval_worklist_complete = false;
  require(!validate_torus_interval_row_worklist(production, spec, mutated),
          "D2a stale completion accepted");

  PatchComplex empty_interval = production;
  empty_interval.tori.front().exposed_intervals.clear();
  const auto empty_fallback =
      build_torus_interval_row_worklist(empty_interval, spec);
  const auto empty_full = build_patch_ray_row_worklist(empty_interval, spec);
  require(empty_fallback.status == "pass_with_fallback" &&
              empty_fallback.fallback_tori == 1U &&
              validate_torus_interval_row_worklist(empty_interval, spec,
                                                   empty_fallback),
          "D2a empty interval did not fallback");
  PatchRefId first_torus;
  require(PatchRefId::try_make(PatchKind::toroidal, 0U, first_torus),
          "D2a fallback ref construction failed");
  require_same_patch_membership(empty_full, empty_fallback, first_torus,
                                "D2a empty fallback pruned full torus rows");
  PatchComplex malformed_interval = production;
  malformed_interval.tori.front().exposed_intervals.front().begin = 2.0;
  malformed_interval.tori.front().exposed_intervals.front().end = 1.0;
  const auto malformed_fallback =
      build_torus_interval_row_worklist(malformed_interval, spec);
  const auto malformed_full =
      build_patch_ray_row_worklist(malformed_interval, spec);
  require(malformed_fallback.status == "pass_with_fallback" &&
              malformed_fallback.fallback_tori == 1U &&
              validate_torus_interval_row_worklist(malformed_interval, spec,
                                                   malformed_fallback),
          "D2a malformed interval did not fallback");
  require_same_patch_membership(
      malformed_full, malformed_fallback, first_torus,
      "D2a malformed fallback pruned full torus rows");
  PatchComplex nonfinite_interval = production;
  nonfinite_interval.tori.front().exposed_intervals.front().begin =
      std::numeric_limits<double>::quiet_NaN();
  const auto nonfinite_fallback =
      build_torus_interval_row_worklist(nonfinite_interval, spec);
  require(nonfinite_fallback.status == "pass_with_fallback" &&
              nonfinite_fallback.fallback_tori == 1U &&
              !nonfinite_fallback.complex_fingerprint.valid &&
              validate_torus_interval_row_worklist(nonfinite_interval, spec,
                                                   nonfinite_fallback),
          "D2a nonfinite interval did not fallback");
  PatchComplex nonfinite_reference = nonfinite_interval;
  nonfinite_reference.tori.front().exposed_intervals.front() = {
      0.0, 2.0 * std::numbers::pi};
  const auto nonfinite_full =
      build_patch_ray_row_worklist(nonfinite_reference, spec);
  require_same_patch_membership(
      nonfinite_full, nonfinite_fallback, first_torus,
      "D2a nonfinite fallback pruned full torus rows");

  PatchComplex invalid_frame = production;
  invalid_frame.tori.front().basis_u.x += 1.0e-4;
  const auto frame_fallback =
      build_torus_interval_row_worklist(invalid_frame, spec);
  const auto frame_full = build_patch_ray_row_worklist(invalid_frame, spec);
  require(frame_fallback.status == "pass_with_fallback" &&
              frame_fallback.fallback_tori == 1U &&
              frame_fallback.torus_fallback_reasons.front() ==
                  TorusIntervalFallbackReason::invalid_frame &&
              validate_torus_interval_row_worklist(invalid_frame, spec,
                                                   frame_fallback),
          "D2a invalid torus frame did not fallback");
  require_same_patch_membership(
      frame_full, frame_fallback, first_torus,
      "D2a invalid-frame fallback pruned full torus rows");

  TorusIntervalRowWorklistBuildOptions limited_options;
  limited_options.raw_candidate_pair_limit = 0U;
  const auto limited = build_torus_interval_row_worklist(
      production, spec, nullptr, limited_options);
  require(limited.status == "candidate_limit" &&
              !limited.interval_worklist_complete &&
              limited.axes[0].offsets.empty(),
          "D2a candidate limit did not fail closed");

  const PatchComplex one_pair = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  const PatchComplex two_pairs =
      build_patch_complex(std::vector<Atom>{{-6.9, 0.0, 0.0, 0.0, 1.0},
                                            {-3.1, 0.0, 0.0, 0.0, 1.0},
                                            {3.1, 0.0, 0.0, 0.0, 1.0},
                                            {6.9, 0.0, 0.0, 0.0, 1.0}},
                          1.0, nullptr, false);
  const PatchMappingGridSpec memory_spec{
      {-10.0, -3.0, -3.0}, 0.5, 41U, 13U, 13U};
  const auto one = build_torus_interval_row_worklist(one_pair, memory_spec);
  const auto two = build_torus_interval_row_worklist(two_pairs, memory_spec);
  require(
      validate_torus_interval_row_worklist(one_pair, memory_spec, one) &&
          validate_torus_interval_row_worklist(two_pairs, memory_spec, two) &&
          two.torus_interval_bounds.size() ==
              2U * one.torus_interval_bounds.size() &&
          two.candidate_pairs == 2U * one.candidate_pairs &&
          two.accounted_peak_temporary_bytes <=
              2U * one.accounted_peak_temporary_bytes &&
          one.retained_bytes == exact_retained_bytes(one) &&
          two.retained_bytes == exact_retained_bytes(two),
      "D2a retained/accounted memory is not linear");
}

} // namespace

int main() {
  try {
    check_extrema_and_boundary_intervals();
    check_fixtures_determinism_and_coverage();
    check_fallback_mutations_and_memory();
    std::cout << "tabi_ses_torus_interval_bounds_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_torus_interval_bounds_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
