#include "tabi/tabi_ses_patch_ray.h"
#include "tabi/tabi_ses_torus_domain_bounds.h"
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
#include <span>
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
          "D2b production complex is incomplete");
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
    const double rx = c * x - s * y;
    const double ry = s * x + c * y;
    atom.x = translation.x + rx;
    atom.y = translation.y + cx * ry - sx * z;
    atom.z = translation.z + sx * ry + cx * z;
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

double normalize_angle(double angle) {
  const double period = 2.0 * std::numbers::pi;
  angle = std::fmod(angle, period);
  return angle < 0.0 ? angle + period : angle;
}

double ccw_distance(double from, double to) {
  const double period = 2.0 * std::numbers::pi;
  double delta = normalize_angle(to) - normalize_angle(from);
  if (delta < 0.0)
    delta += period;
  return delta;
}

void check_directed_arc_envelope() {
  const double period = 2.0 * std::numbers::pi;
  const std::array<std::pair<double, double>, 8> fixed = {
      std::pair{0.0, 0.0},
      std::pair{0.0, std::nextafter(0.0, 1.0)},
      std::pair{0.0, std::numbers::pi},
      std::pair{0.0, std::nextafter(std::numbers::pi, 0.0)},
      std::pair{0.0, std::nextafter(std::numbers::pi, period)},
      std::pair{6.1, 0.2},
      std::pair{0.2, 6.1},
      std::pair{-0.3, period + 0.7}};
  const auto check = [&](double endpoint_a, double endpoint_b) {
    TorusBetaMinorArc arc;
    require(directed_minor_arc_from_endpoints_for_audit(endpoint_a, endpoint_b,
                                                        arc),
            "D2b directed arc construction failed");
    const double forward = ccw_distance(endpoint_a, endpoint_b);
    const double expected_start =
        normalize_angle(forward <= std::numbers::pi ? endpoint_a : endpoint_b);
    const double expected_span =
        forward <= std::numbers::pi ? forward : period - forward;
    require(std::bit_cast<std::uint64_t>(arc.start) ==
                    std::bit_cast<std::uint64_t>(expected_start) &&
                std::bit_cast<std::uint64_t>(arc.span) ==
                    std::bit_cast<std::uint64_t>(expected_span),
            "D2b directed arc branch drifted");
    const std::array<double, 8> probes = {
        arc.start,
        std::nextafter(arc.start, -std::numeric_limits<double>::infinity()),
        normalize_angle(arc.start + 0.5 * arc.span),
        normalize_angle(arc.start + arc.span),
        normalize_angle(arc.start + arc.span + kAnalyticSesAngleTolerance),
        normalize_angle(
            std::nextafter(arc.start + arc.span + kAnalyticSesAngleTolerance,
                           std::numeric_limits<double>::infinity())),
        normalize_angle(arc.start - 0.5),
        normalize_angle(arc.start + arc.span + 0.5)};
    for (const double probe : probes) {
      const bool expected_envelope = ccw_distance(arc.start, probe) <=
                                     arc.span + kAnalyticSesAngleTolerance;
      require(angle_on_torus_beta_minor_arc_for_audit(probe, arc) ==
                  expected_envelope,
              "D2b directed audit envelope drifted");
      const bool on_forward = ccw_distance(endpoint_a, probe) <=
                              forward + kAnalyticSesAngleTolerance;
      const bool expected_production =
          forward <= std::numbers::pi ? on_forward : !on_forward;
      const bool production =
          angle_on_torus_beta_production_minor_arc_for_audit(probe, endpoint_a,
                                                             endpoint_b);
      require(production == expected_production,
              "D2b production minor-arc predicate changed");
      require(!production ||
                  angle_on_torus_beta_minor_arc_for_audit(probe, arc),
              "D2b audit envelope omitted a production-accepted beta");
    }
    if (forward > std::numbers::pi) {
      const double a = normalize_angle(endpoint_a);
      const double b = normalize_angle(endpoint_b);
      const double just_after_b =
          normalize_angle(b + 0.5 * kAnalyticSesAngleTolerance);
      require(!angle_on_torus_beta_production_minor_arc_for_audit(a, endpoint_a,
                                                                  endpoint_b) &&
                  !angle_on_torus_beta_production_minor_arc_for_audit(
                      b, endpoint_a, endpoint_b) &&
                  !angle_on_torus_beta_production_minor_arc_for_audit(
                      just_after_b, endpoint_a, endpoint_b),
              "D2b reverse-branch production endpoint semantics changed");
      require(angle_on_torus_beta_minor_arc_for_audit(a, arc) &&
                  angle_on_torus_beta_minor_arc_for_audit(b, arc),
              "D2b conservative reverse envelope did not close endpoints");
    }
  };
  for (const auto &[a, b] : fixed)
    check(a, b);
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  for (std::size_t sample = 0U; sample < 512U; ++sample) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    const double a = period * static_cast<double>(state >> 11U) /
                     static_cast<double>(std::uint64_t{1} << 53U);
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    const double b = period * static_cast<double>(state >> 11U) /
                     static_cast<double>(std::uint64_t{1} << 53U);
    check(a, b);
  }
}

void check_beta_extrema_boundaries() {
  const double period = 2.0 * std::numbers::pi;
  const std::array<TorusBetaMinorArc, 6> arcs = {
      TorusBetaMinorArc{0.0, 0.0},
      TorusBetaMinorArc{0.0, std::nextafter(0.0, 1.0)},
      TorusBetaMinorArc{0.37, std::nextafter(std::numbers::pi, 0.0)},
      TorusBetaMinorArc{0.37, std::numbers::pi},
      TorusBetaMinorArc{std::nextafter(period, 0.0), 0.25},
      TorusBetaMinorArc{6.1, 0.4}};
  for (const TorusBetaMinorArc &arc : arcs) {
    double cos_minimum = 0.0;
    double cos_maximum = 0.0;
    double sin_minimum = 0.0;
    double sin_maximum = 0.0;
    require(torus_beta_trigonometric_ranges_for_audit(
                arc, cos_minimum, cos_maximum, sin_minimum, sin_maximum),
            "D2b beta extrema construction failed");
    constexpr std::size_t samples = 4097U;
    const long double accepted_span =
        static_cast<long double>(arc.span) + kAnalyticSesAngleTolerance;
    const auto check = [&](long double angle) {
      const long double cosine = std::cos(angle);
      const long double sine = std::sin(angle);
      require(cosine >= static_cast<long double>(cos_minimum) &&
                  cosine <= static_cast<long double>(cos_maximum) &&
                  sine >= static_cast<long double>(sin_minimum) &&
                  sine <= static_cast<long double>(sin_maximum),
              "D2b beta long-double sample escaped extrema range");
    };
    for (std::size_t sample = 0U; sample < samples; ++sample) {
      const long double angle =
          static_cast<long double>(arc.start) +
          accepted_span * static_cast<long double>(sample) / (samples - 1U);
      check(angle);
    }
    check(arc.start);
    check(static_cast<long double>(arc.start) + arc.span);
    check(static_cast<long double>(arc.start) + arc.span +
          kAnalyticSesAngleTolerance);
    require(
        !angle_on_torus_beta_minor_arc_for_audit(
            std::nextafter(arc.start, -std::numeric_limits<double>::infinity()),
            arc),
        "D2b beta extrema fixture expanded before start");
  }
}

void require_dense_domain_enclosure(const PatchComplex &complex,
                                    std::size_t torus, std::size_t interval,
                                    double spacing) {
  PatchBounds bounds;
  require(torus_domain_interval_bound_for_audit(complex, torus, interval,
                                                spacing, bounds),
          "D2b domain bound construction failed");
  TorusBetaMinorArc beta;
  require(torus_beta_minor_arc_for_audit(complex, torus, beta),
          "D2b beta arc construction failed");
  const ToroidalPatch &patch = complex.tori[torus];
  const AngularInterval &theta_range = patch.exposed_intervals[interval];
  constexpr std::size_t theta_samples = 129U;
  constexpr std::size_t beta_samples = 129U;
  const long double period = 2.0L * std::numbers::pi_v<long double>;
  const long double theta_begin =
      static_cast<long double>(theta_range.begin) - kAnalyticSesAngleTolerance;
  const long double theta_end =
      static_cast<long double>(theta_range.end) + kAnalyticSesAngleTolerance;
  std::array<std::array<long double, 2>, 3> theta_segments{};
  std::size_t theta_segment_count = 0U;
  theta_segments[theta_segment_count++] = {std::max(0.0L, theta_begin),
                                           std::min(period, theta_end)};
  if (theta_begin < 0.0L)
    theta_segments[theta_segment_count++] = {period + theta_begin, period};
  if (theta_end > period)
    theta_segments[theta_segment_count++] = {0.0L, theta_end - period};

  const long double beta_start = beta.start;
  const long double beta_end =
      beta_start + beta.span + kAnalyticSesAngleTolerance;
  std::array<std::array<long double, 2>, 2> beta_segments{};
  std::size_t beta_segment_count = 0U;
  if (beta_end <= period) {
    beta_segments[beta_segment_count++] = {beta_start, beta_end};
  } else {
    beta_segments[beta_segment_count++] = {beta_start, period};
    beta_segments[beta_segment_count++] = {0.0L, beta_end - period};
  }

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
  const auto check_point = [&](long double theta, long double beta_value) {
    const long double ct = std::cos(theta);
    const long double st = std::sin(theta);
    const long double cb = std::cos(beta_value);
    const long double sb = std::sin(beta_value);
    for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
      const long double radial = u[coordinate] * ct + v[coordinate] * st;
      const long double value =
          center[coordinate] +
          (static_cast<long double>(patch.circle_radius) +
           static_cast<long double>(complex.probe_radius) * cb) *
              radial +
          static_cast<long double>(complex.probe_radius) * sb *
              axis[coordinate];
      require(value >= static_cast<long double>(lower[coordinate]) &&
                  value <= static_cast<long double>(upper[coordinate]),
              "D2b long-double sample escaped domain bound");
    }
  };
  for (std::size_t ts = 0U; ts < theta_segment_count; ++ts) {
    for (std::size_t bs = 0U; bs < beta_segment_count; ++bs) {
      for (std::size_t it = 0U; it < theta_samples; ++it) {
        const long double theta =
            theta_segments[ts][0] +
            static_cast<long double>(it) / (theta_samples - 1U) *
                (theta_segments[ts][1] - theta_segments[ts][0]);
        for (std::size_t ib = 0U; ib < beta_samples; ++ib) {
          const long double beta_value =
              beta_segments[bs][0] +
              static_cast<long double>(ib) / (beta_samples - 1U) *
                  (beta_segments[bs][1] - beta_segments[bs][0]);
          check_point(theta, beta_value);
        }
      }
    }
  }
  const double terminal = normalize_angle(beta.start + beta.span);
  require(
      !angle_on_torus_beta_minor_arc_for_audit(
          std::nextafter(beta.start, -std::numeric_limits<double>::infinity()),
          beta),
      "D2b beta predicate expanded before directed start");
  for (const double boundary :
       {beta.start, terminal,
        normalize_angle(beta.start + beta.span + kAnalyticSesAngleTolerance)}) {
    require(angle_on_torus_beta_minor_arc_for_audit(boundary, beta),
            "D2b beta boundary predicate rejected accepted point");
  }
}

void require_all_domain_enclosures(const PatchComplex &complex,
                                   double spacing) {
  for (std::size_t torus = 0U; torus < complex.tori.size(); ++torus) {
    for (std::size_t interval = 0U;
         interval < complex.tori[torus].exposed_intervals.size(); ++interval)
      require_dense_domain_enclosure(complex, torus, interval, spacing);
  }
}

struct Coverage {
  std::size_t full_candidates{0U};
  std::size_t theta_candidates{0U};
  std::size_t domain_candidates{0U};
  std::size_t required{0U};
  std::size_t full_false_positives{0U};
  std::size_t theta_false_positives{0U};
  std::size_t domain_false_positives{0U};
  std::size_t false_negatives{0U};
  double oracle_ms{0.0};
};

void accumulate(Coverage &target, const Coverage &source) {
  target.full_candidates += source.full_candidates;
  target.theta_candidates += source.theta_candidates;
  target.domain_candidates += source.domain_candidates;
  target.required += source.required;
  target.full_false_positives += source.full_false_positives;
  target.theta_false_positives += source.theta_false_positives;
  target.domain_false_positives += source.domain_false_positives;
  target.false_negatives += source.false_negatives;
  target.oracle_ms += source.oracle_ms;
}

Coverage require_coverage(const PatchComplex &complex,
                          const PatchMappingGridSpec &spec,
                          const PreparedPatchRayRowWorklist &full,
                          const PreparedTorusIntervalRowWorklist &theta,
                          const PreparedTorusDomainRowWorklist &domain) {
  require(validate_patch_ray_row_worklist(complex, spec, full) &&
              validate_torus_interval_row_worklist(complex, spec, theta) &&
              validate_torus_domain_row_worklist(complex, spec, domain),
          "D2b coverage inputs failed validation");
  const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
  Coverage result;
  const auto started = std::chrono::steady_clock::now();
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    for (std::size_t row = 0U; row < full.axes[axis_index].row_count(); ++row) {
      const auto full_row = full.axes[axis_index].row(row);
      const auto theta_row = theta.axes[axis_index].row(row);
      const auto domain_row = domain.axes[axis_index].row(row);
      require(std::includes(full_row.begin(), full_row.end(), theta_row.begin(),
                            theta_row.end()) &&
                  std::includes(theta_row.begin(), theta_row.end(),
                                domain_row.begin(), domain_row.end()),
              "D2b subset chain D2b subset D2a subset D1 failed");
      for (const PatchRefId ref : theta_row) {
        if (ref.kind() != PatchKind::toroidal) {
          require(std::binary_search(domain_row.begin(), domain_row.end(), ref),
                  "D2b changed contact/reentrant membership");
        }
      }
      DecodedPatchRayRow decoded;
      require(decode_patch_ray_row(spec, axis, row, decoded),
              "D2b row decode failed");
      for (std::size_t torus = 0U; torus < complex.tori.size(); ++torus) {
        PatchRefId ref;
        require(PatchRefId::try_make(PatchKind::toroidal, torus, ref),
                "D2b torus ref overflow");
        const bool in_full =
            std::binary_search(full_row.begin(), full_row.end(), ref);
        const bool in_theta =
            std::binary_search(theta_row.begin(), theta_row.end(), ref);
        const bool in_domain =
            std::binary_search(domain_row.begin(), domain_row.end(), ref);
        const PatchAxisRoots roots =
            intersect_patch_axis(complex, probes, PatchKind::toroidal, torus,
                                 axis, decoded.u, decoded.v);
        const bool required = roots.count != 0U || !roots.certified;
        result.full_candidates += in_full;
        result.theta_candidates += in_theta;
        result.domain_candidates += in_domain;
        result.required += required;
        result.full_false_positives += in_full && !required;
        result.theta_false_positives += in_theta && !required;
        result.domain_false_positives += in_domain && !required;
        result.false_negatives += required && !in_domain;
      }
    }
  }
  result.oracle_ms = std::chrono::duration<double, std::milli>(
                         std::chrono::steady_clock::now() - started)
                         .count();
  require(result.false_negatives == 0U,
          "D2b domain worklist missed torus root/uncertainty pair");
  return result;
}

std::uint64_t checksum(const PreparedTorusDomainRowWorklist &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  for (const auto &axis : prepared.axes) {
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const PatchRefId ref : axis.candidates)
      mix(ref.packed());
  }
  for (const std::size_t offset : prepared.torus_domain_offsets)
    mix(offset);
  for (const TorusDomainIntervalBound &bound : prepared.torus_domain_bounds) {
    mix(bound.ref.packed());
    mix(bound.interval_index);
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.x));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.y));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.lower.z));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.x));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.y));
    mix(std::bit_cast<std::uint64_t>(bound.bounds.upper.z));
  }
  for (const TorusDomainFallbackReason reason : prepared.torus_fallback_reasons)
    mix(static_cast<std::uint64_t>(reason));
  mix(prepared.total_rows);
  mix(prepared.theta_candidate_pairs);
  mix(prepared.candidate_pairs);
  mix(prepared.nonempty_rows);
  mix(prepared.maximum_row_occupancy);
  mix(prepared.domain_rectangles);
  mix(prepared.beta_fallback_to_theta_tori);
  mix(prepared.full_fallback_tori);
  for (const std::size_t count : prepared.fallback_reason_counts)
    mix(count);
  mix(prepared.domain_worklist_complete);
  for (const char value : prepared.status)
    mix(static_cast<unsigned char>(value));
  return hash;
}

std::size_t exact_retained(const PreparedTorusDomainRowWorklist &prepared) {
  std::size_t bytes = 0U;
  for (const auto &axis : prepared.axes) {
    bytes += axis.offsets.capacity() * sizeof(std::size_t);
    bytes += axis.candidates.capacity() * sizeof(PatchRefId);
  }
  bytes += prepared.torus_domain_offsets.capacity() * sizeof(std::size_t);
  bytes += prepared.torus_domain_bounds.capacity() *
           sizeof(TorusDomainIntervalBound);
  bytes += prepared.torus_fallback_reasons.capacity() *
           sizeof(TorusDomainFallbackReason);
  return bytes;
}

void check_fixtures_and_telemetry() {
  PatchComplex production = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(production, 0.5);
  TitanExecutor one(1U);
  TitanExecutor eight(8U);
  const auto full = build_patch_ray_row_worklist(production, spec, &one);
  const auto theta = build_torus_interval_row_worklist(production, spec, &one);
  const auto domain = build_torus_domain_row_worklist(production, spec, &one);
  const auto threaded =
      build_torus_domain_row_worklist(production, spec, &eight);
  TorusDomainRowWorklistBuildOptions reverse_options;
  reverse_options.reverse_torus_schedule_for_audit = true;
  const auto reversed = build_torus_domain_row_worklist(
      production, spec, &eight, reverse_options);
  require(validate_torus_domain_row_worklist(production, spec, domain) &&
              validate_torus_domain_row_worklist(production, spec, threaded) &&
              validate_torus_domain_row_worklist(production, spec, reversed) &&
              checksum(domain) == checksum(threaded) &&
              checksum(domain) == checksum(reversed),
          "D2b nt1/8/reverse determinism failed");
  require(threaded.bounds_workers_used > 1U,
          "D2b threaded bounds path did not use multiple workers");
  Coverage production_coverage =
      require_coverage(production, spec, full, theta, domain);
  Coverage aggregate = production_coverage;
  require_all_domain_enclosures(production, spec.grid_spacing);

  const std::vector<Atom> base_atoms = production.atoms;
  for (std::size_t fixture = 0U; fixture < 4U; ++fixture) {
    const Vec3 translation =
        fixture == 3U ? Vec3{1.0e6, -2.0e6, 0.5e6}
                      : Vec3{4.0 * fixture, -3.0 * fixture, 1.5 * fixture};
    const PatchComplex rotated = build_patch_complex(
        rotate_translate_atoms(base_atoms, 0.271 * (fixture + 1U), translation),
        1.0, nullptr, false);
    require(rotated.status == "pass" && !rotated.tori.empty(),
            "D2b rotated fixture is incomplete");
    const PatchMappingGridSpec rotated_spec = enclosing_grid(rotated, 0.5);
    const auto rotated_full =
        build_patch_ray_row_worklist(rotated, rotated_spec);
    const auto rotated_theta =
        build_torus_interval_row_worklist(rotated, rotated_spec);
    const auto rotated_domain =
        build_torus_domain_row_worklist(rotated, rotated_spec);
    require_all_domain_enclosures(rotated, rotated_spec.grid_spacing);
    accumulate(aggregate, require_coverage(rotated, rotated_spec, rotated_full,
                                           rotated_theta, rotated_domain));
  }

  const PatchComplex tetra =
      build_patch_complex(regular_tetrahedron_atoms(3.2), 1.0, nullptr, false);
  require(tetra.status == "pass" && tetra.tori.size() >= 4U,
          "D2b tetra fixture is incomplete");
  const PatchMappingGridSpec tetra_spec = enclosing_grid(tetra, 0.5);
  const auto tetra_full = build_patch_ray_row_worklist(tetra, tetra_spec);
  const auto tetra_theta = build_torus_interval_row_worklist(tetra, tetra_spec);
  const auto tetra_domain = build_torus_domain_row_worklist(tetra, tetra_spec);
  require_all_domain_enclosures(tetra, tetra_spec.grid_spacing);
  accumulate(aggregate, require_coverage(tetra, tetra_spec, tetra_full,
                                         tetra_theta, tetra_domain));

  const PatchComplex spindle = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  require(spindle.status == "pass" && spindle.tori.size() == 1U &&
              spindle.tori[0].circle_radius < spindle.probe_radius,
          "D2b spindle fixture failed");
  const PatchMappingGridSpec spindle_spec = enclosing_grid(spindle, 0.25);
  const auto spindle_full = build_patch_ray_row_worklist(spindle, spindle_spec);
  const auto spindle_theta =
      build_torus_interval_row_worklist(spindle, spindle_spec);
  const auto spindle_domain =
      build_torus_domain_row_worklist(spindle, spindle_spec);
  require_all_domain_enclosures(spindle, spindle_spec.grid_spacing);
  accumulate(aggregate, require_coverage(spindle, spindle_spec, spindle_full,
                                         spindle_theta, spindle_domain));

  const PatchComplex unequal =
      build_patch_complex(std::vector<Atom>{{-1.55, 0.0, 0.0, 0.0, 0.8},
                                            {1.55, 0.0, 0.0, 0.0, 1.35}},
                          1.0, nullptr, false);
  require(unequal.status == "pass" && unequal.tori.size() == 1U,
          "D2b unequal-radii fixture failed");
  const PatchMappingGridSpec unequal_spec = enclosing_grid(unequal, 0.25);
  const auto unequal_full = build_patch_ray_row_worklist(unequal, unequal_spec);
  const auto unequal_theta =
      build_torus_interval_row_worklist(unequal, unequal_spec);
  const auto unequal_domain =
      build_torus_domain_row_worklist(unequal, unequal_spec);
  require_all_domain_enclosures(unequal, unequal_spec.grid_spacing);
  accumulate(aggregate, require_coverage(unequal, unequal_spec, unequal_full,
                                         unequal_theta, unequal_domain));
  require(aggregate.false_negatives == 0U,
          "D2b aggregate coverage has a false negative");

  for (int warmup = 0; warmup < 2; ++warmup)
    require(checksum(build_torus_domain_row_worklist(production, spec)) ==
                checksum(domain),
            "D2b warmup drifted");
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto timed = build_torus_domain_row_worklist(production, spec);
    require(checksum(timed) == checksum(domain), "D2b timed build drifted");
    sample = timed.build_ms;
  }
  auto sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  const auto reduction = [](std::size_t baseline, std::size_t reduced) {
    return baseline == 0U ? 0.0
                          : 100.0 * static_cast<double>(baseline - reduced) /
                                static_cast<double>(baseline);
  };
  std::cout
      << "d2b_checksum=" << checksum(domain)
      << " d2b_production_d1_candidates=" << production_coverage.full_candidates
      << " d2b_production_d2a_candidates="
      << production_coverage.theta_candidates
      << " d2b_production_d2b_candidates="
      << production_coverage.domain_candidates
      << " d2b_production_required=" << production_coverage.required
      << " d2b_production_d1_false_positives="
      << production_coverage.full_false_positives
      << " d2b_production_d2a_false_positives="
      << production_coverage.theta_false_positives
      << " d2b_production_d2b_false_positives="
      << production_coverage.domain_false_positives
      << " d2b_production_false_negatives="
      << production_coverage.false_negatives
      << " d2b_production_incremental_reduction_percent="
      << reduction(production_coverage.theta_candidates,
                   production_coverage.domain_candidates)
      << " d2b_production_total_reduction_percent="
      << reduction(production_coverage.full_candidates,
                   production_coverage.domain_candidates)
      << " d2b_aggregate_fixture_count=8"
      << " d2b_aggregate_d1_candidates=" << aggregate.full_candidates
      << " d2b_aggregate_d2a_candidates=" << aggregate.theta_candidates
      << " d2b_aggregate_d2b_candidates=" << aggregate.domain_candidates
      << " d2b_aggregate_required=" << aggregate.required
      << " d2b_aggregate_d1_false_positives=" << aggregate.full_false_positives
      << " d2b_aggregate_d2a_false_positives="
      << aggregate.theta_false_positives
      << " d2b_aggregate_d2b_false_positives="
      << aggregate.domain_false_positives
      << " d2b_aggregate_false_negatives=" << aggregate.false_negatives
      << " d2b_aggregate_incremental_reduction_percent="
      << reduction(aggregate.theta_candidates, aggregate.domain_candidates)
      << " d2b_aggregate_total_reduction_percent="
      << reduction(aggregate.full_candidates, aggregate.domain_candidates)
      << " d2b_domain_rectangles=" << domain.domain_rectangles
      << " d2b_beta_to_theta_fallback=" << domain.beta_fallback_to_theta_tori
      << " d2b_full_fallback=" << domain.full_fallback_tori
      << " d2b_rectangle_bytes="
      << domain.torus_domain_bounds.capacity() *
             sizeof(TorusDomainIntervalBound)
      << " d2b_retained_bytes=" << domain.retained_bytes
      << " d2b_peak_bytes=" << domain.accounted_peak_temporary_bytes
      << " d2b_oracle_ms=" << production_coverage.oracle_ms
      << " d2b_aggregate_oracle_ms=" << aggregate.oracle_ms
      << " d2b_bounds_workers_nt8=" << threaded.bounds_workers_used
      << " d2b_samples_ms=" << samples[0] << ',' << samples[1] << ','
      << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
      << samples[5] << ',' << samples[6] << " d2b_median_ms=" << sorted[3]
      << '\n';
}

void check_fallback_mutation_and_memory() {
  const PatchComplex production = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(production, 0.5);
  const auto valid = build_torus_domain_row_worklist(production, spec);
  require(validate_torus_domain_row_worklist(production, spec, valid) &&
              valid.retained_bytes == exact_retained(valid),
          "D2b mutation baseline failed");
  auto mutated = valid;
  ++mutated.complex_fingerprint.first;
  require(!validate_torus_domain_row_worklist(production, spec, mutated),
          "D2b stale fingerprint accepted");
  mutated = valid;
  ++mutated.lattice_identity.scalar_bits[0];
  require(!validate_torus_domain_row_worklist(production, spec, mutated),
          "D2b stale lattice accepted");
  mutated = valid;
  mutated.axes[0].offsets.back() = 0U;
  require(!validate_torus_domain_row_worklist(production, spec, mutated),
          "D2b malformed CSR accepted");
  mutated = valid;
  ++mutated.candidate_pairs;
  require(!validate_torus_domain_row_worklist(production, spec, mutated),
          "D2b stale counter accepted");
  mutated = valid;
  mutated.domain_worklist_complete = false;
  require(!validate_torus_domain_row_worklist(production, spec, mutated),
          "D2b stale completion accepted");

  PatchComplex bad_beta = production;
  bad_beta.tori.front().expanded_radius_i += 1.0e-4;
  const auto theta = build_torus_interval_row_worklist(bad_beta, spec);
  const auto fallback = build_torus_domain_row_worklist(bad_beta, spec);
  require(fallback.status == "pass_with_theta_fallback" &&
              fallback.beta_fallback_to_theta_tori >= 1U &&
              validate_torus_domain_row_worklist(bad_beta, spec, fallback),
          "D2b invalid beta did not fallback to theta");
  PatchRefId first_torus;
  require(PatchRefId::try_make(PatchKind::toroidal, 0U, first_torus),
          "D2b fallback torus ref failed");
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    for (std::size_t row = 0U; row < theta.axes[axis].row_count(); ++row) {
      const auto theta_row = theta.axes[axis].row(row);
      const auto fallback_row = fallback.axes[axis].row(row);
      require(
          std::binary_search(theta_row.begin(), theta_row.end(), first_torus) ==
              std::binary_search(fallback_row.begin(), fallback_row.end(),
                                 first_torus),
          "D2b beta fallback changed theta torus membership");
    }
  }

  PatchComplex bad_frame = production;
  bad_frame.tori.front().basis_u.x += 1.0e-4;
  const auto theta_full = build_torus_interval_row_worklist(bad_frame, spec);
  const auto full_fallback = build_torus_domain_row_worklist(bad_frame, spec);
  require(
      full_fallback.status == "pass_with_full_fallback" &&
          full_fallback.full_fallback_tori >= 1U &&
          validate_torus_domain_row_worklist(bad_frame, spec, full_fallback),
      "D2b invalid frame did not preserve D2a full fallback");
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    for (std::size_t row = 0U; row < theta_full.axes[axis].row_count(); ++row) {
      const auto theta_row = theta_full.axes[axis].row(row);
      const auto fallback_row = full_fallback.axes[axis].row(row);
      require(
          std::binary_search(theta_row.begin(), theta_row.end(), first_torus) ==
              std::binary_search(fallback_row.begin(), fallback_row.end(),
                                 first_torus),
          "D2b full fallback changed D2a torus membership");
    }
  }

  TorusDomainRowWorklistBuildOptions limited;
  limited.candidate_pair_limit = 0U;
  const auto rejected =
      build_torus_domain_row_worklist(production, spec, nullptr, limited);
  require(rejected.status == "candidate_limit",
          "D2b candidate limit did not fail closed");

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
  const auto one = build_torus_domain_row_worklist(one_pair, memory_spec);
  const auto two = build_torus_domain_row_worklist(two_pairs, memory_spec);
  require(validate_torus_domain_row_worklist(one_pair, memory_spec, one) &&
              validate_torus_domain_row_worklist(two_pairs, memory_spec, two) &&
              two.torus_domain_bounds.size() ==
                  2U * one.torus_domain_bounds.size() &&
              two.candidate_pairs == 2U * one.candidate_pairs &&
              two.accounted_peak_temporary_bytes <=
                  2U * one.accounted_peak_temporary_bytes &&
              one.retained_bytes == exact_retained(one) &&
              two.retained_bytes == exact_retained(two),
          "D2b retained/accounted memory is not linear");
}

} // namespace

int main() {
  try {
    check_directed_arc_envelope();
    check_beta_extrema_boundaries();
    check_fixtures_and_telemetry();
    check_fallback_mutation_and_memory();
    std::cout << "tabi_ses_torus_domain_bounds_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_torus_domain_bounds_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
