#include "tabi/tabi_ses_patch_ray.h"
#include "tabi/tabi_ses_tiled_row_worklist.h"
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
  PatchComplex complex =
      build_patch_complex(std::vector<Atom>{{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                            {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                            {0.0, root3, 0.0, 0.0, 1.0}},
                          1.0, nullptr, false);
  require(complex.status == "pass" && !complex.contacts.empty() &&
              !complex.tori.empty() && !complex.reentrant.empty(),
          "D2c production complex is incomplete");
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

void require_same_rows(const PreparedTorusDomainRowWorklist &baseline,
                       const PreparedTiledPatchRayRowWorklist &tiled,
                       const char *message) {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(baseline.axes[axis].axis == tiled.axes[axis].axis &&
                baseline.axes[axis].u_count == tiled.axes[axis].u_count &&
                baseline.axes[axis].v_count == tiled.axes[axis].v_count &&
                baseline.axes[axis].offsets == tiled.axes[axis].offsets &&
                baseline.axes[axis].candidates == tiled.axes[axis].candidates,
            message);
  }
}

void require_same_tiled_rows(const PreparedTiledPatchRayRowWorklist &left,
                             const PreparedTiledPatchRayRowWorklist &right,
                             const char *message) {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(left.axes[axis].offsets == right.axes[axis].offsets &&
                left.axes[axis].candidates == right.axes[axis].candidates,
            message);
  }
}

struct Coverage {
  std::size_t candidates{0U};
  std::size_t required{0U};
  std::size_t false_negatives{0U};
};

Coverage require_fixture_parity(const PatchComplex &complex,
                                const PatchMappingGridSpec &spec,
                                TitanExecutor *executor = nullptr) {
  const auto baseline =
      build_torus_domain_row_worklist(complex, spec, executor);
  const auto tiled =
      build_tiled_patch_ray_row_worklist(complex, spec, executor);
  require(validate_torus_domain_row_worklist(complex, spec, baseline) &&
              validate_tiled_patch_ray_row_worklist(complex, spec, tiled),
          "D2c fixture validation failed");
  require_same_rows(baseline, tiled, "D2c row CSR differs from D2b baseline");
  const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
  Coverage coverage;
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    coverage.candidates += static_cast<std::size_t>(std::count_if(
        tiled.axes[axis_index].candidates.begin(),
        tiled.axes[axis_index].candidates.end(),
        [](PatchRefId ref) { return ref.kind() == PatchKind::toroidal; }));
    for (std::size_t row = 0U; row < tiled.axes[axis_index].row_count();
         ++row) {
      DecodedPatchRayRow decoded;
      require(decode_patch_ray_row(spec, axis, row, decoded),
              "D2c row decode failed");
      const auto candidates = tiled.axes[axis_index].row(row);
      for (std::size_t torus = 0U; torus < complex.tori.size(); ++torus) {
        PatchRefId ref;
        require(PatchRefId::try_make(PatchKind::toroidal, torus, ref),
                "D2c torus ref overflow");
        const PatchAxisRoots roots =
            intersect_patch_axis(complex, probes, PatchKind::toroidal, torus,
                                 axis, decoded.u, decoded.v);
        const bool required = roots.count != 0U || !roots.certified;
        coverage.required += required;
        coverage.false_negatives +=
            required &&
            !std::binary_search(candidates.begin(), candidates.end(), ref);
      }
    }
  }
  require(coverage.false_negatives == 0U,
          "D2c direct intersection audit found a false negative");
  return coverage;
}

std::uint64_t checksum(const PreparedTiledPatchRayRowWorklist &prepared) {
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
  for (const auto &axis : prepared.axes) {
    mix(static_cast<std::uint64_t>(axis.axis));
    mix(axis.u_count);
    mix(axis.v_count);
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const PatchRefId ref : axis.candidates)
      mix(ref.packed());
  }
  mix(prepared.total_rows);
  mix(prepared.tile_u_rows);
  mix(prepared.tile_v_rows);
  mix(prepared.bound_rectangles);
  mix(prepared.total_tiles);
  mix(prepared.active_tiles);
  mix(prepared.rectangle_tile_pairs);
  mix(prepared.raw_row_pairs);
  mix(prepared.deduplicated_candidate_pairs);
  mix(prepared.maximum_raw_pairs_per_tile);
  mix(prepared.torus_domain_plans);
  mix(prepared.torus_theta_fallback_plans);
  mix(prepared.torus_full_fallback_plans);
  mix(prepared.tiled_worklist_complete);
  for (const char value : prepared.status)
    mix(static_cast<unsigned char>(value));
  return hash;
}

std::size_t exact_retained(const PreparedTiledPatchRayRowWorklist &prepared) {
  std::size_t bytes = 0U;
  for (const auto &axis : prepared.axes) {
    bytes += axis.offsets.capacity() * sizeof(std::size_t);
    bytes += axis.candidates.capacity() * sizeof(PatchRefId);
  }
  return bytes;
}

void check_production_determinism_and_timing() {
  const PatchComplex complex = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(complex, 0.5);
  TitanExecutor one(1U);
  TitanExecutor eight(8U);
  const auto baseline = build_torus_domain_row_worklist(complex, spec, &one);
  const auto serial = build_tiled_patch_ray_row_worklist(complex, spec, &one);
  const auto threaded =
      build_tiled_patch_ray_row_worklist(complex, spec, &eight);
  TiledPatchRayRowWorklistBuildOptions reversed_options;
  reversed_options.reverse_rectangle_schedule_for_audit = true;
  reversed_options.reverse_tile_schedule_for_audit = true;
  const auto reversed = build_tiled_patch_ray_row_worklist(
      complex, spec, &eight, reversed_options);
  require(validate_tiled_patch_ray_row_worklist(complex, spec, serial) &&
              validate_tiled_patch_ray_row_worklist(complex, spec, threaded) &&
              validate_tiled_patch_ray_row_worklist(complex, spec, reversed),
          "D2c production worklist failed validation");
  require_same_rows(baseline, serial, "D2c production baseline parity failed");
  require_same_tiled_rows(serial, threaded, "D2c nt1/8 rows differ");
  require_same_tiled_rows(serial, reversed, "D2c reverse rows differ");
  require(checksum(serial) == checksum(threaded) &&
              checksum(serial) == checksum(reversed),
          "D2c semantic checksum differs across schedules");
  const Coverage production_coverage =
      require_fixture_parity(complex, spec, &one);
  require(production_coverage.candidates == 87U &&
              threaded.scheduled_workers > 1U &&
              threaded.count_workers_used >= 1U &&
              threaded.fill_workers_used >= 1U,
          "D2c production torus-candidate/dynamic-worker invariant failed");
  require(serial.retained_bytes == exact_retained(serial),
          "D2c retained byte accounting failed");

  for (int warmup = 0; warmup < 2; ++warmup)
    require(checksum(build_tiled_patch_ray_row_worklist(complex, spec)) ==
                checksum(serial),
            "D2c warmup drifted");
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto timed = build_tiled_patch_ray_row_worklist(complex, spec);
    require(checksum(timed) == checksum(serial), "D2c timed build drifted");
    sample = timed.build_ms;
  }
  auto sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  std::cout << "d2c_checksum=" << checksum(serial)
            << " d2c_candidates=" << serial.deduplicated_candidate_pairs
            << " d2c_torus_candidates=" << production_coverage.candidates
            << " d2c_raw_row_pairs=" << serial.raw_row_pairs
            << " d2c_rectangles=" << serial.bound_rectangles
            << " d2c_total_tiles=" << serial.total_tiles
            << " d2c_active_tiles=" << serial.active_tiles
            << " d2c_rectangle_tile_pairs=" << serial.rectangle_tile_pairs
            << " d2c_max_raw_pairs_per_tile="
            << serial.maximum_raw_pairs_per_tile
            << " d2c_count_workers_nt8=" << threaded.count_workers_used
            << " d2c_fill_workers_nt8=" << threaded.fill_workers_used
            << " d2c_retained_bytes=" << serial.retained_bytes
            << " d2c_peak_upper_bytes="
            << serial.accounted_peak_temporary_upper_bound_bytes
            << " d2c_observed_scratch_nt8="
            << threaded.observed_worker_scratch_bytes
            << " d2c_samples_ms=" << samples[0] << ',' << samples[1] << ','
            << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
            << samples[5] << ',' << samples[6] << " d2c_median_ms=" << sorted[3]
            << '\n';
}

void check_eight_fixture_parity() {
  PatchComplex production = production_three_atom_complex();
  Coverage aggregate =
      require_fixture_parity(production, enclosing_grid(production, 0.5));
  const std::vector<Atom> base = production.atoms;
  for (std::size_t fixture = 0U; fixture < 4U; ++fixture) {
    const Vec3 translation =
        fixture == 3U ? Vec3{1.0e6, -2.0e6, 0.5e6}
                      : Vec3{4.0 * fixture, -3.0 * fixture, 1.5 * fixture};
    const PatchComplex rotated = build_patch_complex(
        rotate_translate_atoms(base, 0.271 * (fixture + 1U), translation), 1.0,
        nullptr, false);
    const Coverage current =
        require_fixture_parity(rotated, enclosing_grid(rotated, 0.5));
    aggregate.candidates += current.candidates;
    aggregate.required += current.required;
    aggregate.false_negatives += current.false_negatives;
  }
  const PatchComplex tetra =
      build_patch_complex(regular_tetrahedron_atoms(3.2), 1.0, nullptr, false);
  const Coverage tetra_coverage =
      require_fixture_parity(tetra, enclosing_grid(tetra, 0.5));
  aggregate.candidates += tetra_coverage.candidates;
  aggregate.required += tetra_coverage.required;
  aggregate.false_negatives += tetra_coverage.false_negatives;
  const PatchComplex spindle = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  const Coverage spindle_coverage =
      require_fixture_parity(spindle, enclosing_grid(spindle, 0.25));
  aggregate.candidates += spindle_coverage.candidates;
  aggregate.required += spindle_coverage.required;
  aggregate.false_negatives += spindle_coverage.false_negatives;
  const PatchComplex unequal =
      build_patch_complex(std::vector<Atom>{{-1.55, 0.0, 0.0, 0.0, 0.8},
                                            {1.55, 0.0, 0.0, 0.0, 1.35}},
                          1.0, nullptr, false);
  const Coverage unequal_coverage =
      require_fixture_parity(unequal, enclosing_grid(unequal, 0.25));
  aggregate.candidates += unequal_coverage.candidates;
  aggregate.required += unequal_coverage.required;
  aggregate.false_negatives += unequal_coverage.false_negatives;
  require(aggregate.candidates == 782U && aggregate.false_negatives == 0U,
          "D2c aggregate parity/candidate invariant failed");
  std::cout << "d2c_aggregate_fixtures=8 d2c_aggregate_candidates="
            << aggregate.candidates
            << " d2c_aggregate_required=" << aggregate.required
            << " d2c_aggregate_false_negatives=" << aggregate.false_negatives
            << '\n';
}

PreparedTiledPatchRayRowWorklist
require_mutation_parity(const PatchComplex &complex,
                        const PatchMappingGridSpec &spec, const char *message) {
  const auto baseline = build_torus_domain_row_worklist(complex, spec);
  const auto tiled = build_tiled_patch_ray_row_worklist(complex, spec);
  require(validate_torus_domain_row_worklist(complex, spec, baseline) &&
              validate_tiled_patch_ray_row_worklist(complex, spec, tiled),
          message);
  require_same_rows(baseline, tiled, message);
  return tiled;
}

void check_fallbacks_mutations_and_limits() {
  const PatchComplex production = production_three_atom_complex();
  const PatchMappingGridSpec spec = enclosing_grid(production, 0.5);
  PatchComplex beta = production;
  beta.tori.front().expanded_radius_i += 1.0e-4;
  require(require_mutation_parity(beta, spec, "D2c beta-to-theta parity failed")
                  .torus_theta_fallback_plans >= 1U,
          "D2c beta mutation did not record theta fallback");
  PatchComplex frame = production;
  frame.tori.front().basis_u.x += 1.0e-4;
  require(require_mutation_parity(frame, spec,
                                  "D2c invalid-frame full parity failed")
                  .torus_full_fallback_plans >= 1U,
          "D2c invalid frame did not record full fallback");
  PatchComplex empty = production;
  empty.tori.front().exposed_intervals.clear();
  require(
      require_mutation_parity(empty, spec, "D2c empty-interval parity failed")
              .torus_full_fallback_plans >= 1U,
      "D2c empty interval did not record full fallback");
  PatchComplex invalid = production;
  invalid.tori.front().exposed_intervals.front() = {2.0, 1.0};
  require(require_mutation_parity(invalid, spec,
                                  "D2c invalid-interval parity failed")
                  .torus_full_fallback_plans >= 1U,
          "D2c invalid interval did not record full fallback");
  PatchComplex nonfinite = production;
  nonfinite.tori.front().exposed_intervals.front().begin =
      std::numeric_limits<double>::quiet_NaN();
  require(require_mutation_parity(nonfinite, spec,
                                  "D2c nonfinite-interval parity failed")
                  .torus_full_fallback_plans >= 1U,
          "D2c nonfinite interval did not record full fallback");

  const auto valid = build_tiled_patch_ray_row_worklist(production, spec);
  auto mutated = valid;
  ++mutated.complex_fingerprint.first;
  require(!validate_tiled_patch_ray_row_worklist(production, spec, mutated),
          "D2c stale fingerprint accepted");
  mutated = valid;
  ++mutated.lattice_identity.scalar_bits[0];
  require(!validate_tiled_patch_ray_row_worklist(production, spec, mutated),
          "D2c stale lattice accepted");
  mutated = valid;
  mutated.axes[0].offsets.back() = 0U;
  require(!validate_tiled_patch_ray_row_worklist(production, spec, mutated),
          "D2c malformed CSR accepted");
  mutated = valid;
  ++mutated.deduplicated_candidate_pairs;
  require(!validate_tiled_patch_ray_row_worklist(production, spec, mutated),
          "D2c stale counter accepted");
  mutated = valid;
  mutated.tiled_worklist_complete = false;
  require(!validate_tiled_patch_ray_row_worklist(production, spec, mutated),
          "D2c stale completion accepted");

  TiledPatchRayRowWorklistBuildOptions options;
  options.raw_row_pair_limit = 0U;
  require(build_tiled_patch_ray_row_worklist(production, spec, nullptr, options)
                  .status == "raw_row_pair_limit",
          "D2c raw-pair limit did not fail closed");
  options = {};
  options.rectangle_tile_pair_limit = 0U;
  require(build_tiled_patch_ray_row_worklist(production, spec, nullptr, options)
                  .status == "rectangle_tile_pair_limit",
          "D2c tile-pair limit did not fail closed");
  options = {};
  options.final_candidate_pair_limit = 0U;
  require(build_tiled_patch_ray_row_worklist(production, spec, nullptr, options)
                  .status == "final_candidate_limit",
          "D2c final-pair limit did not fail closed");
  options = {};
  options.per_tile_scratch_pair_limit = 0U;
  require(build_tiled_patch_ray_row_worklist(production, spec, nullptr, options)
                  .status == "per_tile_scratch_limit",
          "D2c scratch limit did not fail closed");
  options = {};
  options.tile_u_rows = 0U;
  require(build_tiled_patch_ray_row_worklist(production, spec, nullptr, options)
                  .status == "invalid_input",
          "D2c zero tile dimension did not fail closed");
  options = {};
  options.tile_u_rows = std::numeric_limits<std::size_t>::max();
  options.tile_v_rows = std::numeric_limits<std::size_t>::max();
  const auto maximum_tile =
      build_tiled_patch_ray_row_worklist(production, spec, nullptr, options);
  require(validate_tiled_patch_ray_row_worklist(production, spec, maximum_tile),
          "D2c maximum tile dimensions overflowed or failed validation");
  require_same_rows(build_torus_domain_row_worklist(production, spec),
                    maximum_tile,
                    "D2c maximum tile dimensions changed membership");
}

void check_memory_scaling_and_peak() {
  const PatchComplex one = build_patch_complex(
      std::vector<Atom>{{-1.9, 0.0, 0.0, 0.0, 1.0}, {1.9, 0.0, 0.0, 0.0, 1.0}},
      1.0, nullptr, false);
  const PatchComplex two =
      build_patch_complex(std::vector<Atom>{{-6.9, 0.0, 0.0, 0.0, 1.0},
                                            {-3.1, 0.0, 0.0, 0.0, 1.0},
                                            {3.1, 0.0, 0.0, 0.0, 1.0},
                                            {6.9, 0.0, 0.0, 0.0, 1.0}},
                          1.0, nullptr, false);
  const PatchMappingGridSpec spec{{-10.0, -3.0, -3.0}, 0.25, 81U, 25U, 25U};
  TitanExecutor serial_executor(1U);
  const auto tiled_one =
      build_tiled_patch_ray_row_worklist(one, spec, &serial_executor);
  const auto tiled_two =
      build_tiled_patch_ray_row_worklist(two, spec, &serial_executor);
  const auto baseline_two =
      build_torus_domain_row_worklist(two, spec, &serial_executor);
  const bool peak_at_most_twice =
      tiled_one.accounted_peak_temporary_upper_bound_bytes <=
          std::numeric_limits<std::size_t>::max() / 2U &&
      tiled_two.accounted_peak_temporary_upper_bound_bytes <=
          tiled_one.accounted_peak_temporary_upper_bound_bytes * 2U;
  require_same_rows(baseline_two, tiled_two,
                    "D2c memory fixture row parity failed");
  require(validate_tiled_patch_ray_row_worklist(one, spec, tiled_one) &&
              validate_tiled_patch_ray_row_worklist(two, spec, tiled_two) &&
              validate_torus_domain_row_worklist(two, spec, baseline_two) &&
              tiled_one.retained_bytes == exact_retained(tiled_one) &&
              tiled_two.retained_bytes == exact_retained(tiled_two) &&
              peak_at_most_twice &&
              tiled_two.accounted_peak_temporary_upper_bound_bytes <
                  baseline_two.accounted_peak_temporary_bytes &&
              tiled_two.maximum_raw_pairs_per_tile < tiled_two.raw_row_pairs,
          "D2c tiled memory/peak scaling gate failed");
  std::cout << "d2c_memory_raw_pairs=" << tiled_two.raw_row_pairs
            << " d2c_memory_max_tile_pairs="
            << tiled_two.maximum_raw_pairs_per_tile << " d2c_memory_tiled_peak="
            << tiled_two.accounted_peak_temporary_upper_bound_bytes
            << " d2c_memory_d2b_peak="
            << baseline_two.accounted_peak_temporary_bytes << '\n';
}

} // namespace

int main() {
  try {
    check_production_determinism_and_timing();
    check_eight_fixture_parity();
    check_fallbacks_mutations_and_limits();
    check_memory_scaling_and_peak();
    std::cout << "tabi_ses_tiled_row_worklist_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_tiled_row_worklist_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
