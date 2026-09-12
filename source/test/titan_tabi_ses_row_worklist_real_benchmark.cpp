#include "tabi/tabi_ses_patch_intersections.h"
#include "tabi/tabi_ses_patch_mapping.h"
#include "tabi/tabi_ses_patch_row_worklist.h"
#include "tabi/tabi_ses_quartic.h"
#include "tabi/tabi_ses_tiled_row_worklist.h"
#include "tabi/tabi_ses_torus_domain_bounds.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

constexpr double kProbeRadius = 1.4;
constexpr std::size_t kMaximumTotalRows = 200000000U;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

[[nodiscard]] bool checked_add(std::size_t &value,
                               std::size_t increment) noexcept {
  if (value > std::numeric_limits<std::size_t>::max() - increment)
    return false;
  value += increment;
  return true;
}

[[nodiscard]] bool checked_multiply(std::size_t left, std::size_t right,
                                    std::size_t &result) noexcept {
  result = 0U;
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  result = left * right;
  return true;
}

std::string trimmed(std::string value) {
  const auto first = value.find_first_not_of(' ');
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(' ');
  return value.substr(first, last - first + 1U);
}

// This is intentionally identical to the deterministic TITAN-owned parser in
// titan_tabi_ses_patch_ray_benchmark.cpp. PDB radii are geometry-only proxies.
double pdb_radius(std::string element) {
  element = trimmed(std::move(element));
  if (element == "H")
    return 1.20;
  if (element == "N")
    return 1.55;
  if (element == "O")
    return 1.52;
  if (element == "S" || element == "P")
    return 1.80;
  return 1.70;
}

std::vector<Atom> read_atoms(const std::string &path, std::size_t max_atoms) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open input: " + path);
  std::vector<Atom> atoms;
  if (max_atoms != 0U)
    atoms.reserve(max_atoms);
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM")))
      continue;
    Atom atom;
    const bool fixed_pdb = line.size() >= 54U && !path.ends_with(".pqr") &&
                           !path.ends_with(".PQR");
    if (fixed_pdb) {
      atom.x = std::stod(line.substr(30U, 8U));
      atom.y = std::stod(line.substr(38U, 8U));
      atom.z = std::stod(line.substr(46U, 8U));
      atom.q = 0.0;
      atom.radius =
          pdb_radius(line.size() >= 78U ? line.substr(76U, 2U) : std::string{});
    } else {
      std::istringstream stream(line);
      std::vector<std::string> fields;
      std::string field;
      while (stream >> field)
        fields.push_back(field);
      if (fields.size() < 9U)
        throw std::runtime_error("malformed selected atom record");
      const std::size_t n = fields.size();
      atom.x = std::stod(fields[n - 5U]);
      atom.y = std::stod(fields[n - 4U]);
      atom.z = std::stod(fields[n - 3U]);
      atom.q = std::stod(fields[n - 2U]);
      atom.radius = std::stod(fields[n - 1U]);
    }
    if (!std::isfinite(atom.x) || !std::isfinite(atom.y) ||
        !std::isfinite(atom.z) || !std::isfinite(atom.q) ||
        !std::isfinite(atom.radius) || atom.radius <= 0.0)
      throw std::runtime_error("invalid selected atom geometry or parameters");
    atoms.push_back(atom);
    if (max_atoms != 0U && atoms.size() >= max_atoms)
      break;
  }
  if (atoms.empty())
    throw std::runtime_error("input contains no atoms: " + path);
  return atoms;
}

PatchMappingGridSpec surface_grid(std::span<const Atom> atoms, double spacing) {
  require(!atoms.empty() && std::isfinite(spacing) && spacing > 0.0,
          "invalid grid input");
  std::array<double, 3> lower = {std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::infinity()};
  std::array<double, 3> upper = {-std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity()};
  for (const Atom &atom : atoms) {
    const std::array<double, 3> center = {atom.x, atom.y, atom.z};
    const double radius = atom.radius + kProbeRadius;
    require(std::isfinite(radius) && radius > 0.0, "invalid atom radius");
    for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
      require(std::isfinite(center[coordinate]), "nonfinite atom coordinate");
      lower[coordinate] =
          std::min(lower[coordinate], center[coordinate] - radius);
      upper[coordinate] =
          std::max(upper[coordinate], center[coordinate] + radius);
    }
  }
  std::array<std::size_t, 3> counts{};
  for (std::size_t coordinate = 0U; coordinate < 3U; ++coordinate) {
    lower[coordinate] -= 2.0 * spacing;
    upper[coordinate] += 2.0 * spacing;
    const double intervals =
        std::ceil((upper[coordinate] - lower[coordinate]) / spacing);
    require(std::isfinite(intervals) && intervals >= 0.0 &&
                intervals < static_cast<double>(
                                std::numeric_limits<std::size_t>::max()),
            "grid count overflow");
    counts[coordinate] = static_cast<std::size_t>(intervals);
    require(checked_add(counts[coordinate], 1U), "grid count overflow");
  }
  std::size_t xy = 0U;
  std::size_t xz = 0U;
  std::size_t yz = 0U;
  std::size_t rows = 0U;
  require(checked_multiply(counts[0], counts[1], xy) &&
              checked_multiply(counts[0], counts[2], xz) &&
              checked_multiply(counts[1], counts[2], yz) &&
              checked_add(rows, xy) && checked_add(rows, xz) &&
              checked_add(rows, yz) && rows <= kMaximumTotalRows,
          "grid row cap exceeded");
  return {
      {lower[0], lower[1], lower[2]}, spacing, counts[0], counts[1], counts[2]};
}

std::size_t total_rows(const PatchMappingGridSpec &spec) {
  std::size_t xy = 0U;
  std::size_t xz = 0U;
  std::size_t yz = 0U;
  std::size_t rows = 0U;
  require(checked_multiply(spec.nx, spec.ny, xy) &&
              checked_multiply(spec.nx, spec.nz, xz) &&
              checked_multiply(spec.ny, spec.nz, yz) && checked_add(rows, xy) &&
              checked_add(rows, xz) && checked_add(rows, yz),
          "total row overflow");
  return rows;
}

template <typename Prepared>
void require_same_rows(const Prepared &left, const Prepared &right,
                       const char *message) {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(left.axes[axis].axis == right.axes[axis].axis &&
                left.axes[axis].u_count == right.axes[axis].u_count &&
                left.axes[axis].v_count == right.axes[axis].v_count &&
                left.axes[axis].offsets == right.axes[axis].offsets &&
                left.axes[axis].candidates == right.axes[axis].candidates,
            message);
  }
}

template <typename Left, typename Right>
void require_same_rows_cross(const Left &left, const Right &right,
                             const char *message) {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(left.axes[axis].axis == right.axes[axis].axis &&
                left.axes[axis].u_count == right.axes[axis].u_count &&
                left.axes[axis].v_count == right.axes[axis].v_count &&
                left.axes[axis].offsets == right.axes[axis].offsets &&
                left.axes[axis].candidates == right.axes[axis].candidates,
            message);
  }
}

std::uint64_t checksum(const PreparedTiledPatchRayRowWorklist &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  for (const std::uint64_t bits : prepared.lattice_identity.scalar_bits)
    mix(bits);
  for (const std::size_t count : prepared.lattice_identity.counts)
    mix(count);
  for (const PatchRayAxisRowWorklist &axis : prepared.axes) {
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const PatchRefId ref : axis.candidates)
      mix(ref.packed());
  }
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
  return hash;
}

double median(std::vector<double> values) {
  require(!values.empty(), "empty timing samples");
  std::sort(values.begin(), values.end());
  return values[values.size() / 2U];
}

struct TiledTiming {
  double build_ms{0.0};
  double bounds_ms{0.0};
  double tile_ms{0.0};
  double count_ms{0.0};
  double prefix_ms{0.0};
  double fill_ms{0.0};
};

TiledTiming benchmark_tiled(std::size_t samples, const PatchComplex &complex,
                            const PatchMappingGridSpec &spec,
                            TitanExecutor *executor,
                            const PreparedTiledPatchRayRowWorklist &canonical) {
  TiledTiming result;
  {
    const auto warmup =
        build_tiled_patch_ray_row_worklist(complex, spec, executor);
    require(validate_tiled_patch_ray_row_worklist(complex, spec, warmup),
            "D2c warmup validation failed");
    require_same_rows(canonical, warmup, "D2c warmup CSR mismatch");
    require(checksum(canonical) == checksum(warmup),
            "D2c warmup checksum mismatch");
  }
  std::vector<double> build;
  std::vector<double> bounds;
  std::vector<double> tile;
  std::vector<double> count;
  std::vector<double> prefix;
  std::vector<double> fill;
  build.reserve(samples);
  bounds.reserve(samples);
  tile.reserve(samples);
  count.reserve(samples);
  prefix.reserve(samples);
  fill.reserve(samples);
  for (std::size_t sample = 0U; sample < samples; ++sample) {
    auto current = build_tiled_patch_ray_row_worklist(complex, spec, executor);
    require(validate_tiled_patch_ray_row_worklist(complex, spec, current),
            "timed D2c validation failed");
    require_same_rows(canonical, current, "timed D2c CSR mismatch");
    require(checksum(canonical) == checksum(current),
            "timed D2c checksum mismatch");
    build.push_back(current.build_ms);
    bounds.push_back(current.bounds_plan_ms);
    tile.push_back(current.tile_plan_ms);
    count.push_back(current.count_ms);
    prefix.push_back(current.prefix_ms);
    fill.push_back(current.fill_ms);
  }
  result.build_ms = median(std::move(build));
  result.bounds_ms = median(std::move(bounds));
  result.tile_ms = median(std::move(tile));
  result.count_ms = median(std::move(count));
  result.prefix_ms = median(std::move(prefix));
  result.fill_ms = median(std::move(fill));
  return result;
}

struct RemovedAudit {
  std::size_t removed{0U};
  std::size_t audited{0U};
  std::size_t required_removed_observed{0U};
  std::size_t roots_count_nonzero{0U};
  std::size_t roots_uncertified{0U};
  std::size_t roots_count_nonzero_and_uncertified{0U};
  struct FirstUnsafeWitness {
    bool present{false};
    RayAxis axis{RayAxis::x};
    std::size_t row{0U};
    DecodedPatchRayRow decoded;
    PatchKind kind{PatchKind::contact};
    std::size_t patch_index{0U};
    PatchAxisRoots roots;
  } first_unsafe;
  bool complete{true};
  double elapsed_ms{0.0};
};

[[nodiscard]] const char *patch_kind_name(PatchKind kind) noexcept {
  switch (kind) {
  case PatchKind::contact:
    return "contact";
  case PatchKind::toroidal:
    return "toroidal";
  case PatchKind::reentrant:
    return "reentrant";
  }
  return "invalid";
}

[[nodiscard]] const char *
quartic_failure_name(QuarticDiagnosticFailure failure) noexcept {
  switch (failure) {
  case QuarticDiagnosticFailure::none:
    return "none";
  case QuarticDiagnosticFailure::invalid_or_degenerate_coefficients:
    return "invalid_or_degenerate_coefficients";
  case QuarticDiagnosticFailure::derivative_cubic_discriminant_unresolved:
    return "derivative_cubic_discriminant_unresolved";
  case QuarticDiagnosticFailure::derivative_residual_uncertified:
    return "derivative_residual_uncertified";
  case QuarticDiagnosticFailure::critical_boundary_sign_unresolved:
    return "critical_boundary_sign_unresolved";
  case QuarticDiagnosticFailure::repeated_root_localization_unavailable:
    return "repeated_root_localization_unavailable";
  case QuarticDiagnosticFailure::repeated_root_localization_mismatch:
    return "repeated_root_localization_mismatch";
  case QuarticDiagnosticFailure::boundary_nonfinite:
    return "boundary_nonfinite";
  case QuarticDiagnosticFailure::boundary_near_zero:
    return "boundary_near_zero";
  case QuarticDiagnosticFailure::root_capacity_or_nonfinite:
    return "root_capacity_or_nonfinite";
  case QuarticDiagnosticFailure::root_residual_failure:
    return "root_residual_failure";
  }
  return "invalid";
}

[[nodiscard]] const char *
quartic_refinement_name(QuarticRefinementTermination termination) noexcept {
  switch (termination) {
  case QuarticRefinementTermination::not_bracketed:
    return "not_bracketed";
  case QuarticRefinementTermination::exact_or_residual:
    return "exact_or_residual";
  case QuarticRefinementTermination::width_floor_one:
    return "width_floor_one";
  case QuarticRefinementTermination::width_scaled:
    return "width_scaled";
  case QuarticRefinementTermination::iteration_limit:
    return "iteration_limit";
  }
  return "invalid";
}

[[nodiscard]] const char *
torus_filter_name(TorusRootFilterReason reason) noexcept {
  switch (reason) {
  case TorusRootFilterReason::accepted:
    return "accepted";
  case TorusRootFilterReason::radial_axis_singular:
    return "radial_axis_singular";
  case TorusRootFilterReason::theta_interval:
    return "theta_interval";
  case TorusRootFilterReason::tube_consistency:
    return "tube_consistency";
  case TorusRootFilterReason::beta_minor_arc:
    return "beta_minor_arc";
  case TorusRootFilterReason::reentrant_coverage:
    return "reentrant_coverage";
  case TorusRootFilterReason::tangential:
    return "tangential";
  }
  return "invalid";
}

[[nodiscard]] bool same_patch_axis_roots(const PatchAxisRoots &left,
                                         const PatchAxisRoots &right) noexcept {
  if (left.count != right.count || left.certified != right.certified)
    return false;
  for (std::size_t index = 0U; index < left.coordinates.size(); ++index) {
    if (std::bit_cast<std::uint64_t>(left.coordinates[index]) !=
            std::bit_cast<std::uint64_t>(right.coordinates[index]) ||
        left.primitive_simplicity_certified[index] !=
            right.primitive_simplicity_certified[index])
      return false;
  }
  return true;
}

[[nodiscard]] std::uint64_t double_bits(double value) noexcept {
  return std::bit_cast<std::uint64_t>(value);
}

struct OrthogonalBound {
  double u_lower{0.0};
  double u_upper{0.0};
  double v_lower{0.0};
  double v_upper{0.0};
};

[[nodiscard]] OrthogonalBound orthogonal_bound(const PatchBounds &bounds,
                                               RayAxis axis) noexcept {
  switch (axis) {
  case RayAxis::x:
    return {bounds.lower.y, bounds.upper.y, bounds.lower.z, bounds.upper.z};
  case RayAxis::y:
    return {bounds.lower.x, bounds.upper.x, bounds.lower.z, bounds.upper.z};
  case RayAxis::z:
    return {bounds.lower.x, bounds.upper.x, bounds.lower.y, bounds.upper.y};
  }
  return {};
}

[[nodiscard]] double signed_interval_gap(double value, double lower,
                                         double upper) noexcept {
  if (value < lower)
    return value - lower;
  if (value > upper)
    return value - upper;
  return 0.0;
}

[[nodiscard]] bool row_inside_bound(const DecodedPatchRayRow &row,
                                    const PatchBounds &bounds) noexcept {
  const OrthogonalBound projected = orthogonal_bound(bounds, row.axis);
  return row.u >= projected.u_lower && row.u <= projected.u_upper &&
         row.v >= projected.v_lower && row.v <= projected.v_upper;
}

void print_bound_diagnostic(std::ostream &output, std::string_view label,
                            const DecodedPatchRayRow &row,
                            const PatchBounds &bounds) {
  const OrthogonalBound projected = orthogonal_bound(bounds, row.axis);
  const double u_gap =
      signed_interval_gap(row.u, projected.u_lower, projected.u_upper);
  const double v_gap =
      signed_interval_gap(row.v, projected.v_lower, projected.v_upper);
  output << ' ' << label << "_u_lower=" << projected.u_lower << ' ' << label
         << "_u_upper=" << projected.u_upper << ' ' << label
         << "_v_lower=" << projected.v_lower << ' ' << label
         << "_v_upper=" << projected.v_upper << ' ' << label
         << "_u_gap=" << u_gap << ' ' << label << "_v_gap=" << v_gap << ' '
         << label << "_contains=" << (row_inside_bound(row, bounds) ? 1 : 0);
}

RemovedAudit audit_removed_candidates(
    const PatchComplex &complex, const PatchMappingGridSpec &spec,
    const ReentrantProbeIndex &probes, const PreparedPatchRayRowWorklist &full,
    const PreparedTiledPatchRayRowWorklist &compact, double budget_seconds) {
  using Clock = std::chrono::steady_clock;
  const auto started = Clock::now();
  RemovedAudit result;
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    for (std::size_t row = 0U; row < full.axes[axis_index].row_count(); ++row) {
      const auto full_refs = full.axes[axis_index].row(row);
      const auto compact_refs = compact.axes[axis_index].row(row);
      DecodedPatchRayRow decoded;
      bool decoded_row = false;
      const auto audit_removed_ref = [&](PatchRefId ref) {
        require(ref.valid(), "removed audit invalid patch ref");
        require(checked_add(result.removed, 1U), "removed counter overflow");
        if (!result.complete)
          return;
        const double elapsed_seconds =
            std::chrono::duration<double>(Clock::now() - started).count();
        if (elapsed_seconds > budget_seconds) {
          result.complete = false;
          return;
        }
        if (!decoded_row) {
          require(decode_patch_ray_row(spec, axis, row, decoded),
                  "removed audit row decode failed");
          decoded_row = true;
        }
        const PatchAxisRoots roots =
            intersect_patch_axis(complex, probes, ref.kind(), ref.index(), axis,
                                 decoded.u, decoded.v);
        require(checked_add(result.audited, 1U),
                "removed audit counter overflow");
        if (roots.count != 0U) {
          require(checked_add(result.roots_count_nonzero, 1U),
                  "nonzero-root counter overflow");
        }
        if (!roots.certified) {
          require(checked_add(result.roots_uncertified, 1U),
                  "uncertified-root counter overflow");
        }
        if (roots.count != 0U && !roots.certified) {
          require(checked_add(result.roots_count_nonzero_and_uncertified, 1U),
                  "overlap-root counter overflow");
        }
        if (roots.count != 0U || !roots.certified) {
          require(checked_add(result.required_removed_observed, 1U),
                  "required removed counter overflow");
          if (!result.first_unsafe.present) {
            result.first_unsafe.present = true;
            result.first_unsafe.axis = axis;
            result.first_unsafe.row = row;
            result.first_unsafe.decoded = decoded;
            result.first_unsafe.kind = ref.kind();
            result.first_unsafe.patch_index = ref.index();
            result.first_unsafe.roots = roots;
          }
        }
      };
      std::size_t full_position = 0U;
      std::size_t compact_position = 0U;
      while (full_position < full_refs.size()) {
        if (compact_position == compact_refs.size() ||
            full_refs[full_position] < compact_refs[compact_position]) {
          audit_removed_ref(full_refs[full_position]);
          ++full_position;
          continue;
        }
        require(!(compact_refs[compact_position] < full_refs[full_position]),
                "D2c row is not a subset of D1");
        ++full_position;
        ++compact_position;
      }
      require(compact_position == compact_refs.size(),
              "D2c row is not a subset of D1");
    }
  }
  result.elapsed_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - started).count();
  return result;
}

void emit_first_unsafe_diagnostic(const PatchComplex &complex,
                                  const PatchMappingGridSpec &spec,
                                  const ReentrantProbeIndex &probes,
                                  const RemovedAudit &audit) {
  if (!audit.first_unsafe.present)
    return;
  const RemovedAudit::FirstUnsafeWitness &witness = audit.first_unsafe;
  std::cerr << std::setprecision(17) << "first_unsafe_removed axis="
            << static_cast<unsigned int>(witness.axis) << " row=" << witness.row
            << " u=" << witness.decoded.u << " v=" << witness.decoded.v
            << " u_bits=0x" << std::hex << double_bits(witness.decoded.u)
            << " v_bits=0x" << double_bits(witness.decoded.v) << std::dec
            << " kind=" << patch_kind_name(witness.kind)
            << " patch=" << witness.patch_index
            << " roots_count=" << witness.roots.count
            << " roots_certified=" << (witness.roots.certified ? 1 : 0) << '\n';
  for (std::size_t root = 0U; root < witness.roots.count; ++root) {
    std::cerr << std::setprecision(17) << "first_unsafe_root ordinal=" << root
              << " coordinate=" << witness.roots.coordinates[root]
              << " coordinate_bits=0x" << std::hex
              << double_bits(witness.roots.coordinates[root]) << std::dec
              << " primitive_simplicity="
              << static_cast<unsigned int>(
                     witness.roots.primitive_simplicity_certified[root])
              << '\n';
  }
  if (witness.kind != PatchKind::toroidal) {
    std::cerr << "first_unsafe_bound_audit=stopped_non_toroidal\n";
    return;
  }
  require(witness.patch_index < complex.tori.size(),
          "unsafe torus witness index out of range");
  QuarticDiagnostics quartic_diagnostics;
  TorusPatchAxisDiagnostics filter_diagnostics;
  const PatchAxisRoots diagnostic_roots =
      intersect_patch_axis_with_quartic_diagnostics(
          complex, probes, witness.kind, witness.patch_index, witness.axis,
          witness.decoded.u, witness.decoded.v, &quartic_diagnostics,
          &filter_diagnostics);
  require(same_patch_axis_roots(witness.roots, diagnostic_roots),
          "diagnostic quartic path changed patch-axis roots");
  std::cerr << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "first_unsafe_quartic failure_mask="
            << quartic_diagnostics.failure_mask << " failure_mask_hex=0x"
            << std::hex << quartic_diagnostics.failure_mask << std::dec
            << " first_failure="
            << quartic_failure_name(quartic_diagnostics.first_failure)
            << " coefficient_scale=" << quartic_diagnostics.coefficient_scale
            << " root_bound=" << quartic_diagnostics.root_bound
            << " output_root_count=" << quartic_diagnostics.output_root_count
            << " output_max_residual="
            << quartic_diagnostics.output_max_normalized_residual << '\n';
  for (std::size_t index = 0U;
       index < quartic_diagnostics.input_coefficients.size(); ++index) {
    std::cerr << std::setprecision(17)
              << "first_unsafe_quartic_coefficient index=" << index
              << " input=" << quartic_diagnostics.input_coefficients[index]
              << " input_bits=0x" << std::hex
              << double_bits(quartic_diagnostics.input_coefficients[index])
              << std::dec
              << std::setprecision(
                     std::numeric_limits<long double>::max_digits10)
              << " scaled=" << quartic_diagnostics.scaled_coefficients[index]
              << '\n';
  }
  for (std::size_t index = 0U; index < quartic_diagnostics.output_root_count;
       ++index) {
    std::cerr
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << "first_unsafe_quartic_root index=" << index
        << " refined=" << quartic_diagnostics.refined_roots[index]
        << std::setprecision(17)
        << " output=" << quartic_diagnostics.output_roots[index]
        << " output_bits=0x" << std::hex
        << double_bits(quartic_diagnostics.output_roots[index]) << std::dec
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << " pre_cast_residual="
        << quartic_diagnostics.pre_cast_normalized_residuals[index]
        << " post_cast_residual="
        << quartic_diagnostics.post_cast_normalized_residuals[index]
        << " derivative_magnitude="
        << quartic_diagnostics.derivative_magnitudes[index]
        << " condition_proxy="
        << quartic_diagnostics.root_condition_proxies[index] << " bracketed="
        << static_cast<unsigned int>(quartic_diagnostics.bracketed_roots[index])
        << " repeated="
        << static_cast<unsigned int>(quartic_diagnostics.repeated_roots[index])
        << " bracket_lower=" << quartic_diagnostics.bracket_lowers[index]
        << " bracket_upper=" << quartic_diagnostics.bracket_uppers[index]
        << " bracket_lower_value="
        << quartic_diagnostics.bracket_lower_values[index]
        << " bracket_upper_value="
        << quartic_diagnostics.bracket_upper_values[index]
        << " final_bracket_lower="
        << quartic_diagnostics.final_bracket_lowers[index]
        << " final_bracket_upper="
        << quartic_diagnostics.final_bracket_uppers[index]
        << " final_bracket_lower_value="
        << quartic_diagnostics.final_bracket_lower_values[index]
        << " final_bracket_upper_value="
        << quartic_diagnostics.final_bracket_upper_values[index]
        << " final_bracket_width="
        << quartic_diagnostics.final_bracket_widths[index]
        << " termination_residual="
        << quartic_diagnostics.termination_normalized_residuals[index]
        << " termination="
        << quartic_refinement_name(
               quartic_diagnostics.refinement_terminations[index])
        << '\n';
  }
  for (std::size_t index = 0U; index < filter_diagnostics.root_count; ++index) {
    const TorusRootFilterDiagnostic &filter = filter_diagnostics.roots[index];
    std::cerr << std::setprecision(17)
              << "first_unsafe_torus_root_filter index=" << index
              << " reason=" << torus_filter_name(filter.reason)
              << " quartic_root=" << filter.quartic_root
              << " coordinate=" << filter.coordinate << " point=("
              << filter.point.x << ',' << filter.point.y << ','
              << filter.point.z << ") theta=" << filter.theta
              << " beta=" << filter.beta
              << " radial_length=" << filter.radial_length
              << " tube_length=" << filter.tube_length
              << " tube_residual=" << filter.tube_residual
              << " directional_normal=" << filter.directional_normal << '\n';
  }
  std::cerr << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "first_unsafe_quartic_cubic a=" << quartic_diagnostics.cubic_a
            << " b=" << quartic_diagnostics.cubic_b
            << " c=" << quartic_diagnostics.cubic_c
            << " p=" << quartic_diagnostics.cubic_p
            << " q=" << quartic_diagnostics.cubic_q
            << " discriminant=" << quartic_diagnostics.cubic_discriminant
            << " discriminant_tolerance="
            << quartic_diagnostics.cubic_discriminant_tolerance
            << " critical_count=" << quartic_diagnostics.critical_count << '\n';
  for (std::size_t index = 0U; index < quartic_diagnostics.critical_count;
       ++index) {
    std::cerr << "first_unsafe_quartic_critical index=" << index
              << " point=" << quartic_diagnostics.critical_points[index]
              << " value=" << quartic_diagnostics.critical_values[index]
              << " normalized_residual="
              << quartic_diagnostics.critical_normalized_residuals[index]
              << " derivative_normalized_residual="
              << quartic_diagnostics.derivative_normalized_residuals[index]
              << '\n';
  }
  std::cerr << "first_unsafe_quartic_repeated attempted="
            << (quartic_diagnostics.repeated_root_localization_attempted ? 1
                                                                         : 0)
            << " available="
            << (quartic_diagnostics.repeated_root_localization_available ? 1
                                                                         : 0)
            << " valid="
            << (quartic_diagnostics.repeated_root_localization_valid ? 1 : 0)
            << " count=" << quartic_diagnostics.repeated_root_count << '\n';
  std::cerr << "first_unsafe_quartic_boundary index=";
  if (quartic_diagnostics.first_failure_boundary_index ==
      std::numeric_limits<std::size_t>::max()) {
    std::cerr << "NA point=NA value=NA normalized_residual=NA critical=NA\n";
  } else {
    std::cerr << quartic_diagnostics.first_failure_boundary_index
              << " point=" << quartic_diagnostics.first_failure_point
              << " value=" << quartic_diagnostics.first_failure_value
              << " normalized_residual="
              << quartic_diagnostics.first_failure_normalized_residual
              << " critical="
              << (quartic_diagnostics.first_failure_is_critical ? 1 : 0)
              << '\n';
  }
  const ToroidalPatch &torus = complex.tori[witness.patch_index];
  TorusBetaMinorArc beta_arc;
  const bool beta_arc_valid =
      torus_beta_minor_arc_for_audit(complex, witness.patch_index, beta_arc);
  const PatchBounds outward = outward_padded_patch_bounds(
      complex, PatchKind::toroidal, witness.patch_index);
  PatchBounds d1_full;
  const bool d1_valid = patch_ray_full_padded_bound_for_audit(
      complex, PatchKind::toroidal, witness.patch_index, spec.grid_spacing,
      d1_full);
  std::cerr << std::setprecision(17)
            << "first_unsafe_torus atom_i=" << torus.atom_i
            << " atom_j=" << torus.atom_j
            << " exposed_intervals=" << torus.exposed_intervals.size()
            << " beta_arc_valid=" << (beta_arc_valid ? 1 : 0)
            << " beta_start=" << beta_arc.start
            << " beta_span=" << beta_arc.span << " beta_start_bits=0x"
            << std::hex << double_bits(beta_arc.start) << " beta_span_bits=0x"
            << double_bits(beta_arc.span) << std::dec << '\n';
  std::cerr << "first_unsafe_outward_bounds lower_x=" << outward.lower.x
            << " lower_y=" << outward.lower.y << " lower_z=" << outward.lower.z
            << " upper_x=" << outward.upper.x << " upper_y=" << outward.upper.y
            << " upper_z=" << outward.upper.z;
  print_bound_diagnostic(std::cerr, "outward", witness.decoded, outward);
  std::cerr << '\n';
  std::cerr << "first_unsafe_d1 valid=" << (d1_valid ? 1 : 0);
  if (d1_valid)
    print_bound_diagnostic(std::cerr, "d1", witness.decoded, d1_full);
  std::cerr << '\n';

  bool any_d2a_contains = false;
  bool any_d2b_contains = false;
  std::size_t relevant_interval = std::numeric_limits<std::size_t>::max();
  double relevant_metric = std::numeric_limits<double>::infinity();
  for (std::size_t interval_index = 0U;
       interval_index < torus.exposed_intervals.size(); ++interval_index) {
    PatchBounds d2a;
    PatchBounds d2b;
    const bool d2a_valid = torus_exposed_interval_bound_for_audit(
        complex, witness.patch_index, interval_index, spec.grid_spacing, d2a);
    const bool d2b_valid = torus_domain_interval_bound_for_audit(
        complex, witness.patch_index, interval_index, spec.grid_spacing, d2b);
    const bool d2a_contains =
        d2a_valid && row_inside_bound(witness.decoded, d2a);
    const bool d2b_contains =
        d2b_valid && row_inside_bound(witness.decoded, d2b);
    any_d2a_contains = any_d2a_contains || d2a_contains;
    any_d2b_contains = any_d2b_contains || d2b_contains;
    if (d2b_valid) {
      const OrthogonalBound projected = orthogonal_bound(d2b, witness.axis);
      const double u_gap = signed_interval_gap(
          witness.decoded.u, projected.u_lower, projected.u_upper);
      const double v_gap = signed_interval_gap(
          witness.decoded.v, projected.v_lower, projected.v_upper);
      const double metric = std::max(std::abs(u_gap), std::abs(v_gap));
      if (metric < relevant_metric) {
        relevant_metric = metric;
        relevant_interval = interval_index;
      }
    }
    const AngularInterval &interval = torus.exposed_intervals[interval_index];
    std::cerr << std::setprecision(17)
              << "first_unsafe_interval index=" << interval_index
              << " theta_begin=" << interval.begin
              << " theta_end=" << interval.end << " theta_begin_bits=0x"
              << std::hex << double_bits(interval.begin) << " theta_end_bits=0x"
              << double_bits(interval.end) << std::dec
              << " d2a_valid=" << (d2a_valid ? 1 : 0)
              << " d2b_valid=" << (d2b_valid ? 1 : 0);
    if (d2a_valid)
      print_bound_diagnostic(std::cerr, "d2a", witness.decoded, d2a);
    if (d2b_valid)
      print_bound_diagnostic(std::cerr, "d2b", witness.decoded, d2b);
    std::cerr << '\n';
  }
  std::cerr << std::setprecision(17) << "first_unsafe_membership d1="
            << (d1_valid && row_inside_bound(witness.decoded, d1_full) ? 1 : 0)
            << " any_d2a=" << (any_d2a_contains ? 1 : 0)
            << " any_d2b=" << (any_d2b_contains ? 1 : 0)
            << " relevant_interval=";
  if (relevant_interval == std::numeric_limits<std::size_t>::max())
    std::cerr << "NA relevant_metric=NA\n";
  else
    std::cerr << relevant_interval << " relevant_metric=" << relevant_metric
              << '\n';
}

std::string file_name(const std::string &path) {
  const std::size_t slash = path.find_last_of("/\\");
  return slash == std::string::npos ? path : path.substr(slash + 1U);
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 9) {
      std::cerr << "usage: benchmark <input> <case_label> <spacing> "
                   "<max_atoms> <radii_provenance> <sha256> <samples> "
                   "<audit_budget_seconds>\n";
      return EXIT_FAILURE;
    }
    const std::string input_path = argv[1];
    const std::string case_label = argv[2];
    const double spacing = std::stod(argv[3]);
    const std::size_t max_atoms =
        static_cast<std::size_t>(std::stoull(argv[4]));
    const std::string radii_provenance = argv[5];
    const std::string input_sha256 = argv[6];
    const std::size_t samples = static_cast<std::size_t>(std::stoull(argv[7]));
    const double audit_budget_seconds = std::stod(argv[8]);
    require(samples >= 5U && samples <= 7U &&
                std::isfinite(audit_budget_seconds) &&
                audit_budget_seconds > 0.0,
            "invalid benchmark options");

    using Clock = std::chrono::steady_clock;
    const auto case_started = Clock::now();
    const std::vector<Atom> atoms = read_atoms(input_path, max_atoms);
    TitanExecutor one(1U);
    TitanExecutor eight(8U);
    const auto patch_started = Clock::now();
    const PatchComplex complex =
        build_patch_complex(atoms, kProbeRadius, &eight, false);
    const double patch_build_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - patch_started)
            .count();
    require(complex.status == "pass", "patch complex build failed");
    const PatchMappingGridSpec spec = surface_grid(atoms, spacing);

    auto d1 = build_patch_ray_row_worklist(complex, spec, &eight);
    require(validate_patch_ray_row_worklist(complex, spec, d1),
            "D1 validation failed");
    const double d1_build_ms = d1.build_ms;
    const auto d2c = build_tiled_patch_ray_row_worklist(complex, spec, &eight);
    require(validate_tiled_patch_ray_row_worklist(complex, spec, d2c),
            "canonical D2c validation failed");
    double d2b_build_ms = 0.0;
    std::size_t d2b_candidates = 0U;
    std::size_t d2b_retained_bytes = 0U;
    std::size_t d2b_peak_bytes = 0U;
    {
      const auto d2b = build_torus_domain_row_worklist(complex, spec, &eight);
      require(validate_torus_domain_row_worklist(complex, spec, d2b),
              "D2b validation failed");
      require_same_rows_cross(d2b, d2c, "D2b/D2c CSR mismatch");
      d2b_build_ms = d2b.build_ms;
      d2b_candidates = d2b.candidate_pairs;
      d2b_retained_bytes = d2b.retained_bytes;
      d2b_peak_bytes = d2b.accounted_peak_temporary_bytes;
    }
    {
      const auto d2c_nt1 =
          build_tiled_patch_ray_row_worklist(complex, spec, &one);
      require(validate_tiled_patch_ray_row_worklist(complex, spec, d2c_nt1),
              "D2c nt1 validation failed");
      require_same_rows(d2c, d2c_nt1, "D2c nt1/8 CSR mismatch");
      require(checksum(d2c) == checksum(d2c_nt1),
              "D2c nt1/8 checksum mismatch");
    }
    {
      TiledPatchRayRowWorklistBuildOptions reverse_options;
      reverse_options.reverse_rectangle_schedule_for_audit = true;
      reverse_options.reverse_tile_schedule_for_audit = true;
      const auto d2c_reverse = build_tiled_patch_ray_row_worklist(
          complex, spec, &eight, reverse_options);
      require(validate_tiled_patch_ray_row_worklist(complex, spec, d2c_reverse),
              "D2c reverse validation failed");
      require_same_rows(d2c, d2c_reverse, "D2c reverse CSR mismatch");
      require(checksum(d2c) == checksum(d2c_reverse),
              "D2c reverse checksum mismatch");
    }
    const TiledTiming timing_one =
        benchmark_tiled(samples, complex, spec, &one, d2c);
    const TiledTiming timing_eight =
        benchmark_tiled(samples, complex, spec, &eight, d2c);

    const auto probe_index_started = Clock::now();
    const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
    const double probe_index_ms = std::chrono::duration<double, std::milli>(
                                      Clock::now() - probe_index_started)
                                      .count();
    const RemovedAudit removed = audit_removed_candidates(
        complex, spec, probes, d1, d2c, audit_budget_seconds);
    emit_first_unsafe_diagnostic(complex, spec, probes, removed);
    const std::string status =
        removed.complete && removed.required_removed_observed != 0U
            ? "false_negative"
        : !removed.complete ? "partial_audit_budget"
                            : "pass";
    const bool d1_denominator_zero = d1.candidate_pairs == 0U;
    const double reduction_percent =
        d1_denominator_zero
            ? 0.0
            : 100.0 *
                  (1.0 - static_cast<double>(d2c.deduplicated_candidate_pairs) /
                             static_cast<double>(d1.candidate_pairs));
    const double case_wall_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - case_started)
            .count();

    std::cout
        << std::setprecision(12)
        << "case_label,input_file,input_sha256_provenance,radii_provenance,"
           "spacing,max_atoms,atom_selection,atoms,"
           "contacts,tori,reentrant,patch_build_ms,nx,ny,nz,total_rows,"
           "samples,d1_candidates,d2b_candidates,d2c_candidates,"
           "d2c_raw_pairs,candidate_reduction_percent,d1_denominator_zero,"
           "d1_build_ms,"
           "d2b_build_ms,"
           "d2c_nt1_median_ms,d2c_nt8_median_ms,d2c_nt8_bounds_median_ms,"
           "d2c_nt8_tile_median_ms,d2c_nt8_count_median_ms,"
           "d2c_nt8_prefix_median_ms,d2c_nt8_fill_median_ms,d2c_rectangles,"
           "d2c_total_tiles,d2c_active_tiles,d2c_rectangle_tile_pairs,"
           "d2c_maximum_raw_pairs_per_tile,d2c_torus_domain_plans,"
           "d2c_torus_theta_fallback_plans,d2c_torus_full_fallback_plans,"
           "d2c_status,d2c_nt8_requested_workers,d2c_nt8_scheduled_workers,"
           "d2c_nt8_count_workers_used,d2c_nt8_fill_workers_used,"
           "d1_retained_bytes,"
           "d1_accounted_peak_bytes,d2b_retained_bytes,"
           "d2b_accounted_peak_bytes,d2c_retained_bytes,"
           "d2c_accounted_peak_bytes,d2c_checksum,d2b_d2c_csr_equal,"
           "d2c_nt_parity,d2c_reverse_parity,removed_total,"
           "removed_audited,audit_complete,required_removed,false_negatives,"
           "unsafe_roots_count_nonzero,unsafe_roots_uncertified,"
           "unsafe_roots_overlap,first_unsafe_axis,first_unsafe_row,"
           "first_unsafe_u,first_unsafe_v,first_unsafe_kind,"
           "first_unsafe_patch,first_unsafe_root_count,"
           "first_unsafe_roots_certified,first_unsafe_root_coordinates,"
           "first_unsafe_root_simplicity,probe_index_ms,"
           "removed_audit_ms,audit_claim,implementation_hypothesis,"
           "case_wall_ms,status\n"
        << case_label << ',' << file_name(input_path) << ',' << input_sha256
        << ',' << radii_provenance << ',' << spacing << ',' << max_atoms << ','
        << (max_atoms == 0U ? "full" : "prefix_cap") << ',' << atoms.size()
        << ',' << complex.contacts.size() << ',' << complex.tori.size() << ','
        << complex.reentrant.size() << ',' << patch_build_ms << ',' << spec.nx
        << ',' << spec.ny << ',' << spec.nz << ',' << total_rows(spec) << ','
        << samples << ',' << d1.candidate_pairs << ',' << d2b_candidates << ','
        << d2c.deduplicated_candidate_pairs << ',' << d2c.raw_row_pairs << ',';
    if (d1_denominator_zero)
      std::cout << "NA";
    else
      std::cout << reduction_percent;
    std::cout << ',' << (d1_denominator_zero ? 1 : 0) << ',' << d1_build_ms
              << ',' << d2b_build_ms << ',' << timing_one.build_ms << ','
              << timing_eight.build_ms << ',' << timing_eight.bounds_ms << ','
              << timing_eight.tile_ms << ',' << timing_eight.count_ms << ','
              << timing_eight.prefix_ms << ',' << timing_eight.fill_ms << ','
              << d2c.bound_rectangles << ',' << d2c.total_tiles << ','
              << d2c.active_tiles << ',' << d2c.rectangle_tile_pairs << ','
              << d2c.maximum_raw_pairs_per_tile << ',' << d2c.torus_domain_plans
              << ',' << d2c.torus_theta_fallback_plans << ','
              << d2c.torus_full_fallback_plans << ',' << d2c.status << ','
              << d2c.requested_workers << ',' << d2c.scheduled_workers << ','
              << d2c.count_workers_used << ',' << d2c.fill_workers_used << ','
              << d1.retained_bytes << ',' << d1.accounted_peak_temporary_bytes
              << ',' << d2b_retained_bytes << ',' << d2b_peak_bytes << ','
              << d2c.retained_bytes << ','
              << d2c.accounted_peak_temporary_upper_bound_bytes << ','
              << checksum(d2c) << ",1,1,1," << removed.removed << ','
              << removed.audited << ',' << (removed.complete ? 1 : 0) << ',';
    if (removed.complete) {
      std::cout << removed.required_removed_observed << ','
                << removed.required_removed_observed;
    } else {
      std::cout << "NA,NA";
    }
    std::cout << ',' << removed.roots_count_nonzero << ','
              << removed.roots_uncertified << ','
              << removed.roots_count_nonzero_and_uncertified << ',';
    if (removed.first_unsafe.present) {
      const RemovedAudit::FirstUnsafeWitness &witness = removed.first_unsafe;
      std::ostringstream coordinates;
      std::ostringstream simplicity;
      coordinates << std::setprecision(17);
      for (std::size_t root = 0U; root < witness.roots.count; ++root) {
        if (root != 0U) {
          coordinates << '|';
          simplicity << '|';
        }
        coordinates << witness.roots.coordinates[root];
        simplicity << static_cast<unsigned int>(
            witness.roots.primitive_simplicity_certified[root]);
      }
      std::cout << static_cast<unsigned int>(witness.axis) << ',' << witness.row
                << ',' << std::setprecision(17) << witness.decoded.u << ','
                << witness.decoded.v << ',' << patch_kind_name(witness.kind)
                << ',' << witness.patch_index << ',' << witness.roots.count
                << ',' << (witness.roots.certified ? 1 : 0) << ','
                << coordinates.str() << ',' << simplicity.str();
    } else {
      std::cout << "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
    }
    std::cout << ',' << probe_index_ms << ',' << removed.elapsed_ms << ','
              << (removed.complete ? "complete" : "unavailable") << ','
              << "theta_then_domain_recompute_baseline" << ',' << case_wall_ms
              << ',' << status << '\n';
    return status == "pass" ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_row_worklist_real_benchmark=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
