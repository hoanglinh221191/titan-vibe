#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_ray.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

std::vector<Atom> translate(std::vector<Atom> atoms, double x, double y,
                            double z) {
  for (Atom &atom : atoms) {
    atom.x += x;
    atom.y += y;
    atom.z += z;
  }
  return atoms;
}

std::size_t total_odd(const PatchRayGrid &grid) {
  std::size_t result = 0;
  for (const auto &axis : grid.axes)
    result += axis.odd_rows;
  return result;
}

std::uint64_t grid_checksum(const PatchRayGrid &grid) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(grid.raw_events);
  mix(grid.duplicate_events);
  for (const AxisRayCrossings &axis : grid.axes) {
    mix(static_cast<std::uint64_t>(axis.axis));
    mix(axis.u_count);
    mix(axis.v_count);
    mix(axis.odd_rows);
    for (const std::size_t offset : axis.offsets)
      mix(offset);
    for (const double coordinate : axis.coordinates)
      mix(std::bit_cast<std::uint64_t>(coordinate));
  }
  return hash;
}

void require_same_grid_contract(const PatchRayGrid &expected,
                                const PatchRayGrid &actual,
                                const char *message) {
  require(expected.status == actual.status, message);
  require(expected.raw_events == actual.raw_events, message);
  require(expected.duplicate_events == actual.duplicate_events, message);
  require(expected.uncertified_intersection_queries ==
              actual.uncertified_intersection_queries,
          message);
  require(expected.uncertified_rows_repaired ==
              actual.uncertified_rows_repaired,
          message);
  require(expected.unresolved_uncertified_rows ==
              actual.unresolved_uncertified_rows,
          message);
  require(expected.degenerate_rows_repaired ==
              actual.degenerate_rows_repaired,
          message);
  for (std::size_t axis = 0U; axis < expected.axes.size(); ++axis) {
    require(expected.axes[axis].axis == actual.axes[axis].axis, message);
    require(expected.axes[axis].u_count == actual.axes[axis].u_count,
            message);
    require(expected.axes[axis].v_count == actual.axes[axis].v_count,
            message);
    require(expected.axes[axis].odd_rows == actual.axes[axis].odd_rows,
            message);
    require(expected.axes[axis].offsets == actual.axes[axis].offsets,
            message);
    require(expected.axes[axis].coordinates ==
                actual.axes[axis].coordinates,
            message);
  }
  require(grid_checksum(expected) == grid_checksum(actual), message);
}

void require_same_provenance(const PatchRayProvenance &expected,
                             const PatchRayProvenance &actual,
                             const char *message) {
  require(expected.enabled == actual.enabled, message);
  require(expected.certified == actual.certified, message);
  require(expected.status == actual.status, message);
  require(expected.lattice_identity == actual.lattice_identity, message);
  require(expected.unique_crossings == actual.unique_crossings, message);
  require(expected.seam_groups == actual.seam_groups, message);
  require(expected.refs_total == actual.refs_total, message);
  require(expected.max_refs_per_crossing ==
              actual.max_refs_per_crossing,
          message);
  require(expected.roots_total == actual.roots_total, message);
  require(expected.repeated_root_groups ==
              actual.repeated_root_groups,
          message);
  require(expected.same_primitive_near_pair_unclassified ==
              actual.same_primitive_near_pair_unclassified,
          message);
  require(expected.retained_bytes == actual.retained_bytes, message);
  require(expected.coordinate_merge_tolerance ==
              actual.coordinate_merge_tolerance,
          message);
  require(expected.near_pair_tolerance ==
              actual.near_pair_tolerance,
          message);
  for (std::size_t axis = 0U; axis < expected.axes.size(); ++axis) {
    require(expected.axes[axis].canonical_refs ==
                actual.axes[axis].canonical_refs,
            message);
    require(expected.axes[axis].seam_groups ==
                actual.axes[axis].seam_groups,
            message);
    require(expected.axes[axis].seam_refs ==
                actual.axes[axis].seam_refs,
            message);
    require(expected.axes[axis].seam_ref_multiplicities ==
                actual.axes[axis].seam_ref_multiplicities,
            message);
    require(expected.axes[axis].single_ref_multiplicities ==
                actual.axes[axis].single_ref_multiplicities,
            message);
  }
}

PatchRefId make_patch_ref(PatchKind kind, std::size_t index) {
  PatchRefId ref;
  require(PatchRefId::try_make(kind, index, ref) && ref.valid(),
          "failed to pack audit patch ref");
  return ref;
}

PatchRayGrid make_synthetic_audit_grid(
    const std::vector<double> &coordinates,
    const std::vector<std::vector<PatchRefId>> &refs,
    double near_pair_tolerance = 1.0e-8) {
  require(coordinates.size() == refs.size(),
          "synthetic audit coordinate/ref mismatch");
  PatchRayGrid grid;
  for (std::size_t axis = 0U; axis < grid.axes.size(); ++axis) {
    grid.axes[axis].axis = static_cast<RayAxis>(axis);
    grid.axes[axis].u_count = 1U;
    grid.axes[axis].v_count = 1U;
    grid.axes[axis].offsets = {
        0U, axis == 0U ? coordinates.size() : 0U};
  }
  grid.axes[0].coordinates = coordinates;
  grid.status = "pass";
  grid.provenance.enabled = true;
  grid.provenance.certified = true;
  grid.provenance.status = "pass";
  grid.provenance.coordinate_merge_tolerance = 1.0e-10;
  grid.provenance.near_pair_tolerance = near_pair_tolerance;
  grid.provenance.unique_crossings = coordinates.size();
  grid.provenance.roots_total = 0U;
  for (const std::vector<PatchRefId> &group : refs) {
    require(!group.empty() &&
                grid.provenance.axes[0].append(group),
            "synthetic audit provenance append failed");
    grid.provenance.refs_total += group.size();
    grid.provenance.roots_total += group.size();
    grid.provenance.max_refs_per_crossing =
        std::max(grid.provenance.max_refs_per_crossing,
                 group.size());
    grid.provenance.seam_groups +=
        static_cast<std::size_t>(group.size() > 1U);
  }
  return grid;
}

} // namespace

int main() {
  try {
    const double root3 = std::sqrt(3.0);
    const std::vector<Atom> atoms = {
        {-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
        {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
        {0.0, root3, 0.0, 0.0, 1.0}};
    const PatchComplex complex = build_patch_complex(atoms, 1.0);
    require(complex.status == "pass", "patch complex build failed");
    TitanExecutor executor(4);
    // Avoid deliberately routing grid rays through atom centres or patch
    // vertices; those are tangent diagnostics rather than parity crossings.
    const Vec3 lower{-5.137, -5.071, -5.193};
    const PatchRayGrid grid =
        build_patch_ray_grid(complex, lower, 0.5, 21, 21, 21, &executor);
    require(grid.status == "pass", "patch-major grid build failed");
    require(grid.raw_events > 0 && grid.retained_bytes > 0,
            "patch-major grid is empty");
    if (total_odd(grid) != 0) {
      std::cerr << "odd_rows x=" << grid.axes[0].odd_rows
                << " y=" << grid.axes[1].odd_rows
                << " z=" << grid.axes[2].odd_rows
                << " raw=" << grid.raw_events
                << " duplicates=" << grid.duplicate_events << '\n';
    }
    require(total_odd(grid) == 0, "closed triple has odd ray rows");

    PatchRayBuildOptions provenance_options;
    provenance_options.retain_provenance = true;
    TitanExecutor executor_one(1);
    TitanExecutor executor_eight(8);
    const PatchRayGrid provenance_one = build_patch_ray_grid(
        complex, lower, 0.5, 21, 21, 21, &executor_one,
        provenance_options);
    const PatchRayGrid provenance_eight = build_patch_ray_grid(
        complex, lower, 0.5, 21, 21, 21, &executor_eight,
        provenance_options);
    require_same_grid_contract(
        grid, provenance_one,
        "provenance-on changed patch-ray coordinates/status");
    require_same_grid_contract(
        provenance_one, provenance_eight,
        "provenance grid is not deterministic at nt=1/8");
    require(provenance_one.provenance.enabled &&
                provenance_one.provenance.certified &&
                provenance_one.provenance.status == "pass" &&
                !grid.provenance.lattice_identity.valid &&
                provenance_one.provenance.lattice_identity ==
                    make_patch_ray_lattice_identity(lower, 0.5, 21U, 21U,
                                                    21U),
            "closed triple provenance was not certified");
    require_same_provenance(
        provenance_one.provenance, provenance_eight.provenance,
        "provenance payload is not deterministic at nt=1/8");
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      require(provenance_one.provenance.axes[axis].validate(
                  provenance_one.axes[axis].coordinates.size()),
              "provenance CSR validation failed");
      for (std::size_t crossing = 0U;
           crossing <
           provenance_one.provenance.axes[axis].canonical_refs.size();
           ++crossing) {
        const auto refs =
            provenance_one.provenance.axes[axis].refs(crossing);
        require(!refs.empty() &&
                    std::is_sorted(refs.begin(), refs.end()),
                "crossing refs are not sorted/nonempty");
      }
    }

    PatchRefId contact_zero;
    PatchRefId contact_one;
    require(PatchRefId::try_make(PatchKind::contact, 0U, contact_zero) &&
                PatchRefId::try_make(PatchKind::contact, 1U, contact_one),
            "failed to create test patch refs");
    PatchRefId overflow_ref;
    require(!PatchRefId::try_make(
                PatchKind::contact,
                static_cast<std::size_t>(PatchRefId::index_mask) + 1U,
                overflow_ref) &&
                !overflow_ref.valid(),
            "PatchRefId index overflow did not fail closed");

    const LocalPatchProjectionIndex projection_index =
        build_local_patch_projection_index(
            complex, LocalPatchProjectionIndexBuildOptions{true});
    const AnalyticPatchAdjacencyOracle oracle =
        build_analytic_patch_adjacency_oracle(complex,
                                              projection_index);
    require(oracle.certified() && oracle.status() == "pass",
            "patch-ray audit oracle was not certified");
    const PatchRayProvenanceAudit real_audit =
        audit_patch_ray_provenance(provenance_one, oracle);
    require(real_audit.status != "not_run" &&
                real_audit.seam_groups_audited ==
                    provenance_one.provenance.seam_groups &&
                !real_audit.crossing_location_certified &&
                !real_audit.active_cell_mapping_certified &&
                !real_audit.mesh_mapping_certified &&
                real_audit.active_cell_mapping_status ==
                    "active_cell_mapping_uncertified" &&
                real_audit.mesh_mapping_status ==
                    "mesh_mapping_uncertified" &&
                real_audit.oracle_build_ms == oracle.build_ms() &&
                real_audit.oracle_retained_bytes ==
                    oracle.retained_payload_bytes() &&
                real_audit.audit_ms >= 0.0,
            "real provenance audit overclaimed mapping/location");

    const std::vector<Atom> pair_atoms = {
        {-1.5, 0.0, 0.0, 0.0, 1.0},
        {1.5, 0.0, 0.0, 0.0, 1.0}};
    const PatchComplex pair_complex =
        build_patch_complex(pair_atoms, 1.0, nullptr, false);
    const LocalPatchProjectionIndex pair_projection =
        build_local_patch_projection_index(
            pair_complex, LocalPatchProjectionIndexBuildOptions{true});
    const AnalyticPatchAdjacencyOracle pair_oracle =
        build_analytic_patch_adjacency_oracle(pair_complex,
                                              pair_projection);
    const PatchRefId pair_torus =
        make_patch_ref(PatchKind::toroidal, 0U);
    require(pair_oracle.certified(),
            "synthetic audit pair oracle was not certified");

    const PatchRayGrid curve_grid = make_synthetic_audit_grid(
        {0.0}, {{contact_zero, pair_torus}});
    const PatchRayProvenanceAudit curve_audit =
        audit_patch_ray_provenance(curve_grid, pair_oracle);
    require(curve_audit.relation_audit_complete &&
                !curve_audit.crossing_location_certified &&
                curve_audit.curve_adjacent_groups == 1U &&
                curve_audit.relation_counts[static_cast<std::size_t>(
                    PatchRelation::curve_adjacent)] == 1U &&
                curve_audit.status ==
                    "crossing_location_uncertified",
            "curve relation audit overclaimed seam location");

    const PatchRefId triple_reentrant =
        make_patch_ref(PatchKind::reentrant, 0U);
    const std::uint32_t triple_atom =
        static_cast<std::uint32_t>(
            complex.reentrant[0].atoms[0]);
    const PatchRefId triple_contact =
        make_patch_ref(PatchKind::contact, triple_atom);
    const auto &triple_seams = oracle.seam_certificate();
    const PatchRefId triple_torus_first =
        triple_seams.reentrant_edges[0].torus;
    const PatchRefId triple_torus_second =
        triple_seams.reentrant_edges[2].torus;
    const PatchRayGrid point_grid = make_synthetic_audit_grid(
        {0.0}, {{triple_contact, triple_reentrant}});
    const PatchRayProvenanceAudit point_audit =
        audit_patch_ray_provenance(point_grid, oracle);
    require(point_audit.relation_audit_complete &&
                point_audit.point_only_groups == 1U &&
                !point_audit.crossing_location_certified &&
                point_audit.status ==
                    "crossing_location_uncertified",
            "point-only audit authorized an unlocated crossing");

    const PatchRayGrid multi_junction_grid =
        make_synthetic_audit_grid(
            {0.0},
            {{triple_torus_first, triple_torus_second}});
    const PatchRayProvenanceAudit multi_junction_audit =
        audit_patch_ray_provenance(multi_junction_grid, oracle);
    require(!multi_junction_audit.relation_audit_complete &&
                multi_junction_audit.uncertified_crossing_groups == 1U &&
                multi_junction_audit.status ==
                    "uncertified_relation",
            "multi-junction torus pair was classified as point-only");

    const PatchRayGrid separated_grid =
        make_synthetic_audit_grid(
            {0.0}, {{contact_zero, contact_one}});
    const PatchRayProvenanceAudit separated_audit =
        audit_patch_ray_provenance(separated_grid, pair_oracle);
    require(separated_audit.relation_audit_complete &&
                separated_audit.incompatible_crossing_groups == 1U &&
                separated_audit.status == "incompatible_crossing",
            "separated crossing group was not rejected");

    const PatchRefId absent_reentrant =
        make_patch_ref(PatchKind::reentrant, 0U);
    const PatchRayGrid mixed_near_grid =
        make_synthetic_audit_grid(
            {0.0, 1.0e-12},
            {{contact_zero, absent_reentrant}, {contact_one}},
            1.0e-8);
    const PatchRayProvenanceAudit mixed_near_audit =
        audit_patch_ray_provenance(mixed_near_grid, pair_oracle);
    require(mixed_near_audit.incompatible_near_pairs == 1U &&
                mixed_near_audit.uncertified_near_pairs == 1U &&
                !mixed_near_audit.relation_audit_complete,
            "multi-ref near pair lost separated/uncertified evidence");

    PatchRayGrid merged_same_primitive =
        make_synthetic_audit_grid({0.0}, {{contact_zero}});
    merged_same_primitive.provenance.axes[0]
        .single_ref_multiplicities.push_back({0U, 2U});
    merged_same_primitive.provenance.roots_total = 2U;
    merged_same_primitive.provenance.repeated_root_groups = 1U;
    merged_same_primitive.provenance
        .same_primitive_near_pair_unclassified = 1U;
    const PatchRayProvenanceAudit merged_same_audit =
        audit_patch_ray_provenance(merged_same_primitive,
                                   pair_oracle);
    require(merged_same_audit.same_primitive_merged_groups == 1U &&
                merged_same_audit
                        .same_primitive_adjacent_near_pairs ==
                    0U &&
                merged_same_audit
                        .same_primitive_near_pair_unclassified ==
                    1U &&
                merged_same_audit.status ==
                    "same_primitive_unclassified",
            "merged same-primitive roots were mislabeled");

    PatchRayGrid seam_multiplicity =
        make_synthetic_audit_grid(
            {0.0}, {{contact_zero, pair_torus}});
    AxisPatchProvenance seam_multiplicity_axis;
    const std::array<PatchRefId, 2> seam_refs = {
        contact_zero, pair_torus};
    const std::array<std::uint32_t, 2> seam_multiplicities = {
        2U, 1U};
    require(seam_multiplicity_axis.append(
                seam_refs, seam_multiplicities),
            "seam multiplicity fixture setup failed");
    seam_multiplicity.provenance.axes[0] =
        std::move(seam_multiplicity_axis);
    seam_multiplicity.provenance.roots_total = 3U;
    seam_multiplicity.provenance.repeated_root_groups = 1U;
    seam_multiplicity.provenance
        .same_primitive_near_pair_unclassified = 1U;
    const PatchRayProvenanceAudit seam_multiplicity_audit =
        audit_patch_ray_provenance(seam_multiplicity,
                                   pair_oracle);
    require(seam_multiplicity_audit
                    .same_primitive_merged_groups ==
                1U &&
                seam_multiplicity_audit
                        .same_primitive_adjacent_near_pairs ==
                    0U &&
                seam_multiplicity_audit
                        .same_primitive_near_pair_unclassified ==
                    1U &&
                seam_multiplicity_audit.status ==
                    "same_primitive_unclassified",
            "seam-group repeated root was treated as certified");

    const PatchRayGrid adjacent_same_primitive =
        make_synthetic_audit_grid(
            {0.0, 1.0e-12},
            {{contact_zero}, {contact_zero}}, 1.0e-8);
    const PatchRayProvenanceAudit adjacent_same_audit =
        audit_patch_ray_provenance(adjacent_same_primitive,
                                   pair_oracle);
    require(adjacent_same_audit.same_primitive_merged_groups == 0U &&
                adjacent_same_audit
                        .same_primitive_adjacent_near_pairs ==
                    1U &&
                adjacent_same_audit
                        .same_primitive_near_pair_unclassified ==
                    1U,
            "adjacent same-primitive roots were mislabeled");

    PatchRayGrid stale_unique = curve_grid;
    ++stale_unique.provenance.unique_crossings;
    require(audit_patch_ray_provenance(stale_unique, pair_oracle)
                    .status == "invalid_provenance",
            "stale unique-crossing telemetry was accepted");
    PatchRayGrid stale_seams = curve_grid;
    ++stale_seams.provenance.seam_groups;
    require(audit_patch_ray_provenance(stale_seams, pair_oracle)
                    .status == "invalid_provenance",
            "stale seam-group telemetry was accepted");
    PatchRayGrid nonfinite_coordinate = curve_grid;
    nonfinite_coordinate.axes[0].coordinates[0] =
        std::numeric_limits<double>::quiet_NaN();
    require(audit_patch_ray_provenance(nonfinite_coordinate,
                                       pair_oracle)
                    .status == "invalid_provenance",
            "nonfinite crossing coordinate was accepted");
    PatchRayGrid equal_coordinates = adjacent_same_primitive;
    equal_coordinates.axes[0].coordinates[1] = 0.0;
    require(audit_patch_ray_provenance(equal_coordinates,
                                       pair_oracle)
                    .status == "invalid_provenance",
            "nonstrict crossing coordinates were accepted");
    PatchRayGrid invalid_merge_tolerance = curve_grid;
    invalid_merge_tolerance.provenance.coordinate_merge_tolerance =
        std::numeric_limits<double>::quiet_NaN();
    require(audit_patch_ray_provenance(
                invalid_merge_tolerance, pair_oracle)
                    .status == "invalid_provenance",
            "nonfinite merge tolerance was accepted");

    AxisPatchProvenance multiplicity_fixture;
    const std::array<PatchRefId, 1> one_ref = {contact_zero};
    const std::array<std::uint32_t, 1> two_roots = {2U};
    require(multiplicity_fixture.append(one_ref, two_roots) &&
                multiplicity_fixture.validate(1U) &&
                multiplicity_fixture.ref_multiplicity(0U, 0U) == 2U,
            "same-ref root multiplicity was not retained");
    std::size_t multiplicity_count = 0U;
    require(multiplicity_fixture.root_event_count(
                0U, multiplicity_count) &&
                multiplicity_count == 2U,
            "same-ref root multiplicity count is inconsistent");

    AxisPatchProvenance transactional_fixture;
    require(transactional_fixture.append(one_ref),
            "transactional fixture setup failed");
    const AxisPatchProvenance transactional_before =
        transactional_fixture;
    const std::array<PatchRefId, 2> two_refs = {
        contact_zero, contact_one};
    const std::array<std::uint32_t, 2> unit_multiplicity = {1U, 1U};
    require(!transactional_fixture.append_with_limits_for_test(
                two_refs, unit_multiplicity, 8U, 1U, 8U),
            "seam-ref overflow test unexpectedly succeeded");
    require(transactional_fixture.canonical_refs ==
                transactional_before.canonical_refs &&
                transactional_fixture.seam_groups ==
                    transactional_before.seam_groups &&
                transactional_fixture.seam_refs ==
                    transactional_before.seam_refs &&
                transactional_fixture.seam_ref_multiplicities ==
                    transactional_before.seam_ref_multiplicities &&
                transactional_fixture.single_ref_multiplicities ==
                    transactional_before.single_ref_multiplicities &&
                transactional_fixture.validate(1U),
            "failed append left a partial provenance object");

    AxisPatchProvenance zero_multiplicity_fixture =
        transactional_fixture;
    const AxisPatchProvenance zero_multiplicity_before =
        zero_multiplicity_fixture;
    const std::array<std::uint32_t, 1> zero_multiplicity = {0U};
    require(!zero_multiplicity_fixture.append(
                one_ref, zero_multiplicity) &&
                zero_multiplicity_fixture.canonical_refs ==
                    zero_multiplicity_before.canonical_refs &&
                zero_multiplicity_fixture.seam_groups ==
                    zero_multiplicity_before.seam_groups &&
                zero_multiplicity_fixture.seam_refs ==
                    zero_multiplicity_before.seam_refs &&
                zero_multiplicity_fixture.seam_ref_multiplicities ==
                    zero_multiplicity_before
                        .seam_ref_multiplicities &&
                zero_multiplicity_fixture.single_ref_multiplicities ==
                    zero_multiplicity_before
                        .single_ref_multiplicities &&
                zero_multiplicity_fixture.validate(1U),
            "zero root multiplicity did not fail transactionally");
    std::uint32_t sticky_multiplicity =
        std::numeric_limits<std::uint32_t>::max();
    patch_ray_increment_root_multiplicity_for_test(
        sticky_multiplicity);
    patch_ray_increment_root_multiplicity_for_test(
        sticky_multiplicity);
    require(sticky_multiplicity == 0U,
            "root multiplicity overflow sentinel was not sticky");

    AxisPatchProvenance corrupted_fixture = transactional_fixture;
    corrupted_fixture.seam_ref_multiplicities.push_back(1U);
    require(!corrupted_fixture.validate(1U),
            "corrupt multiplicity CSR did not fail validation");

    PatchRayBuildOptions allocation_failure_options =
        provenance_options;
    allocation_failure_options.simulate_allocation_failure = true;
    const PatchRayGrid allocation_failure = build_patch_ray_grid(
        complex, lower, 0.5, 21, 21, 21, &executor_one,
        allocation_failure_options);
    require(allocation_failure.status == "alloc_fail" &&
                allocation_failure.provenance.enabled &&
                !allocation_failure.provenance.certified,
            "allocation failure did not return a fail-closed audit grid");
    for (const AxisPatchProvenance &axis :
         allocation_failure.provenance.axes) {
      require(axis.canonical_refs.empty() && axis.seam_groups.empty() &&
                  axis.seam_refs.empty() &&
                  axis.seam_ref_multiplicities.empty() &&
                  axis.single_ref_multiplicities.empty(),
              "allocation failure leaked partial provenance");
    }

    const PatchRayGrid row_overflow = build_patch_ray_grid(
        complex, lower, 0.5, std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(), 2U, &executor_one,
        provenance_options);
    require(row_overflow.status == "ray_grid_too_large" &&
                !row_overflow.provenance.certified,
            "row-count overflow did not fail closed");
    const PatchRayGrid invalid_spacing = build_patch_ray_grid(
        complex, lower, std::numeric_limits<double>::quiet_NaN(), 21, 21,
        21, &executor_one, provenance_options);
    require(invalid_spacing.status == "invalid_input" &&
                !invalid_spacing.provenance.certified,
            "nonfinite spacing did not fail closed");

    const std::vector<Atom> shifted_atoms =
        translate(atoms, 1.0e6, -17.0, 29.0);
    const PatchComplex shifted = build_patch_complex(shifted_atoms, 1.0);
    const Vec3 shifted_lower{lower.x + 1.0e6, lower.y - 17.0,
                             lower.z + 29.0};
    const PatchRayGrid shifted_grid = build_patch_ray_grid(
        shifted, shifted_lower, 0.5, 21, 21, 21, &executor);
    require(shifted_grid.status == "pass",
            "translated patch-major grid build failed");
    require(total_odd(shifted_grid) == 0,
            "translation introduced odd ray rows");
    for (std::size_t axis = 0; axis < 3; ++axis) {
      require(grid.axes[axis].coordinates.size() ==
                  shifted_grid.axes[axis].coordinates.size(),
              "translation changed crossing count");
      const double shift = axis == 0 ? 1.0e6 : (axis == 1 ? -17.0 : 29.0);
      for (std::size_t i = 0; i < grid.axes[axis].coordinates.size(); ++i) {
        require(std::abs(shifted_grid.axes[axis].coordinates[i] - shift -
                         grid.axes[axis].coordinates[i]) < 2.0e-8,
                "translation changed patch-major crossing");
      }
    }

    const std::vector<Atom> residue = {
        {-6.406, 5.469, -3.259, -0.3000, 1.8500},
        {-5.676, 6.671, -2.760, 0.2100, 2.2750},
        {-4.542, 6.236, -1.802, 0.5100, 2.0000},
        {-3.896, 5.233, -2.028, -0.5100, 1.7000},
        {-6.754, 7.501, -2.044, -0.1800, 2.1750},
        {-7.649, 6.581, -1.205, -0.1400, 2.1750},
        {-9.351, 6.657, -1.829, -0.0900, 2.0000},
        {-9.441, 4.973, -2.494, -0.2200, 2.0600},
        {-7.379, 5.677, -3.311, 0.3300, 0.2245},
        {-5.260, 7.122, -3.545, 0.1000, 1.3200}};
    const double residue_probe = 1.4;
    const double residue_spacing = 0.75;
    Vec3 residue_lower{1.0e9, 1.0e9, 1.0e9};
    Vec3 residue_upper{-1.0e9, -1.0e9, -1.0e9};
    for (const Atom &atom : residue) {
      const double radius = atom.radius + residue_probe;
      residue_lower.x = std::min(residue_lower.x, atom.x - radius);
      residue_lower.y = std::min(residue_lower.y, atom.y - radius);
      residue_lower.z = std::min(residue_lower.z, atom.z - radius);
      residue_upper.x = std::max(residue_upper.x, atom.x + radius);
      residue_upper.y = std::max(residue_upper.y, atom.y + radius);
      residue_upper.z = std::max(residue_upper.z, atom.z + radius);
    }
    residue_lower.x -= 2.0 * residue_spacing;
    residue_lower.y -= 2.0 * residue_spacing;
    residue_lower.z -= 2.0 * residue_spacing;
    residue_upper.x += 2.0 * residue_spacing;
    residue_upper.y += 2.0 * residue_spacing;
    residue_upper.z += 2.0 * residue_spacing;
    const auto count = [&](double lo, double hi) {
      return static_cast<std::size_t>(
                 std::ceil((hi - lo) / residue_spacing)) +
             1U;
    };
    const std::size_t residue_nx = count(residue_lower.x, residue_upper.x);
    const std::size_t residue_ny = count(residue_lower.y, residue_upper.y);
    const std::size_t residue_nz = count(residue_lower.z, residue_upper.z);
    const PatchComplex residue_complex =
        build_patch_complex(residue, residue_probe);
    const PatchRayGrid residue_grid = build_patch_ray_grid(
        residue_complex, residue_lower, residue_spacing, residue_nx,
        residue_ny, residue_nz, &executor);
    require(residue_complex.status == "pass",
            "unequal-radius residue patch complex failed");
    require(residue_grid.status == "pass",
            "unequal-radius residue patch-ray grid failed");
    std::size_t legacy_x_odd = 0U;
    std::size_t x_parity_disagreements = 0U;
    for (std::size_t iz = 0; iz < residue_nz; ++iz) {
      const double z = residue_lower.z + residue_spacing * iz;
      for (std::size_t iy = 0; iy < residue_ny; ++iy) {
        const double y = residue_lower.y + residue_spacing * iy;
        const auto legacy = intersect_ray_x(residue_complex, y, z);
        legacy_x_odd += legacy.size() & 1U;
        const auto compact = residue_grid.axes[0].row(iy + residue_ny * iz);
        x_parity_disagreements +=
            static_cast<std::size_t>((legacy.size() & 1U) !=
                                     (compact.size() & 1U));
      }
    }
    require(total_odd(residue_grid) == 0U,
            "unequal-radius torus selected a solvent-facing major arc");
    require(legacy_x_odd == 0U,
            "legacy unequal-radius residue has odd x-ray rows");
    require(x_parity_disagreements == 0U,
            "compact and legacy x-ray parity disagree");

    // A near-tangent unequal-radius pair produces a spindle torus.  The
    // quartic implicitization also contains its self-intersecting inner sheet;
    // accepting that algebraic branch leaves odd ray rows on real proteins.
    const std::vector<Atom> spindle_atoms = {
        {-0.493, 5.059, 7.359, 0.0, 1.320},
        {4.041, 1.752, 9.679, 0.0, 2.275}};
    const PatchComplex spindle =
        build_patch_complex(spindle_atoms, residue_probe);
    require(spindle.status == "pass" && spindle.tori.size() == 1U,
            "spindle patch complex build failed");
    require(spindle.tori[0].circle_radius < residue_probe,
            "spindle fixture no longer exercises a self-intersecting torus");
    const PatchRayGrid spindle_grid = build_patch_ray_grid(
        spindle, {-3.713, -2.423, 4.139}, 0.25, 49, 44, 40, &executor);
    require(spindle_grid.status == "pass", "spindle ray grid failed");
    require(total_odd(spindle_grid) == 0U,
            "spindle algebraic inner sheet created odd ray rows");

    const std::vector<Atom> single_atom = {
        {0.0, 0.0, 0.0, 0.0, 1.0}};
    PatchComplex duplicated_contact =
        build_patch_complex(single_atom, 1.0);
    require(duplicated_contact.status == "pass" &&
                duplicated_contact.contacts.size() == 1U,
            "single-contact provenance fixture build failed");
    duplicated_contact.contacts.push_back(
        duplicated_contact.contacts.front());
    const Vec3 duplicate_lower{-1.371, -1.293, -1.417};
    const PatchRayGrid duplicate_off = build_patch_ray_grid(
        duplicated_contact, duplicate_lower, 0.4, 9, 9, 9,
        &executor_one);
    const PatchRayGrid duplicate_one = build_patch_ray_grid(
        duplicated_contact, duplicate_lower, 0.4, 9, 9, 9,
        &executor_one, provenance_options);
    const PatchRayGrid duplicate_eight = build_patch_ray_grid(
        duplicated_contact, duplicate_lower, 0.4, 9, 9, 9,
        &executor_eight, provenance_options);
    require_same_grid_contract(
        duplicate_off, duplicate_one,
        "cross-patch provenance merge changed coordinates");
    require_same_grid_contract(
        duplicate_one, duplicate_eight,
        "cross-patch provenance merge is not nt1/nt8 deterministic");
    require_same_provenance(
        duplicate_one.provenance, duplicate_eight.provenance,
        "cross-patch provenance payload is not nt1/nt8 deterministic");
    require(duplicate_one.provenance.certified &&
                duplicate_one.provenance.unique_crossings > 0U &&
                duplicate_one.provenance.seam_groups ==
                    duplicate_one.provenance.unique_crossings &&
                duplicate_one.provenance.refs_total ==
                    2U * duplicate_one.provenance.unique_crossings,
            "cross-patch near groups were not retained sparsely");
    for (const AxisPatchProvenance &axis :
         duplicate_one.provenance.axes) {
      for (std::size_t crossing = 0U;
           crossing < axis.canonical_refs.size(); ++crossing) {
        const auto refs = axis.refs(crossing);
        require(refs.size() == 2U && refs[0] == contact_zero &&
                    refs[1] == contact_one &&
                    axis.ref_multiplicity(crossing, 0U) == 1U &&
                    axis.ref_multiplicity(crossing, 1U) == 1U,
                "cross-patch seam group lost identity/multiplicity");
      }
    }

    PatchRayGrid synthetic_odd_off;
    synthetic_odd_off.axes[0].axis = RayAxis::x;
    synthetic_odd_off.axes[0].u_count = 1U;
    synthetic_odd_off.axes[0].v_count = 1U;
    synthetic_odd_off.axes[0].offsets = {0U, 1U};
    synthetic_odd_off.axes[0].coordinates = {0.0};
    synthetic_odd_off.axes[0].odd_rows = 1U;
    synthetic_odd_off.raw_events = 1U;
    synthetic_odd_off.status = "pass";
    PatchRayGrid synthetic_odd_on = synthetic_odd_off;
    synthetic_odd_on.provenance.enabled = true;
    synthetic_odd_on.provenance.status = "building";
    require(synthetic_odd_on.provenance.axes[0].append(one_ref),
            "synthetic odd-row provenance setup failed");
    patch_ray_finalize_provenance_for_test(synthetic_odd_on);
    require_same_grid_contract(
        synthetic_odd_off, synthetic_odd_on,
        "odd-row finalization changed production grid contract");
    require(synthetic_odd_on.status == "pass" &&
                !synthetic_odd_on.provenance.certified &&
                synthetic_odd_on.provenance.status ==
                    "unresolved_odd_rows",
            "odd-row finalization overclaimed provenance");

    bool exercised_perturbed_repair = false;
    std::size_t perturbed_repair_rows = 0U;
    for (const ReentrantPatch &patch : complex.reentrant) {
      for (const Vec3 &contact_direction : patch.contact_directions) {
        const Vec3 seam_point{
            patch.probe_center.x +
                complex.probe_radius * contact_direction.x,
            patch.probe_center.y +
                complex.probe_radius * contact_direction.y,
            patch.probe_center.z +
                complex.probe_radius * contact_direction.z};
        const Vec3 aligned_lower{seam_point.x - 3.0,
                                 seam_point.y - 3.0,
                                 seam_point.z - 3.0};
        const PatchRayGrid aligned_off = build_patch_ray_grid(
            complex, aligned_lower, 0.5, 13, 13, 13,
            &executor_one);
        if (aligned_off.degenerate_rows_repaired == 0U)
          continue;
        const PatchRayGrid aligned_one = build_patch_ray_grid(
            complex, aligned_lower, 0.5, 13, 13, 13,
            &executor_one, provenance_options);
        const PatchRayGrid aligned_eight = build_patch_ray_grid(
            complex, aligned_lower, 0.5, 13, 13, 13,
            &executor_eight, provenance_options);
        require_same_grid_contract(
            aligned_off, aligned_one,
            "perturbed-row provenance changed coordinates/status");
        require_same_grid_contract(
            aligned_one, aligned_eight,
            "perturbed-row repair is not nt1/nt8 deterministic");
        require_same_provenance(
            aligned_one.provenance, aligned_eight.provenance,
            "perturbed-row provenance is not nt1/nt8 deterministic");
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          require(aligned_one.provenance.axes[axis].validate(
                      aligned_one.axes[axis].coordinates.size()),
                  "perturbed-row provenance CSR is invalid");
        }
        if (total_odd(aligned_off) == 0U &&
            aligned_off.unresolved_uncertified_rows == 0U) {
          require(aligned_one.provenance.certified &&
                      aligned_one.provenance.status == "pass",
                  "resolved perturbed row did not retain a certificate");
          exercised_perturbed_repair = true;
          perturbed_repair_rows =
              aligned_off.degenerate_rows_repaired;
          break;
        }
      }
      if (exercised_perturbed_repair)
        break;
    }
    require(exercised_perturbed_repair,
            "no deterministic perturbed-row repair fixture was found");

    std::cout << "tabi_ses_patch_ray_smoke=pass events=" << grid.raw_events
              << " duplicates=" << grid.duplicate_events
              << " checksum=" << grid_checksum(provenance_one)
              << " provenance_bytes="
              << provenance_one.provenance.retained_bytes
              << " perturbed_rows=" << perturbed_repair_rows
              << " synthetic_odd_status="
              << synthetic_odd_on.provenance.status
              << " retained_bytes=" << grid.retained_bytes << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_ray_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
