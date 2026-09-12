#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_intersections.h"
#include "tabi/tabi_ses_patch_mapping.h"
#include "tabi/tabi_ses_patch_ray.h"
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
#include <memory>
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

PatchMappingGridSpec sphere_grid(const Atom &atom, double probe_radius,
                                 double spacing) {
  const double half_extent = atom.radius + probe_radius + 2.0 * spacing;
  const Vec3 lower{atom.x - half_extent, atom.y - half_extent,
                   atom.z - half_extent};
  const auto count = [&](double extent) {
    return static_cast<std::size_t>(std::ceil((2.0 * extent) / spacing)) + 1U;
  };
  const std::size_t n = count(half_extent);
  return {lower, spacing, n, n, n};
}

std::uint64_t grid_checksum(const PatchRayGrid &grid) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(grid.raw_events);
  mix(grid.duplicate_events);
  mix(grid.uncertified_intersection_queries);
  mix(grid.uncertified_rows_repaired);
  mix(grid.unresolved_uncertified_rows);
  mix(grid.degenerate_rows_repaired);
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

void require_same_grid(const PatchRayGrid &left, const PatchRayGrid &right) {
  require(left.status == right.status && left.raw_events == right.raw_events &&
              left.duplicate_events == right.duplicate_events &&
              left.uncertified_intersection_queries ==
                  right.uncertified_intersection_queries &&
              left.uncertified_rows_repaired ==
                  right.uncertified_rows_repaired &&
              left.unresolved_uncertified_rows ==
                  right.unresolved_uncertified_rows &&
              left.degenerate_rows_repaired == right.degenerate_rows_repaired,
          "patch mapping changed scalar ray contract");
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    require(left.axes[axis].axis == right.axes[axis].axis &&
                left.axes[axis].u_count == right.axes[axis].u_count &&
                left.axes[axis].v_count == right.axes[axis].v_count &&
                left.axes[axis].offsets == right.axes[axis].offsets &&
                left.axes[axis].coordinates == right.axes[axis].coordinates &&
                left.axes[axis].odd_rows == right.axes[axis].odd_rows,
            "patch mapping changed ray CSR contract");
  }
  require(grid_checksum(left) == grid_checksum(right),
          "patch mapping changed ray checksum");
}

std::uint64_t witness_checksum(const PreparedPatchCrossingWitnesses &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  for (const std::size_t offset : prepared.axis_offsets)
    mix(offset);
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  mix(static_cast<std::uint64_t>(prepared.complex_fingerprint.valid));
  for (const std::uint64_t bits : prepared.lattice_identity.scalar_bits)
    mix(bits);
  for (const std::size_t count : prepared.lattice_identity.counts)
    mix(count);
  mix(static_cast<std::uint64_t>(prepared.lattice_identity.valid));
  mix(prepared.refs_audited);
  mix(prepared.root_events_audited);
  mix(prepared.identity_certified_crossings);
  mix(prepared.root_reproduced_crossings);
  mix(prepared.bounded_patch_location_certified_crossings);
  mix(prepared.transversality_certified_crossings);
  mix(prepared.multi_ref_crossings);
  mix(prepared.multi_ref_relation_uncertified_crossings);
  mix(prepared.primitive_simplicity_uncertified_crossings);
  mix(prepared.repeated_root_ambiguities);
  mix(prepared.intersection_uncertified);
  mix(prepared.roots_not_reproduced);
  mix(prepared.root_multiplicity_mismatches);
  mix(prepared.projector_uncertified);
  mix(prepared.patch_residual_failures);
  mix(prepared.nonfinite_bounded_residuals);
  for (const CertifiedPatchCrossingWitness &witness : prepared.witnesses) {
    mix(witness.address.packed_route());
    mix(witness.address.crossing());
    mix(std::bit_cast<std::uint64_t>(witness.point.x));
    mix(std::bit_cast<std::uint64_t>(witness.point.y));
    mix(std::bit_cast<std::uint64_t>(witness.point.z));
    mix(std::bit_cast<std::uint64_t>(witness.coordinate));
    mix(std::bit_cast<std::uint64_t>(witness.root_match_tolerance));
    mix(std::bit_cast<std::uint64_t>(witness.bounded_patch_tolerance));
    mix(std::bit_cast<std::uint64_t>(witness.maximum_root_residual));
    mix(std::bit_cast<std::uint64_t>(witness.maximum_patch_distance));
    mix(witness.ref_count);
    mix(witness.total_root_multiplicity);
    mix(witness.maximum_ref_multiplicity);
    mix(static_cast<std::uint64_t>(witness.identity_status));
    mix(static_cast<std::uint64_t>(witness.root_status));
    mix(static_cast<std::uint64_t>(witness.location_status));
    mix(static_cast<std::uint64_t>(witness.transversality_status));
    mix(static_cast<std::uint64_t>(witness.orientation_sign_certified));
    mix(static_cast<std::uint64_t>(witness.primitive_simplicity_certified));
  }
  return hash;
}

void check_row_decoding() {
  const PatchMappingGridSpec spec{{10.0, 20.0, 30.0}, 0.5, 5U, 7U, 9U};
  std::array<std::size_t, 3> rows{};
  std::size_t vertices = 0U;
  require(validate_patch_mapping_grid_spec(spec, &rows, &vertices) &&
              rows == std::array<std::size_t, 3>{63U, 45U, 35U} &&
              vertices == 315U,
          "mapping grid validation failed");

  DecodedPatchRayRow decoded;
  require(decode_patch_ray_row(spec, RayAxis::x, 3U + 7U * 4U, decoded) &&
              decoded.u_index == 3U && decoded.v_index == 4U &&
              decoded.u == 21.5 && decoded.v == 32.0,
          "x-row decoding failed");
  PatchCrossingAddress address;
  Vec3 point;
  require(
      PatchCrossingAddress::try_make(RayAxis::x, decoded.row, 5U, address) &&
          reconstruct_patch_crossing_point(spec, address, 12.0, point) &&
          point.x == 12.0 && point.y == 21.5 && point.z == 32.0,
      "x crossing reconstruction failed");

  require(decode_patch_ray_row(spec, RayAxis::y, 2U + 5U * 4U, decoded) &&
              decoded.u == 11.0 && decoded.v == 32.0 &&
              PatchCrossingAddress::try_make(RayAxis::y, decoded.row, 6U,
                                             address) &&
              reconstruct_patch_crossing_point(spec, address, 22.25, point) &&
              point.x == 11.0 && point.y == 22.25 && point.z == 32.0,
          "y-row decoding failed");
  require(decode_patch_ray_row(spec, RayAxis::z, 2U + 5U * 3U, decoded) &&
              decoded.u == 11.0 && decoded.v == 21.5 &&
              PatchCrossingAddress::try_make(RayAxis::z, decoded.row, 7U,
                                             address) &&
              reconstruct_patch_crossing_point(spec, address, 31.75, point) &&
              point.x == 11.0 && point.y == 21.5 && point.z == 31.75,
          "z-row decoding failed");
}

void check_overflow_guards() {
  PatchCrossingAddress address;
  require(!PatchCrossingAddress::try_make(static_cast<RayAxis>(3U), 0U, 0U,
                                          address) &&
              !address.valid(),
          "invalid axis produced a crossing address");
  if constexpr (sizeof(std::size_t) > sizeof(std::uint32_t)) {
    require(!PatchCrossingAddress::try_make(
                RayAxis::x, 0U,
                static_cast<std::size_t>(
                    std::numeric_limits<std::uint32_t>::max()) +
                    1U,
                address),
            "crossing index overflow was accepted");
  }
  if (std::numeric_limits<std::size_t>::max() >
      PatchCrossingAddress::row_mask) {
    require(!PatchCrossingAddress::try_make(
                RayAxis::x,
                static_cast<std::size_t>(PatchCrossingAddress::row_mask + 1U),
                0U, address),
            "row index overflow was accepted");
  }
  const PatchMappingGridSpec overflow{
      {0.0, 0.0, 0.0}, 1.0, std::numeric_limits<std::size_t>::max(), 2U, 2U};
  require(!validate_patch_mapping_grid_spec(overflow),
          "overflow grid was accepted");
}

void check_tolerance_boundaries() {
  DecodedPatchRayRow row;
  row.u = 2.0;
  row.v = -3.0;
  const double merge_dominated =
      patch_crossing_root_match_tolerance(1.0e-6, row, 4.0, 0.5);
  require(merge_dominated == 1.0e-6,
          "merge tolerance did not dominate the reproduction band");
  const double just_below = std::nextafter(merge_dominated, 0.0);
  const double just_above =
      std::nextafter(merge_dominated, std::numeric_limits<double>::infinity());
  require(
      patch_crossing_residual_within_tolerance(just_below, merge_dominated) &&
          patch_crossing_residual_within_tolerance(merge_dominated,
                                                   merge_dominated) &&
          !patch_crossing_residual_within_tolerance(just_above,
                                                    merge_dominated),
      "root reproduction tolerance boundary is not closed/exact");

  row.u = 1000000000.0;
  row.v = -2000000000.0;
  const double forward_error_dominated =
      patch_crossing_root_match_tolerance(0.0, row, 3000000000.0, 0.125);
  const double expected =
      2048.0 * std::numeric_limits<double>::epsilon() * 3000000000.0;
  require(forward_error_dominated == expected &&
              forward_error_dominated > merge_dominated,
          "coordinate-scale forward-error band was not selected");
}

struct SphereFixture {
  std::unique_ptr<PatchComplex> complex;
  LocalPatchProjectionIndex index;
  PatchMappingGridSpec spec;
  PatchRayGrid off;
  PatchRayGrid on;
  PreparedPatchCrossingWitnesses prepared;
};

LocalPatchProjectionIndex
build_audit_projection_index(const PatchComplex &complex) {
  LocalPatchProjectionIndexBuildOptions options;
  options.retain_complex_fingerprint = true;
  return build_local_patch_projection_index(complex, options);
}

SphereFixture build_sphere_fixture(const Atom &atom, double probe_radius,
                                   double spacing, std::size_t threads) {
  SphereFixture fixture;
  fixture.complex = std::make_unique<PatchComplex>(build_patch_complex(
      std::vector<Atom>{atom}, probe_radius, nullptr, false));
  require(fixture.complex->status == "pass" &&
              fixture.complex->contacts.size() == 1U,
          "single-sphere patch fixture failed");
  fixture.index = build_audit_projection_index(*fixture.complex);
  require(fixture.index.bounded_domains_certified &&
              fixture.index.bvh_invariants_certified &&
              fixture.index.source_complex_fingerprint.valid,
          "single-sphere projection index failed");
  fixture.spec = sphere_grid(atom, probe_radius, spacing);
  TitanExecutor executor(threads);
  fixture.off = build_patch_ray_grid(
      *fixture.complex, fixture.spec.lower, fixture.spec.grid_spacing,
      fixture.spec.nx, fixture.spec.ny, fixture.spec.nz, &executor);
  PatchRayBuildOptions options;
  options.retain_provenance = true;
  fixture.on = build_patch_ray_grid(
      *fixture.complex, fixture.spec.lower, fixture.spec.grid_spacing,
      fixture.spec.nx, fixture.spec.ny, fixture.spec.nz, &executor, options);
  require_same_grid(fixture.off, fixture.on);
  require(!fixture.off.provenance.lattice_identity.valid &&
              fixture.on.provenance.lattice_identity.valid,
          "provenance-off/on lattice binding policy changed");
  fixture.prepared = build_certified_patch_crossing_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec);
  return fixture;
}

void check_simple_sphere_and_determinism() {
  const Atom atom{0.25, -0.5, 0.75, 0.0, 1.2};
  SphereFixture serial = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  SphereFixture parallel = build_sphere_fixture(atom, 1.4, 0.5, 8U);
  require(serial.prepared.status == "pass" &&
              serial.prepared.identity_certificate_complete &&
              serial.prepared.root_reproduction_complete &&
              serial.prepared.location_certificate_complete &&
              serial.prepared.transversality_certificate_complete &&
              !serial.prepared.witnesses.empty() &&
              serial.prepared.bounded_patch_location_certified_crossings ==
                  serial.prepared.witnesses.size(),
          "simple sphere crossing mapping was not certified");
  for (const CertifiedPatchCrossingWitness &witness :
       serial.prepared.witnesses) {
    require(witness.identity_certified() && witness.root_reproduced() &&
                witness.location_certified() &&
                witness.transversality_certified() &&
                witness.primitive_simplicity_certified &&
                !witness.orientation_sign_certified &&
                witness.maximum_root_residual <= witness.root_match_tolerance &&
                witness.maximum_patch_distance <=
                    witness.bounded_patch_tolerance &&
                witness.total_root_multiplicity == 1U &&
                witness.maximum_ref_multiplicity == 1U,
            "simple sphere witness contract failed");
  }
  require_same_grid(serial.on, parallel.on);
  require(witness_checksum(serial.prepared) ==
                  witness_checksum(parallel.prepared) &&
              serial.prepared.refs_audited == parallel.prepared.refs_audited &&
              serial.prepared.transversality_certified_crossings ==
                  parallel.prepared.transversality_certified_crossings,
          "nt=1/8 witness mapping was not deterministic");
}

void check_malformed_and_repeated_fail_closed() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture fixture = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  std::size_t axis_index = 0U;
  while (axis_index < 3U && fixture.on.axes[axis_index].coordinates.empty())
    ++axis_index;
  require(axis_index < 3U, "sphere fixture has no crossing to corrupt");

  PatchRayGrid malformed = fixture.on;
  malformed.axes[axis_index].offsets.front() = 1U;
  const PreparedPatchCrossingWitnesses malformed_result =
      build_certified_patch_crossing_witnesses(*fixture.complex, fixture.index,
                                               malformed, fixture.spec);
  require(malformed_result.status == "invalid_grid_contract" &&
              malformed_result.witnesses.empty() &&
              !malformed_result.location_certificate_complete,
          "malformed CSR did not fail closed");

  PatchRayGrid invalid_ref = fixture.on;
  invalid_ref.provenance.axes[axis_index].canonical_refs.front() = {};
  const PreparedPatchCrossingWitnesses invalid_ref_result =
      build_certified_patch_crossing_witnesses(*fixture.complex, fixture.index,
                                               invalid_ref, fixture.spec);
  require(invalid_ref_result.status == "invalid_grid_contract" &&
              invalid_ref_result.witnesses.empty(),
          "invalid patch ref did not fail closed");

  PatchRayGrid repeated = fixture.on;
  repeated.provenance.axes[axis_index].single_ref_multiplicities.push_back(
      {0U, 2U});
  ++repeated.provenance.roots_total;
  ++repeated.provenance.repeated_root_groups;
  ++repeated.provenance.same_primitive_near_pair_unclassified;
  const PreparedPatchCrossingWitnesses repeated_result =
      build_certified_patch_crossing_witnesses(*fixture.complex, fixture.index,
                                               repeated, fixture.spec);
  require(repeated_result.status == "root_reproduction_uncertified" &&
              repeated_result.identity_certificate_complete &&
              !repeated_result.root_reproduction_complete &&
              !repeated_result.location_certificate_complete &&
              !repeated_result.transversality_certificate_complete &&
              repeated_result.repeated_root_ambiguities == 1U &&
              repeated_result.root_multiplicity_mismatches == 1U &&
              repeated_result.root_events_audited ==
                  repeated.provenance.roots_total &&
              repeated_result.witnesses.front().total_root_multiplicity == 2U &&
              repeated_result.witnesses.front().maximum_ref_multiplicity ==
                  2U &&
              repeated_result.witnesses.front().identity_status ==
                  CrossingIdentityCertificateStatus::valid_single_ref &&
              repeated_result.witnesses.front().root_status ==
                  CrossingRootReproductionStatus::root_multiplicity_mismatch &&
              repeated_result.witnesses.front().location_status ==
                  CrossingBoundedLocationStatus::root_unreproduced &&
              repeated_result.witnesses.front().transversality_status ==
                  CrossingTransversalityStatus::location_uncertified,
          "fabricated repeated-root provenance did not fail closed");

  PatchRayGrid repeated_three = repeated;
  repeated_three.provenance.axes[axis_index]
      .single_ref_multiplicities.front()
      .multiplicity = 3U;
  ++repeated_three.provenance.roots_total;
  const PreparedPatchCrossingWitnesses repeated_three_result =
      build_certified_patch_crossing_witnesses(*fixture.complex, fixture.index,
                                               repeated_three, fixture.spec);
  require(
      repeated_three_result.status == "root_reproduction_uncertified" &&
          repeated_three_result.witnesses.front().total_root_multiplicity ==
              3U &&
          repeated_three_result.witnesses.front().maximum_ref_multiplicity ==
              3U &&
          witness_checksum(repeated_result) !=
              witness_checksum(repeated_three_result),
      "exact root multiplicity was not retained in the witness checksum");

  const CrossingMultiplicityClassification reproduced_simple =
      classify_patch_crossing_multiplicity(1U, 1U);
  const CrossingMultiplicityClassification reproduced_repeated =
      classify_patch_crossing_multiplicity(2U, 2U);
  const CrossingMultiplicityClassification multiplicity_mismatch =
      classify_patch_crossing_multiplicity(2U, 1U);
  const CrossingTransversalityStatus unresolved_multi_ref =
      classify_patch_crossing_group_transversality(2U, 1U, true, true);
  const CrossingTransversalityStatus unresolved_primitive_simplicity =
      classify_patch_crossing_group_transversality(1U, 1U, true, false);
  require(
      reproduced_simple.root_status ==
              CrossingRootReproductionStatus::reproduced &&
          reproduced_simple.transversality_if_location_certified ==
              CrossingTransversalityStatus::
                  primitive_filtered_transversal_unoriented &&
          reproduced_repeated.root_status ==
              CrossingRootReproductionStatus::reproduced &&
          reproduced_repeated.transversality_if_location_certified ==
              CrossingTransversalityStatus::repeated_root_ambiguous &&
          multiplicity_mismatch.root_status ==
              CrossingRootReproductionStatus::root_multiplicity_mismatch &&
          multiplicity_mismatch.transversality_if_location_certified ==
              CrossingTransversalityStatus::location_uncertified &&
          unresolved_multi_ref ==
              CrossingTransversalityStatus::multi_ref_relation_uncertified &&
          unresolved_primitive_simplicity ==
              CrossingTransversalityStatus::primitive_simplicity_uncertified,
      "root reproduction and transversality classifications were conflated");

  PatchRayGrid stale_summary = repeated;
  ++stale_summary.provenance.roots_total;
  const PreparedPatchCrossingWitnesses stale_summary_result =
      build_certified_patch_crossing_witnesses(*fixture.complex, fixture.index,
                                               stale_summary, fixture.spec);
  require(stale_summary_result.status == "invalid_provenance_summary" &&
              stale_summary_result.witnesses.empty() &&
              stale_summary_result.root_events_audited == 0U,
          "mutated provenance summary was silently accepted");
}

void check_fingerprint_and_bounded_projection_fail_closed() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  PatchComplex nongeometric_mutation =
      build_patch_complex(std::vector<Atom>{atom}, 1.4, nullptr, false);
  const PatchComplexFingerprint geometry_only_fingerprint =
      fingerprint_patch_complex(nongeometric_mutation);
  nongeometric_mutation.atoms[0].q += 3.0;
  ++nongeometric_mutation.candidate_pairs;
  ++nongeometric_mutation.candidate_triples;
  ++nongeometric_mutation.buried_pair_circles;
  ++nongeometric_mutation.rejected_probe_centres;
  require(fingerprint_patch_complex(nongeometric_mutation) ==
              geometry_only_fingerprint,
          "nongeometric diagnostics changed the intersection fingerprint");

  SphereFixture stale = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  const PatchComplexFingerprint original =
      stale.on.provenance.complex_fingerprint;
  stale.complex->contacts[0].center.x += 1.0e-6;
  require(fingerprint_patch_complex(*stale.complex) != original,
          "same-address complex mutation did not change fingerprint");
  const PreparedPatchCrossingWitnesses stale_result =
      build_certified_patch_crossing_witnesses(*stale.complex, stale.index,
                                               stale.on, stale.spec);
  require(stale_result.status == "stale_grid_complex" &&
              stale_result.witnesses.empty() &&
              stale_result.refs_audited == 0U &&
              !stale_result.identity_certificate_complete &&
              !stale_result.location_certificate_complete,
          "stale same-address grid/complex binding did not fail closed");

  SphereFixture stale_index = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  const LocalPatchProjectionIndex pre_mutation_index = stale_index.index;
  stale_index.complex->atoms[0].x += 1.0e-4;
  stale_index.complex->contacts[0].center.x += 1.0e-4;
  TitanExecutor executor(1U);
  PatchRayBuildOptions provenance_options;
  provenance_options.retain_provenance = true;
  const PatchRayGrid fresh_grid = build_patch_ray_grid(
      *stale_index.complex, stale_index.spec.lower,
      stale_index.spec.grid_spacing, stale_index.spec.nx, stale_index.spec.ny,
      stale_index.spec.nz, &executor, provenance_options);
  require(fresh_grid.provenance.certified &&
              fresh_grid.provenance.complex_fingerprint ==
                  fingerprint_patch_complex(*stale_index.complex),
          "fresh post-mutation provenance grid was not certified");
  const PreparedPatchCrossingWitnesses stale_index_result =
      build_certified_patch_crossing_witnesses(*stale_index.complex,
                                               pre_mutation_index, fresh_grid,
                                               stale_index.spec);
  require(stale_index_result.status == "invalid_projection_index" &&
              stale_index_result.witnesses.empty(),
          "fresh grid accepted a stale same-address projection index");

  SphereFixture unbound_fixture = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  const LocalPatchProjectionIndex unbound_index =
      build_local_patch_projection_index(*unbound_fixture.complex);
  require(!unbound_index.source_complex_fingerprint.valid &&
              unbound_index.bounded_domains_certified &&
              unbound_index.bvh_invariants_certified,
          "default projection index unexpectedly retained an audit binding");
  const PreparedPatchCrossingWitnesses unbound_result =
      build_certified_patch_crossing_witnesses(
          *unbound_fixture.complex, unbound_index, unbound_fixture.on,
          unbound_fixture.spec);
  require(unbound_result.status == "invalid_projection_index" &&
              unbound_result.witnesses.empty(),
          "mapping accepted a production-default unbound projection index");

  SphereFixture stale_lattice = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  PatchMappingGridSpec stale_lower = stale_lattice.spec;
  stale_lower.lower.x = std::nextafter(stale_lower.lower.x,
                                       std::numeric_limits<double>::infinity());
  const PreparedPatchCrossingWitnesses stale_lower_result =
      build_certified_patch_crossing_witnesses(*stale_lattice.complex,
                                               stale_lattice.index,
                                               stale_lattice.on, stale_lower);
  PatchMappingGridSpec stale_spacing = stale_lattice.spec;
  stale_spacing.grid_spacing = std::nextafter(
      stale_spacing.grid_spacing, std::numeric_limits<double>::infinity());
  const PreparedPatchCrossingWitnesses stale_spacing_result =
      build_certified_patch_crossing_witnesses(*stale_lattice.complex,
                                               stale_lattice.index,
                                               stale_lattice.on, stale_spacing);
  require(stale_lower_result.status == "stale_grid_lattice" &&
              stale_spacing_result.status == "stale_grid_lattice" &&
              stale_lower_result.witnesses.empty() &&
              stale_spacing_result.witnesses.empty(),
          "stale lower/spacing lattice identity was accepted");

  PatchRayGrid out_of_bounds = stale_lattice.on;
  std::size_t out_axis = 0U;
  while (out_axis < out_of_bounds.axes.size() &&
         out_of_bounds.axes[out_axis].coordinates.empty())
    ++out_axis;
  require(out_axis < out_of_bounds.axes.size(),
          "out-of-bounds fixture has no crossing");
  const double axis_lower = out_axis == 0U   ? stale_lattice.spec.lower.x
                            : out_axis == 1U ? stale_lattice.spec.lower.y
                                             : stale_lattice.spec.lower.z;
  out_of_bounds.axes[out_axis].coordinates.front() = axis_lower - 100.0;
  const PreparedPatchCrossingWitnesses out_of_bounds_result =
      build_certified_patch_crossing_witnesses(
          *stale_lattice.complex, stale_lattice.index, out_of_bounds,
          stale_lattice.spec);
  require(out_of_bounds_result.status == "crossing_out_of_axis_bounds" &&
              out_of_bounds_result.witnesses.empty(),
          "crossing outside its lattice axis was accepted");

  SphereFixture invalid_projector = build_sphere_fixture(atom, 1.4, 0.5, 1U);
  require(!invalid_projector.index.contact_domain_empty.empty(),
          "contact empty-domain fixture is missing");
  invalid_projector.index.contact_domain_empty[0] = 1U;
  const PreparedPatchCrossingWitnesses projector_result =
      build_certified_patch_crossing_witnesses(
          *invalid_projector.complex, invalid_projector.index,
          invalid_projector.on, invalid_projector.spec);
  require(projector_result.status == "location_uncertified" &&
              projector_result.identity_certificate_complete &&
              projector_result.root_reproduction_complete &&
              !projector_result.location_certificate_complete &&
              !projector_result.transversality_certificate_complete &&
              projector_result.projector_uncertified > 0U,
          "bounded-patch projector failure was silently certified");
}

void check_primitive_filtered_transversality_contract() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  const PatchComplex complex =
      build_patch_complex(std::vector<Atom>{atom}, 1.4, nullptr, false);
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  const PatchAxisRoots tangent =
      intersect_patch_axis(complex, index.reentrant_probes, PatchKind::contact,
                           0U, RayAxis::x, atom.y + atom.radius, atom.z);
  const PatchAxisRoots retained_tangent =
      intersect_patch_axis(complex, index.reentrant_probes, PatchKind::contact,
                           0U, RayAxis::x, atom.y + atom.radius, atom.z, true);
  const PatchAxisRoots one_ulp_inside = intersect_patch_axis(
      complex, index.reentrant_probes, PatchKind::contact, 0U, RayAxis::x,
      std::nextafter(atom.y + atom.radius, atom.y), atom.z);
  const PatchAxisRoots clearly_transversal =
      intersect_patch_axis(complex, index.reentrant_probes, PatchKind::contact,
                           0U, RayAxis::x, atom.y, atom.z);
  require(tangent.certified && tangent.count == 0U &&
              retained_tangent.certified && retained_tangent.count == 1U &&
              retained_tangent.primitive_simplicity_certified[0] == 0U &&
              one_ulp_inside.certified && one_ulp_inside.count == 1U,
          "sphere tangent/collapsed-pair coordinate contract changed");
  require(one_ulp_inside.primitive_simplicity_certified[0] == 0U &&
              clearly_transversal.certified &&
              clearly_transversal.count == 2U &&
              clearly_transversal.primitive_simplicity_certified[0] != 0U &&
              clearly_transversal.primitive_simplicity_certified[1] != 0U,
          "sphere primitive-simplicity metadata is not fail closed");

  const std::vector<Atom> paired_atoms = {{-3.0, 0.0, 0.0, 0.0, 1.0},
                                          {3.0, 0.0, 0.0, 0.0, 1.0}};
  const PatchComplex paired_complex =
      build_patch_complex(paired_atoms, 1.4, nullptr, false);
  const LocalPatchProjectionIndex paired_index =
      build_audit_projection_index(paired_complex);
  const double paired_u = std::nextafter(1.0, 0.0);
  const PatchAxisRoots paired_left =
      intersect_patch_axis(paired_complex, paired_index.reentrant_probes,
                           PatchKind::contact, 0U, RayAxis::x, paired_u, 0.0);
  const PatchAxisRoots paired_right =
      intersect_patch_axis(paired_complex, paired_index.reentrant_probes,
                           PatchKind::contact, 1U, RayAxis::x, paired_u, 0.0);
  require(paired_complex.status == "pass" &&
              paired_complex.contacts.size() == 2U &&
              paired_complex.tori.empty() && paired_left.count == 1U &&
              paired_right.count == 1U &&
              paired_left.primitive_simplicity_certified[0] == 0U &&
              paired_right.primitive_simplicity_certified[0] == 0U,
          "paired collapsed-root fixture failed");
  const PatchMappingGridSpec spec{{-5.0, paired_u, 0.0}, 10.0, 2U, 1U, 1U};
  PatchRayGrid collapsed_grid;
  collapsed_grid.status = "pass";
  collapsed_grid.raw_events = 2U;
  collapsed_grid.axes[0] = {
      RayAxis::x,
      1U,
      1U,
      {0U, 2U},
      {paired_left.coordinates[0], paired_right.coordinates[0]},
      0U};
  collapsed_grid.axes[1] = {RayAxis::y, 2U, 1U, {0U, 0U, 0U}, {}, 0U};
  collapsed_grid.axes[2] = {RayAxis::z, 2U, 1U, {0U, 0U, 0U}, {}, 0U};
  collapsed_grid.provenance.enabled = true;
  collapsed_grid.provenance.certified = true;
  collapsed_grid.provenance.status = "pass";
  collapsed_grid.provenance.complex_fingerprint =
      fingerprint_patch_complex(paired_complex);
  collapsed_grid.provenance.lattice_identity = make_patch_ray_lattice_identity(
      spec.lower, spec.grid_spacing, spec.nx, spec.ny, spec.nz);
  collapsed_grid.provenance.coordinate_merge_tolerance = 0.0;
  collapsed_grid.provenance.unique_crossings = 2U;
  collapsed_grid.provenance.refs_total = 2U;
  collapsed_grid.provenance.max_refs_per_crossing = 1U;
  collapsed_grid.provenance.roots_total = 2U;
  PatchRefId left_ref;
  PatchRefId right_ref;
  require(PatchRefId::try_make(PatchKind::contact, 0U, left_ref) &&
              PatchRefId::try_make(PatchKind::contact, 1U, right_ref),
          "paired contact ref packing failed");
  const std::array<PatchRefId, 1> left_refs{left_ref};
  const std::array<PatchRefId, 1> right_refs{right_ref};
  require(collapsed_grid.provenance.axes[0].append(left_refs) &&
              collapsed_grid.provenance.axes[0].append(right_refs),
          "collapsed-pair provenance fixture failed");
  const PreparedPatchCrossingWitnesses collapsed_result =
      build_certified_patch_crossing_witnesses(paired_complex, paired_index,
                                               collapsed_grid, spec);
  require(
      collapsed_result.status == "primitive_simplicity_uncertified" &&
          collapsed_result.identity_certificate_complete &&
          collapsed_result.root_reproduction_complete &&
          collapsed_result.location_certificate_complete &&
          !collapsed_result.transversality_certificate_complete &&
          collapsed_result.primitive_simplicity_uncertified_crossings == 2U &&
          collapsed_result.witnesses.size() == 2U &&
          !collapsed_result.witnesses[0].primitive_simplicity_certified &&
          !collapsed_result.witnesses[1].primitive_simplicity_certified &&
          collapsed_result.witnesses[0].transversality_status ==
              CrossingTransversalityStatus::primitive_simplicity_uncertified,
      "collapsed near-tangent pair was certified as a parity crossing");
}

void check_translation_scale_tolerance() {
  const Atom translated{1000000.25, -2000000.5, 3000000.75, 0.0, 1.1};
  SphereFixture fixture = build_sphere_fixture(translated, 1.4, 0.125, 8U);
  require(fixture.prepared.status == "pass" &&
              fixture.prepared.location_certificate_complete &&
              fixture.prepared.transversality_certificate_complete &&
              fixture.prepared.bounded_patch_location_certified_crossings ==
                  fixture.prepared.witnesses.size(),
          "translated scale-aware root certificate failed");
  double maximum_residual = 0.0;
  double maximum_tolerance = 0.0;
  for (const CertifiedPatchCrossingWitness &witness :
       fixture.prepared.witnesses) {
    maximum_residual =
        std::max(maximum_residual, witness.maximum_root_residual);
    maximum_tolerance =
        std::max(maximum_tolerance, witness.root_match_tolerance);
  }
  require(std::isfinite(maximum_residual) &&
              maximum_residual <= maximum_tolerance,
          "translated root residual exceeded its certificate tolerance");
}

std::size_t grid_coordinate_index(double lower, double spacing,
                                  double coordinate) {
  const double scaled = (coordinate - lower) / spacing;
  const auto rounded = static_cast<std::size_t>(std::llround(scaled));
  require(
      std::abs((lower + spacing * static_cast<double>(rounded)) - coordinate) <=
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(coordinate)),
      "requested coordinate is not on the fixture lattice");
  return rounded;
}

PatchGridVertexMagnitudeRequest make_vertex_request(
    const PatchMappingGridSpec &spec, double x, double y, double z,
    PatchGridTopologyState state = PatchGridTopologyState::inside,
    PatchGridTopologySource source =
        PatchGridTopologySource::inside_atom_union) {
  PatchGridVertexMagnitudeRequest request;
  request.ix = grid_coordinate_index(spec.lower.x, spec.grid_spacing, x);
  request.iy = grid_coordinate_index(spec.lower.y, spec.grid_spacing, y);
  request.iz = grid_coordinate_index(spec.lower.z, spec.grid_spacing, z);
  request.topology_state = state;
  request.topology_source = source;
  return request;
}

double reference_production_raw_magnitude(
    const PatchRayGrid &grid, const PatchMappingGridSpec &spec,
    const PatchGridVertexMagnitudeRequest &request) {
  const std::array<std::size_t, 3> rows = {request.iy + spec.ny * request.iz,
                                           request.ix + spec.nx * request.iz,
                                           request.ix + spec.nx * request.iy};
  const std::array<double, 3> coordinates = {
      spec.lower.x + spec.grid_spacing * static_cast<double>(request.ix),
      spec.lower.y + spec.grid_spacing * static_cast<double>(request.iy),
      spec.lower.z + spec.grid_spacing * static_cast<double>(request.iz)};
  double distance = request.sparse_magnitude_present
                        ? static_cast<double>(request.sparse_magnitude)
                        : spec.grid_spacing;
  const auto nearest = [](std::span<const double> crossings, double coordinate,
                          double fallback) {
    const auto next =
        std::lower_bound(crossings.begin(), crossings.end(), coordinate);
    double value = fallback;
    if (next != crossings.end())
      value = std::min(value, std::abs(*next - coordinate));
    if (next != crossings.begin())
      value = std::min(value, std::abs(*std::prev(next) - coordinate));
    return value;
  };
  distance = nearest(grid.axes[0].row(rows[0]), coordinates[0], distance);
  if (!request.sparse_magnitude_present) {
    distance = nearest(grid.axes[1].row(rows[1]), coordinates[1], distance);
    distance = nearest(grid.axes[2].row(rows[2]), coordinates[2], distance);
  }
  return distance;
}

std::span<const PatchGridMagnitudeSourceCandidate>
vertex_candidates(const PreparedPatchGridVertexMagnitudeWitnesses &prepared,
                  std::size_t vertex) {
  require(vertex + 1U < prepared.candidate_offsets.size(),
          "candidate CSR vertex is out of range");
  return {prepared.candidates.data() + prepared.candidate_offsets[vertex],
          prepared.candidate_offsets[vertex + 1U] -
              prepared.candidate_offsets[vertex]};
}

std::size_t count_candidate_kind(
    std::span<const PatchGridMagnitudeSourceCandidate> candidates,
    PatchGridMagnitudeSourceKind kind) {
  return static_cast<std::size_t>(
      std::count_if(candidates.begin(), candidates.end(),
                    [&](const PatchGridMagnitudeSourceCandidate &candidate) {
                      return candidate.kind == kind;
                    }));
}

std::uint64_t magnitude_witness_checksum(
    const PreparedPatchGridVertexMagnitudeWitnesses &prepared) {
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
  for (const PatchGridVertexMagnitudeWitness &vertex : prepared.vertices) {
    mix(vertex.linear_index);
    mix(vertex.ix);
    mix(vertex.iy);
    mix(vertex.iz);
    mix(std::bit_cast<std::uint64_t>(vertex.raw_magnitude));
    mix(std::bit_cast<std::uint64_t>(vertex.clamped_magnitude));
    mix(std::bit_cast<std::uint64_t>(vertex.double_arithmetic_tie_tolerance));
    mix(std::bit_cast<std::uint64_t>(vertex.sparse_float_tie_tolerance));
    mix(std::bit_cast<std::uint64_t>(vertex.source_tie_tolerance));
    mix(static_cast<std::uint64_t>(vertex.topology_state));
    mix(static_cast<std::uint64_t>(vertex.topology_source));
    mix(std::bit_cast<std::uint32_t>(vertex.sparse_magnitude));
    mix(static_cast<std::uint64_t>(vertex.topology_evidence_status));
    mix(static_cast<std::uint64_t>(vertex.sparse_magnitude_present));
    mix(static_cast<std::uint64_t>(vertex.sparse_suppressed_yz));
    mix(static_cast<std::uint64_t>(vertex.cavity_fill));
    mix(static_cast<std::uint64_t>(vertex.floor_clamp_applied));
    mix(static_cast<std::uint64_t>(vertex.topology_analytic_evidence));
    mix(static_cast<std::uint64_t>(vertex.topology_exact_vertex_certified));
    mix(static_cast<std::uint64_t>(vertex.raw_source_certificate_complete));
    mix(static_cast<std::uint64_t>(
        vertex.final_scalar_source_certificate_complete));
    mix(static_cast<std::uint64_t>(vertex.surface_location_certified));
  }
  for (const std::size_t offset : prepared.candidate_offsets)
    mix(offset);
  for (const PatchGridMagnitudeSourceCandidate &candidate :
       prepared.candidates) {
    mix(static_cast<std::uint64_t>(candidate.kind));
    mix(static_cast<std::uint64_t>(candidate.axis));
    mix(candidate.crossing.packed_route());
    mix(candidate.crossing.crossing());
    mix(std::bit_cast<std::uint64_t>(candidate.distance));
    mix(candidate.refs_begin);
    mix(candidate.refs_end);
    mix(static_cast<std::uint64_t>(candidate.evidence_status));
    mix(static_cast<std::uint64_t>(candidate.raw_exact_owner));
    mix(static_cast<std::uint64_t>(candidate.raw_tolerance_tied));
    mix(static_cast<std::uint64_t>(candidate.final_scalar_exact_owner));
  }
  for (const PatchRefId ref : prepared.candidate_refs)
    mix(ref.packed());
  mix(prepared.analytic_candidates);
  mix(prepared.certified_analytic_candidates);
  mix(prepared.uncertified_analytic_candidates);
  mix(prepared.sparse_candidates);
  mix(prepared.fallback_candidates);
  mix(prepared.clamp_floor_candidates);
  mix(prepared.raw_exact_owner_candidates);
  mix(prepared.tolerance_tied_candidates);
  mix(prepared.final_scalar_exact_owner_candidates);
  mix(prepared.raw_source_certified_vertices);
  mix(prepared.final_scalar_source_certified_vertices);
  mix(prepared.floor_clamped_vertices);
  mix(prepared.sparse_suppressed_yz_vertices);
  mix(prepared.topology_analytic_evidence_vertices);
  mix(prepared.topology_exact_vertex_certified_vertices);
  mix(prepared.topology_nonanalytic_vertices);
  mix(static_cast<std::uint64_t>(prepared.raw_source_certificate_complete));
  mix(static_cast<std::uint64_t>(
      prepared.final_scalar_source_certificate_complete));
  mix(static_cast<std::uint64_t>(
      prepared.topology_source_classification_complete));
  mix(static_cast<std::uint64_t>(
      prepared.topology_exact_vertex_certificate_complete));
  mix(static_cast<std::uint64_t>(prepared.surface_location_certified));
  return hash;
}

std::uint64_t
cell_lineage_checksum(const PreparedPatchGridCellLineage &prepared) {
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
  for (const PatchGridCanonicalCellLineage &cell : prepared.cells) {
    mix(cell.address.linear());
    mix(cell.ix);
    mix(cell.iy);
    mix(cell.iz);
    mix(cell.inside_mask);
    for (const std::size_t vertex : cell.corner_vertex_linear_indices)
      mix(vertex);
    mix(cell.topology_sign_evidence_complete);
  }
  for (const std::size_t offset : prepared.cell_tetra_offsets)
    mix(offset);
  for (const PatchGridCanonicalTetraLineage &tetra : prepared.tetrahedra) {
    mix(tetra.cell_position);
    mix(tetra.tetra_index);
    mix(tetra.inside_mask);
    for (const std::uint8_t corner : tetra.cell_corners)
      mix(corner);
    for (const std::size_t vertex : tetra.vertex_linear_indices)
      mix(vertex);
    mix(tetra.topology_sign_evidence_complete);
  }
  for (const std::size_t offset : prepared.tetra_edge_offsets)
    mix(offset);
  for (const PatchGridSignChangingEdgeLineage &edge : prepared.edges) {
    mix(edge.tetra_position);
    mix(edge.edge_ordinal);
    mix(edge.first_cell_corner);
    mix(edge.second_cell_corner);
    mix(edge.first_vertex_linear_index);
    mix(edge.second_vertex_linear_index);
    mix(static_cast<std::uint64_t>(edge.cartesian_axis));
    for (const auto status : edge.endpoint_scalar_status)
      mix(static_cast<std::uint64_t>(status));
    for (const auto status : edge.endpoint_topology_status)
      mix(static_cast<std::uint64_t>(status));
    mix(static_cast<std::uint64_t>(edge.scalar_status));
    mix(static_cast<std::uint64_t>(edge.location_status));
    mix(static_cast<std::uint64_t>(edge.topology_status));
    mix(edge.location_crossing.packed_route());
    mix(edge.location_crossing.crossing());
    mix(std::bit_cast<std::uint64_t>(edge.interpolated_point.x));
    mix(std::bit_cast<std::uint64_t>(edge.interpolated_point.y));
    mix(std::bit_cast<std::uint64_t>(edge.interpolated_point.z));
    mix(std::bit_cast<std::uint64_t>(edge.interpolation_parameter));
    mix(std::bit_cast<std::uint64_t>(edge.interpolation_residual));
    mix(std::bit_cast<std::uint64_t>(edge.interpolation_tolerance));
    mix(edge.production_interpolation_complete);
    mix(edge.analytic_location_residual_available);
    mix(edge.cartesian_axis_edge);
    mix(edge.scalar_lineage_complete);
    mix(edge.patch_candidate_collection_complete);
    mix(edge.patch_relation_classification_complete);
    mix(edge.surface_location_certified);
    mix(edge.topology_sign_evidence_complete);
  }
  for (const std::size_t offset : prepared.edge_endpoint_patch_offsets)
    mix(offset);
  for (const PatchRefId ref : prepared.endpoint_patch_refs)
    mix(ref.packed());
  for (const std::size_t offset : prepared.edge_relation_offsets)
    mix(offset);
  for (const PatchGridEdgePatchRelationEvidence &relation :
       prepared.relations) {
    mix(relation.first.packed());
    mix(relation.second.packed());
    mix(static_cast<std::uint64_t>(relation.relation));
  }
  for (const std::size_t offset : prepared.edge_location_ref_offsets)
    mix(offset);
  for (const PatchRefId ref : prepared.location_patch_refs)
    mix(ref.packed());
  for (const std::size_t count : prepared.relation_counts)
    mix(count);
  mix(prepared.axis_edges);
  mix(prepared.diagonal_edges);
  mix(prepared.scalar_lineage_certified_edges);
  mix(prepared.surface_location_certified_edges);
  mix(prepared.topology_sign_certified_edges);
  mix(static_cast<std::uint64_t>(prepared.relation_oracle_binding));
  mix(prepared.scalar_lineage_complete);
  mix(prepared.patch_candidate_collection_complete);
  mix(prepared.patch_relation_classification_complete);
  mix(prepared.surface_location_certificate_complete);
  mix(prepared.topology_sign_certificate_complete);
  return hash;
}

std::uint64_t triangle_lineage_checksum(
    const PreparedPatchGridTriangleOccurrenceLineage &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  mix(static_cast<std::uint64_t>(prepared.complex_fingerprint.valid));
  for (const std::uint64_t bits : prepared.lattice_identity.scalar_bits)
    mix(bits);
  for (const std::size_t count : prepared.lattice_identity.counts)
    mix(count);
  mix(static_cast<std::uint64_t>(prepared.lattice_identity.valid));
  for (const std::size_t value : prepared.edge_occurrence_to_canonical_vertex)
    mix(value);
  for (const PatchGridCanonicalInterpolationVertexLineage &vertex :
       prepared.canonical_vertices) {
    mix(vertex.key.first_vertex_linear_index());
    mix(vertex.key.second_vertex_linear_index());
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.x));
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.y));
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.z));
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_parameter));
    mix(vertex.exact_location_crossing.packed_route());
    mix(vertex.exact_location_crossing.crossing());
    mix(static_cast<std::uint64_t>(vertex.canonical_interpolation_complete));
    mix(static_cast<std::uint64_t>(
        vertex.canonical_edge_occurrence_geometry_bitwise_consistent));
    mix(static_cast<std::uint64_t>(
        vertex.directed_production_geometry_bitwise_certified));
    mix(static_cast<std::uint64_t>(vertex.patch_candidate_collection_complete));
    mix(static_cast<std::uint64_t>(
        vertex.patch_relation_classification_complete));
    mix(static_cast<std::uint64_t>(vertex.exact_analytic_location_certified));
    mix(static_cast<std::uint64_t>(
        vertex.candidate_refs_present_without_exact_location));
  }
  for (const std::size_t value : prepared.canonical_vertex_occurrence_offsets)
    mix(value);
  for (const std::size_t value : prepared.canonical_vertex_edge_occurrences)
    mix(value);
  for (const std::size_t value :
       prepared.canonical_vertex_candidate_ref_offsets)
    mix(value);
  for (const PatchRefId ref : prepared.canonical_vertex_candidate_refs)
    mix(ref.packed());
  for (const std::size_t value : prepared.canonical_vertex_relation_offsets)
    mix(value);
  for (const PatchGridCanonicalVertexPatchRelationEvidence &relation :
       prepared.canonical_vertex_relations) {
    mix(relation.first.packed());
    mix(relation.second.packed());
    mix(static_cast<std::uint64_t>(relation.relation));
  }
  for (const std::size_t value :
       prepared.canonical_vertex_exact_location_ref_offsets)
    mix(value);
  for (const PatchRefId ref : prepared.canonical_vertex_exact_location_refs)
    mix(ref.packed());
  for (const PatchGridPreOrientationTriangleOccurrenceLineage &triangle :
       prepared.triangles) {
    mix(triangle.tetra_position);
    mix(triangle.triangle_ordinal);
    for (const std::size_t value : triangle.source_edge_occurrence_positions)
      mix(value);
    for (const std::size_t value : triangle.canonical_vertex_positions)
      mix(value);
    mix(std::bit_cast<std::uint64_t>(triangle.canonical_twice_area));
    mix(triangle.single_common_candidate_ref.packed());
    mix(static_cast<std::uint64_t>(triangle.orientation_stage));
    mix(static_cast<std::uint64_t>(triangle.finite_geometry_observed));
    mix(static_cast<std::uint64_t>(triangle.computed_zero_area_observed));
    mix(static_cast<std::uint64_t>(triangle.patch_candidate_summary_complete));
    mix(static_cast<std::uint64_t>(
        triangle.single_common_candidate_ref_available));
    mix(static_cast<std::uint64_t>(
        triangle.vertex_triplet_patch_identity_certified));
    mix(static_cast<std::uint64_t>(
        triangle.planar_face_surface_location_certified));
    mix(static_cast<std::uint64_t>(triangle.final_mesh_orientation_certified));
  }
  for (const std::size_t value : prepared.triangle_candidate_ref_offsets)
    mix(value);
  for (const PatchRefId ref : prepared.triangle_candidate_refs)
    mix(ref.packed());
  for (const std::size_t value : prepared.relation_counts)
    mix(value);
  mix(prepared.interpolation_occurrences);
  mix(prepared.canonical_interpolation_vertices);
  mix(prepared.bitwise_interpolation_conflicts);
  mix(prepared.exact_location_vertices);
  mix(prepared.candidate_refs_without_exact_location_vertices);
  mix(prepared.triangle_occurrences);
  mix(prepared.finite_triangle_occurrences);
  mix(prepared.computed_zero_area_triangle_occurrences);
  mix(prepared.nonfinite_triangle_occurrences);
  mix(prepared.single_common_candidate_triangles);
  mix(prepared.vertex_triplet_patch_identity_certified_triangles);
  mix(static_cast<std::uint64_t>(prepared.interpolation_identity_complete));
  mix(static_cast<std::uint64_t>(
      prepared.patch_candidate_attribution_complete));
  mix(static_cast<std::uint64_t>(
      prepared.triangle_occurrence_lineage_complete));
  mix(static_cast<std::uint64_t>(
      prepared.vertex_triplet_patch_identity_classification_complete));
  mix(static_cast<std::uint64_t>(
      prepared.planar_face_surface_location_certified));
  mix(static_cast<std::uint64_t>(prepared.final_mesh_orientation_certified));
  for (const unsigned char character : prepared.status)
    mix(character);
  return hash;
}

struct HardCodedTriangleOracle {
  std::size_t count{0U};
  std::array<PatchGridPreOrientationTriangleEdges, 2> triangles{};
};

HardCodedTriangleOracle hard_coded_triangle_oracle(std::uint8_t mask) {
  using Edge = PatchGridLocalTetraEdge;
  HardCodedTriangleOracle result;
  switch (mask) {
  case 0U:
  case 15U:
    break;
  case 1U:
    result = {1U, {{{{Edge{0U, 1U}, Edge{0U, 2U}, Edge{0U, 3U}}}, {}}}};
    break;
  case 2U:
    result = {1U, {{{{Edge{0U, 1U}, Edge{1U, 2U}, Edge{1U, 3U}}}, {}}}};
    break;
  case 4U:
    result = {1U, {{{{Edge{0U, 2U}, Edge{1U, 2U}, Edge{2U, 3U}}}, {}}}};
    break;
  case 8U:
    result = {1U, {{{{Edge{0U, 3U}, Edge{1U, 3U}, Edge{2U, 3U}}}, {}}}};
    break;
  case 3U:
    result = {2U,
              {{{{Edge{0U, 2U}, Edge{0U, 3U}, Edge{1U, 3U}}},
                {{Edge{0U, 2U}, Edge{1U, 3U}, Edge{1U, 2U}}}}}};
    break;
  case 5U:
    result = {2U,
              {{{{Edge{0U, 1U}, Edge{0U, 3U}, Edge{2U, 3U}}},
                {{Edge{0U, 1U}, Edge{2U, 3U}, Edge{1U, 2U}}}}}};
    break;
  case 6U:
    result = {2U,
              {{{{Edge{0U, 1U}, Edge{1U, 3U}, Edge{2U, 3U}}},
                {{Edge{0U, 1U}, Edge{2U, 3U}, Edge{0U, 2U}}}}}};
    break;
  case 9U:
    result = {2U,
              {{{{Edge{0U, 1U}, Edge{0U, 2U}, Edge{2U, 3U}}},
                {{Edge{0U, 1U}, Edge{2U, 3U}, Edge{1U, 3U}}}}}};
    break;
  case 10U:
    result = {2U,
              {{{{Edge{0U, 1U}, Edge{1U, 2U}, Edge{2U, 3U}}},
                {{Edge{0U, 1U}, Edge{2U, 3U}, Edge{0U, 3U}}}}}};
    break;
  case 12U:
    result = {2U,
              {{{{Edge{0U, 2U}, Edge{1U, 2U}, Edge{1U, 3U}}},
                {{Edge{0U, 2U}, Edge{1U, 3U}, Edge{0U, 3U}}}}}};
    break;
  case 7U:
    result = {1U, {{{{Edge{0U, 3U}, Edge{2U, 3U}, Edge{1U, 3U}}}, {}}}};
    break;
  case 11U:
    result = {1U, {{{{Edge{0U, 2U}, Edge{2U, 3U}, Edge{1U, 2U}}}, {}}}};
    break;
  case 13U:
    result = {1U, {{{{Edge{0U, 1U}, Edge{1U, 3U}, Edge{1U, 2U}}}, {}}}};
    break;
  case 14U:
    result = {1U, {{{{Edge{0U, 1U}, Edge{0U, 3U}, Edge{0U, 2U}}}, {}}}};
    break;
  default:
    result.count = std::numeric_limits<std::size_t>::max();
    break;
  }
  return result;
}

void check_c2a_key_and_triangle_oracle() {
  PatchGridCanonicalInterpolationVertexKey key;
  require(PatchGridCanonicalInterpolationVertexKey::try_make(7U, 3U, key) &&
              key.valid() && key.first_vertex_linear_index() == 3U &&
              key.second_vertex_linear_index() == 7U,
          "C2a reverse edge key was not canonicalized");
  PatchGridCanonicalInterpolationVertexKey forward;
  require(PatchGridCanonicalInterpolationVertexKey::try_make(3U, 7U, forward) &&
              forward == key,
          "C2a forward/reverse keys differ");
  require(!PatchGridCanonicalInterpolationVertexKey::try_make(3U, 3U, key) &&
              !key.valid(),
          "C2a equal endpoint key was accepted or not reset");
  const std::size_t invalid = std::numeric_limits<std::size_t>::max();
  require(
      !PatchGridCanonicalInterpolationVertexKey::try_make(invalid, 3U, key) &&
          !key.valid() &&
          !PatchGridCanonicalInterpolationVertexKey::try_make(3U, invalid,
                                                              key) &&
          !key.valid(),
      "C2a invalid endpoint key was accepted or not reset");

  for (std::uint8_t mask = 0U; mask < 16U; ++mask) {
    const HardCodedTriangleOracle expected = hard_coded_triangle_oracle(mask);
    std::array<PatchGridPreOrientationTriangleEdges, 2> actual{};
    const std::size_t count =
        enumerate_patch_grid_pre_orientation_triangles(mask, actual);
    require(count == expected.count, "C2a tetra mask triangle count changed");
    for (std::size_t triangle = 0U; triangle < count; ++triangle)
      require(actual[triangle] == expected.triangles[triangle],
              "C2a tetra mask pre-orientation order changed");
  }
  std::array<PatchGridPreOrientationTriangleEdges, 2> invalid_output{};
  invalid_output[0].edges = {PatchGridLocalTetraEdge{1U, 2U},
                             PatchGridLocalTetraEdge{1U, 3U},
                             PatchGridLocalTetraEdge{2U, 3U}};
  require(enumerate_patch_grid_pre_orientation_triangles(
              0x10U, invalid_output) == invalid &&
              invalid_output ==
                  std::array<PatchGridPreOrientationTriangleEdges, 2>{},
          "C2a high mask bits did not fail with zeroed output");
}

std::array<PatchGridVertexMagnitudeRequest, 8>
exact_cell_vertex_requests(const PatchMappingGridSpec &spec,
                           const PatchGridActiveCellRequest &cell) {
  // Synthetic C1 input: these requests exercise already-classified topology
  // source metadata. This fixture does not independently re-prove the
  // atom-union or probe-expanded-union geometry behind that metadata.
  std::array<PatchGridVertexMagnitudeRequest, 8> requests{};
  for (std::size_t corner = 0U; corner < requests.size(); ++corner) {
    const auto &offset = patch_grid_production_corners[corner];
    PatchGridVertexMagnitudeRequest &request = requests[corner];
    request.ix = cell.ix + offset[0];
    request.iy = cell.iy + offset[1];
    request.iz = cell.iz + offset[2];
    std::size_t linear = 0U;
    require(try_patch_grid_vertex_linear_index(spec, request.ix, request.iy,
                                               request.iz, linear),
            "cell corner request is out of range");
    const bool inside = (cell.inside_mask & (std::uint8_t{1} << corner)) != 0U;
    request.topology_state = inside ? PatchGridTopologyState::inside
                                    : PatchGridTopologyState::outside;
    request.topology_source =
        inside ? PatchGridTopologySource::inside_atom_union
               : PatchGridTopologySource::outside_probe_expanded_union;
  }
  return requests;
}

const PatchGridSignChangingEdgeLineage *
find_cell_corner_edge(const PreparedPatchGridCellLineage &prepared,
                      std::uint8_t first, std::uint8_t second) {
  if (second < first)
    std::swap(first, second);
  const auto found =
      std::find_if(prepared.edges.begin(), prepared.edges.end(),
                   [&](const PatchGridSignChangingEdgeLineage &edge) {
                     std::uint8_t edge_first = edge.first_cell_corner;
                     std::uint8_t edge_second = edge.second_cell_corner;
                     if (edge_second < edge_first)
                       std::swap(edge_first, edge_second);
                     return edge_first == first && edge_second == second;
                   });
  return found == prepared.edges.end() ? nullptr : &*found;
}

std::size_t
find_tetra_local_edge_occurrence(const PreparedPatchGridCellLineage &prepared,
                                 std::size_t tetra_position,
                                 PatchGridLocalTetraEdge local_edge) {
  require(tetra_position < prepared.tetrahedra.size() &&
              tetra_position + 1U < prepared.tetra_edge_offsets.size() &&
              local_edge.first < 4U && local_edge.second < 4U,
          "invalid tetra edge lookup request");
  const auto &tetra = prepared.tetrahedra[tetra_position];
  std::uint8_t first = tetra.cell_corners[local_edge.first];
  std::uint8_t second = tetra.cell_corners[local_edge.second];
  if (second < first)
    std::swap(first, second);
  for (std::size_t edge = prepared.tetra_edge_offsets[tetra_position];
       edge < prepared.tetra_edge_offsets[tetra_position + 1U]; ++edge) {
    std::uint8_t actual_first = prepared.edges[edge].first_cell_corner;
    std::uint8_t actual_second = prepared.edges[edge].second_cell_corner;
    if (actual_second < actual_first)
      std::swap(actual_first, actual_second);
    if (actual_first == first && actual_second == second)
      return edge;
  }
  return std::numeric_limits<std::size_t>::max();
}

double canonical_triangle_twice_area(const Vec3 &a, const Vec3 &b,
                                     const Vec3 &c) {
  const double abx = b.x - a.x;
  const double aby = b.y - a.y;
  const double abz = b.z - a.z;
  const double acx = c.x - a.x;
  const double acy = c.y - a.y;
  const double acz = c.z - a.z;
  const double cx = aby * acz - abz * acy;
  const double cy = abz * acx - abx * acz;
  const double cz = abx * acy - aby * acx;
  return std::hypot(cx, cy, cz);
}

void check_c2a_runtime_happy_path() {
  SphereFixture fixture =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 1U);
  PatchGridActiveCellRequest cell;
  cell.ix = grid_coordinate_index(fixture.spec.lower.x,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iy = grid_coordinate_index(fixture.spec.lower.y,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iz = grid_coordinate_index(fixture.spec.lower.z,
                                  fixture.spec.grid_spacing, 0.0);
  cell.inside_mask = 1U;
  const auto requests = exact_cell_vertex_requests(fixture.spec, cell);
  const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, requests);
  const auto cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, std::span(&cell, 1U));
  const auto triangles = build_patch_grid_pre_orientation_triangle_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, cells);

  require((triangles.status == "pass" ||
           triangles.status == "pass_with_computed_zero_canonical_area" ||
           triangles.status == "pass_with_candidate_attribution_uncertified") &&
              triangles.interpolation_occurrences == cells.edges.size() &&
              triangles.edge_occurrence_to_canonical_vertex.size() ==
                  cells.edges.size() &&
              triangles.canonical_interpolation_vertices ==
                  triangles.canonical_vertices.size() &&
              triangles.canonical_vertex_occurrence_offsets.size() ==
                  triangles.canonical_vertices.size() + 1U &&
              triangles.canonical_vertex_candidate_ref_offsets.size() ==
                  triangles.canonical_vertices.size() + 1U &&
              triangles.canonical_vertex_relation_offsets.size() ==
                  triangles.canonical_vertices.size() + 1U &&
              triangles.canonical_vertex_exact_location_ref_offsets.size() ==
                  triangles.canonical_vertices.size() + 1U &&
              triangles.triangle_candidate_ref_offsets.size() ==
                  triangles.triangles.size() + 1U &&
              triangles.interpolation_identity_complete &&
              triangles.triangle_occurrence_lineage_complete &&
              triangles.vertex_triplet_patch_identity_classification_complete,
          "C2a public builder did not produce complete canonical lineage");
  require(!triangles.canonical_vertices.empty() &&
              triangles.canonical_vertices.size() < cells.edges.size(),
          "C2a fixture did not exercise shared-edge canonical deduplication");

  const auto monotone = [](const std::vector<std::size_t> &offsets,
                           std::size_t sentinel) {
    return !offsets.empty() && offsets.front() == 0U &&
           offsets.back() == sentinel &&
           std::is_sorted(offsets.begin(), offsets.end());
  };
  require(monotone(triangles.canonical_vertex_occurrence_offsets,
                   triangles.canonical_vertex_edge_occurrences.size()) &&
              monotone(triangles.canonical_vertex_candidate_ref_offsets,
                       triangles.canonical_vertex_candidate_refs.size()) &&
              monotone(triangles.canonical_vertex_relation_offsets,
                       triangles.canonical_vertex_relations.size()) &&
              monotone(triangles.canonical_vertex_exact_location_ref_offsets,
                       triangles.canonical_vertex_exact_location_refs.size()) &&
              monotone(triangles.triangle_candidate_ref_offsets,
                       triangles.triangle_candidate_refs.size()),
          "C2a CSR sentinel/monotonicity contract failed");

  for (std::size_t canonical = 0U;
       canonical < triangles.canonical_vertices.size(); ++canonical) {
    const auto &vertex = triangles.canonical_vertices[canonical];
    require(vertex.key.valid() && vertex.canonical_interpolation_complete &&
                vertex.canonical_edge_occurrence_geometry_bitwise_consistent &&
                !vertex.directed_production_geometry_bitwise_certified &&
                (canonical == 0U ||
                 triangles.canonical_vertices[canonical - 1U].key < vertex.key),
            "C2a canonical vertex ordering/certificate boundary failed");
    const std::size_t occurrence_begin =
        triangles.canonical_vertex_occurrence_offsets[canonical];
    const std::size_t occurrence_end =
        triangles.canonical_vertex_occurrence_offsets[canonical + 1U];
    require(occurrence_begin < occurrence_end,
            "C2a canonical vertex lost all occurrences");
    for (std::size_t position = occurrence_begin; position < occurrence_end;
         ++position) {
      const std::size_t edge_position =
          triangles.canonical_vertex_edge_occurrences[position];
      require(
          edge_position < cells.edges.size() &&
              triangles.edge_occurrence_to_canonical_vertex[edge_position] ==
                  canonical,
          "C2a occurrence-to-canonical reverse CSR is inconsistent");
      PatchGridCanonicalInterpolationVertexKey expected_key;
      const auto &edge = cells.edges[edge_position];
      require(PatchGridCanonicalInterpolationVertexKey::try_make(
                  edge.first_vertex_linear_index,
                  edge.second_vertex_linear_index, expected_key) &&
                  expected_key == vertex.key &&
                  std::bit_cast<std::uint64_t>(edge.interpolation_parameter) ==
                      std::bit_cast<std::uint64_t>(
                          vertex.canonical_interpolation_parameter) &&
                  std::bit_cast<std::uint64_t>(edge.interpolated_point.x) ==
                      std::bit_cast<std::uint64_t>(
                          vertex.canonical_interpolation_point.x) &&
                  std::bit_cast<std::uint64_t>(edge.interpolated_point.y) ==
                      std::bit_cast<std::uint64_t>(
                          vertex.canonical_interpolation_point.y) &&
                  std::bit_cast<std::uint64_t>(edge.interpolated_point.z) ==
                      std::bit_cast<std::uint64_t>(
                          vertex.canonical_interpolation_point.z),
              "C2a shared occurrence changed canonical interpolation bits");
    }
    const std::size_t candidate_begin =
        triangles.canonical_vertex_candidate_ref_offsets[canonical];
    const std::size_t candidate_end =
        triangles.canonical_vertex_candidate_ref_offsets[canonical + 1U];
    require(
        candidate_end <= triangles.canonical_vertex_candidate_refs.size() &&
            std::is_sorted(triangles.canonical_vertex_candidate_refs.begin() +
                               candidate_begin,
                           triangles.canonical_vertex_candidate_refs.begin() +
                               candidate_end) &&
            std::adjacent_find(
                triangles.canonical_vertex_candidate_refs.begin() +
                    candidate_begin,
                triangles.canonical_vertex_candidate_refs.begin() +
                    candidate_end) ==
                triangles.canonical_vertex_candidate_refs.begin() +
                    candidate_end,
        "C2a canonical candidate refs are not sorted unique");
    const std::size_t relation_begin =
        triangles.canonical_vertex_relation_offsets[canonical];
    const std::size_t relation_end =
        triangles.canonical_vertex_relation_offsets[canonical + 1U];
    const auto relation_less =
        [](const PatchGridCanonicalVertexPatchRelationEvidence &left,
           const PatchGridCanonicalVertexPatchRelationEvidence &right) {
          if (left.first != right.first)
            return left.first < right.first;
          if (left.second != right.second)
            return left.second < right.second;
          return static_cast<std::uint8_t>(left.relation) <
                 static_cast<std::uint8_t>(right.relation);
        };
    require(
        relation_end <= triangles.canonical_vertex_relations.size() &&
            std::is_sorted(
                triangles.canonical_vertex_relations.begin() + relation_begin,
                triangles.canonical_vertex_relations.begin() + relation_end,
                relation_less) &&
            std::adjacent_find(
                triangles.canonical_vertex_relations.begin() + relation_begin,
                triangles.canonical_vertex_relations.begin() + relation_end,
                [](const PatchGridCanonicalVertexPatchRelationEvidence &left,
                   const PatchGridCanonicalVertexPatchRelationEvidence &right) {
                  return left.first == right.first &&
                         left.second == right.second &&
                         left.relation == right.relation;
                }) ==
                triangles.canonical_vertex_relations.begin() + relation_end,
        "C2a canonical relation evidence is not normalized");
    const std::size_t exact_begin =
        triangles.canonical_vertex_exact_location_ref_offsets[canonical];
    const std::size_t exact_end =
        triangles.canonical_vertex_exact_location_ref_offsets[canonical + 1U];
    require(exact_end <=
                    triangles.canonical_vertex_exact_location_refs.size() &&
                std::is_sorted(
                    triangles.canonical_vertex_exact_location_refs.begin() +
                        exact_begin,
                    triangles.canonical_vertex_exact_location_refs.begin() +
                        exact_end) &&
                std::includes(
                    triangles.canonical_vertex_candidate_refs.begin() +
                        candidate_begin,
                    triangles.canonical_vertex_candidate_refs.begin() +
                        candidate_end,
                    triangles.canonical_vertex_exact_location_refs.begin() +
                        exact_begin,
                    triangles.canonical_vertex_exact_location_refs.begin() +
                        exact_end) &&
                vertex.exact_analytic_location_certified ==
                    (exact_end > exact_begin) &&
                vertex.candidate_refs_present_without_exact_location ==
                    (candidate_end > candidate_begin &&
                     !vertex.exact_analytic_location_certified),
            "C2a exact-location refs were mixed with candidate-only refs");
  }

  std::size_t expected_triangle_count = 0U;
  for (std::size_t tetra_position = 0U;
       tetra_position < cells.tetrahedra.size(); ++tetra_position) {
    const HardCodedTriangleOracle expected = hard_coded_triangle_oracle(
        cells.tetrahedra[tetra_position].inside_mask);
    require(expected.count != std::numeric_limits<std::size_t>::max(),
            "C2a happy fixture has invalid tetra mask");
    for (std::size_t ordinal = 0U; ordinal < expected.count; ++ordinal) {
      require(expected_triangle_count < triangles.triangles.size(),
              "C2a omitted a hard-coded triangle occurrence");
      const auto &actual = triangles.triangles[expected_triangle_count++];
      require(actual.tetra_position == tetra_position &&
                  actual.triangle_ordinal == ordinal,
              "C2a triangle occurrence ordering changed");
      for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
        const std::size_t edge_position = find_tetra_local_edge_occurrence(
            cells, tetra_position, expected.triangles[ordinal].edges[vertex]);
        require(
            edge_position != std::numeric_limits<std::size_t>::max() &&
                actual.source_edge_occurrence_positions[vertex] ==
                    edge_position &&
                actual.canonical_vertex_positions[vertex] ==
                    triangles
                        .edge_occurrence_to_canonical_vertex[edge_position],
            "C2a triangle source/canonical order differs from literal oracle");
      }
      const Vec3 &a =
          triangles.canonical_vertices[actual.canonical_vertex_positions[0]]
              .canonical_interpolation_point;
      const Vec3 &b =
          triangles.canonical_vertices[actual.canonical_vertex_positions[1]]
              .canonical_interpolation_point;
      const Vec3 &c =
          triangles.canonical_vertices[actual.canonical_vertex_positions[2]]
              .canonical_interpolation_point;
      require(
          std::bit_cast<std::uint64_t>(actual.canonical_twice_area) ==
                  std::bit_cast<std::uint64_t>(
                      canonical_triangle_twice_area(a, b, c)) &&
              actual.orientation_stage ==
                  PatchGridTriangleOrientationStage::pre_atom_sdf_orientation &&
              !actual.planar_face_surface_location_certified &&
              !actual.final_mesh_orientation_certified,
          "C2a canonical triangle geometry/orientation boundary changed");
      const std::size_t candidate_begin =
          triangles
              .triangle_candidate_ref_offsets[expected_triangle_count - 1U];
      const std::size_t candidate_end =
          triangles.triangle_candidate_ref_offsets[expected_triangle_count];
      require(
          candidate_end <= triangles.triangle_candidate_refs.size() &&
              std::is_sorted(
                  triangles.triangle_candidate_refs.begin() + candidate_begin,
                  triangles.triangle_candidate_refs.begin() + candidate_end) &&
              actual.single_common_candidate_ref_available ==
                  (actual.patch_candidate_summary_complete &&
                   candidate_end == candidate_begin + 1U) &&
              (!actual.single_common_candidate_ref_available ||
               actual.single_common_candidate_ref ==
                   triangles.triangle_candidate_refs[candidate_begin]),
          "C2a triangle candidate intersection summary is inconsistent");
    }
  }
  require(expected_triangle_count == triangles.triangles.size() &&
              triangles.triangle_occurrences == triangles.triangles.size() &&
              triangles.finite_triangle_occurrences +
                      triangles.nonfinite_triangle_occurrences ==
                  triangles.triangle_occurrences &&
              triangles.computed_zero_area_triangle_occurrences <=
                  triangles.finite_triangle_occurrences &&
              !triangles.planar_face_surface_location_certified &&
              !triangles.final_mesh_orientation_certified,
          "C2a triangle counters or global claim boundary failed");

  const PatchGridSignChangingEdgeLineage *diagonal =
      find_cell_corner_edge(cells, 0U, 5U);
  require(diagonal != nullptr && !diagonal->surface_location_certified,
          "C2a candidate-only diagonal fixture disappeared");
  const std::size_t diagonal_edge =
      static_cast<std::size_t>(diagonal - cells.edges.data());
  const std::size_t diagonal_vertex =
      triangles.edge_occurrence_to_canonical_vertex[diagonal_edge];
  require(diagonal_vertex < triangles.canonical_vertices.size() &&
              !triangles.canonical_vertices[diagonal_vertex]
                   .exact_analytic_location_certified &&
              triangles.canonical_vertices[diagonal_vertex]
                  .candidate_refs_present_without_exact_location,
          "C2a diagonal candidate evidence upgraded exact location");
  require(triangles.single_common_candidate_triangles > 0U &&
              std::all_of(
                  triangles.triangles.begin(), triangles.triangles.end(),
                  [](const PatchGridPreOrientationTriangleOccurrenceLineage
                         &triangle) {
                    return !triangle.planar_face_surface_location_certified &&
                           !triangle.final_mesh_orientation_certified;
                  }),
          "C2a common candidate incorrectly certified a planar face");

  std::cout << "slice_c2a_happy_status=" << triangles.status
            << " slice_c2a_checksum=" << triangle_lineage_checksum(triangles)
            << " slice_c2a_vertices=" << triangles.canonical_vertices.size()
            << " slice_c2a_occurrences=" << triangles.interpolation_occurrences
            << " slice_c2a_triangles=" << triangles.triangle_occurrences
            << " slice_c2a_exact_vertices=" << triangles.exact_location_vertices
            << " slice_c2a_candidate_only_vertices="
            << triangles.candidate_refs_without_exact_location_vertices
            << " slice_c2a_triplet_certified="
            << triangles.vertex_triplet_patch_identity_certified_triangles
            << '\n';
}

void require_c2a_failure(
    const PreparedPatchGridTriangleOccurrenceLineage &result,
    std::span<const std::string> accepted_statuses, const char *message) {
  require(std::find(accepted_statuses.begin(), accepted_statuses.end(),
                    result.status) != accepted_statuses.end() &&
              result.edge_occurrence_to_canonical_vertex.empty() &&
              result.canonical_vertices.empty() &&
              result.canonical_vertex_occurrence_offsets.empty() &&
              result.canonical_vertex_edge_occurrences.empty() &&
              result.canonical_vertex_candidate_ref_offsets.empty() &&
              result.canonical_vertex_candidate_refs.empty() &&
              result.canonical_vertex_relation_offsets.empty() &&
              result.canonical_vertex_relations.empty() &&
              result.canonical_vertex_exact_location_ref_offsets.empty() &&
              result.canonical_vertex_exact_location_refs.empty() &&
              result.triangles.empty() &&
              result.triangle_candidate_ref_offsets.empty() &&
              result.triangle_candidate_refs.empty() &&
              !result.interpolation_identity_complete &&
              !result.patch_candidate_attribution_complete &&
              !result.triangle_occurrence_lineage_complete &&
              !result.vertex_triplet_patch_identity_classification_complete &&
              !result.planar_face_surface_location_certified &&
              !result.final_mesh_orientation_certified,
          message);
}

void check_c2a_public_builder_fail_closed() {
  SphereFixture fixture =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 1U);
  PatchGridActiveCellRequest cell;
  cell.ix = grid_coordinate_index(fixture.spec.lower.x,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iy = grid_coordinate_index(fixture.spec.lower.y,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iz = grid_coordinate_index(fixture.spec.lower.z,
                                  fixture.spec.grid_spacing, 0.0);
  cell.inside_mask = 1U;
  const auto requests = exact_cell_vertex_requests(fixture.spec, cell);
  const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, requests);
  const auto cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, std::span(&cell, 1U));
  const auto baseline = build_patch_grid_pre_orientation_triangle_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, cells);
  require(baseline.status == "pass" && !baseline.triangles.empty(),
          "C2a mutation fixture did not start from a pass payload");

  const std::array<std::string, 1> invalid_cell = {"invalid_cell_payload"};
  const std::array<std::string, 1> invalid_vertex = {"invalid_vertex_payload"};
  const std::array<std::string, 3> canonical_conflict = {
      "invalid_cell_payload", "canonical_geometry_conflict",
      "canonical_evidence_conflict"};
  const auto build =
      [&](const PatchRayGrid &grid, const PatchMappingGridSpec &spec,
          const PreparedPatchCrossingWitnesses &crossings,
          const PreparedPatchGridVertexMagnitudeWitnesses &vertex_payload,
          const PreparedPatchGridCellLineage &cell_payload) {
        return build_patch_grid_pre_orientation_triangle_lineage(
            *fixture.complex, fixture.index, grid, spec, crossings,
            vertex_payload, cell_payload);
      };
  const auto fail_cell = [&](const PreparedPatchGridCellLineage &payload,
                             const char *message) {
    require_c2a_failure(
        build(fixture.on, fixture.spec, fixture.prepared, vertices, payload),
        invalid_cell, message);
  };
  const auto fail_vertex =
      [&](const PreparedPatchGridVertexMagnitudeWitnesses &payload,
          const char *message) {
        require_c2a_failure(
            build(fixture.on, fixture.spec, fixture.prepared, payload, cells),
            invalid_vertex, message);
      };

  PreparedPatchGridCellLineage swapped_tetra = cells;
  std::swap(swapped_tetra.tetrahedra[0], swapped_tetra.tetrahedra[1]);
  fail_cell(swapped_tetra, "C2a accepted swapped tetra order");
  PreparedPatchGridCellLineage bad_cell_sentinel = cells;
  ++bad_cell_sentinel.cell_tetra_offsets.back();
  fail_cell(bad_cell_sentinel, "C2a accepted malformed cell CSR sentinel");
  PreparedPatchGridCellLineage bad_edge_monotonic = cells;
  bad_edge_monotonic.tetra_edge_offsets[1] =
      bad_edge_monotonic.tetra_edge_offsets[2] + 1U;
  fail_cell(bad_edge_monotonic, "C2a accepted nonmonotone tetra edge CSR");
  PreparedPatchGridCellLineage bad_relation_range = cells;
  ++bad_relation_range.edge_relation_offsets.back();
  fail_cell(bad_relation_range, "C2a accepted out-of-range relation CSR");
  PreparedPatchGridCellLineage bad_location_range = cells;
  ++bad_location_range.edge_location_ref_offsets.back();
  fail_cell(bad_location_range, "C2a accepted out-of-range exact-location CSR");

  require(!cells.endpoint_patch_refs.empty() && !cells.relations.empty() &&
              !cells.location_patch_refs.empty(),
          "C2a mutation fixture lacks evidence payload");
  PreparedPatchGridCellLineage bad_endpoint_ref = cells;
  bad_endpoint_ref.endpoint_patch_refs.front() = {};
  fail_cell(bad_endpoint_ref, "C2a accepted altered endpoint ref");
  PreparedPatchGridCellLineage bad_relation_ref = cells;
  bad_relation_ref.relations.front().first = {};
  fail_cell(bad_relation_ref, "C2a accepted altered relation ref/pair");
  PreparedPatchGridCellLineage bad_relation_enum = cells;
  bad_relation_enum.relations.front().relation =
      static_cast<PatchRelation>(std::numeric_limits<std::uint8_t>::max());
  fail_cell(bad_relation_enum, "C2a accepted invalid relation enum");
  PreparedPatchGridCellLineage bad_relation_count = cells;
  ++bad_relation_count.relation_counts[0];
  fail_cell(bad_relation_count, "C2a accepted stale relation counters");
  PreparedPatchGridCellLineage bad_status = cells;
  bad_status.status = "pass";
  fail_cell(bad_status, "C2a accepted stale C1 status summary");
  PreparedPatchGridCellLineage bad_relation_summary = cells;
  bad_relation_summary.patch_relation_classification_complete =
      !bad_relation_summary.patch_relation_classification_complete;
  fail_cell(bad_relation_summary, "C2a accepted stale relation summary");

  SphereFixture duplicated;
  duplicated.complex = std::make_unique<PatchComplex>(build_patch_complex(
      std::vector<Atom>{Atom{0.0, 0.0, 0.0, 0.0, 1.0}}, 1.0, nullptr, false));
  duplicated.complex->contacts.push_back(duplicated.complex->contacts.front());
  duplicated.index = build_audit_projection_index(*duplicated.complex);
  duplicated.spec = sphere_grid(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0);
  TitanExecutor duplicate_executor(1U);
  PatchRayBuildOptions duplicate_options;
  duplicate_options.retain_provenance = true;
  duplicated.on = build_patch_ray_grid(
      *duplicated.complex, duplicated.spec.lower, duplicated.spec.grid_spacing,
      duplicated.spec.nx, duplicated.spec.ny, duplicated.spec.nz,
      &duplicate_executor, duplicate_options);
  duplicated.prepared = build_certified_patch_crossing_witnesses(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec);
  const auto duplicated_requests =
      exact_cell_vertex_requests(duplicated.spec, cell);
  const auto duplicated_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
      duplicated.prepared, duplicated_requests);
  const auto duplicated_cells = build_patch_grid_cell_lineage(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
      duplicated.prepared, duplicated_vertices, std::span(&cell, 1U));
  const auto duplicated_triangles =
      build_patch_grid_pre_orientation_triangle_lineage(
          *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
          duplicated.prepared, duplicated_vertices, duplicated_cells);
  require(duplicated_triangles.status.rfind("pass", 0U) == 0U,
          "C2a valid multi-ref relation fixture did not build");
  PreparedPatchGridCellLineage bad_relation_pair_order = duplicated_cells;
  bool swapped_relation_pair = false;
  for (std::size_t edge = 0U; edge < bad_relation_pair_order.edges.size();
       ++edge) {
    const std::size_t begin =
        bad_relation_pair_order.edge_relation_offsets[edge];
    const std::size_t end =
        bad_relation_pair_order.edge_relation_offsets[edge + 1U];
    if (end - begin > 1U) {
      std::swap(bad_relation_pair_order.relations[begin],
                bad_relation_pair_order.relations[begin + 1U]);
      swapped_relation_pair = true;
      break;
    }
  }
  require(swapped_relation_pair,
          "C2a multi-ref fixture lacks a relation Cartesian product");
  const auto bad_pair_result =
      build_patch_grid_pre_orientation_triangle_lineage(
          *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
          duplicated.prepared, duplicated_vertices, bad_relation_pair_order);
  require_c2a_failure(bad_pair_result, invalid_cell,
                      "C2a accepted reordered valid relation pair");

  PreparedPatchGridCellLineage bad_scalar_endpoint = cells;
  bad_scalar_endpoint.edges.front().endpoint_scalar_status[0] =
      PatchGridEdgeScalarLineageStatus::clamp_floor;
  fail_cell(bad_scalar_endpoint, "C2a accepted altered endpoint scalar status");
  PreparedPatchGridCellLineage bad_scalar_aggregate = cells;
  bad_scalar_aggregate.edges.front().scalar_lineage_complete =
      !bad_scalar_aggregate.edges.front().scalar_lineage_complete;
  fail_cell(bad_scalar_aggregate, "C2a accepted altered scalar aggregate");
  PreparedPatchGridCellLineage bad_topology_endpoint = cells;
  bad_topology_endpoint.edges.front().endpoint_topology_status[0] =
      PatchGridTopologyEvidenceStatus::unknown;
  fail_cell(bad_topology_endpoint,
            "C2a accepted altered endpoint topology status");
  PreparedPatchGridCellLineage bad_topology_aggregate = cells;
  bad_topology_aggregate.edges.front().topology_sign_evidence_complete = false;
  fail_cell(bad_topology_aggregate, "C2a accepted altered topology aggregate");

  PreparedPatchGridCellLineage bad_location_ref = cells;
  bad_location_ref.location_patch_refs.front() = {};
  fail_cell(bad_location_ref, "C2a accepted altered exact-location ref");
  const auto certified_edge =
      std::find_if(cells.edges.begin(), cells.edges.end(),
                   [](const PatchGridSignChangingEdgeLineage &edge) {
                     return edge.surface_location_certified;
                   });
  require(certified_edge != cells.edges.end(),
          "C2a mutation fixture lacks exact-location edge");
  const std::size_t certified_position =
      static_cast<std::size_t>(certified_edge - cells.edges.begin());
  PreparedPatchGridCellLineage bad_crossing = cells;
  PatchCrossingAddress altered_crossing;
  require(PatchCrossingAddress::try_make(
              certified_edge->location_crossing.axis(),
              certified_edge->location_crossing.row(),
              certified_edge->location_crossing.crossing() + 1U,
              altered_crossing),
          "C2a could not form altered crossing address");
  bad_crossing.edges[certified_position].location_crossing = altered_crossing;
  fail_cell(bad_crossing, "C2a accepted altered exact crossing address");
  PreparedPatchGridCellLineage bad_residual = cells;
  bad_residual.edges[certified_position].interpolation_residual =
      std::nextafter(certified_edge->interpolation_residual,
                     std::numeric_limits<double>::infinity());
  fail_cell(bad_residual, "C2a accepted altered interpolation residual");
  PreparedPatchGridCellLineage bad_tolerance = cells;
  bad_tolerance.edges[certified_position].interpolation_tolerance =
      std::nextafter(certified_edge->interpolation_tolerance,
                     std::numeric_limits<double>::infinity());
  fail_cell(bad_tolerance, "C2a accepted altered interpolation tolerance");

  PreparedPatchGridVertexMagnitudeWitnesses bad_owner = vertices;
  const auto candidate = std::find_if(
      bad_owner.candidates.begin(), bad_owner.candidates.end(),
      [&](const PatchGridMagnitudeSourceCandidate &item) {
        return item.crossing == certified_edge->location_crossing &&
               item.raw_exact_owner && item.final_scalar_exact_owner;
      });
  require(candidate != bad_owner.candidates.end(),
          "C2a mutation fixture lacks eligible analytic owner");
  candidate->raw_tolerance_tied = true;
  require_c2a_failure(
      build(fixture.on, fixture.spec, fixture.prepared, bad_owner, cells),
      invalid_cell, "C2a accepted altered eligible-owner/tie evidence");
  PreparedPatchGridVertexMagnitudeWitnesses bad_sparse = vertices;
  bad_sparse.vertices.front().sparse_magnitude = 1.0F;
  fail_vertex(bad_sparse, "C2a accepted mutated absent sparse payload");
  PreparedPatchGridVertexMagnitudeWitnesses bad_nonanalytic_refs = vertices;
  const auto nonanalytic = std::find_if(
      bad_nonanalytic_refs.candidates.begin(),
      bad_nonanalytic_refs.candidates.end(),
      [](const PatchGridMagnitudeSourceCandidate &item) {
        return item.kind != PatchGridMagnitudeSourceKind::analytic_crossing;
      });
  require(nonanalytic != bad_nonanalytic_refs.candidates.end() &&
              !bad_nonanalytic_refs.candidate_refs.empty(),
          "C2a mutation fixture lacks nonanalytic source");
  nonanalytic->refs_begin = 1U;
  nonanalytic->refs_end = 1U;
  fail_vertex(bad_nonanalytic_refs,
              "C2a accepted nonanalytic source with nonzero ref range");
  PreparedPatchGridVertexMagnitudeWitnesses bad_candidate_sentinel = vertices;
  ++bad_candidate_sentinel.candidate_offsets.back();
  fail_vertex(bad_candidate_sentinel,
              "C2a accepted malformed Slice-B CSR sentinel");
  PreparedPatchGridVertexMagnitudeWitnesses bad_candidate_monotonic = vertices;
  require(bad_candidate_monotonic.candidate_offsets.size() > 2U,
          "C2a mutation fixture lacks Slice-B CSR rows");
  bad_candidate_monotonic.candidate_offsets[1] =
      bad_candidate_monotonic.candidate_offsets[2] + 1U;
  fail_vertex(bad_candidate_monotonic, "C2a accepted nonmonotone Slice-B CSR");
  PreparedPatchGridVertexMagnitudeWitnesses bad_candidate_range = vertices;
  const auto analytic_candidate = std::find_if(
      bad_candidate_range.candidates.begin(),
      bad_candidate_range.candidates.end(),
      [](const PatchGridMagnitudeSourceCandidate &item) {
        return item.kind == PatchGridMagnitudeSourceKind::analytic_crossing;
      });
  require(analytic_candidate != bad_candidate_range.candidates.end(),
          "C2a mutation fixture lacks analytic Slice-B candidate");
  ++analytic_candidate->refs_end;
  fail_vertex(bad_candidate_range,
              "C2a accepted out-of-range Slice-B candidate refs");

  PreparedPatchGridVertexMagnitudeWitnesses stale_vertex_fingerprint = vertices;
  ++stale_vertex_fingerprint.complex_fingerprint.first;
  fail_vertex(stale_vertex_fingerprint,
              "C2a accepted stale Slice-B fingerprint");
  PreparedPatchGridVertexMagnitudeWitnesses stale_vertex_lattice = vertices;
  ++stale_vertex_lattice.lattice_identity.scalar_bits[0];
  fail_vertex(stale_vertex_lattice, "C2a accepted stale Slice-B lattice");
  PreparedPatchCrossingWitnesses stale_crossings = fixture.prepared;
  ++stale_crossings.complex_fingerprint.first;
  require_c2a_failure(
      build(fixture.on, fixture.spec, stale_crossings, vertices, cells),
      invalid_vertex, "C2a accepted stale Slice-A complex fingerprint");
  PreparedPatchCrossingWitnesses stale_crossing_lattice = fixture.prepared;
  ++stale_crossing_lattice.lattice_identity.scalar_bits[0];
  require_c2a_failure(
      build(fixture.on, fixture.spec, stale_crossing_lattice, vertices, cells),
      invalid_vertex, "C2a accepted stale Slice-A lattice");
  PreparedPatchCrossingWitnesses bad_crossing_point = fixture.prepared;
  bad_crossing_point.witnesses.front().point.x =
      std::nextafter(bad_crossing_point.witnesses.front().point.x,
                     std::numeric_limits<double>::infinity());
  require_c2a_failure(
      build(fixture.on, fixture.spec, bad_crossing_point, vertices, cells),
      invalid_vertex, "C2a accepted mutated Slice-A crossing point");
  PreparedPatchCrossingWitnesses bad_crossing_multiplicity = fixture.prepared;
  ++bad_crossing_multiplicity.witnesses.front().total_root_multiplicity;
  require_c2a_failure(build(fixture.on, fixture.spec, bad_crossing_multiplicity,
                            vertices, cells),
                      invalid_vertex,
                      "C2a accepted mutated Slice-A root multiplicity");
  PatchRayGrid stale_grid = fixture.on;
  ++stale_grid.provenance.lattice_identity.scalar_bits[0];
  require_c2a_failure(
      build(stale_grid, fixture.spec, fixture.prepared, vertices, cells),
      invalid_vertex, "C2a accepted stale ray-grid lattice");
  PatchRayGrid stale_grid_complex = fixture.on;
  ++stale_grid_complex.provenance.complex_fingerprint.first;
  require_c2a_failure(build(stale_grid_complex, fixture.spec, fixture.prepared,
                            vertices, cells),
                      invalid_vertex,
                      "C2a accepted stale ray-grid complex fingerprint");
  PatchMappingGridSpec stale_spec = fixture.spec;
  stale_spec.lower.x = std::nextafter(stale_spec.lower.x,
                                      std::numeric_limits<double>::infinity());
  require_c2a_failure(
      build(fixture.on, stale_spec, fixture.prepared, vertices, cells),
      invalid_vertex, "C2a accepted stale mapping grid spec");

  PreparedPatchGridCellLineage geometry_conflict = cells;
  std::size_t duplicate_first = std::numeric_limits<std::size_t>::max();
  std::size_t duplicate_second = std::numeric_limits<std::size_t>::max();
  for (std::size_t first = 0U; first < baseline.canonical_vertices.size();
       ++first) {
    const std::size_t begin =
        baseline.canonical_vertex_occurrence_offsets[first];
    const std::size_t end =
        baseline.canonical_vertex_occurrence_offsets[first + 1U];
    if (end - begin > 1U) {
      duplicate_first = baseline.canonical_vertex_edge_occurrences[begin];
      duplicate_second = baseline.canonical_vertex_edge_occurrences[begin + 1U];
      break;
    }
  }
  require(duplicate_first != std::numeric_limits<std::size_t>::max() &&
              duplicate_second < geometry_conflict.edges.size(),
          "C2a mutation fixture lacks shared canonical edge");
  geometry_conflict.edges[duplicate_second].interpolated_point.x =
      std::nextafter(
          geometry_conflict.edges[duplicate_second].interpolated_point.x,
          std::numeric_limits<double>::infinity());
  require_c2a_failure(build(fixture.on, fixture.spec, fixture.prepared,
                            vertices, geometry_conflict),
                      canonical_conflict,
                      "C2a accepted same-key geometry mutation");
  PreparedPatchGridCellLineage evidence_conflict = cells;
  evidence_conflict.edges[duplicate_second]
      .patch_relation_classification_complete =
      !evidence_conflict.edges[duplicate_second]
           .patch_relation_classification_complete;
  require_c2a_failure(build(fixture.on, fixture.spec, fixture.prepared,
                            vertices, evidence_conflict),
                      canonical_conflict,
                      "C2a accepted same-key evidence mutation");

  const auto empty_cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices,
      std::span<const PatchGridActiveCellRequest>{});
  const auto empty =
      build(fixture.on, fixture.spec, fixture.prepared, vertices, empty_cells);
  require(empty.status == "empty_input" && empty.triangles.empty() &&
              !empty.interpolation_identity_complete &&
              !empty.patch_candidate_attribution_complete &&
              !empty.triangle_occurrence_lineage_complete &&
              !empty.vertex_triplet_patch_identity_classification_complete &&
              !empty.planar_face_surface_location_certified &&
              !empty.final_mesh_orientation_certified,
          "C2a empty input became a vacuous certificate");
  PreparedPatchGridCellLineage bad_empty_binding = empty_cells;
  bad_empty_binding.complex_fingerprint = fixture.prepared.complex_fingerprint;
  require_c2a_failure(build(fixture.on, fixture.spec, fixture.prepared,
                            vertices, bad_empty_binding),
                      invalid_cell,
                      "C2a accepted mutated empty complex binding");
  PreparedPatchGridCellLineage bad_empty_oracle = empty_cells;
  bad_empty_oracle.oracle_accounted_setup_bytes = 1U;
  require_c2a_failure(build(fixture.on, fixture.spec, fixture.prepared,
                            vertices, bad_empty_oracle),
                      invalid_cell,
                      "C2a accepted mutated empty oracle telemetry");
  std::cout << "slice_c2a_mutation_suite=pass"
            << " slice_c2a_empty_status=" << empty.status << '\n';
}

void check_c2a_all_axis_candidate_boundary() {
  SphereFixture fixture =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 0.2, 1U);
  PatchGridActiveCellRequest cell{8U, 9U, 14U, 0x4eU};
  const auto requests = exact_cell_vertex_requests(fixture.spec, cell);
  const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, requests);
  const auto cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, std::span(&cell, 1U));
  const auto triangles = build_patch_grid_pre_orientation_triangle_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, cells);
  require(cells.cells.size() == 1U && cells.tetrahedra.size() == 6U &&
              cells.tetrahedra[0].inside_mask == 12U &&
              !triangles.triangles.empty(),
          "C2a two-inside physical sphere fixture did not reach tetra zero");
  const auto triangle = std::find_if(
      triangles.triangles.begin(), triangles.triangles.end(),
      [](const PatchGridPreOrientationTriangleOccurrenceLineage &item) {
        return item.tetra_position == 0U && item.triangle_ordinal == 0U;
      });
  require(triangle != triangles.triangles.end() &&
              triangle->patch_candidate_summary_complete &&
              triangle->single_common_candidate_ref_available &&
              !triangle->planar_face_surface_location_certified &&
              !triangle->final_mesh_orientation_certified,
          "C2a two-inside all-axis triangle lost candidate-only evidence");
  PatchRefId exact_ref;
  bool have_ref = false;
  std::size_t exact_vertices = 0U;
  std::size_t no_common_axis_vertices = 0U;
  for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
    const std::size_t canonical = triangle->canonical_vertex_positions[vertex];
    const std::size_t source_edge =
        triangle->source_edge_occurrence_positions[vertex];
    require(canonical < triangles.canonical_vertices.size(),
            "C2a triplet canonical vertex is out of range");
    require(source_edge < cells.edges.size() &&
                cells.edges[source_edge].cartesian_axis_edge,
            "C2a all-axis fixture lost a Cartesian source edge");
    if (cells.edges[source_edge].location_status ==
        PatchGridEdgeSurfaceLocationStatus::no_common_axis_crossing)
      ++no_common_axis_vertices;
    const std::size_t begin =
        triangles.canonical_vertex_exact_location_ref_offsets[canonical];
    const std::size_t end =
        triangles.canonical_vertex_exact_location_ref_offsets[canonical + 1U];
    if (triangles.canonical_vertices[canonical]
            .exact_analytic_location_certified) {
      ++exact_vertices;
      require(end == begin + 1U &&
                  end <= triangles.canonical_vertex_exact_location_refs.size(),
              "C2a triplet exact-location ref is not singleton");
      const PatchRefId ref =
          triangles.canonical_vertex_exact_location_refs[begin];
      if (!have_ref) {
        exact_ref = ref;
        have_ref = true;
      } else {
        require(ref == exact_ref,
                "C2a exact triplet vertices disagree on PatchRefId");
      }
    } else {
      require(begin == end,
              "C2a uncertified triplet vertex retained exact-location refs");
    }
    const auto &key = triangles.canonical_vertices[canonical].key;
    const auto find_corner = [&](std::size_t linear) {
      const auto found = std::find(
          cells.cells.front().corner_vertex_linear_indices.begin(),
          cells.cells.front().corner_vertex_linear_indices.end(), linear);
      require(found != cells.cells.front().corner_vertex_linear_indices.end(),
              "C2a triplet key is not a cell edge");
      return static_cast<std::size_t>(
          found - cells.cells.front().corner_vertex_linear_indices.begin());
    };
    const auto &first_corner = patch_grid_production_corners[find_corner(
        key.first_vertex_linear_index())];
    const auto &second_corner = patch_grid_production_corners[find_corner(
        key.second_vertex_linear_index())];
    const unsigned manhattan =
        static_cast<unsigned>(std::abs(static_cast<int>(first_corner[0]) -
                                       static_cast<int>(second_corner[0]))) +
        static_cast<unsigned>(std::abs(static_cast<int>(first_corner[1]) -
                                       static_cast<int>(second_corner[1]))) +
        static_cast<unsigned>(std::abs(static_cast<int>(first_corner[2]) -
                                       static_cast<int>(second_corner[2])));
    require(manhattan == 1U,
            "C2a all-axis candidate fixture contains a tetra diagonal");
  }
  require(!triangle->vertex_triplet_patch_identity_certified && have_ref &&
              exact_vertices == 1U && no_common_axis_vertices == 2U &&
              triangles.vertex_triplet_patch_identity_certified_triangles == 0U,
          "C2a all-axis candidate boundary was incorrectly upgraded");
  std::cout << "slice_c2a_triplet_mask="
            << static_cast<unsigned>(cell.inside_mask)
            << " slice_c2a_triplet_available="
            << triangle->vertex_triplet_patch_identity_certified
            << " slice_c2a_triplet_exact_vertices=" << exact_vertices
            << " slice_c2a_triplet_no_common=" << no_common_axis_vertices
            << " slice_c2a_triplet_count="
            << triangles.vertex_triplet_patch_identity_certified_triangles
            << " slice_c2a_triplet_limitation="
               "all_axis_not_exact_under_magnitude_owner"
            << '\n';
}

struct C2aPipeline {
  PreparedPatchGridVertexMagnitudeWitnesses vertices;
  PreparedPatchGridCellLineage cells;
  PreparedPatchGridTriangleOccurrenceLineage triangles;
};

std::size_t exact_c2a_retained_bytes(
    const PreparedPatchGridTriangleOccurrenceLineage &prepared) {
  return prepared.edge_occurrence_to_canonical_vertex.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertices.capacity() *
             sizeof(PatchGridCanonicalInterpolationVertexLineage) +
         prepared.canonical_vertex_occurrence_offsets.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertex_edge_occurrences.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertex_candidate_ref_offsets.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertex_candidate_refs.capacity() *
             sizeof(PatchRefId) +
         prepared.canonical_vertex_relation_offsets.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertex_relations.capacity() *
             sizeof(PatchGridCanonicalVertexPatchRelationEvidence) +
         prepared.canonical_vertex_exact_location_ref_offsets.capacity() *
             sizeof(std::size_t) +
         prepared.canonical_vertex_exact_location_refs.capacity() *
             sizeof(PatchRefId) +
         prepared.triangles.capacity() *
             sizeof(PatchGridPreOrientationTriangleOccurrenceLineage) +
         prepared.triangle_candidate_ref_offsets.capacity() *
             sizeof(std::size_t) +
         prepared.triangle_candidate_refs.capacity() * sizeof(PatchRefId);
}

C2aPipeline
build_c2a_pipeline(SphereFixture &fixture,
                   std::span<const PatchGridActiveCellRequest> active_cells) {
  std::vector<PatchGridVertexMagnitudeRequest> requests;
  requests.reserve(active_cells.size() * 8U);
  for (const PatchGridActiveCellRequest &cell : active_cells) {
    const auto corners = exact_cell_vertex_requests(fixture.spec, cell);
    requests.insert(requests.end(), corners.begin(), corners.end());
  }
  C2aPipeline result;
  result.vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, requests);
  result.cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, result.vertices, active_cells);
  result.triangles = build_patch_grid_pre_orientation_triangle_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, result.vertices, result.cells);
  return result;
}

void check_c2a_determinism_timing_and_memory() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture serial = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  SphereFixture parallel = build_sphere_fixture(atom, 1.0, 2.0, 8U);
  PatchGridActiveCellRequest centered;
  centered.ix =
      grid_coordinate_index(serial.spec.lower.x, serial.spec.grid_spacing, 0.0);
  centered.iy =
      grid_coordinate_index(serial.spec.lower.y, serial.spec.grid_spacing, 0.0);
  centered.iz =
      grid_coordinate_index(serial.spec.lower.z, serial.spec.grid_spacing, 0.0);
  centered.inside_mask = 1U;
  C2aPipeline serial_center =
      build_c2a_pipeline(serial, std::span(&centered, 1U));
  C2aPipeline parallel_center =
      build_c2a_pipeline(parallel, std::span(&centered, 1U));
  constexpr std::uint64_t expected_c2a_checksum = 5292789444643985863ULL;
  require(grid_checksum(serial.on) == grid_checksum(parallel.on) &&
              witness_checksum(serial.prepared) ==
                  witness_checksum(parallel.prepared) &&
              triangle_lineage_checksum(serial_center.triangles) ==
                  triangle_lineage_checksum(parallel_center.triangles) &&
              triangle_lineage_checksum(serial_center.triangles) ==
                  expected_c2a_checksum &&
              cell_lineage_checksum(serial_center.cells) ==
                  cell_lineage_checksum(parallel_center.cells) &&
              magnitude_witness_checksum(serial_center.vertices) ==
                  magnitude_witness_checksum(parallel_center.vertices),
          "C2a nt=1/8 semantic payload is not deterministic");

  PatchGridActiveCellRequest first{1U, 1U, 1U, 1U};
  PatchGridActiveCellRequest second{4U, 4U, 4U, 1U};
  const std::array<PatchGridActiveCellRequest, 1> one_cell = {first};
  const std::array<PatchGridActiveCellRequest, 2> two_cells = {first, second};
  const std::array<PatchGridActiveCellRequest, 2> reversed_cells = {second,
                                                                    first};
  C2aPipeline one = build_c2a_pipeline(serial, one_cell);
  C2aPipeline two = build_c2a_pipeline(serial, two_cells);
  C2aPipeline reversed = build_c2a_pipeline(serial, reversed_cells);
  require(one.triangles.status.rfind("pass", 0U) == 0U &&
              two.triangles.status.rfind("pass", 0U) == 0U &&
              magnitude_witness_checksum(two.vertices) ==
                  magnitude_witness_checksum(reversed.vertices) &&
              cell_lineage_checksum(two.cells) ==
                  cell_lineage_checksum(reversed.cells) &&
              triangle_lineage_checksum(two.triangles) ==
                  triangle_lineage_checksum(reversed.triangles),
          "C2a unsorted cells changed canonical semantic output");

  for (std::size_t warmup = 0U; warmup < 2U; ++warmup) {
    const auto warmed = build_patch_grid_pre_orientation_triangle_lineage(
        *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
        two.vertices, two.cells);
    require(triangle_lineage_checksum(warmed) ==
                triangle_lineage_checksum(two.triangles),
            "C2a warmup changed semantic checksum");
  }
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto measured = build_patch_grid_pre_orientation_triangle_lineage(
        *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
        two.vertices, two.cells);
    require(triangle_lineage_checksum(measured) ==
                triangle_lineage_checksum(two.triangles),
            "C2a timed sample changed semantic checksum");
    sample = measured.build_ms;
  }
  std::array<double, 7> sorted_samples = samples;
  std::sort(sorted_samples.begin(), sorted_samples.end());

  const double retained_normalized_ratio =
      (static_cast<double>(two.triangles.retained_bytes) /
       static_cast<double>(two.triangles.interpolation_occurrences)) /
      (static_cast<double>(one.triangles.retained_bytes) /
       static_cast<double>(one.triangles.interpolation_occurrences));
  const double peak_normalized_ratio =
      (static_cast<double>(two.triangles.accounted_peak_temporary_bytes) /
       static_cast<double>(two.triangles.interpolation_occurrences)) /
      (static_cast<double>(one.triangles.accounted_peak_temporary_bytes) /
       static_cast<double>(one.triangles.interpolation_occurrences));
  require(one.triangles.retained_bytes > 0U &&
              two.triangles.retained_bytes > one.triangles.retained_bytes &&
              one.triangles.retained_bytes ==
                  exact_c2a_retained_bytes(one.triangles) &&
              two.triangles.retained_bytes ==
                  exact_c2a_retained_bytes(two.triangles) &&
              one.triangles.accounted_peak_temporary_bytes >=
                  one.triangles.retained_bytes &&
              two.triangles.accounted_peak_temporary_bytes >=
                  two.triangles.retained_bytes &&
              retained_normalized_ratio >= 0.5 &&
              retained_normalized_ratio <= 1.5 &&
              peak_normalized_ratio >= 0.5 && peak_normalized_ratio <= 1.5 &&
              two.triangles.canonical_vertex_candidate_refs.size() <=
                  2U * one.triangles.canonical_vertex_candidate_refs.size() &&
              two.triangles.canonical_vertex_relations.size() <=
                  2U * one.triangles.canonical_vertex_relations.size(),
          "C2a retained/accounted memory or evidence growth is superlinear");

  const auto scaled_surface_corner = [](double scale) {
    SphereFixture fixture =
        build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, scale}, scale, scale, 1U);
    PatchGridActiveCellRequest surface_corner;
    surface_corner.ix = grid_coordinate_index(fixture.spec.lower.x,
                                              fixture.spec.grid_spacing, scale);
    surface_corner.iy = grid_coordinate_index(fixture.spec.lower.y,
                                              fixture.spec.grid_spacing, 0.0);
    surface_corner.iz = grid_coordinate_index(fixture.spec.lower.z,
                                              fixture.spec.grid_spacing, 0.0);
    surface_corner.inside_mask = 1U;
    return build_c2a_pipeline(fixture, std::span(&surface_corner, 1U))
        .triangles;
  };
  const auto small_geometry = scaled_surface_corner(0.125);
  const auto unit_geometry = scaled_surface_corner(1.0);
  const auto large_geometry = scaled_surface_corner(8.0);
  require(
      small_geometry.status == "pass" && unit_geometry.status == "pass" &&
          large_geometry.status == "pass" &&
          small_geometry.triangles.size() == unit_geometry.triangles.size() &&
          unit_geometry.triangles.size() == large_geometry.triangles.size() &&
          small_geometry.computed_zero_area_triangle_occurrences == 0U &&
          unit_geometry.computed_zero_area_triangle_occurrences == 0U &&
          large_geometry.computed_zero_area_triangle_occurrences == 0U &&
          small_geometry.nonfinite_triangle_occurrences == 0U &&
          unit_geometry.nonfinite_triangle_occurrences == 0U &&
          large_geometry.nonfinite_triangle_occurrences == 0U,
      "C2a computed-zero observation changed under uniform scaling");
  for (std::size_t triangle = 0U; triangle < unit_geometry.triangles.size();
       ++triangle) {
    require(
        small_geometry.triangles[triangle].canonical_twice_area > 0.0 &&
            unit_geometry.triangles[triangle].canonical_twice_area > 0.0 &&
            large_geometry.triangles[triangle].canonical_twice_area > 0.0 &&
            small_geometry.triangles[triangle].computed_zero_area_observed ==
                unit_geometry.triangles[triangle].computed_zero_area_observed &&
            unit_geometry.triangles[triangle].computed_zero_area_observed ==
                large_geometry.triangles[triangle]
                    .computed_zero_area_observed &&
            unit_geometry.triangles[triangle].computed_zero_area_observed ==
                (unit_geometry.triangles[triangle].canonical_twice_area == 0.0),
        "C2a computed-zero flag is not its documented exact observation");
  }

  std::cout << "slice_c2a_checksum="
            << triangle_lineage_checksum(serial_center.triangles)
            << " slice_c2a_samples_ms=" << samples[0] << ',' << samples[1]
            << ',' << samples[2] << ',' << samples[3] << ',' << samples[4]
            << ',' << samples[5] << ',' << samples[6]
            << " slice_c2a_median_ms=" << sorted_samples[3]
            << " slice_c2a_retained_one=" << one.triangles.retained_bytes
            << " slice_c2a_retained_two=" << two.triangles.retained_bytes
            << " slice_c2a_peak_one="
            << one.triangles.accounted_peak_temporary_bytes
            << " slice_c2a_peak_two="
            << two.triangles.accounted_peak_temporary_bytes
            << " slice_c2a_retained_normalized=" << retained_normalized_ratio
            << " slice_c2a_peak_normalized=" << peak_normalized_ratio
            << " slice_c2a_computed_zero="
            << unit_geometry.computed_zero_area_triangle_occurrences << '\n';
}

void check_c1_tetra_enumeration_and_axis_location() {
  static constexpr std::array<PatchGridLocalTetraEdge, 6> pairs = {
      PatchGridLocalTetraEdge{0U, 1U}, PatchGridLocalTetraEdge{0U, 2U},
      PatchGridLocalTetraEdge{0U, 3U}, PatchGridLocalTetraEdge{1U, 2U},
      PatchGridLocalTetraEdge{1U, 3U}, PatchGridLocalTetraEdge{2U, 3U}};
  for (std::uint8_t mask = 0U; mask < 16U; ++mask) {
    std::array<PatchGridLocalTetraEdge, 6> actual{};
    const std::size_t count =
        enumerate_patch_grid_tetra_sign_changing_edges(mask, actual);
    std::size_t expected = 0U;
    for (const PatchGridLocalTetraEdge pair : pairs) {
      const bool a = (mask & (std::uint8_t{1} << pair.first)) != 0U;
      const bool b = (mask & (std::uint8_t{1} << pair.second)) != 0U;
      if (a != b)
        require(actual[expected++] == pair,
                "tetra edge order changed from fixed pair order");
    }
    require(count == expected,
            "tetra sign-changing edge enumeration count is wrong");
  }
  std::array<PatchGridLocalTetraEdge, 6> invalid_edges{};
  require(enumerate_patch_grid_tetra_sign_changing_edges(
              0x10U, invalid_edges) == std::numeric_limits<std::size_t>::max(),
          "tetra high mask bits were not rejected");

  SphereFixture fixture =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 1U);
  PatchGridActiveCellRequest cell;
  cell.ix = grid_coordinate_index(fixture.spec.lower.x,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iy = grid_coordinate_index(fixture.spec.lower.y,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iz = grid_coordinate_index(fixture.spec.lower.z,
                                  fixture.spec.grid_spacing, 0.0);
  cell.inside_mask = 1U;
  const auto requests = exact_cell_vertex_requests(fixture.spec, cell);
  const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, requests);
  const auto lineage = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, vertices, std::span(&cell, 1U));
  require(
      lineage.status == "pass_with_uncertified_edges" &&
          lineage.cells.size() == 1U &&
          lineage.cell_tetra_offsets == std::vector<std::size_t>{0U, 6U} &&
          lineage.tetrahedra.size() == 6U &&
          lineage.tetra_edge_offsets.size() == 7U &&
          lineage.edges.size() == 18U &&
          lineage.edge_endpoint_patch_offsets.size() ==
              2U * lineage.edges.size() + 1U &&
          lineage.edge_relation_offsets.size() == lineage.edges.size() + 1U &&
          lineage.edge_location_ref_offsets.size() ==
              lineage.edges.size() + 1U &&
          std::all_of(lineage.edges.begin(), lineage.edges.end(),
                      [](const PatchGridSignChangingEdgeLineage &edge) {
                        return edge.production_interpolation_complete;
                      }) &&
          lineage.relation_oracle_binding ==
              PatchGridRelationOracleBindingStatus::built_once_from_bound_batch,
      "C1 canonical cell/tetra/edge CSR contract failed");
  const PatchGridSignChangingEdgeLineage *axis_edge =
      find_cell_corner_edge(lineage, 0U, 1U);
  require(axis_edge != nullptr && axis_edge->cartesian_axis_edge &&
              axis_edge->cartesian_axis == RayAxis::x &&
              axis_edge->production_interpolation_complete &&
              axis_edge->analytic_location_residual_available &&
              axis_edge->scalar_lineage_complete &&
              axis_edge->surface_location_certified &&
              axis_edge->location_status ==
                  PatchGridEdgeSurfaceLocationStatus::certified_axis_crossing &&
              axis_edge->interpolated_point.x == 1.0 &&
              axis_edge->interpolated_point.y == 0.0 &&
              axis_edge->interpolated_point.z == 0.0 &&
              axis_edge->interpolation_residual <=
                  axis_edge->interpolation_tolerance,
          "aligned sphere axis edge did not certify exact interpolation");
  const std::size_t axis_edge_position =
      static_cast<std::size_t>(axis_edge - lineage.edges.data());
  const std::size_t location_begin =
      lineage.edge_location_ref_offsets[axis_edge_position];
  const std::size_t location_end =
      lineage.edge_location_ref_offsets[axis_edge_position + 1U];
  const std::size_t relation_begin =
      lineage.edge_relation_offsets[axis_edge_position];
  const std::size_t relation_end =
      lineage.edge_relation_offsets[axis_edge_position + 1U];
  const std::span<const PatchRefId> exact_crossing_refs =
      lineage.location_patch_refs.empty()
          ? std::span<const PatchRefId>{}
          : fixture.on.provenance
                .axes[static_cast<std::size_t>(
                    axis_edge->location_crossing.axis())]
                .refs(axis_edge->location_crossing.crossing());
  require(
      axis_edge->location_crossing.valid() &&
          axis_edge->location_crossing.axis() == RayAxis::x &&
          location_end == location_begin + 1U &&
          exact_crossing_refs.size() == 1U &&
          lineage.location_patch_refs[location_begin] ==
              exact_crossing_refs.front() &&
          relation_end > relation_begin &&
          std::all_of(lineage.relations.begin() + relation_begin,
                      lineage.relations.begin() + relation_end,
                      [](const PatchGridEdgePatchRelationEvidence &relation) {
                        return relation.relation == PatchRelation::same_surface;
                      }),
      "axis certificate lost its exact ref or bound same-surface relation");
  const PatchGridSignChangingEdgeLineage *diagonal =
      find_cell_corner_edge(lineage, 0U, 5U);
  require(diagonal != nullptr && !diagonal->cartesian_axis_edge &&
              diagonal->production_interpolation_complete &&
              !diagonal->analytic_location_residual_available &&
              !diagonal->surface_location_certified &&
              diagonal->location_status ==
                  PatchGridEdgeSurfaceLocationStatus::
                      internal_diagonal_location_uncertified,
          "tetra diagonal was incorrectly location-certified");
  for (std::size_t tetra = 0U; tetra < lineage.tetrahedra.size(); ++tetra) {
    const std::size_t begin_edge = lineage.tetra_edge_offsets[tetra];
    const std::size_t end_edge = lineage.tetra_edge_offsets[tetra + 1U];
    for (std::size_t edge = begin_edge; edge < end_edge; ++edge)
      require(lineage.edges[edge].edge_ordinal == edge - begin_edge,
              "edge ordinal is not filtered-list local order");
  }
}

void check_c1_semantic_downgrades_and_fail_closed() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture coarse = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  PatchGridActiveCellRequest cell;
  cell.ix =
      grid_coordinate_index(coarse.spec.lower.x, coarse.spec.grid_spacing, 0.0);
  cell.iy =
      grid_coordinate_index(coarse.spec.lower.y, coarse.spec.grid_spacing, 0.0);
  cell.iz =
      grid_coordinate_index(coarse.spec.lower.z, coarse.spec.grid_spacing, 0.0);
  cell.inside_mask = 1U;
  const auto exact_requests = exact_cell_vertex_requests(coarse.spec, cell);
  const auto exact_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_requests);
  const auto exact_lineage = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span(&cell, 1U));
  const PatchGridSignChangingEdgeLineage *exact_edge =
      find_cell_corner_edge(exact_lineage, 0U, 1U);
  require(exact_edge != nullptr &&
              exact_edge->topology_sign_evidence_complete &&
              exact_edge->topology_status == PatchGridEdgeTopologyStatus::exact,
          "exact union topology evidence did not remain independent");

  auto sparse_requests = exact_requests;
  sparse_requests[0].sparse_magnitude_present = true;
  sparse_requests[0].sparse_magnitude = 0.5F;
  const auto sparse_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      sparse_requests);
  const auto sparse_lineage = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      sparse_vertices, std::span(&cell, 1U));
  const auto *sparse_edge = find_cell_corner_edge(sparse_lineage, 0U, 1U);
  require(sparse_edge != nullptr &&
              sparse_edge->production_interpolation_complete &&
              !sparse_edge->analytic_location_residual_available &&
              sparse_edge->endpoint_scalar_status[0] ==
                  PatchGridEdgeScalarLineageStatus::sparse_sample &&
              !sparse_edge->scalar_lineage_complete &&
              sparse_edge->location_status ==
                  PatchGridEdgeSurfaceLocationStatus::
                      endpoint_scalar_lineage_uncertified,
          "sparse endpoint downgrade was lost");

  const auto *fallback_edge = find_cell_corner_edge(exact_lineage, 0U, 6U);
  require(fallback_edge != nullptr &&
              (fallback_edge->endpoint_scalar_status[0] ==
                   PatchGridEdgeScalarLineageStatus::grid_spacing_fallback ||
               fallback_edge->endpoint_scalar_status[1] ==
                   PatchGridEdgeScalarLineageStatus::grid_spacing_fallback),
          "grid-spacing fallback endpoint was not explicit");

  SphereFixture fine = build_sphere_fixture(atom, 1.5, 0.5, 1U);
  PatchGridActiveCellRequest fine_cell;
  fine_cell.ix =
      grid_coordinate_index(fine.spec.lower.x, fine.spec.grid_spacing, 0.5);
  fine_cell.iy =
      grid_coordinate_index(fine.spec.lower.y, fine.spec.grid_spacing, 0.0);
  fine_cell.iz =
      grid_coordinate_index(fine.spec.lower.z, fine.spec.grid_spacing, 0.0);
  fine_cell.inside_mask = 1U;
  auto fine_requests = exact_cell_vertex_requests(fine.spec, fine_cell);
  const auto fine_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      fine_requests);
  const auto fine_lineage = build_patch_grid_cell_lineage(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      fine_vertices, std::span(&fine_cell, 1U));
  const auto *floor_edge = find_cell_corner_edge(fine_lineage, 0U, 1U);
  require(floor_edge != nullptr &&
              (floor_edge->endpoint_scalar_status[0] ==
                   PatchGridEdgeScalarLineageStatus::clamp_floor ||
               floor_edge->endpoint_scalar_status[1] ==
                   PatchGridEdgeScalarLineageStatus::clamp_floor),
          "clamp-floor endpoint downgrade was lost");
  const auto *same_patch_diagonal = find_cell_corner_edge(fine_lineage, 0U, 5U);
  require(same_patch_diagonal != nullptr &&
              !same_patch_diagonal->surface_location_certified &&
              same_patch_diagonal->location_status ==
                  PatchGridEdgeSurfaceLocationStatus::
                      internal_diagonal_location_uncertified,
          "same-surface diagonal was upgraded to a location certificate");
  const std::size_t diagonal_position =
      static_cast<std::size_t>(same_patch_diagonal - fine_lineage.edges.data());
  const std::size_t first_ref_begin =
      fine_lineage.edge_endpoint_patch_offsets[2U * diagonal_position];
  const std::size_t first_ref_end =
      fine_lineage.edge_endpoint_patch_offsets[2U * diagonal_position + 1U];
  const std::size_t second_ref_end =
      fine_lineage.edge_endpoint_patch_offsets[2U * diagonal_position + 2U];
  require(first_ref_end > first_ref_begin && second_ref_end > first_ref_end &&
              fine_lineage.endpoint_patch_refs[first_ref_begin] ==
                  fine_lineage.endpoint_patch_refs[first_ref_end],
          "diagonal negative fixture did not retain the same patch candidate");
  const std::size_t diagonal_relation_begin =
      fine_lineage.edge_relation_offsets[diagonal_position];
  const std::size_t diagonal_relation_end =
      fine_lineage.edge_relation_offsets[diagonal_position + 1U];
  require(
      same_patch_diagonal->patch_relation_classification_complete &&
          diagonal_relation_end > diagonal_relation_begin &&
          std::any_of(fine_lineage.relations.begin() + diagonal_relation_begin,
                      fine_lineage.relations.begin() + diagonal_relation_end,
                      [](const PatchGridEdgePatchRelationEvidence &relation) {
                        return relation.relation == PatchRelation::same_surface;
                      }) &&
          !same_patch_diagonal->surface_location_certified,
      "same-surface relation evidence authorized a diagonal location");

  const auto find_fine_vertex =
      [&](std::size_t linear_index) -> const PatchGridVertexMagnitudeWitness & {
    const auto found = std::find_if(
        fine_vertices.vertices.begin(), fine_vertices.vertices.end(),
        [&](const PatchGridVertexMagnitudeWitness &vertex) {
          return vertex.linear_index == linear_index;
        });
    require(found != fine_vertices.vertices.end(),
            "diagonal interpolation endpoint was not retained");
    return *found;
  };
  const PatchGridVertexMagnitudeWitness &diagonal_first =
      find_fine_vertex(same_patch_diagonal->first_vertex_linear_index);
  const PatchGridVertexMagnitudeWitness &diagonal_second =
      find_fine_vertex(same_patch_diagonal->second_vertex_linear_index);
  const Vec3 diagonal_first_point{
      fine.spec.lower.x +
          fine.spec.grid_spacing * static_cast<double>(diagonal_first.ix),
      fine.spec.lower.y +
          fine.spec.grid_spacing * static_cast<double>(diagonal_first.iy),
      fine.spec.lower.z +
          fine.spec.grid_spacing * static_cast<double>(diagonal_first.iz)};
  const Vec3 diagonal_second_point{
      fine.spec.lower.x +
          fine.spec.grid_spacing * static_cast<double>(diagonal_second.ix),
      fine.spec.lower.y +
          fine.spec.grid_spacing * static_cast<double>(diagonal_second.iy),
      fine.spec.lower.z +
          fine.spec.grid_spacing * static_cast<double>(diagonal_second.iz)};
  const bool diagonal_first_inside =
      (fine_cell.inside_mask &
       (std::uint8_t{1} << same_patch_diagonal->first_cell_corner)) != 0U;
  const bool diagonal_second_inside =
      (fine_cell.inside_mask &
       (std::uint8_t{1} << same_patch_diagonal->second_cell_corner)) != 0U;
  const double diagonal_fa = diagonal_first_inside
                                 ? -diagonal_first.clamped_magnitude
                                 : diagonal_first.clamped_magnitude;
  const double diagonal_fb = diagonal_second_inside
                                 ? -diagonal_second.clamped_magnitude
                                 : diagonal_second.clamped_magnitude;
  const double diagonal_denominator = diagonal_fa - diagonal_fb;
  const double diagonal_t =
      std::abs(diagonal_denominator) > 1.0e-300
          ? std::clamp(diagonal_fa / diagonal_denominator, 0.0, 1.0)
          : 0.5;
  const Vec3 diagonal_point{
      diagonal_first_point.x +
          diagonal_t * (diagonal_second_point.x - diagonal_first_point.x),
      diagonal_first_point.y +
          diagonal_t * (diagonal_second_point.y - diagonal_first_point.y),
      diagonal_first_point.z +
          diagonal_t * (diagonal_second_point.z - diagonal_first_point.z)};
  require(same_patch_diagonal->production_interpolation_complete &&
              !same_patch_diagonal->analytic_location_residual_available &&
              std::bit_cast<std::uint64_t>(
                  same_patch_diagonal->interpolation_parameter) ==
                  std::bit_cast<std::uint64_t>(diagonal_t) &&
              std::bit_cast<std::uint64_t>(
                  same_patch_diagonal->interpolated_point.x) ==
                  std::bit_cast<std::uint64_t>(diagonal_point.x) &&
              std::bit_cast<std::uint64_t>(
                  same_patch_diagonal->interpolated_point.y) ==
                  std::bit_cast<std::uint64_t>(diagonal_point.y) &&
              std::bit_cast<std::uint64_t>(
                  same_patch_diagonal->interpolated_point.z) ==
                  std::bit_cast<std::uint64_t>(diagonal_point.z),
          "diagonal did not retain bitwise production interpolation geometry");

  auto tolerance_requests = fine_requests;
  tolerance_requests[0].sparse_magnitude_present = true;
  tolerance_requests[0].sparse_magnitude =
      std::nextafter(0.5F, std::numeric_limits<float>::infinity());
  const auto tolerance_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      tolerance_requests);
  const auto tolerance_lineage = build_patch_grid_cell_lineage(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      tolerance_vertices, std::span(&fine_cell, 1U));
  const auto *tolerance_edge = find_cell_corner_edge(tolerance_lineage, 0U, 5U);
  require(tolerance_edge != nullptr &&
              tolerance_edge->endpoint_scalar_status[0] ==
                  PatchGridEdgeScalarLineageStatus::tolerance_tied_source,
          "float-quantized tolerance tie was not retained on the endpoint");

  const std::size_t center_position = static_cast<std::size_t>(
      std::find_if(exact_vertices.vertices.begin(),
                   exact_vertices.vertices.end(),
                   [&](const PatchGridVertexMagnitudeWitness &vertex) {
                     return vertex.ix == cell.ix && vertex.iy == cell.iy &&
                            vertex.iz == cell.iz;
                   }) -
      exact_vertices.vertices.begin());
  require(center_position < exact_vertices.vertices.size(),
          "center vertex was not retained");
  const auto center_candidates =
      vertex_candidates(exact_vertices, center_position);
  const auto target_candidate = std::find_if(
      center_candidates.begin(), center_candidates.end(),
      [&](const PatchGridMagnitudeSourceCandidate &candidate) {
        return candidate.kind ==
                   PatchGridMagnitudeSourceKind::analytic_crossing &&
               candidate.axis == RayAxis::x &&
               coarse.on.axes[0].coordinates[candidate.crossing.crossing()] >
                   0.0;
      });
  require(target_candidate != center_candidates.end(),
          "aligned target crossing was not found");
  const PatchCrossingAddress target = target_candidate->crossing;

  const auto require_mutated_status =
      [&](PreparedPatchCrossingWitnesses prepared,
          PatchGridEdgeScalarLineageStatus expected, const char *message) {
        const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
            *coarse.complex, coarse.index, coarse.on, coarse.spec, prepared,
            exact_requests);
        const auto lineage = build_patch_grid_cell_lineage(
            *coarse.complex, coarse.index, coarse.on, coarse.spec, prepared,
            vertices, std::span(&cell, 1U));
        const auto *edge = find_cell_corner_edge(lineage, 0U, 1U);
        require(edge != nullptr &&
                    (edge->endpoint_scalar_status[0] == expected ||
                     edge->endpoint_scalar_status[1] == expected),
                message);
      };
  const auto mutate_target = [&](auto mutate) {
    PreparedPatchCrossingWitnesses prepared = coarse.prepared;
    const auto found =
        std::find_if(prepared.witnesses.begin(), prepared.witnesses.end(),
                     [&](const CertifiedPatchCrossingWitness &witness) {
                       return witness.address == target;
                     });
    require(found != prepared.witnesses.end(),
            "target witness disappeared from mutation fixture");
    mutate(*found);
    return prepared;
  };
  PreparedPatchCrossingWitnesses missing = coarse.prepared;
  missing.witnesses.clear();
  missing.axis_offsets = {};
  require_mutated_status(
      missing, PatchGridEdgeScalarLineageStatus::missing_crossing_witness,
      "missing crossing reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.identity_status = CrossingIdentityCertificateStatus::invalid;
      }),
      PatchGridEdgeScalarLineageStatus::identity_uncertified,
      "identity reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.root_status =
            CrossingRootReproductionStatus::root_not_reproduced;
        witness.location_status =
            CrossingBoundedLocationStatus::root_unreproduced;
      }),
      PatchGridEdgeScalarLineageStatus::root_uncertified,
      "root reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.location_status =
            CrossingBoundedLocationStatus::patch_residual_failed;
      }),
      PatchGridEdgeScalarLineageStatus::bounded_location_uncertified,
      "bounded-location reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.transversality_status =
            CrossingTransversalityStatus::repeated_root_ambiguous;
      }),
      PatchGridEdgeScalarLineageStatus::repeated_root_ambiguous,
      "repeated-root reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.transversality_status =
            CrossingTransversalityStatus::primitive_simplicity_uncertified;
      }),
      PatchGridEdgeScalarLineageStatus::primitive_simplicity_uncertified,
      "primitive-simplicity reason was lost in C1");
  require_mutated_status(
      mutate_target([](CertifiedPatchCrossingWitness &witness) {
        witness.transversality_status =
            CrossingTransversalityStatus::location_uncertified;
      }),
      PatchGridEdgeScalarLineageStatus::transversality_uncertified,
      "transversality reason was lost in C1");

  SphereFixture duplicated;
  duplicated.complex = std::make_unique<PatchComplex>(
      build_patch_complex(std::vector<Atom>{atom}, 1.0, nullptr, false));
  duplicated.complex->contacts.push_back(duplicated.complex->contacts.front());
  duplicated.index = build_audit_projection_index(*duplicated.complex);
  duplicated.spec = sphere_grid(atom, 1.0, 2.0);
  TitanExecutor executor(1U);
  PatchRayBuildOptions options;
  options.retain_provenance = true;
  duplicated.on = build_patch_ray_grid(
      *duplicated.complex, duplicated.spec.lower, duplicated.spec.grid_spacing,
      duplicated.spec.nx, duplicated.spec.ny, duplicated.spec.nz, &executor,
      options);
  duplicated.prepared = build_certified_patch_crossing_witnesses(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec);
  const auto duplicated_requests =
      exact_cell_vertex_requests(duplicated.spec, cell);
  const auto duplicated_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
      duplicated.prepared, duplicated_requests);
  const auto duplicated_lineage = build_patch_grid_cell_lineage(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
      duplicated.prepared, duplicated_vertices, std::span(&cell, 1U));
  const auto *multi_ref_edge =
      find_cell_corner_edge(duplicated_lineage, 0U, 1U);
  require(multi_ref_edge != nullptr &&
              (multi_ref_edge->endpoint_scalar_status[0] ==
                   PatchGridEdgeScalarLineageStatus::
                       multi_ref_relation_uncertified ||
               multi_ref_edge->endpoint_scalar_status[1] ==
                   PatchGridEdgeScalarLineageStatus::
                       multi_ref_relation_uncertified),
          "multi-ref endpoint reason was lost in C1");

  struct TopologyCase {
    PatchGridTopologyState state;
    PatchGridTopologySource source;
    bool cavity;
    PatchGridTopologyEvidenceStatus endpoint;
    PatchGridEdgeTopologyStatus aggregate;
  };
  const std::array<TopologyCase, 5> topology_cases = {
      TopologyCase{PatchGridTopologyState::inside,
                   PatchGridTopologySource::direct_unanimous_parity, false,
                   PatchGridTopologyEvidenceStatus::
                       direct_unanimous_parity_not_yet_audited,
                   PatchGridEdgeTopologyStatus::direct_parity_not_yet_audited},
      TopologyCase{
          PatchGridTopologyState::inside,
          PatchGridTopologySource::analytic_local_repair, false,
          PatchGridTopologyEvidenceStatus::analytic_perturbed_point_parity,
          PatchGridEdgeTopologyStatus::analytic_local_repair_perturbed},
      TopologyCase{
          PatchGridTopologyState::inside,
          PatchGridTopologySource::morphology_fallback, false,
          PatchGridTopologyEvidenceStatus::morphology_fallback_nonanalytic,
          PatchGridEdgeTopologyStatus::morphology_fallback_nonanalytic},
      TopologyCase{PatchGridTopologyState::inside,
                   PatchGridTopologySource::cavity_fill, true,
                   PatchGridTopologyEvidenceStatus::cavity_fill_nonanalytic,
                   PatchGridEdgeTopologyStatus::cavity_fill_nonanalytic},
      TopologyCase{PatchGridTopologyState::unresolved,
                   PatchGridTopologySource::unknown, false,
                   PatchGridTopologyEvidenceStatus::unknown,
                   PatchGridEdgeTopologyStatus::topology_evidence_missing}};
  for (const TopologyCase &test : topology_cases) {
    auto requests = exact_requests;
    requests[0].topology_state = test.state;
    requests[0].topology_source = test.source;
    requests[0].cavity_fill = test.cavity;
    const auto vertices = build_patch_grid_vertex_magnitude_witnesses(
        *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
        requests);
    const auto lineage = build_patch_grid_cell_lineage(
        *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
        vertices, std::span(&cell, 1U));
    const auto *edge = find_cell_corner_edge(lineage, 0U, 1U);
    require(edge != nullptr &&
                edge->endpoint_topology_status[0] == test.endpoint &&
                edge->topology_status == test.aggregate &&
                !edge->topology_sign_evidence_complete,
            "non-exact topology evidence was upgraded or collapsed");
  }

  auto opposite_requests = exact_requests;
  opposite_requests[0].topology_state = PatchGridTopologyState::outside;
  opposite_requests[0].topology_source =
      PatchGridTopologySource::outside_probe_expanded_union;
  const auto opposite_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      opposite_requests);
  const auto opposite = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      opposite_vertices, std::span(&cell, 1U));
  require(opposite.status == "topology_sign_mismatch" &&
              opposite.cells.empty() && opposite.edges.empty(),
          "explicit topology state opposite to inside_mask was accepted");

  const auto empty = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span<const PatchGridActiveCellRequest>{});
  require(empty.status == "empty_active_cells" && empty.cells.empty() &&
              !empty.scalar_lineage_complete &&
              !empty.surface_location_certificate_complete &&
              !empty.topology_sign_certificate_complete,
          "empty C1 input became a vacuous certificate");
  const std::array<PatchGridActiveCellRequest, 2> duplicate_cells = {cell,
                                                                     cell};
  const auto duplicate_result = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, duplicate_cells);
  PatchGridActiveCellRequest inactive = cell;
  inactive.inside_mask = 0U;
  const auto inactive_result = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span(&inactive, 1U));
  PatchGridActiveCellRequest full = cell;
  full.inside_mask = 0xffU;
  const auto full_result = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span(&full, 1U));
  PatchGridActiveCellRequest out_of_range = cell;
  out_of_range.ix = coarse.spec.nx;
  const auto out_of_range_result = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span(&out_of_range, 1U));
  require(duplicate_result.status == "duplicate_active_cell" &&
              inactive_result.status == "invalid_active_cell" &&
              full_result.status == "invalid_active_cell" &&
              out_of_range_result.status == "invalid_active_cell" &&
              duplicate_result.cells.empty() && inactive_result.cells.empty() &&
              full_result.cells.empty() && out_of_range_result.cells.empty(),
          "invalid/duplicate active cell did not fail closed");

  std::array<PatchGridVertexMagnitudeRequest, 7> missing_requests{};
  std::copy_n(exact_requests.begin(), missing_requests.size(),
              missing_requests.begin());
  const auto missing_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      missing_requests);
  const auto missing_corner = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      missing_vertices, std::span(&cell, 1U));
  require(missing_corner.status == "missing_corner_vertex" &&
              missing_corner.cells.empty(),
          "missing cell corner did not fail closed");

  PreparedPatchGridVertexMagnitudeWitnesses stale_fingerprint = exact_vertices;
  ++stale_fingerprint.complex_fingerprint.first;
  PreparedPatchGridVertexMagnitudeWitnesses stale_lattice = exact_vertices;
  ++stale_lattice.lattice_identity.scalar_bits[0];
  PreparedPatchGridVertexMagnitudeWitnesses malformed_offsets = exact_vertices;
  ++malformed_offsets.candidate_offsets.back();
  PreparedPatchGridVertexMagnitudeWitnesses malformed_candidate =
      exact_vertices;
  require(!malformed_candidate.candidates.empty(),
          "malformed candidate fixture has no candidate");
  ++malformed_candidate.candidates.front().refs_end;
  PreparedPatchGridVertexMagnitudeWitnesses sparse_mutation = exact_vertices;
  sparse_mutation.vertices.front().sparse_magnitude = 1.0F;
  const auto require_invalid_vertex_payload =
      [&](const PreparedPatchGridVertexMagnitudeWitnesses &payload,
          const char *message) {
        const auto failed = build_patch_grid_cell_lineage(
            *coarse.complex, coarse.index, coarse.on, coarse.spec,
            coarse.prepared, payload, std::span(&cell, 1U));
        require(failed.status == "invalid_vertex_payload" &&
                    failed.cells.empty() && failed.edges.empty(),
                message);
      };
  require_invalid_vertex_payload(stale_fingerprint,
                                 "stale vertex fingerprint was accepted");
  require_invalid_vertex_payload(stale_lattice,
                                 "stale vertex lattice was accepted");
  require_invalid_vertex_payload(malformed_offsets,
                                 "malformed vertex CSR was accepted");
  require_invalid_vertex_payload(malformed_candidate,
                                 "malformed candidate ref range was accepted");
  require_invalid_vertex_payload(sparse_mutation,
                                 "exact sparse payload mutation was accepted");

  const auto require_invalid_crossing_payload =
      [&](PreparedPatchCrossingWitnesses prepared, const char *message) {
        const auto failed = build_patch_grid_cell_lineage(
            *coarse.complex, coarse.index, coarse.on, coarse.spec, prepared,
            exact_vertices, std::span(&cell, 1U));
        require(failed.status == "invalid_vertex_payload" &&
                    failed.cells.empty() && failed.edges.empty(),
                message);
      };
  PreparedPatchCrossingWitnesses stale_crossing_fingerprint = coarse.prepared;
  ++stale_crossing_fingerprint.complex_fingerprint.first;
  require_invalid_crossing_payload(stale_crossing_fingerprint,
                                   "stale Slice-A fingerprint was accepted");
  PreparedPatchCrossingWitnesses stale_crossing_lattice = coarse.prepared;
  ++stale_crossing_lattice.lattice_identity.scalar_bits[0];
  require_invalid_crossing_payload(stale_crossing_lattice,
                                   "stale Slice-A lattice was accepted");
  PreparedPatchCrossingWitnesses malformed_crossing_point = coarse.prepared;
  require(!malformed_crossing_point.witnesses.empty(),
          "malformed Slice-A point fixture has no witness");
  malformed_crossing_point.witnesses.front().point.x =
      std::nextafter(malformed_crossing_point.witnesses.front().point.x,
                     std::numeric_limits<double>::infinity());
  require_invalid_crossing_payload(malformed_crossing_point,
                                   "malformed Slice-A point was accepted");
  PreparedPatchCrossingWitnesses malformed_crossing_multiplicity =
      coarse.prepared;
  ++malformed_crossing_multiplicity.witnesses.front().total_root_multiplicity;
  require_invalid_crossing_payload(
      malformed_crossing_multiplicity,
      "malformed Slice-A ref/multiplicity payload was accepted");

  PatchGridActiveCellRequest overflow = cell;
  overflow.ix = std::numeric_limits<std::size_t>::max();
  const auto overflow_result = build_patch_grid_cell_lineage(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      exact_vertices, std::span(&overflow, 1U));
  require(overflow_result.status == "invalid_active_cell" &&
              overflow_result.cells.empty(),
          "overflow active cell did not fail closed");
}

void check_c1_determinism_and_linear_memory() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture serial = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  SphereFixture parallel = build_sphere_fixture(atom, 1.0, 2.0, 8U);
  PatchGridActiveCellRequest centered;
  centered.ix =
      grid_coordinate_index(serial.spec.lower.x, serial.spec.grid_spacing, 0.0);
  centered.iy =
      grid_coordinate_index(serial.spec.lower.y, serial.spec.grid_spacing, 0.0);
  centered.iz =
      grid_coordinate_index(serial.spec.lower.z, serial.spec.grid_spacing, 0.0);
  centered.inside_mask = 1U;
  const auto requests = exact_cell_vertex_requests(serial.spec, centered);
  const auto serial_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      requests);
  const auto parallel_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *parallel.complex, parallel.index, parallel.on, parallel.spec,
      parallel.prepared, requests);
  const auto serial_lineage = build_patch_grid_cell_lineage(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      serial_vertices, std::span(&centered, 1U));
  const auto parallel_lineage = build_patch_grid_cell_lineage(
      *parallel.complex, parallel.index, parallel.on, parallel.spec,
      parallel.prepared, parallel_vertices, std::span(&centered, 1U));
  require(
      cell_lineage_checksum(serial_lineage) ==
              cell_lineage_checksum(parallel_lineage) &&
          serial_lineage.relation_oracle_binding ==
              PatchGridRelationOracleBindingStatus::
                  built_once_from_bound_batch &&
          parallel_lineage.relation_oracle_binding ==
              PatchGridRelationOracleBindingStatus::built_once_from_bound_batch,
      "nt=1/8 C1 lineage is not deterministic or oracle-bound once");

  PatchGridActiveCellRequest first{1U, 1U, 1U, 1U};
  PatchGridActiveCellRequest second{4U, 4U, 4U, 1U};
  std::vector<PatchGridVertexMagnitudeRequest> two_cell_requests;
  const auto first_requests = exact_cell_vertex_requests(serial.spec, first);
  const auto second_requests = exact_cell_vertex_requests(serial.spec, second);
  two_cell_requests.insert(two_cell_requests.end(), first_requests.begin(),
                           first_requests.end());
  two_cell_requests.insert(two_cell_requests.end(), second_requests.begin(),
                           second_requests.end());
  const auto two_cell_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      two_cell_requests);
  const std::array<PatchGridActiveCellRequest, 2> sorted_cells = {first,
                                                                  second};
  const std::array<PatchGridActiveCellRequest, 2> unsorted_cells = {second,
                                                                    first};
  const auto sorted_lineage = build_patch_grid_cell_lineage(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      two_cell_vertices, sorted_cells);
  const auto unsorted_lineage = build_patch_grid_cell_lineage(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      two_cell_vertices, unsorted_cells);
  const auto one_lineage = build_patch_grid_cell_lineage(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      two_cell_vertices, std::span(&first, 1U));
  for (std::size_t warmup = 0U; warmup < 2U; ++warmup) {
    const auto warmed = build_patch_grid_cell_lineage(
        *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
        two_cell_vertices, sorted_cells);
    require(cell_lineage_checksum(warmed) ==
                cell_lineage_checksum(sorted_lineage),
            "C1 warmup changed deterministic lineage");
  }
  std::array<double, 7> c1_samples{};
  for (double &sample : c1_samples) {
    const auto measured = build_patch_grid_cell_lineage(
        *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
        two_cell_vertices, sorted_cells);
    require(cell_lineage_checksum(measured) ==
                cell_lineage_checksum(sorted_lineage),
            "C1 timed sample changed deterministic lineage");
    sample = measured.build_ms;
  }
  std::array<double, 7> sorted_samples = c1_samples;
  std::sort(sorted_samples.begin(), sorted_samples.end());
  const double normalized_ratio =
      (static_cast<double>(sorted_lineage.retained_bytes) / 2.0) /
      static_cast<double>(one_lineage.retained_bytes);
  require(cell_lineage_checksum(sorted_lineage) ==
                  cell_lineage_checksum(unsorted_lineage) &&
              sorted_lineage.cells.size() == 2U &&
              sorted_lineage.cells[0].address <
                  sorted_lineage.cells[1].address &&
              one_lineage.retained_bytes > 0U &&
              sorted_lineage.retained_bytes > one_lineage.retained_bytes &&
              normalized_ratio >= 0.5 && normalized_ratio <= 1.5 &&
              sorted_lineage.accounted_peak_temporary_bytes >=
                  sorted_lineage.retained_bytes &&
              sorted_lineage.oracle_accounted_setup_bytes > 0U,
          "C1 canonical ordering or retained memory is not linear/bounded");
  std::cout << "slice_c1_checksum=" << cell_lineage_checksum(serial_lineage)
            << " slice_c1_edges=" << serial_lineage.edges.size()
            << " slice_c1_location_certified="
            << serial_lineage.surface_location_certified_edges
            << " slice_c1_retained_one=" << one_lineage.retained_bytes
            << " slice_c1_retained_two=" << sorted_lineage.retained_bytes
            << " slice_c1_normalized_ratio=" << normalized_ratio
            << " slice_c1_oracle_bytes="
            << sorted_lineage.oracle_accounted_setup_bytes
            << " slice_c1_samples_ms=" << c1_samples[0] << ',' << c1_samples[1]
            << ',' << c1_samples[2] << ',' << c1_samples[3] << ','
            << c1_samples[4] << ',' << c1_samples[5] << ',' << c1_samples[6]
            << " slice_c1_median_ms=" << sorted_samples[3] << '\n';
}

void require_scalar_parity(
    const PatchRayGrid &grid, const PatchMappingGridSpec &spec,
    std::span<const PatchGridVertexMagnitudeRequest> requests,
    const PreparedPatchGridVertexMagnitudeWitnesses &prepared) {
  require(prepared.vertices.size() == requests.size(),
          "magnitude batch lost a requested vertex");
  for (const PatchGridVertexMagnitudeWitness &vertex : prepared.vertices) {
    const auto found = std::find_if(
        requests.begin(), requests.end(),
        [&](const PatchGridVertexMagnitudeRequest &request) {
          std::size_t linear = 0U;
          return try_patch_grid_vertex_linear_index(
                     spec, request.ix, request.iy, request.iz, linear) &&
                 linear == vertex.linear_index;
        });
    require(found != requests.end(), "sorted vertex lost its source request");
    const double raw = reference_production_raw_magnitude(grid, spec, *found);
    const double clamped = std::max(raw, 1.0e-9 * spec.grid_spacing);
    require(std::bit_cast<std::uint64_t>(raw) ==
                    std::bit_cast<std::uint64_t>(vertex.raw_magnitude) &&
                std::bit_cast<std::uint64_t>(clamped) ==
                    std::bit_cast<std::uint64_t>(vertex.clamped_magnitude),
            "Slice B changed the production magnitude scalar");
  }
}

void check_topology_source_classification() {
  const auto inside = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::inside,
      PatchGridTopologySource::inside_atom_union, false);
  const auto outside = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::outside,
      PatchGridTopologySource::outside_probe_expanded_union, false);
  const auto direct = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::inside,
      PatchGridTopologySource::direct_unanimous_parity, false);
  const auto repaired = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::inside,
      PatchGridTopologySource::analytic_local_repair, false);
  const auto morphology = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::inside,
      PatchGridTopologySource::morphology_fallback, false);
  const auto cavity = classify_patch_grid_topology_evidence(
      PatchGridTopologyState::inside, PatchGridTopologySource::cavity_fill,
      true);
  require(inside.exact_vertex_certified && inside.analytic_evidence &&
              outside.exact_vertex_certified && outside.analytic_evidence &&
              direct.status == PatchGridTopologyEvidenceStatus::
                                   direct_unanimous_parity_not_yet_audited &&
              !direct.analytic_evidence && !direct.exact_vertex_certified &&
              repaired.status == PatchGridTopologyEvidenceStatus::
                                     analytic_perturbed_point_parity &&
              repaired.analytic_evidence && !repaired.exact_vertex_certified &&
              morphology.status == PatchGridTopologyEvidenceStatus::
                                       morphology_fallback_nonanalytic &&
              !morphology.analytic_evidence &&
              cavity.status ==
                  PatchGridTopologyEvidenceStatus::cavity_fill_nonanalytic &&
              !cavity.analytic_evidence,
          "topology source classification over-certified a branch");
}

void check_magnitude_scalar_and_source_selection() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture fine = build_sphere_fixture(atom, 1.5, 0.5, 1U);
  require(fine.prepared.status == "pass",
          "fine aligned sphere fixture did not certify crossings");

  PatchGridVertexMagnitudeRequest sparse_x_wins =
      make_vertex_request(fine.spec, 0.5, 0.0, 0.0);
  sparse_x_wins.sparse_magnitude_present = true;
  sparse_x_wins.sparse_magnitude = 0.75F;
  const auto x_wins = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      sparse_x_wins);
  require(
      x_wins.status == "pass" && x_wins.vertices.size() == 1U &&
          x_wins.vertices[0].raw_magnitude == 0.5 &&
          std::bit_cast<std::uint32_t>(x_wins.vertices[0].sparse_magnitude) ==
              std::bit_cast<std::uint32_t>(0.75F) &&
          x_wins.vertices[0].sparse_suppressed_yz &&
          x_wins.raw_source_certificate_complete &&
          x_wins.final_scalar_source_certificate_complete &&
          count_candidate_kind(
              vertex_candidates(x_wins, 0U),
              PatchGridMagnitudeSourceKind::analytic_crossing) == 1U &&
          x_wins.sparse_candidates == 0U &&
          !x_wins.vertices[0].surface_location_certified &&
          !x_wins.surface_location_certified,
      "X did not beat a present sparse magnitude exactly");
  require_scalar_parity(fine.on, fine.spec, std::span(&sparse_x_wins, 1U),
                        x_wins);

  PatchGridVertexMagnitudeRequest sparse_suppresses_yz =
      make_vertex_request(fine.spec, 0.0, 0.5, 0.0);
  sparse_suppresses_yz.sparse_magnitude_present = true;
  sparse_suppresses_yz.sparse_magnitude = 0.75F;
  const auto yz_suppressed = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      sparse_suppresses_yz);
  require(yz_suppressed.vertices[0].raw_magnitude == 0.75 &&
              yz_suppressed.candidates.size() == 1U &&
              yz_suppressed.candidates[0].kind ==
                  PatchGridMagnitudeSourceKind::sparse_sample &&
              yz_suppressed.status == "no_usable_analytic_evidence" &&
              yz_suppressed.analytic_candidates == 0U,
          "Y/Z were not suppressed by a present sparse magnitude");
  require_scalar_parity(fine.on, fine.spec,
                        std::span(&sparse_suppresses_yz, 1U), yz_suppressed);

  SphereFixture xyz = build_sphere_fixture(atom, 1.25, 0.75, 1U);
  PatchGridVertexMagnitudeRequest y_wins =
      make_vertex_request(xyz.spec, 0.0, 0.75, 0.0);
  const auto xyz_result = select_patch_grid_vertex_magnitude_sources(
      *xyz.complex, xyz.index, xyz.on, xyz.spec, xyz.prepared, y_wins);
  require(xyz_result.status == "pass" &&
              xyz_result.vertices[0].raw_magnitude == 0.25 &&
              xyz_result.candidates.size() == 1U &&
              xyz_result.candidates[0].kind ==
                  PatchGridMagnitudeSourceKind::analytic_crossing &&
              xyz_result.candidates[0].axis == RayAxis::y,
          "sparse-absent selector did not consider all X/Y/Z rows");
  require_scalar_parity(xyz.on, xyz.spec, std::span(&y_wins, 1U), xyz_result);

  PatchGridVertexMagnitudeRequest fallback = make_vertex_request(
      fine.spec, fine.spec.lower.x, fine.spec.lower.y, fine.spec.lower.z,
      PatchGridTopologyState::outside,
      PatchGridTopologySource::outside_probe_expanded_union);
  const auto fallback_result = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared, fallback);
  require(fallback_result.vertices[0].raw_magnitude == fine.spec.grid_spacing &&
              fallback_result.candidates.size() == 1U &&
              fallback_result.candidates[0].kind ==
                  PatchGridMagnitudeSourceKind::grid_spacing_fallback &&
              fallback_result.status == "no_usable_analytic_evidence",
          "implicit grid-spacing fallback was not retained explicitly");

  PatchGridVertexMagnitudeRequest exact_sparse_tie = sparse_x_wins;
  exact_sparse_tie.sparse_magnitude = 0.5F;
  const auto exact_tie = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      exact_sparse_tie);
  require(exact_tie.candidates.size() == 2U &&
              exact_tie.sparse_candidates == 1U &&
              exact_tie.analytic_candidates == 1U &&
              exact_tie.raw_exact_owner_candidates == 2U &&
              exact_tie.final_scalar_exact_owner_candidates == 2U &&
              !exact_tie.raw_source_certificate_complete &&
              exact_tie.status == "raw_source_uncertified",
          "exact sparse/X tie was not retained conservatively");

  PatchGridVertexMagnitudeRequest quantized_sparse_tie = sparse_x_wins;
  quantized_sparse_tie.sparse_magnitude =
      std::nextafter(0.5F, std::numeric_limits<float>::infinity());
  const auto quantized_tie = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      quantized_sparse_tie);
  require(quantized_tie.vertices[0].raw_magnitude == 0.5 &&
              quantized_tie.vertices[0].sparse_float_tie_tolerance > 0.0 &&
              quantized_tie.sparse_candidates == 1U &&
              quantized_tie.tolerance_tied_candidates >= 1U &&
              !quantized_tie.raw_source_certificate_complete,
          "float quantization band lost a plausible sparse/X tie");
  const auto quantized_candidates = vertex_candidates(quantized_tie, 0U);
  const auto quantized_sparse = std::find_if(
      quantized_candidates.begin(), quantized_candidates.end(),
      [](const PatchGridMagnitudeSourceCandidate &candidate) {
        return candidate.kind == PatchGridMagnitudeSourceKind::sparse_sample;
      });
  require(
      quantized_sparse != quantized_candidates.end() &&
          quantized_sparse->raw_tolerance_tied &&
          !quantized_sparse->raw_exact_owner &&
          !quantized_sparse->final_scalar_exact_owner &&
          quantized_tie.raw_exact_owner_candidates == 1U &&
          quantized_tie.final_scalar_exact_owner_candidates == 1U &&
          std::count_if(
              quantized_candidates.begin(), quantized_candidates.end(),
              [](const PatchGridMagnitudeSourceCandidate &candidate) {
                return candidate.kind ==
                           PatchGridMagnitudeSourceKind::analytic_crossing &&
                       candidate.raw_exact_owner &&
                       candidate.final_scalar_exact_owner;
              }) == 1,
      "tolerance-only sparse source was overstated as an exact owner");

  PatchGridVertexMagnitudeRequest floor_request =
      make_vertex_request(fine.spec, 1.0, 0.0, 0.0);
  const auto floor_result = select_patch_grid_vertex_magnitude_sources(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      floor_request);
  require(floor_result.vertices[0].raw_magnitude == 0.0 &&
              std::bit_cast<std::uint32_t>(
                  floor_result.vertices[0].sparse_magnitude) ==
                  std::bit_cast<std::uint32_t>(0.0F) &&
              floor_result.vertices[0].clamped_magnitude ==
                  1.0e-9 * fine.spec.grid_spacing &&
              floor_result.vertices[0].floor_clamp_applied &&
              floor_result.raw_source_certificate_complete &&
              !floor_result.final_scalar_source_certificate_complete &&
              floor_result.clamp_floor_candidates == 1U &&
              floor_result.raw_exact_owner_candidates == 1U &&
              floor_result.final_scalar_exact_owner_candidates == 1U &&
              count_candidate_kind(vertex_candidates(floor_result, 0U),
                                   PatchGridMagnitudeSourceKind::clamp_floor) ==
                  1U &&
              floor_result.status == "final_scalar_source_uncertified",
          "floor clamp was not represented as a nonanalytic final owner");
  require_scalar_parity(fine.on, fine.spec, std::span(&floor_request, 1U),
                        floor_result);

  const auto empty = build_patch_grid_vertex_magnitude_witnesses(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared,
      std::span<const PatchGridVertexMagnitudeRequest>{});
  require(empty.status == "no_vertices" && empty.vertices.empty() &&
              empty.candidate_offsets.empty() && empty.candidates.empty() &&
              empty.candidate_refs.empty() &&
              !empty.raw_source_certificate_complete &&
              !empty.final_scalar_source_certificate_complete &&
              !empty.topology_source_classification_complete &&
              !empty.topology_exact_vertex_certificate_complete &&
              !empty.surface_location_certified,
          "empty magnitude request became a vacuous certificate");

  std::array<PatchGridVertexMagnitudeRequest, 3> unsorted = {
      floor_request, sparse_x_wins, fallback};
  const auto sorted = build_patch_grid_vertex_magnitude_witnesses(
      *fine.complex, fine.index, fine.on, fine.spec, fine.prepared, unsorted);
  require(sorted.vertices.size() == unsorted.size() &&
              sorted.candidate_offsets.size() == unsorted.size() + 1U &&
              sorted.candidate_offsets.front() == 0U &&
              sorted.candidate_offsets.back() == sorted.candidates.size() &&
              sorted.vertices[0].linear_index <
                  sorted.vertices[1].linear_index &&
              sorted.vertices[1].linear_index < sorted.vertices[2].linear_index,
          "unsorted vertex batch did not produce canonical CSR order");
  require_scalar_parity(fine.on, fine.spec, unsorted, sorted);
}

void check_tied_crossings_cavity_and_full_refs() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture coarse = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  PatchGridVertexMagnitudeRequest center =
      make_vertex_request(coarse.spec, 0.0, 0.0, 0.0);
  const auto tied = select_patch_grid_vertex_magnitude_sources(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      center);
  require(tied.status == "pass" && tied.candidates.size() == 6U &&
              tied.candidate_refs.size() == 6U &&
              tied.certified_analytic_candidates == 6U &&
              tied.fallback_candidates == 0U &&
              tied.raw_exact_owner_candidates == 6U &&
              tied.raw_source_certificate_complete &&
              tied.final_scalar_source_certificate_complete,
          "within-row/cross-axis ties did not retain all analytic sources");
  for (const PatchGridMagnitudeSourceCandidate &candidate : tied.candidates) {
    require(candidate.refs_end - candidate.refs_begin == 1U &&
                candidate.source_evidence_certified() &&
                candidate.raw_exact_owner && candidate.final_scalar_exact_owner,
            "tied crossing lost its complete PatchRef group");
  }

  PatchGridVertexMagnitudeRequest cavity = center;
  cavity.cavity_fill = true;
  cavity.topology_source = PatchGridTopologySource::cavity_fill;
  const auto cavity_result = select_patch_grid_vertex_magnitude_sources(
      *coarse.complex, coarse.index, coarse.on, coarse.spec, coarse.prepared,
      cavity);
  require(cavity_result.analytic_candidates == 6U &&
              cavity_result.raw_source_certificate_complete &&
              cavity_result.topology_nonanalytic_vertices == 1U &&
              !cavity_result.topology_exact_vertex_certificate_complete &&
              cavity_result.status == "topology_exact_vertex_uncertified",
          "cavity-fill without sparse magnitude did not use X/Y/Z exactly");

  PatchGridVertexMagnitudeRequest malformed_cavity = cavity;
  malformed_cavity.sparse_magnitude_present = true;
  malformed_cavity.sparse_magnitude = 1.0F;
  const auto malformed_cavity_result =
      select_patch_grid_vertex_magnitude_sources(
          *coarse.complex, coarse.index, coarse.on, coarse.spec,
          coarse.prepared, malformed_cavity);
  require(malformed_cavity_result.status == "invalid_vertex_request" &&
              malformed_cavity_result.vertices.empty(),
          "cavity-fill plus sparse magnitude did not fail closed");

  SphereFixture duplicated;
  duplicated.complex = std::make_unique<PatchComplex>(
      build_patch_complex(std::vector<Atom>{atom}, 1.0, nullptr, false));
  duplicated.complex->contacts.push_back(duplicated.complex->contacts.front());
  duplicated.index = build_audit_projection_index(*duplicated.complex);
  duplicated.spec = sphere_grid(atom, 1.0, 2.0);
  TitanExecutor executor(1U);
  PatchRayBuildOptions options;
  options.retain_provenance = true;
  duplicated.on = build_patch_ray_grid(
      *duplicated.complex, duplicated.spec.lower, duplicated.spec.grid_spacing,
      duplicated.spec.nx, duplicated.spec.ny, duplicated.spec.nz, &executor,
      options);
  duplicated.prepared = build_certified_patch_crossing_witnesses(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec);
  require(duplicated.index.bounded_domains_certified &&
              duplicated.on.provenance.certified &&
              duplicated.on.provenance.seam_groups > 0U &&
              duplicated.prepared.status == "multi_ref_relation_uncertified",
          "duplicated-contact multi-ref fixture failed");
  const auto multi_ref = select_patch_grid_vertex_magnitude_sources(
      *duplicated.complex, duplicated.index, duplicated.on, duplicated.spec,
      duplicated.prepared, center);
  require(multi_ref.analytic_candidates == 6U &&
              multi_ref.uncertified_analytic_candidates == 6U &&
              multi_ref.candidate_refs.size() == 12U &&
              !multi_ref.raw_source_certificate_complete &&
              multi_ref.status == "no_usable_analytic_evidence",
          "multi-ref groups were truncated or incorrectly upgraded");
  for (const PatchGridMagnitudeSourceCandidate &candidate :
       multi_ref.candidates) {
    require(candidate.refs_end - candidate.refs_begin == 2U &&
                candidate.evidence_status == PatchGridMagnitudeEvidenceStatus::
                                                 multi_ref_relation_uncertified,
            "multi-ref candidate did not retain both PatchRefIds");
  }
}

void check_magnitude_fail_closed_and_explicit_uncertainty() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture fixture = build_sphere_fixture(atom, 1.5, 0.5, 1U);
  PatchGridVertexMagnitudeRequest request =
      make_vertex_request(fixture.spec, 0.5, 0.0, 0.0);
  request.sparse_magnitude_present = true;
  request.sparse_magnitude = 0.75F;
  const auto certified = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, request);
  require(certified.candidates.size() == 1U &&
              certified.candidates[0].source_evidence_certified(),
          "uncertainty fixture lacks a certified crossing");
  const PatchCrossingAddress target = certified.candidates[0].crossing;

  PreparedPatchCrossingWitnesses missing = fixture.prepared;
  missing.witnesses.clear();
  missing.axis_offsets = {};
  const auto missing_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, missing,
      request);
  require(missing_result.candidates[0].evidence_status ==
                  PatchGridMagnitudeEvidenceStatus::missing_crossing_witness &&
              !missing_result.raw_source_certificate_complete,
          "missing prepared crossing was not retained explicitly");

  const auto mutate_target = [&](CrossingRootReproductionStatus root,
                                 CrossingBoundedLocationStatus location,
                                 CrossingTransversalityStatus transversality) {
    PreparedPatchCrossingWitnesses changed = fixture.prepared;
    const auto found =
        std::find_if(changed.witnesses.begin(), changed.witnesses.end(),
                     [&](const CertifiedPatchCrossingWitness &witness) {
                       return witness.address == target;
                     });
    require(found != changed.witnesses.end(),
            "target crossing witness was not found");
    found->root_status = root;
    found->location_status = location;
    found->transversality_status = transversality;
    return changed;
  };
  const PreparedPatchCrossingWitnesses root_bad =
      mutate_target(CrossingRootReproductionStatus::root_not_reproduced,
                    CrossingBoundedLocationStatus::root_unreproduced,
                    CrossingTransversalityStatus::location_uncertified);
  const PreparedPatchCrossingWitnesses location_bad =
      mutate_target(CrossingRootReproductionStatus::reproduced,
                    CrossingBoundedLocationStatus::patch_residual_failed,
                    CrossingTransversalityStatus::location_uncertified);
  const PreparedPatchCrossingWitnesses repeated =
      mutate_target(CrossingRootReproductionStatus::reproduced,
                    CrossingBoundedLocationStatus::certified,
                    CrossingTransversalityStatus::repeated_root_ambiguous);
  const PreparedPatchCrossingWitnesses non_simple = mutate_target(
      CrossingRootReproductionStatus::reproduced,
      CrossingBoundedLocationStatus::certified,
      CrossingTransversalityStatus::primitive_simplicity_uncertified);
  PreparedPatchCrossingWitnesses identity_bad = fixture.prepared;
  const auto identity_target =
      std::find_if(identity_bad.witnesses.begin(), identity_bad.witnesses.end(),
                   [&](const CertifiedPatchCrossingWitness &witness) {
                     return witness.address == target;
                   });
  require(identity_target != identity_bad.witnesses.end(),
          "identity fixture lost selected witness");
  identity_target->identity_status = CrossingIdentityCertificateStatus::invalid;
  const auto root_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, root_bad,
      request);
  const auto location_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, location_bad,
      request);
  const auto repeated_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, repeated,
      request);
  const auto non_simple_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, non_simple,
      request);
  const auto identity_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, identity_bad,
      request);
  require(
      root_result.candidates[0].evidence_status ==
              PatchGridMagnitudeEvidenceStatus::root_uncertified &&
          location_result.candidates[0].evidence_status ==
              PatchGridMagnitudeEvidenceStatus::bounded_location_uncertified &&
          repeated_result.candidates[0].evidence_status ==
              PatchGridMagnitudeEvidenceStatus::repeated_root_ambiguous &&
          non_simple_result.candidates[0].evidence_status ==
              PatchGridMagnitudeEvidenceStatus::
                  primitive_simplicity_uncertified &&
          identity_result.candidates[0].evidence_status ==
              PatchGridMagnitudeEvidenceStatus::identity_uncertified &&
          identity_result.status == "no_usable_analytic_evidence" &&
          !root_result.raw_source_certificate_complete &&
          !location_result.raw_source_certificate_complete &&
          !repeated_result.raw_source_certificate_complete &&
          !non_simple_result.raw_source_certificate_complete &&
          !identity_result.raw_source_certificate_complete,
      "uncertified crossing reason was lost or upgraded");

  const std::array<PatchGridVertexMagnitudeRequest, 2> duplicate = {request,
                                                                    request};
  const auto duplicate_result = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, duplicate);
  require(duplicate_result.status == "duplicate_vertex" &&
              duplicate_result.vertices.empty(),
          "duplicate vertex did not fail closed");

  PatchGridVertexMagnitudeRequest out_of_range = request;
  out_of_range.ix = fixture.spec.nx;
  PatchGridVertexMagnitudeRequest nonfinite_sparse = request;
  nonfinite_sparse.sparse_magnitude_present = true;
  nonfinite_sparse.sparse_magnitude = std::numeric_limits<float>::quiet_NaN();
  const auto out_of_range_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, out_of_range);
  const auto nonfinite_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec,
      fixture.prepared, nonfinite_sparse);
  require(out_of_range_result.status == "invalid_vertex_request" &&
              nonfinite_result.status == "invalid_vertex_request",
          "invalid vertex/sparse value did not fail closed");

  PreparedPatchCrossingWitnesses stale_prepared = fixture.prepared;
  ++stale_prepared.complex_fingerprint.first;
  PatchMappingGridSpec stale_spec = fixture.spec;
  stale_spec.lower.x = std::nextafter(stale_spec.lower.x,
                                      std::numeric_limits<double>::infinity());
  LocalPatchProjectionIndex stale_index = fixture.index;
  ++stale_index.source_complex_fingerprint.second;
  const auto stale_prepared_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, fixture.spec, stale_prepared,
      request);
  const auto stale_lattice_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, fixture.on, stale_spec, fixture.prepared,
      request);
  const auto stale_index_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, stale_index, fixture.on, fixture.spec, fixture.prepared,
      request);
  require(stale_prepared_result.status == "stale_or_invalid_binding" &&
              stale_lattice_result.status == "stale_or_invalid_binding" &&
              stale_index_result.status == "stale_or_invalid_binding",
          "stale prepared complex/lattice/index binding was accepted");

  const auto corrupt_target = [&](int kind) {
    PreparedPatchCrossingWitnesses malformed_prepared = fixture.prepared;
    const auto found =
        std::find_if(malformed_prepared.witnesses.begin(),
                     malformed_prepared.witnesses.end(),
                     [&](const CertifiedPatchCrossingWitness &witness) {
                       return witness.address == target;
                     });
    require(found != malformed_prepared.witnesses.end(),
            "malformed prepared fixture lost target");
    if (kind == 0) {
      found->point.x = std::nextafter(found->point.x,
                                      std::numeric_limits<double>::infinity());
    } else if (kind == 1) {
      ++found->ref_count;
    } else {
      ++found->total_root_multiplicity;
    }
    return select_patch_grid_vertex_magnitude_sources(
        *fixture.complex, fixture.index, fixture.on, fixture.spec,
        malformed_prepared, request);
  };
  require(corrupt_target(0).status == "stale_or_invalid_binding" &&
              corrupt_target(1).status == "stale_or_invalid_binding" &&
              corrupt_target(2).status == "stale_or_invalid_binding",
          "malformed prepared point/ref/multiplicity metadata was accepted");

  PatchRayGrid malformed = fixture.on;
  malformed.axes[0].offsets.front() = 1U;
  const auto malformed_result = select_patch_grid_vertex_magnitude_sources(
      *fixture.complex, fixture.index, malformed, fixture.spec,
      fixture.prepared, request);
  require(malformed_result.status == "stale_or_invalid_binding" &&
              malformed_result.vertices.empty(),
          "malformed crossing CSR did not fail closed");

  const PatchMappingGridSpec overflow_spec{
      {0.0, 0.0, 0.0}, 1.0, std::numeric_limits<std::size_t>::max(), 2U, 2U};
  std::size_t linear = 0U;
  require(
      !try_patch_grid_vertex_linear_index(overflow_spec, 0U, 0U, 0U, linear),
      "overflow vertex lattice produced a linear index");
}

void check_magnitude_determinism_translation_and_memory() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  SphereFixture serial = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  SphereFixture parallel = build_sphere_fixture(atom, 1.0, 2.0, 8U);
  const PatchGridVertexMagnitudeRequest center =
      make_vertex_request(serial.spec, 0.0, 0.0, 0.0);
  const auto serial_result = select_patch_grid_vertex_magnitude_sources(
      *serial.complex, serial.index, serial.on, serial.spec, serial.prepared,
      center);
  const auto parallel_result = select_patch_grid_vertex_magnitude_sources(
      *parallel.complex, parallel.index, parallel.on, parallel.spec,
      parallel.prepared, center);
  require(magnitude_witness_checksum(serial_result) ==
                  magnitude_witness_checksum(parallel_result) &&
              serial_result.candidate_refs == parallel_result.candidate_refs,
          "nt=1/8 magnitude witness selection is not deterministic");

  const Atom shifted_atom{1000000.0, -2000000.0, 3000000.0, 0.0, 1.0};
  SphereFixture shifted = build_sphere_fixture(shifted_atom, 1.0, 2.0, 8U);
  const PatchGridVertexMagnitudeRequest shifted_center = make_vertex_request(
      shifted.spec, shifted_atom.x, shifted_atom.y, shifted_atom.z);
  const auto shifted_result = select_patch_grid_vertex_magnitude_sources(
      *shifted.complex, shifted.index, shifted.on, shifted.spec,
      shifted.prepared, shifted_center);
  require(shifted_result.status == "pass" &&
              shifted_result.candidates.size() == 6U &&
              shifted_result.vertices[0].raw_magnitude == 1.0 &&
              shifted_result.vertices[0].source_tie_tolerance >
                  serial_result.vertices[0].source_tie_tolerance,
          "translated scale-aware tie selection failed");

  SphereFixture memory_fixture = build_sphere_fixture(atom, 1.5, 0.5, 1U);
  std::vector<PatchGridVertexMagnitudeRequest> all_requests;
  all_requests.reserve(memory_fixture.spec.nx * memory_fixture.spec.ny *
                       memory_fixture.spec.nz);
  for (std::size_t iz = 0U; iz < memory_fixture.spec.nz; ++iz) {
    for (std::size_t iy = 0U; iy < memory_fixture.spec.ny; ++iy) {
      for (std::size_t ix = 0U; ix < memory_fixture.spec.nx; ++ix) {
        PatchGridVertexMagnitudeRequest request;
        request.ix = ix;
        request.iy = iy;
        request.iz = iz;
        request.topology_state = PatchGridTopologyState::unresolved;
        request.topology_source = PatchGridTopologySource::unknown;
        all_requests.push_back(request);
      }
    }
  }
  const std::span<const PatchGridVertexMagnitudeRequest> half_requests(
      all_requests.data(), all_requests.size() / 2U);
  const auto half = build_patch_grid_vertex_magnitude_witnesses(
      *memory_fixture.complex, memory_fixture.index, memory_fixture.on,
      memory_fixture.spec, memory_fixture.prepared, half_requests);
  for (int warmup = 0; warmup < 2; ++warmup) {
    (void)build_patch_grid_vertex_magnitude_witnesses(
        *memory_fixture.complex, memory_fixture.index, memory_fixture.on,
        memory_fixture.spec, memory_fixture.prepared, all_requests);
  }
  std::array<double, 7> samples{};
  PreparedPatchGridVertexMagnitudeWitnesses full;
  for (double &sample : samples) {
    full = build_patch_grid_vertex_magnitude_witnesses(
        *memory_fixture.complex, memory_fixture.index, memory_fixture.on,
        memory_fixture.spec, memory_fixture.prepared, all_requests);
    sample = full.build_ms;
  }
  std::sort(samples.begin(), samples.end());
  const auto exact_payload = [](const auto &result) {
    return result.vertices.capacity() *
               sizeof(PatchGridVertexMagnitudeWitness) +
           result.candidate_offsets.capacity() * sizeof(std::size_t) +
           result.candidates.capacity() *
               sizeof(PatchGridMagnitudeSourceCandidate) +
           result.candidate_refs.capacity() * sizeof(PatchRefId);
  };
  const auto normalized = [](const auto &result) {
    const std::size_t units =
        result.vertices.size() + result.candidate_offsets.size() +
        result.candidates.size() + result.candidate_refs.size();
    return static_cast<double>(result.retained_bytes) /
           static_cast<double>(std::max<std::size_t>(units, 1U));
  };
  const double normalized_ratio = normalized(full) / normalized(half);
  require(!half.vertices.empty() &&
              full.vertices.size() == all_requests.size() &&
              half.retained_bytes == exact_payload(half) &&
              full.retained_bytes == exact_payload(full) &&
              half.peak_temporary_bytes >= half.retained_bytes &&
              full.peak_temporary_bytes >= full.retained_bytes &&
              normalized_ratio >= 0.5 && normalized_ratio <= 2.0,
          "magnitude witness storage is not exact/linearly bounded");
  std::cout << "slice_b_samples_ms=" << samples[0] << ',' << samples[1] << ','
            << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
            << samples[5] << ',' << samples[6]
            << " slice_b_median_ms=" << samples[3]
            << " slice_b_retained_half=" << half.retained_bytes
            << " slice_b_retained_full=" << full.retained_bytes
            << " slice_b_peak_full=" << full.peak_temporary_bytes
            << " slice_b_normalized_ratio=" << normalized_ratio << '\n';
}

} // namespace

int main() {
  try {
    check_row_decoding();
    check_overflow_guards();
    check_tolerance_boundaries();
    check_simple_sphere_and_determinism();
    check_malformed_and_repeated_fail_closed();
    check_fingerprint_and_bounded_projection_fail_closed();
    check_primitive_filtered_transversality_contract();
    check_translation_scale_tolerance();
    check_topology_source_classification();
    check_magnitude_scalar_and_source_selection();
    check_tied_crossings_cavity_and_full_refs();
    check_magnitude_fail_closed_and_explicit_uncertainty();
    check_magnitude_determinism_translation_and_memory();
    check_c2a_key_and_triangle_oracle();
    check_c1_tetra_enumeration_and_axis_location();
    check_c1_semantic_downgrades_and_fail_closed();
    check_c1_determinism_and_linear_memory();
    check_c2a_runtime_happy_path();
    check_c2a_public_builder_fail_closed();
    check_c2a_all_axis_candidate_boundary();
    check_c2a_determinism_timing_and_memory();
    const SphereFixture summary =
        build_sphere_fixture(Atom{0.25, -0.5, 0.75, 0.0, 1.2}, 1.4, 0.5, 8U);
    SphereFixture magnitude_summary =
        build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 8U);
    const PatchGridVertexMagnitudeRequest center =
        make_vertex_request(magnitude_summary.spec, 0.0, 0.0, 0.0);
    const auto selected = select_patch_grid_vertex_magnitude_sources(
        *magnitude_summary.complex, magnitude_summary.index,
        magnitude_summary.on, magnitude_summary.spec,
        magnitude_summary.prepared, center);
    std::cout << "tabi_ses_patch_mapping_smoke=pass crossings="
              << summary.prepared.witnesses.size()
              << " refs=" << summary.prepared.refs_audited
              << " checksum=" << witness_checksum(summary.prepared)
              << " vertex_candidates=" << selected.candidates.size()
              << " vertex_checksum=" << magnitude_witness_checksum(selected)
              << " vertex_retained_bytes=" << selected.retained_bytes
              << " vertex_peak_bytes=" << selected.peak_temporary_bytes
              << " retained_bytes=" << summary.prepared.retained_bytes
              << " build_ms=" << summary.prepared.build_ms << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_mapping_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
