#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_mapping.h"
#include "tabi/tabi_ses_patch_projection.h"
#include "tabi/tabi_ses_patch_ray.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
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
  const std::size_t count =
      static_cast<std::size_t>(std::ceil(2.0 * half_extent / spacing)) + 1U;
  return {{atom.x - half_extent, atom.y - half_extent, atom.z - half_extent},
          spacing,
          count,
          count,
          count};
}

std::size_t grid_coordinate_index(double lower, double spacing,
                                  double coordinate) {
  const double scaled = (coordinate - lower) / spacing;
  const auto rounded = static_cast<std::size_t>(std::llround(scaled));
  require(
      std::abs(lower + spacing * static_cast<double>(rounded) - coordinate) <=
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(coordinate)),
      "coordinate is not on projection fixture lattice");
  return rounded;
}

struct ProjectionFixture {
  std::unique_ptr<PatchComplex> complex;
  LocalPatchProjectionIndex index;
  PatchMappingGridSpec spec;
  PatchRayGrid grid;
  PreparedPatchCrossingWitnesses crossings;
};

ProjectionFixture build_sphere_fixture(const Atom &atom, double probe_radius,
                                       double spacing, std::size_t threads) {
  ProjectionFixture fixture;
  fixture.complex = std::make_unique<PatchComplex>(build_patch_complex(
      std::vector<Atom>{atom}, probe_radius, nullptr, false));
  require(fixture.complex->status == "pass" &&
              fixture.complex->contacts.size() == 1U,
          "projection sphere patch fixture failed");
  LocalPatchProjectionIndexBuildOptions index_options;
  index_options.retain_complex_fingerprint = true;
  fixture.index =
      build_local_patch_projection_index(*fixture.complex, index_options);
  require(fixture.index.bounded_domains_certified &&
              fixture.index.bvh_invariants_certified &&
              fixture.index.source_complex_fingerprint.valid,
          "projection sphere index failed");
  fixture.spec = sphere_grid(atom, probe_radius, spacing);
  PatchRayBuildOptions ray_options;
  ray_options.retain_provenance = true;
  TitanExecutor executor(threads);
  fixture.grid = build_patch_ray_grid(*fixture.complex, fixture.spec.lower,
                                      fixture.spec.grid_spacing,
                                      fixture.spec.nx, fixture.spec.ny,
                                      fixture.spec.nz, &executor, ray_options);
  require(fixture.grid.status == "pass" && fixture.grid.provenance.certified,
          "projection sphere ray grid failed");
  fixture.crossings = build_certified_patch_crossing_witnesses(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec);
  require(fixture.crossings.complex_fingerprint.valid &&
              fixture.crossings.lattice_identity.valid &&
              !fixture.crossings.witnesses.empty(),
          "projection sphere crossing witness binding failed");
  return fixture;
}

PatchComplex make_torus_reentrant_complex() {
  const double root3 = std::sqrt(3.0);
  const std::vector<Atom> atoms = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                   {0.0, root3, 0.0, 0.0, 1.0}};
  PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  require(complex.status == "pass" && !complex.contacts.empty() &&
              !complex.tori.empty() && !complex.reentrant.empty(),
          "production torus/reentrant complex was not closed");
  return complex;
}

PatchMappingGridSpec enclosing_complex_grid(const PatchComplex &complex,
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

PatchMappingGridSpec seam_aligned_complex_grid(const PatchComplex &complex,
                                               double spacing) {
  PatchMappingGridSpec spec = enclosing_complex_grid(complex, spacing);
  const AnalyticPatchSeamCertificate seams =
      build_analytic_patch_seam_certificate(complex);
  const auto edge =
      std::find_if(seams.reentrant_edges.begin(), seams.reentrant_edges.end(),
                   [](const ReentrantEdgeSeamCertificate &candidate) {
                     return candidate.relation == PatchRelation::curve_adjacent;
                   });
  require(seams.all_seams_certified && edge != seams.reentrant_edges.end() &&
              edge->reentrant.valid(),
          "production complex has no certified torus-reentrant seam");
  const ReentrantPatch &reentrant = complex.reentrant[edge->reentrant.index()];
  const std::size_t first = edge->edge;
  const std::size_t second = (first + 1U) % 3U;
  Vec3 direction{reentrant.contact_directions[first].x +
                     reentrant.contact_directions[second].x,
                 reentrant.contact_directions[first].y +
                     reentrant.contact_directions[second].y,
                 reentrant.contact_directions[first].z +
                     reentrant.contact_directions[second].z};
  const double direction_norm =
      std::hypot(direction.x, direction.y, direction.z);
  require(std::isfinite(direction_norm) && direction_norm > 0.0,
          "certified seam midpoint direction is degenerate");
  direction.x /= direction_norm;
  direction.y /= direction_norm;
  direction.z /= direction_norm;
  const Vec3 seam_point{
      reentrant.probe_center.x + complex.probe_radius * direction.x,
      reentrant.probe_center.y + complex.probe_radius * direction.y,
      reentrant.probe_center.z + complex.probe_radius * direction.z};
  const std::array<double, 3> normal = {
      std::abs(direction.x), std::abs(direction.y), std::abs(direction.z)};
  const std::size_t ray_axis = static_cast<std::size_t>(
      std::max_element(normal.begin(), normal.end()) - normal.begin());
  std::array<double, 3> lower = {spec.lower.x, spec.lower.y, spec.lower.z};
  const std::array<double, 3> seam = {seam_point.x, seam_point.y, seam_point.z};
  std::array<std::size_t, 3> count = {spec.nx, spec.ny, spec.nz};
  for (std::size_t dimension = 0U; dimension < 3U; ++dimension) {
    if (dimension == ray_axis)
      continue;
    const double steps =
        std::ceil((seam[dimension] - lower[dimension]) / spacing) + 1.0;
    require(std::isfinite(steps) && steps > 0.0,
            "certified seam lies outside enclosing grid");
    const double aligned_lower = seam[dimension] - steps * spacing;
    const double old_upper =
        lower[dimension] + spacing * static_cast<double>(count[dimension] - 1U);
    lower[dimension] = aligned_lower;
    count[dimension] = static_cast<std::size_t>(
                           std::ceil((old_upper - aligned_lower) / spacing)) +
                       1U;
  }
  return {
      {lower[0], lower[1], lower[2]}, spacing, count[0], count[1], count[2]};
}

ProjectionFixture build_complex_fixture(PatchComplex complex,
                                        PatchMappingGridSpec spec,
                                        std::size_t threads) {
  ProjectionFixture fixture;
  fixture.complex = std::make_unique<PatchComplex>(std::move(complex));
  LocalPatchProjectionIndexBuildOptions index_options;
  index_options.retain_complex_fingerprint = true;
  fixture.index =
      build_local_patch_projection_index(*fixture.complex, index_options);
  require(fixture.index.bounded_domains_certified &&
              fixture.index.bvh_invariants_certified &&
              fixture.index.source_complex_fingerprint.valid,
          "projection seam index failed");
  fixture.spec = spec;
  PatchRayBuildOptions options;
  options.retain_provenance = true;
  TitanExecutor executor(threads);
  fixture.grid = build_patch_ray_grid(
      *fixture.complex, fixture.spec.lower, fixture.spec.grid_spacing,
      fixture.spec.nx, fixture.spec.ny, fixture.spec.nz, &executor, options);
  if (fixture.grid.status != "pass" || !fixture.grid.provenance.enabled ||
      !fixture.grid.provenance.certified ||
      fixture.grid.provenance.status != "pass") {
    throw std::runtime_error(
        "projection seam ray grid failed status=" + fixture.grid.status +
        " provenance=" + fixture.grid.provenance.status + " odd=" +
        std::to_string(fixture.grid.axes[0].odd_rows +
                       fixture.grid.axes[1].odd_rows +
                       fixture.grid.axes[2].odd_rows) +
        " unresolved=" +
        std::to_string(fixture.grid.unresolved_uncertified_rows));
  }
  fixture.crossings = build_certified_patch_crossing_witnesses(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec);
  require(fixture.crossings.complex_fingerprint.valid &&
              fixture.crossings.lattice_identity.valid,
          "projection seam crossing binding failed");
  return fixture;
}

std::vector<PatchGridVertexMagnitudeRequest>
cell_vertex_requests(const PatchMappingGridSpec &spec,
                     std::span<const PatchGridActiveCellRequest> cells) {
  struct OrderedRequest {
    std::size_t linear{0U};
    PatchGridVertexMagnitudeRequest request;
  };
  std::vector<OrderedRequest> ordered;
  ordered.reserve(cells.size() * 8U);
  for (const PatchGridActiveCellRequest &cell : cells) {
    for (std::size_t corner = 0U; corner < 8U; ++corner) {
      const auto &offset = patch_grid_production_corners[corner];
      PatchGridVertexMagnitudeRequest request;
      request.ix = cell.ix + offset[0];
      request.iy = cell.iy + offset[1];
      request.iz = cell.iz + offset[2];
      const bool inside =
          (cell.inside_mask & (std::uint8_t{1} << corner)) != 0U;
      request.topology_state = inside ? PatchGridTopologyState::inside
                                      : PatchGridTopologyState::outside;
      request.topology_source =
          inside ? PatchGridTopologySource::inside_atom_union
                 : PatchGridTopologySource::outside_probe_expanded_union;
      std::size_t linear = 0U;
      require(try_patch_grid_vertex_linear_index(spec, request.ix, request.iy,
                                                 request.iz, linear),
              "projection cell corner is outside lattice");
      ordered.push_back({linear, request});
    }
  }
  std::sort(ordered.begin(), ordered.end(),
            [](const OrderedRequest &left, const OrderedRequest &right) {
              return left.linear < right.linear;
            });
  std::vector<PatchGridVertexMagnitudeRequest> result;
  result.reserve(ordered.size());
  for (const OrderedRequest &item : ordered) {
    if (!result.empty()) {
      std::size_t previous = 0U;
      require(try_patch_grid_vertex_linear_index(spec, result.back().ix,
                                                 result.back().iy,
                                                 result.back().iz, previous),
              "previous projection request is invalid");
      if (previous == item.linear) {
        require(result.back().topology_state == item.request.topology_state &&
                    result.back().topology_source ==
                        item.request.topology_source,
                "overlapping projection cells disagree on topology state");
        continue;
      }
    }
    result.push_back(item.request);
  }
  return result;
}

struct ProjectionPipeline {
  PreparedPatchGridVertexMagnitudeWitnesses vertices;
  PreparedPatchGridCellLineage cells;
  PreparedPatchGridTriangleOccurrenceLineage triangles;
  PreparedPatchGridProjectionAudit projections;
};

ProjectionPipeline build_projection_pipeline(
    ProjectionFixture &fixture,
    std::span<const PatchGridActiveCellRequest> active_cells) {
  ProjectionPipeline pipeline;
  const auto requests = cell_vertex_requests(fixture.spec, active_cells);
  pipeline.vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, requests);
  pipeline.cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, pipeline.vertices, active_cells);
  pipeline.triangles = build_patch_grid_pre_orientation_triangle_lineage(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, pipeline.vertices, pipeline.cells);
  pipeline.projections = build_patch_grid_projection_audit(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, pipeline.vertices, pipeline.cells, pipeline.triangles);
  return pipeline;
}

bool crossing_cell(const ProjectionFixture &fixture, RayAxis axis,
                   std::size_t crossing_index,
                   PatchGridActiveCellRequest &cell) {
  const AxisRayCrossings &rays =
      fixture.grid.axes[static_cast<std::size_t>(axis)];
  if (crossing_index >= rays.coordinates.size() || rays.offsets.empty())
    return false;
  const auto after = std::upper_bound(rays.offsets.begin(), rays.offsets.end(),
                                      crossing_index);
  if (after == rays.offsets.begin())
    return false;
  const std::size_t row =
      static_cast<std::size_t>(std::prev(after) - rays.offsets.begin());
  if (row + 1U >= rays.offsets.size() || crossing_index < rays.offsets[row] ||
      crossing_index >= rays.offsets[row + 1U])
    return false;
  DecodedPatchRayRow decoded;
  if (!decode_patch_ray_row(fixture.spec, axis, row, decoded))
    return false;
  const double lower = axis == RayAxis::x   ? fixture.spec.lower.x
                       : axis == RayAxis::y ? fixture.spec.lower.y
                                            : fixture.spec.lower.z;
  const std::size_t count = axis == RayAxis::x   ? fixture.spec.nx
                            : axis == RayAxis::y ? fixture.spec.ny
                                                 : fixture.spec.nz;
  const double coordinate = rays.coordinates[crossing_index];
  const double scaled = (coordinate - lower) / fixture.spec.grid_spacing;
  if (!std::isfinite(scaled) || scaled < 0.0)
    return false;
  std::size_t along = static_cast<std::size_t>(std::floor(scaled));
  if (along + 1U >= count) {
    if (along == count - 1U && along != 0U)
      --along;
    else
      return false;
  }
  switch (axis) {
  case RayAxis::x:
    // This is synthetic already-classified C1 topology metadata. The
    // projection smoke does not re-prove atom/probe union topology.
    cell = {along, decoded.u_index, decoded.v_index, 1U};
    break;
  case RayAxis::y:
    cell = {decoded.u_index, along, decoded.v_index, 1U};
    break;
  case RayAxis::z:
    cell = {decoded.u_index, decoded.v_index, along, 1U};
    break;
  }
  return cell.ix + 1U < fixture.spec.nx && cell.iy + 1U < fixture.spec.ny &&
         cell.iz + 1U < fixture.spec.nz;
}

ProjectionPipeline find_patch_kind_pipeline(ProjectionFixture &fixture,
                                            PatchKind kind,
                                            bool require_multi_ref) {
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const RayAxis axis = static_cast<RayAxis>(axis_index);
    const AxisPatchProvenance &provenance =
        fixture.grid.provenance.axes[axis_index];
    for (std::size_t crossing = 0U; crossing < provenance.canonical_refs.size();
         ++crossing) {
      const std::span<const PatchRefId> refs = provenance.refs(crossing);
      if ((require_multi_ref && refs.size() < 2U) ||
          std::none_of(refs.begin(), refs.end(),
                       [kind](PatchRefId ref) { return ref.kind() == kind; }))
        continue;
      PatchGridActiveCellRequest cell;
      if (!crossing_cell(fixture, axis, crossing, cell))
        continue;
      ProjectionPipeline pipeline =
          build_projection_pipeline(fixture, std::span(&cell, 1U));
      const std::size_t slot = static_cast<std::size_t>(kind);
      if (!pipeline.projections.input_binding_certified ||
          slot >= pipeline.projections.vertex_projector_calls_by_kind.size() ||
          pipeline.projections.vertex_projector_calls_by_kind[slot] == 0U)
        continue;
      if (require_multi_ref &&
          std::none_of(
              pipeline.projections.vertices.begin(),
              pipeline.projections.vertices.end(),
              [](const PatchGridCanonicalVertexProjectionAudit &vertex) {
                return vertex.unique_candidate_status ==
                       PatchUniqueCandidateProjectionStatus::
                           multiple_candidate_refs;
              }))
        continue;
      return pipeline;
    }
  }
  throw std::runtime_error(require_multi_ref
                               ? "no certified seam candidate pipeline found"
                               : "no requested patch-kind pipeline found");
}

std::uint64_t
projection_checksum(const PreparedPatchGridProjectionAudit &prepared) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(prepared.complex_fingerprint.first);
  mix(prepared.complex_fingerprint.second);
  mix(prepared.complex_fingerprint.valid);
  for (const std::uint64_t value : prepared.lattice_identity.scalar_bits)
    mix(value);
  for (const std::size_t value : prepared.lattice_identity.counts)
    mix(value);
  mix(prepared.lattice_identity.valid);
  for (const auto &vertex : prepared.vertices) {
    mix(vertex.canonical_vertex_position);
    mix(vertex.key.first_vertex_linear_index());
    mix(vertex.key.second_vertex_linear_index());
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.x));
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.y));
    mix(std::bit_cast<std::uint64_t>(vertex.canonical_interpolation_point.z));
    mix(std::bit_cast<std::uint64_t>(vertex.edge_distance_upper_bound));
    mix(vertex.exact_location_crossing.packed_route());
    mix(vertex.exact_location_crossing.crossing());
    mix(vertex.exact_location_ref.packed());
    mix(vertex.unique_candidate_ref.packed());
    mix(static_cast<std::uint64_t>(vertex.unique_candidate_status));
    mix(static_cast<std::uint64_t>(vertex.exact_location_status));
    mix(vertex.exact_location_preserved);
    mix(vertex.bounded_candidate_projection_complete);
    mix(vertex.unique_candidate_classification_complete);
    mix(vertex.unique_candidate_projection_certified);
  }
  for (const std::size_t offset : prepared.vertex_projection_offsets)
    mix(offset);
  for (const auto &projection : prepared.candidate_projections) {
    mix(projection.ref.packed());
    mix(std::bit_cast<std::uint64_t>(projection.projected_point.x));
    mix(std::bit_cast<std::uint64_t>(projection.projected_point.y));
    mix(std::bit_cast<std::uint64_t>(projection.projected_point.z));
    mix(std::bit_cast<std::uint64_t>(projection.projected_normal.x));
    mix(std::bit_cast<std::uint64_t>(projection.projected_normal.y));
    mix(std::bit_cast<std::uint64_t>(projection.projected_normal.z));
    mix(std::bit_cast<std::uint64_t>(projection.distance));
    mix(static_cast<std::uint64_t>(projection.status));
    mix(projection.bounded_projection_certified);
    mix(projection.normal_certified);
  }
  for (const auto &triangle : prepared.triangles) {
    mix(triangle.triangle_position);
    mix(triangle.common_candidate_ref.packed());
    mix(std::bit_cast<std::uint64_t>(triangle.canonical_centroid.x));
    mix(std::bit_cast<std::uint64_t>(triangle.canonical_centroid.y));
    mix(std::bit_cast<std::uint64_t>(triangle.canonical_centroid.z));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_point.x));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_point.y));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_point.z));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_normal.x));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_normal.y));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projected_normal.z));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_projection_distance));
    mix(std::bit_cast<std::uint64_t>(triangle.centroid_maximum_distance));
    mix(static_cast<std::uint64_t>(triangle.centroid_status));
    mix(std::bit_cast<std::uint64_t>(triangle.orientation_interval.lower));
    mix(std::bit_cast<std::uint64_t>(triangle.orientation_interval.upper));
    mix(static_cast<std::uint64_t>(triangle.orientation_interval.status));
    mix(triangle.orientation_interval.rounding_contract_available);
    mix(triangle.orientation_interval.interval_finite);
    mix(triangle.orientation_interval.sign_certified);
    mix(triangle.exact_location_vertex_count);
    mix(triangle.same_candidate_vertex_projection_certified);
    mix(triangle.centroid_bounded_projection_certified);
    mix(triangle.suggested_swap_preorientation);
    mix(triangle.planar_face_surface_location_certified);
    mix(triangle.final_mesh_orientation_certified);
  }
  for (const std::size_t value : prepared.vertex_projector_calls_by_kind)
    mix(value);
  for (const std::size_t value : prepared.vertex_projector_successes_by_kind)
    mix(value);
  for (const std::size_t value : prepared.vertex_projector_ambiguous_by_kind)
    mix(value);
  mix(prepared.vertex_projector_calls);
  mix(prepared.vertex_projector_successes);
  mix(prepared.vertex_projector_ambiguous);
  mix(prepared.exact_locations_preserved);
  mix(prepared.unique_candidate_vertices);
  mix(prepared.centroid_projector_calls);
  mix(prepared.centroid_projector_successes);
  mix(prepared.orientation_positive);
  mix(prepared.orientation_negative);
  mix(prepared.orientation_ambiguous);
  mix(prepared.input_binding_certified);
  mix(prepared.vertex_projection_classification_complete);
  mix(prepared.triangle_projection_classification_complete);
  mix(prepared.orientation_classification_complete);
  mix(prepared.planar_face_surface_location_certified);
  mix(prepared.final_mesh_orientation_certified);
  for (const unsigned char value : prepared.status)
    mix(value);
  return hash;
}

void check_orientation_interval() {
  require(patch_orientation_rounding_contract_available(),
          "default projection rounding contract is unavailable");
  const Vec3 a{0.0, 0.0, 0.0};
  const Vec3 b{1.0, 0.0, 0.0};
  const Vec3 c{0.0, 1.0, 0.0};
  const auto positive =
      audit_patch_orientation_sign_interval(a, b, c, {0.0, 0.0, 1.0});
  const auto negative =
      audit_patch_orientation_sign_interval(a, b, c, {0.0, 0.0, -1.0});
  const auto zero = audit_patch_orientation_sign_interval(a, b, {2.0, 0.0, 0.0},
                                                          {0.0, 0.0, 1.0});
  require(positive.status == PatchOrientationSignStatus::positive &&
              positive.sign_certified && positive.lower > 0.0 &&
              negative.status == PatchOrientationSignStatus::negative &&
              negative.sign_certified && negative.upper < 0.0 &&
              zero.status ==
                  PatchOrientationSignStatus::ambiguous_zero_interval &&
              !zero.sign_certified,
          "projection orientation interval sign classification failed");

  for (const double scale : {0.125, 1.0, 8.0}) {
    const Vec3 translated{1000000.0, -2000000.0, 3000000.0};
    const Vec3 ta{translated.x, translated.y, translated.z};
    const Vec3 tb{translated.x + scale, translated.y, translated.z};
    const Vec3 tc{translated.x, translated.y + scale, translated.z};
    const auto result =
        audit_patch_orientation_sign_interval(ta, tb, tc, {0.0, 0.0, 1.0});
    require(result.status == PatchOrientationSignStatus::positive &&
                result.sign_certified,
            "projection orientation lost translation/scale sign");
  }

  const int original_rounding = std::fegetround();
  require(original_rounding == FE_TONEAREST && std::fesetround(FE_UPWARD) == 0,
          "projection rounding-mode fixture setup failed");
  const auto wrong_rounding =
      audit_patch_orientation_sign_interval(a, b, c, {0.0, 0.0, 1.0});
  require(std::fesetround(original_rounding) == 0 &&
              wrong_rounding.status ==
                  PatchOrientationSignStatus::ambiguous_rounding_environment &&
              !wrong_rounding.sign_certified,
          "non-nearest rounding mode did not fail closed");
  volatile double denormal = std::numeric_limits<double>::denorm_min();
  volatile double retained = denormal * 1.0;
  require(retained == denormal && retained != 0.0,
          "test environment unexpectedly flushes subnormals");
  const auto nonfinite = audit_patch_orientation_sign_interval(
      a, b, c, {0.0, 0.0, std::numeric_limits<double>::infinity()});
  require(nonfinite.status == PatchOrientationSignStatus::ambiguous_nonfinite,
          "nonfinite normal did not fail closed");
}

PatchGridActiveCellRequest centered_cell(const ProjectionFixture &fixture,
                                         std::uint8_t inside_mask) {
  PatchGridActiveCellRequest cell;
  cell.ix = grid_coordinate_index(fixture.spec.lower.x,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iy = grid_coordinate_index(fixture.spec.lower.y,
                                  fixture.spec.grid_spacing, 0.0);
  cell.iz = grid_coordinate_index(fixture.spec.lower.z,
                                  fixture.spec.grid_spacing, 0.0);
  cell.inside_mask = inside_mask;
  return cell;
}

void check_contact_projection_and_all_axis_boundary() {
  ProjectionFixture coarse =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 1U);
  const PatchGridActiveCellRequest cell = centered_cell(coarse, 1U);
  const ProjectionPipeline pipeline =
      build_projection_pipeline(coarse, std::span(&cell, 1U));
  const auto &audit = pipeline.projections;
  require((audit.status == "pass" ||
           audit.status == "pass_with_ambiguous_projection") &&
              audit.input_binding_certified && !audit.vertices.empty() &&
              audit.vertices.size() ==
                  pipeline.triangles.canonical_vertices.size() &&
              audit.vertex_projection_offsets.size() ==
                  audit.vertices.size() + 1U &&
              audit.vertex_projection_offsets.back() ==
                  audit.candidate_projections.size() &&
              audit.vertex_projector_calls ==
                  audit.candidate_projections.size() &&
              audit.vertex_projector_successes > 0U &&
              audit.vertex_projector_calls_by_kind[0] ==
                  audit.vertex_projector_calls &&
              audit.triangles.size() == pipeline.triangles.triangles.size() &&
              audit.triangle_projection_classification_complete &&
              audit.orientation_classification_complete &&
              !audit.planar_face_surface_location_certified &&
              !audit.final_mesh_orientation_certified,
          "contact C2b-A projection audit failed");
  for (std::size_t vertex = 0U; vertex < audit.vertices.size(); ++vertex) {
    const auto &item = audit.vertices[vertex];
    require(item.exact_location_preserved ==
                    pipeline.triangles.canonical_vertices[vertex]
                        .exact_analytic_location_certified &&
                item.edge_distance_upper_bound >= coarse.spec.grid_spacing,
            "C2b-A exact location or lattice bound drifted");
    if (item.unique_candidate_projection_certified)
      require(item.unique_candidate_status ==
                  PatchUniqueCandidateProjectionStatus::certified_unique,
              "C2b-A unique candidate status mismatch");
  }
  for (const auto &triangle : audit.triangles)
    require(!triangle.planar_face_surface_location_certified &&
                !triangle.final_mesh_orientation_certified,
            "C2b-A promoted planar/final mesh certificate");

  ProjectionFixture fine =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 0.2, 1U);
  const PatchGridActiveCellRequest all_axis{8U, 9U, 14U, 0x4EU};
  const ProjectionPipeline all_axis_pipeline =
      build_projection_pipeline(fine, std::span(&all_axis, 1U));
  require(all_axis_pipeline.triangles.exact_location_vertices == 1U &&
              all_axis_pipeline.projections.exact_locations_preserved == 1U &&
              all_axis_pipeline.projections.vertex_projector_successes > 0U &&
              all_axis_pipeline.projections.centroid_projector_successes > 0U &&
              std::all_of(
                  all_axis_pipeline.projections.triangles.begin(),
                  all_axis_pipeline.projections.triangles.end(),
                  [](const PatchGridTriangleProjectionAudit &triangle) {
                    return !triangle.planar_face_surface_location_certified &&
                           !triangle.final_mesh_orientation_certified;
                  }),
          "C2b-A all-axis candidate boundary changed exact C2a claims");
}

void check_torus_reentrant_and_seam_projection() {
  PatchComplex de_aligned_complex = make_torus_reentrant_complex();
  const PatchMappingGridSpec de_aligned_spec =
      enclosing_complex_grid(de_aligned_complex, 0.5);
  ProjectionFixture fixture =
      build_complex_fixture(std::move(de_aligned_complex), de_aligned_spec, 1U);
  const ProjectionPipeline torus =
      find_patch_kind_pipeline(fixture, PatchKind::toroidal, false);
  const ProjectionPipeline reentrant =
      find_patch_kind_pipeline(fixture, PatchKind::reentrant, false);
  require(torus.projections.vertex_projector_calls_by_kind[1] > 0U &&
              torus.projections.vertex_projector_successes_by_kind[1] > 0U &&
              reentrant.projections.vertex_projector_calls_by_kind[2] > 0U &&
              reentrant.projections.vertex_projector_successes_by_kind[2] > 0U,
          "C2b-B torus/reentrant bounded projection was not exercised");

  PatchComplex seam_complex = make_torus_reentrant_complex();
  const PatchMappingGridSpec seam_spec =
      seam_aligned_complex_grid(seam_complex, 0.5);
  ProjectionFixture seam_fixture =
      build_complex_fixture(std::move(seam_complex), seam_spec, 1U);
  require(seam_fixture.grid.provenance.seam_groups > 0U,
          "production C2b-B fixture retained no true seam group");
  const ProjectionPipeline seam =
      find_patch_kind_pipeline(seam_fixture, PatchKind::toroidal, true);
  std::size_t multi_vertices = 0U;
  bool torus_reentrant_multi_ref_observed = false;
  std::vector<std::size_t> multi_positions;
  for (std::size_t position = 0U; position < seam.projections.vertices.size();
       ++position) {
    const auto &vertex = seam.projections.vertices[position];
    if (vertex.unique_candidate_status ==
        PatchUniqueCandidateProjectionStatus::multiple_candidate_refs) {
      ++multi_vertices;
      multi_positions.push_back(position);
      require(position + 1U < seam.projections.vertex_projection_offsets.size(),
              "C2b-B seam candidate CSR is malformed");
      bool has_torus = false;
      bool has_reentrant = false;
      for (std::size_t projection =
               seam.projections.vertex_projection_offsets[position];
           projection <
           seam.projections.vertex_projection_offsets[position + 1U];
           ++projection) {
        require(projection < seam.projections.candidate_projections.size(),
                "C2b-B seam candidate range is malformed");
        const PatchKind kind =
            seam.projections.candidate_projections[projection].ref.kind();
        has_torus |= kind == PatchKind::toroidal;
        has_reentrant |= kind == PatchKind::reentrant;
      }
      torus_reentrant_multi_ref_observed |= has_torus && has_reentrant;
      require(
          !vertex.unique_candidate_projection_certified &&
              !vertex.exact_location_preserved &&
              vertex.exact_location_status !=
                  PatchExactLocationPreservationStatus::preserved_single_ref,
          "C2b-B seam evidence was promoted to exact/unique identity");
    }
  }
  std::size_t seam_triangles = 0U;
  for (std::size_t triangle_position = 0U;
       triangle_position < seam.triangles.triangles.size();
       ++triangle_position) {
    const auto &source = seam.triangles.triangles[triangle_position];
    const bool contains_multi = std::any_of(
        source.canonical_vertex_positions.begin(),
        source.canonical_vertex_positions.end(),
        [&](std::size_t vertex_position) {
          return std::find(multi_positions.begin(), multi_positions.end(),
                           vertex_position) != multi_positions.end();
        });
    if (!contains_multi)
      continue;
    ++seam_triangles;
    const auto &audit = seam.projections.triangles[triangle_position];
    require(audit.orientation_interval.status ==
                    PatchOrientationSignStatus::ambiguous_candidate_identity &&
                !audit.orientation_interval.sign_certified &&
                !audit.suggested_swap_preorientation,
            "C2b-B multi-ref triangle received an orientation sign");
  }
  require(multi_vertices > 0U && torus_reentrant_multi_ref_observed &&
              seam_triangles > 0U &&
              seam.projections.orientation_ambiguous >= seam_triangles &&
              std::all_of(
                  seam.projections.triangles.begin(),
                  seam.projections.triangles.end(),
                  [](const PatchGridTriangleProjectionAudit &triangle) {
                    return !triangle.planar_face_surface_location_certified &&
                           !triangle.final_mesh_orientation_certified;
                  }),
          "C2b-B seam did not remain orientation/face ambiguous");
  std::cout << "c2b_torus_calls="
            << torus.projections.vertex_projector_calls_by_kind[1]
            << " c2b_reentrant_calls="
            << reentrant.projections.vertex_projector_calls_by_kind[2]
            << " c2b_seam_groups=" << seam_fixture.grid.provenance.seam_groups
            << " c2b_seam_multi_vertices=" << multi_vertices
            << " c2b_seam_ambiguous_triangles=" << seam_triangles << '\n';
}

void require_failed(const PreparedPatchGridProjectionAudit &audit,
                    const char *message) {
  require(
      audit.status != "pass" &&
          audit.status != "pass_with_ambiguous_projection" &&
          audit.vertices.empty() && audit.vertex_projection_offsets.empty() &&
          audit.candidate_projections.empty() && audit.triangles.empty() &&
          audit.vertex_projector_calls == 0U &&
          audit.centroid_projector_calls == 0U &&
          audit.orientation_positive == 0U &&
          audit.orientation_negative == 0U &&
          audit.orientation_ambiguous == 0U && !audit.input_binding_certified &&
          !audit.vertex_projection_classification_complete &&
          !audit.triangle_projection_classification_complete &&
          !audit.orientation_classification_complete &&
          !audit.planar_face_surface_location_certified &&
          !audit.final_mesh_orientation_certified,
      message);
}

void check_fail_closed_mutations_and_empty() {
  ProjectionFixture fixture =
      build_sphere_fixture(Atom{0.0, 0.0, 0.0, 0.0, 1.0}, 1.0, 2.0, 1U);
  const PatchGridActiveCellRequest cell = centered_cell(fixture, 1U);
  const ProjectionPipeline pipeline =
      build_projection_pipeline(fixture, std::span(&cell, 1U));
  const auto run =
      [&](const LocalPatchProjectionIndex &index, const PatchRayGrid &grid,
          const PreparedPatchCrossingWitnesses &crossings,
          const PreparedPatchGridVertexMagnitudeWitnesses &vertices,
          const PreparedPatchGridCellLineage &cells,
          const PreparedPatchGridTriangleOccurrenceLineage &triangles) {
        return build_patch_grid_projection_audit(*fixture.complex, index, grid,
                                                 fixture.spec, crossings,
                                                 vertices, cells, triangles);
      };

  auto bad_triangles = pipeline.triangles;
  ++bad_triangles.complex_fingerprint.first;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a fingerprint was accepted");
  bad_triangles = pipeline.triangles;
  ++bad_triangles.lattice_identity.scalar_bits[3];
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a lattice was accepted");
  bad_triangles = pipeline.triangles;
  bad_triangles.canonical_vertex_candidate_ref_offsets.back() = 0U;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a CSR was accepted");
  bad_triangles = pipeline.triangles;
  require(!bad_triangles.canonical_vertices.empty(),
          "mutation fixture has no C2a vertex");
  bad_triangles.canonical_vertices.front().canonical_interpolation_point.x =
      std::nextafter(bad_triangles.canonical_vertices.front()
                         .canonical_interpolation_point.x,
                     std::numeric_limits<double>::infinity());
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a point was accepted");
  bad_triangles = pipeline.triangles;
  bad_triangles.status = "pass_with_nonfinite_canonical_geometry";
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a status was accepted");
  bad_triangles = pipeline.triangles;
  ++bad_triangles.exact_location_vertices;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, bad_triangles),
                 "mutated C2a counter was accepted");

  auto bad_cells = pipeline.cells;
  require(!bad_cells.endpoint_patch_refs.empty(),
          "mutation fixture has no C1 patch refs");
  bad_cells.endpoint_patch_refs.front() = {};
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, bad_cells, pipeline.triangles),
                 "mutated C1 ref was accepted");
  auto bad_vertices = pipeline.vertices;
  ++bad_vertices.raw_exact_owner_candidates;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-B counter was accepted");
  bad_vertices = pipeline.vertices;
  ++bad_vertices.topology_nonanalytic_vertices;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "second mutated Slice-B counter was accepted");
  bad_vertices = pipeline.vertices;
  bad_vertices.raw_source_certificate_complete =
      !bad_vertices.raw_source_certificate_complete;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-B complete flag was accepted");
  bad_vertices = pipeline.vertices;
  require(!bad_vertices.vertices.empty(),
          "mutation fixture has no Slice-B vertex");
  bad_vertices.vertices.front().topology_analytic_evidence =
      !bad_vertices.vertices.front().topology_analytic_evidence;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-B topology analytic flag was accepted");
  bad_vertices = pipeline.vertices;
  bad_vertices.vertices.front().topology_exact_vertex_certified =
      !bad_vertices.vertices.front().topology_exact_vertex_certified;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-B topology exact flag was accepted");
  bad_vertices = pipeline.vertices;
  bad_vertices.vertices.front().topology_evidence_status =
      PatchGridTopologyEvidenceStatus::unknown;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     bad_vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-B topology status was accepted");
  bad_cells = pipeline.cells;
  ++bad_cells.axis_edges;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, bad_cells, pipeline.triangles),
                 "mutated C1 summary counter was accepted");
  bad_cells = pipeline.cells;
  bad_cells.scalar_lineage_complete = !bad_cells.scalar_lineage_complete;
  require_failed(run(fixture.index, fixture.grid, fixture.crossings,
                     pipeline.vertices, bad_cells, pipeline.triangles),
                 "mutated C1 complete flag was accepted");
  auto bad_crossings = fixture.crossings;
  require(!bad_crossings.witnesses.empty(),
          "mutation fixture has no Slice-A witness");
  ++bad_crossings.root_reproduced_crossings;
  require_failed(run(fixture.index, fixture.grid, bad_crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-A root counter was accepted");
  bad_crossings = fixture.crossings;
  ++bad_crossings.patch_residual_failures;
  require_failed(run(fixture.index, fixture.grid, bad_crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-A residual counter was accepted");
  bad_crossings = fixture.crossings;
  bad_crossings.location_certificate_complete =
      !bad_crossings.location_certificate_complete;
  require_failed(run(fixture.index, fixture.grid, bad_crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-A complete flag was accepted");
  bad_crossings = fixture.crossings;
  bad_crossings.status = "location_uncertified";
  require_failed(run(fixture.index, fixture.grid, bad_crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-A status was accepted");
  bad_crossings = fixture.crossings;
  bad_crossings.witnesses.front().point.x =
      std::nextafter(bad_crossings.witnesses.front().point.x,
                     std::numeric_limits<double>::infinity());
  require_failed(run(fixture.index, fixture.grid, bad_crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated Slice-A point was accepted");
  auto bad_grid = fixture.grid;
  ++bad_grid.provenance.lattice_identity.scalar_bits[0];
  require_failed(run(fixture.index, bad_grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated ray lattice was accepted");
  auto bad_index = fixture.index;
  ++bad_index.source_complex_fingerprint.second;
  require_failed(run(bad_index, fixture.grid, fixture.crossings,
                     pipeline.vertices, pipeline.cells, pipeline.triangles),
                 "mutated projection-index fingerprint was accepted");

  const std::vector<PatchGridActiveCellRequest> no_cells;
  const auto empty_cell_lineage = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, pipeline.vertices, no_cells);
  const auto empty_triangles =
      build_patch_grid_pre_orientation_triangle_lineage(
          *fixture.complex, fixture.index, fixture.grid, fixture.spec,
          fixture.crossings, pipeline.vertices, empty_cell_lineage);
  const auto empty = build_patch_grid_projection_audit(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, pipeline.vertices, empty_cell_lineage,
      empty_triangles);
  require(empty.status == "empty_input" && empty.vertices.empty() &&
              empty.triangles.empty() && !empty.input_binding_certified &&
              !empty.vertex_projection_classification_complete &&
              !empty.triangle_projection_classification_complete &&
              !empty.orientation_classification_complete,
          "empty C2b-A input vacuously certified a claim");

  const std::vector<PatchGridVertexMagnitudeRequest> no_requests;
  const auto exact_empty_vertices = build_patch_grid_vertex_magnitude_witnesses(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, no_requests);
  const auto exact_empty_cells = build_patch_grid_cell_lineage(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, exact_empty_vertices, no_cells);
  const auto exact_empty_triangles =
      build_patch_grid_pre_orientation_triangle_lineage(
          *fixture.complex, fixture.index, fixture.grid, fixture.spec,
          fixture.crossings, exact_empty_vertices, exact_empty_cells);
  require(exact_empty_vertices.status == "no_vertices" &&
              exact_empty_cells.status == "empty_active_cells" &&
              exact_empty_triangles.status == "invalid_vertex_payload",
          "upstream exact-empty compatibility fixture changed");
  const auto exact_empty = build_patch_grid_projection_audit(
      *fixture.complex, fixture.index, fixture.grid, fixture.spec,
      fixture.crossings, exact_empty_vertices, exact_empty_cells,
      exact_empty_triangles);
  require(exact_empty.status == "empty_input" && exact_empty.vertices.empty() &&
              exact_empty.triangles.empty() &&
              !exact_empty.input_binding_certified &&
              !exact_empty.orientation_classification_complete,
          "upstream exact-empty compatibility became a certificate");
}

std::size_t
exact_retained_bytes(const PreparedPatchGridProjectionAudit &audit) {
  return audit.vertices.capacity() *
             sizeof(PatchGridCanonicalVertexProjectionAudit) +
         audit.vertex_projection_offsets.capacity() * sizeof(std::size_t) +
         audit.candidate_projections.capacity() *
             sizeof(PatchGridBoundedCandidateProjection) +
         audit.triangles.capacity() * sizeof(PatchGridTriangleProjectionAudit);
}

void check_determinism_memory_and_timing() {
  const Atom atom{0.0, 0.0, 0.0, 0.0, 1.0};
  ProjectionFixture serial = build_sphere_fixture(atom, 1.0, 2.0, 1U);
  ProjectionFixture parallel = build_sphere_fixture(atom, 1.0, 2.0, 8U);
  const PatchGridActiveCellRequest cell = centered_cell(serial, 1U);
  const ProjectionPipeline serial_pipeline =
      build_projection_pipeline(serial, std::span(&cell, 1U));
  const ProjectionPipeline parallel_pipeline =
      build_projection_pipeline(parallel, std::span(&cell, 1U));
  const std::uint64_t checksum =
      projection_checksum(serial_pipeline.projections);
  require(checksum == projection_checksum(parallel_pipeline.projections),
          "C2b-A nt=1/8 semantic checksum differs");

  PatchGridActiveCellRequest second{0U, 0U, 0U, 1U};
  const std::array<PatchGridActiveCellRequest, 2> two_cells{cell, second};
  const ProjectionPipeline doubled =
      build_projection_pipeline(serial, two_cells);
  const std::array<PatchGridActiveCellRequest, 2> reversed_cells{second, cell};
  const ProjectionPipeline reversed =
      build_projection_pipeline(serial, reversed_cells);
  require(projection_checksum(doubled.projections) ==
              projection_checksum(reversed.projections),
          "C2b active-cell input order changed canonical output");
  require(serial_pipeline.projections.retained_bytes ==
                  exact_retained_bytes(serial_pipeline.projections) &&
              doubled.projections.retained_bytes ==
                  exact_retained_bytes(doubled.projections) &&
              serial_pipeline.projections.accounted_peak_temporary_bytes >=
                  serial_pipeline.projections.retained_bytes &&
              doubled.projections.accounted_peak_temporary_bytes >=
                  doubled.projections.retained_bytes,
          "C2b-A retained/accounted byte contract failed");
  const double retained_normalized =
      static_cast<double>(doubled.projections.retained_bytes) /
      (2.0 * static_cast<double>(serial_pipeline.projections.retained_bytes));
  const double peak_normalized =
      static_cast<double>(doubled.projections.accounted_peak_temporary_bytes) /
      (2.0 * static_cast<double>(
                 serial_pipeline.projections.accounted_peak_temporary_bytes));
  require(retained_normalized > 0.5 && retained_normalized < 1.5 &&
              peak_normalized > 0.5 && peak_normalized < 1.5,
          "C2b-A memory scaling is not linear");

  for (int warmup = 0; warmup < 2; ++warmup) {
    const auto ignored = build_patch_grid_projection_audit(
        *serial.complex, serial.index, serial.grid, serial.spec,
        serial.crossings, serial_pipeline.vertices, serial_pipeline.cells,
        serial_pipeline.triangles);
    require(ignored.input_binding_certified, "C2b-A timing warmup failed");
  }
  std::array<double, 7> samples{};
  for (double &sample : samples) {
    const auto timed = build_patch_grid_projection_audit(
        *serial.complex, serial.index, serial.grid, serial.spec,
        serial.crossings, serial_pipeline.vertices, serial_pipeline.cells,
        serial_pipeline.triangles);
    require(timed.input_binding_certified &&
                projection_checksum(timed) == checksum,
            "C2b-A timed sample changed semantics");
    sample = timed.build_ms;
  }
  std::array<double, 7> sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  std::cout << "c2b_checksum=" << checksum
            << " c2b_vertices=" << serial_pipeline.projections.vertices.size()
            << " c2b_candidate_calls="
            << serial_pipeline.projections.vertex_projector_calls
            << " c2b_centroid_calls="
            << serial_pipeline.projections.centroid_projector_calls
            << " c2b_orientation_positive="
            << serial_pipeline.projections.orientation_positive
            << " c2b_orientation_negative="
            << serial_pipeline.projections.orientation_negative
            << " c2b_orientation_ambiguous="
            << serial_pipeline.projections.orientation_ambiguous
            << " c2b_retained_one="
            << serial_pipeline.projections.retained_bytes
            << " c2b_retained_two=" << doubled.projections.retained_bytes
            << " c2b_peak_one="
            << serial_pipeline.projections.accounted_peak_temporary_bytes
            << " c2b_peak_two="
            << doubled.projections.accounted_peak_temporary_bytes
            << " c2b_retained_normalized=" << retained_normalized
            << " c2b_peak_normalized=" << peak_normalized
            << " c2b_samples_ms=" << samples[0] << ',' << samples[1] << ','
            << samples[2] << ',' << samples[3] << ',' << samples[4] << ','
            << samples[5] << ',' << samples[6] << " c2b_median_ms=" << sorted[3]
            << '\n';
}

void check_rounding_contract_only() {
  const Vec3 a{0.0, 0.0, 0.0};
  const Vec3 b{1.0, 0.0, 0.0};
  const Vec3 c{0.0, 1.0, 0.0};
  const PatchOrientationSignInterval result =
      audit_patch_orientation_sign_interval(a, b, c, {0.0, 0.0, 1.0});
#if defined(__FAST_MATH__)
  require(!patch_orientation_rounding_contract_available() &&
              result.status ==
                  PatchOrientationSignStatus::ambiguous_rounding_environment &&
              !result.sign_certified,
          "fast-math build unexpectedly certified interval orientation");
#else
  require(patch_orientation_rounding_contract_available() &&
              result.status == PatchOrientationSignStatus::positive &&
              result.sign_certified,
          "strict floating-point build lost interval orientation certificate");
#endif
  std::cout << "c2b_rounding_contract_only=pass fast_math="
#if defined(__FAST_MATH__)
            << 1
#else
            << 0
#endif
            << '\n';
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--rounding-contract-only") {
      check_rounding_contract_only();
      return EXIT_SUCCESS;
    }
    check_orientation_interval();
    check_contact_projection_and_all_axis_boundary();
    check_torus_reentrant_and_seam_projection();
    check_fail_closed_mutations_and_empty();
    check_determinism_memory_and_timing();
    std::cout << "tabi_ses_patch_projection_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_projection_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
