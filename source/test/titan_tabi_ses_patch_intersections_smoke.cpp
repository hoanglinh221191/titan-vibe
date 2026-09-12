#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_intersections.h"
#include "tabi/tabi_ses_patch_provenance.h"
#include "tabi/tabi_ses_quartic.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

PatchRefId make_patch_ref(PatchKind kind, std::size_t index) {
  PatchRefId ref;
  require(PatchRefId::try_make(kind, index, ref) && ref.valid(),
          "failed to pack patch ref");
  return ref;
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

Vec3 vadd(const Vec3 &a, const Vec3 &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 vsub(const Vec3 &a, const Vec3 &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 vscale(double value, const Vec3 &a) {
  return {value * a.x, value * a.y, value * a.z};
}

double vdot(const Vec3 &a, const Vec3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 vcross(const Vec3 &a, const Vec3 &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double vnorm(const Vec3 &a) { return std::sqrt(vdot(a, a)); }

Vec3 vnormalized(const Vec3 &a) {
  const double length = vnorm(a);
  require(length > 1.0e-14, "test fixture has zero direction");
  return vscale(1.0 / length, a);
}

Vec3 atom_position(const Atom &atom) { return {atom.x, atom.y, atom.z}; }

ToroidalPatch make_welded_torus(const PatchComplex &complex, std::size_t atom_i,
                                std::size_t atom_j, const Vec3 &probe_center) {
  const Vec3 center_i = atom_position(complex.atoms[atom_i]);
  const Vec3 center_j = atom_position(complex.atoms[atom_j]);
  const double expanded_i = complex.atoms[atom_i].radius + complex.probe_radius;
  const double expanded_j = complex.atoms[atom_j].radius + complex.probe_radius;
  const Vec3 delta = vsub(center_j, center_i);
  const double distance = vnorm(delta);
  const Vec3 axis = vscale(1.0 / distance, delta);
  const double offset = (expanded_i * expanded_i - expanded_j * expanded_j +
                         distance * distance) /
                        (2.0 * distance);
  const Vec3 circle_center = vadd(center_i, vscale(offset, axis));
  const double circle_radius =
      std::sqrt(expanded_i * expanded_i - offset * offset);
  const Vec3 basis_u = vnormalized(vsub(probe_center, circle_center));
  const Vec3 basis_v = vnormalized(vcross(axis, basis_u));
  return {atom_i,     atom_j,     circle_center, axis,
          basis_u,    basis_v,    circle_radius, offset,
          expanded_i, expanded_j, {{0.0, 0.60}}};
}

struct SeamFixture {
  PatchComplex complex;
  Vec3 query;
  Vec3 expected;
};

SeamFixture make_torus_reentrant_seam_fixture() {
  SeamFixture fixture;
  PatchComplex &complex = fixture.complex;
  complex.probe_radius = 1.0;
  complex.status = "pass";
  const double primary_circle_radius = std::sqrt(1.75);
  complex.atoms = {{0.0, 0.0, -1.5, 0.0, 1.0},
                   {0.0, 0.0, 1.5, 0.0, 1.0},
                   {primary_circle_radius, 2.0, 0.0, 0.0, 1.0}};
  const Vec3 reentrant_center{primary_circle_radius, 0.0, 0.0};
  ReentrantPatch reentrant;
  reentrant.atoms = {0U, 1U, 2U};
  reentrant.probe_center = reentrant_center;
  for (std::size_t atom = 0; atom < 3U; ++atom)
    reentrant.contact_directions[atom] =
        vnormalized(vsub(atom_position(complex.atoms[atom]), reentrant_center));
  complex.reentrant.push_back(reentrant);
  complex.tori.push_back(make_welded_torus(complex, 0U, 1U, reentrant_center));
  complex.tori.push_back(make_welded_torus(complex, 1U, 2U, reentrant_center));
  complex.tori.push_back(make_welded_torus(complex, 0U, 2U, reentrant_center));

  const double seam_radius = primary_circle_radius - complex.probe_radius;
  fixture.expected = {seam_radius, 0.0, 0.0};
  constexpr double outside_theta = -0.20;
  const double query_radius = seam_radius + 0.05;
  fixture.query = {query_radius * std::cos(outside_theta),
                   query_radius * std::sin(outside_theta), 0.0};

  // Stale owner has a valid interior contact hit at distance 0.30. The
  // constrained torus/reentrant seam is closer and must replace that hint.
  complex.atoms.push_back({fixture.query.x, fixture.query.y, 0.50, 0.0, 0.20});
  complex.contacts.push_back({3U, atom_position(complex.atoms[3]), 0.20, 1.20});
  complex.atom_neighbors.resize(complex.atoms.size());
  return fixture;
}

void compare_shifted_roots(const PatchComplex &base,
                           const PatchComplex &shifted, double shift) {
  const ReentrantProbeIndex base_index = build_reentrant_probe_index(base);
  const ReentrantProbeIndex shifted_index =
      build_reentrant_probe_index(shifted);
  require(base.tori.size() == shifted.tori.size(),
          "translation changed torus topology");
  for (std::size_t i = 0; i < base.tori.size(); ++i) {
    const double y = base.tori[i].circle_center.y + 0.37;
    const double z = base.tori[i].circle_center.z;
    const PatchAxisRoots a = intersect_patch_axis(
        base, base_index, PatchKind::toroidal, i, RayAxis::x, y, z);
    const PatchAxisRoots b =
        intersect_patch_axis(shifted, shifted_index, PatchKind::toroidal, i,
                             RayAxis::x, y - 17.0, z + 29.0);
    require(a.certified && b.certified, "torus roots were not certified");
    require(a.count == b.count, "translation changed torus root count");
    for (std::size_t root = 0; root < a.count; ++root)
      require(std::abs((b.coordinates[root] - shift) - a.coordinates[root]) <
                  2.0e-9,
              "translation changed local torus root");
  }
}

LocalProjectionHit check_spatial_projection_outside_owner_hint() {
  constexpr std::size_t kAtomCount = 96U;
  std::vector<Atom> atoms;
  atoms.reserve(kAtomCount);
  for (std::size_t atom = 0; atom < kAtomCount; ++atom)
    atoms.push_back({8.0 * static_cast<double>(atom), 0.0, 0.0, 0.0, 1.0});

  const PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  require(complex.status == "pass" && complex.contacts.size() == kAtomCount &&
              complex.tori.empty() && complex.reentrant.empty(),
          "separated-atom BVH fixture failed");
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  require(index.spatial_entries.size() == kAtomCount &&
              !index.spatial_nodes.empty(),
          "spatial patch index was not prepared");

  const Vec3 query{atoms.back().x + atoms.back().radius + 0.05, 0.0, 0.0};
  const LocalProjectionHit hit =
      project_local_patch_closest(complex, index, 0U, query, 0.10);
  const LocalProjectionHit repeated =
      project_local_patch_closest(complex, index, 0U, query, 0.10);
  require(hit.found && hit.kind == PatchKind::contact &&
              hit.patch_index == kAtomCount - 1U &&
              hit.atom_index == kAtomCount - 1U,
          "BVH missed nearest patch outside stale owner hint");
  require(hit.spatial_fallback && !hit.owner_hint_hit && !hit.axis_fallback &&
              hit.bounded_closest_certified,
          "outside-owner result did not come from certified spatial search");
  require(hit.candidate_patches < kAtomCount / 4U &&
              hit.bound_pruned_patches >= kAtomCount / 2U,
          "BVH query disguised a global all-patch scan");
  require(repeated.found && repeated.kind == hit.kind &&
              repeated.patch_index == hit.patch_index &&
              repeated.point.x == hit.point.x &&
              repeated.point.y == hit.point.y &&
              repeated.point.z == hit.point.z &&
              repeated.candidate_patches == hit.candidate_patches &&
              repeated.bound_pruned_patches == hit.bound_pruned_patches &&
              repeated.spatial_nodes_visited == hit.spatial_nodes_visited,
          "BVH projection or decision counters are not deterministic");
  return hit;
}

struct LargeIndexTiming {
  double patch_build_ms{0.0};
  double index_build_ms{0.0};
  std::size_t retained_bytes{0U};
  std::size_t stamp_bytes_per_worker{0U};
};

LargeIndexTiming measure_large_index_setup() {
  constexpr std::size_t kAtomCount = 4096U;
  std::vector<Atom> atoms;
  atoms.reserve(kAtomCount);
  for (std::size_t atom = 0; atom < kAtomCount; ++atom)
    atoms.push_back({8.0 * static_cast<double>(atom), 0.0, 0.0, 0.0, 1.0});
  const auto begin = std::chrono::steady_clock::now();
  const PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  const auto patch_end = std::chrono::steady_clock::now();
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  const auto index_end = std::chrono::steady_clock::now();
  require(complex.status == "pass" &&
              index.spatial_entries.size() == kAtomCount &&
              index.bounded_domains_certified,
          "large deterministic index fixture failed");
  return {
      std::chrono::duration<double, std::milli>(patch_end - begin).count(),
      std::chrono::duration<double, std::milli>(index_end - patch_end).count(),
      index.retained_bytes(), index.stamp_bytes_per_worker()};
}

struct ConstrainedSeamHits {
  LocalProjectionHit torus_reentrant;
  LocalProjectionHit contact_torus;
  LocalProjectionHit reentrant_only;
  std::size_t torus_reentrant_pairs_audited{0U};
  double torus_reentrant_audit_ms{0.0};
};

ConstrainedSeamHits check_constrained_seam_minima() {
  SeamFixture fixture = make_torus_reentrant_seam_fixture();
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(fixture.complex);
  require(index.bounded_domains_certified,
          "welded seam fixture was not runtime-certified");
  const LocalProjectionHit seam = project_local_patch_closest(
      fixture.complex, index, 3U, fixture.query, 0.40);
  const double expected_distance = vnorm(vsub(fixture.query, fixture.expected));
  require(seam.found &&
              (seam.kind == PatchKind::toroidal ||
               seam.kind == PatchKind::reentrant) &&
              vnorm(vsub(seam.point, fixture.expected)) < 2.0e-10 &&
              std::abs(seam.distance - expected_distance) < 2.0e-10,
          "bounded torus-reentrant seam minimum was missed");
  require(seam.spatial_fallback && !seam.owner_hint_hit &&
              !seam.axis_fallback && seam.bounded_closest_certified,
          "torus-reentrant seam result was not BVH-certified");

  // Production pair: the unrestricted radial foot on contact atom 0 lies
  // inside atom 1's expanded sphere. The active half-space minimum is the
  // analytic contact-torus circle.
  const std::vector<Atom> contact_pair = {{-1.5, 0.0, 0.0, 0.0, 1.0},
                                          {1.5, 0.0, 0.0, 0.0, 1.0}};
  PatchComplex contact_complex =
      build_patch_complex(contact_pair, 1.0, nullptr, false);
  require(contact_complex.status == "pass" &&
              contact_complex.contacts.size() == 2U &&
              contact_complex.tori.size() == 1U &&
              contact_complex.atom_neighbors[0].size() == 1U,
          "production contact constraint fixture failed");
  const Vec3 unrestricted_direction{0.90, std::sqrt(1.0 - 0.90 * 0.90), 0.0};
  const Vec3 contact_query = vadd(atom_position(contact_complex.atoms[0]),
                                  vscale(1.05, unrestricted_direction));
  const Vec3 unrestricted_probe = vadd(atom_position(contact_complex.atoms[0]),
                                       vscale(2.0, unrestricted_direction));
  require(vnorm(vsub(unrestricted_probe,
                     atom_position(contact_complex.atoms[1]))) < 2.0,
          "unrestricted contact foot was not blocked");
  constexpr double active_x = 0.75;
  const Vec3 active_direction{active_x, std::sqrt(1.0 - active_x * active_x),
                              0.0};
  const Vec3 contact_expected =
      vadd(atom_position(contact_complex.atoms[0]), active_direction);
  const LocalPatchProjectionIndex contact_index =
      build_local_patch_projection_index(contact_complex);
  const LocalProjectionHit contact_seam = project_patch_bounded_closest(
      contact_complex, contact_index, PatchKind::contact, 0U, contact_query,
      0.50);
  require(contact_index.bounded_domains_certified && contact_seam.found &&
              contact_seam.kind == PatchKind::contact &&
              contact_seam.patch_index == 0U &&
              vnorm(vsub(contact_seam.point, contact_expected)) < 2.0e-10 &&
              contact_seam.bounded_closest_certified,
          "active contact constraint seam minimum was missed");
  const LocalProjectionHit contact_global = project_local_patch_closest(
      contact_complex, contact_index, 0U, contact_query, 0.50);
  require(contact_global.found && contact_global.bounded_closest_certified &&
              contact_global.distance <= contact_seam.distance + 2.0e-12,
          "global pair projection did not remain certified");
  const double expanded_boundary =
      vnorm(vsub(vadd(atom_position(contact_complex.atoms[0]),
                      vscale(2.0, active_direction)),
                 atom_position(contact_complex.atoms[1])));
  require(std::abs(expanded_boundary - 2.0) < 2.0e-12,
          "contact active-set point is not on the blocker seam");

  // Remove the corresponding torus only to exercise the reentrant edge
  // active-set path. Runtime welding must now fail closed even though the
  // geometric boundary point itself remains recoverable.
  PatchComplex reentrant_only = fixture.complex;
  reentrant_only.tori.erase(reentrant_only.tori.begin());
  const LocalPatchProjectionIndex reentrant_index =
      build_local_patch_projection_index(reentrant_only);
  const LocalProjectionHit reentrant_hit = project_local_patch_closest(
      reentrant_only, reentrant_index, 3U, fixture.query, 0.40);
  require(!reentrant_index.bounded_domains_certified && reentrant_hit.found &&
              reentrant_hit.kind == PatchKind::reentrant &&
              vnorm(vsub(reentrant_hit.point, fixture.expected)) < 2.0e-10 &&
              !reentrant_hit.bounded_closest_certified,
          "reentrant seam edge or fail-closed certificate is incorrect");
  return {seam, contact_seam, reentrant_hit,
          index.torus_reentrant_pairs_audited, index.torus_reentrant_audit_ms};
}

void require_cached_matches_exhaustive(const LocalProjectionHit &cached,
                                       const LocalProjectionHit &exhaustive,
                                       const char *message) {
  require(cached.found == exhaustive.found &&
              cached.bounded_closest_certified ==
                  exhaustive.bounded_closest_certified &&
              cached.kind == exhaustive.kind &&
              cached.patch_index == exhaustive.patch_index,
          message);
  if (!cached.found)
    return;
  const double scale = std::max({1.0, cached.distance, exhaustive.distance});
  require(std::abs(cached.distance - exhaustive.distance) <= 2.0e-13 * scale &&
              vnorm(vsub(cached.point, exhaustive.point)) <= 2.0e-13 * scale &&
              vnorm(vsub(cached.normal, exhaustive.normal)) <= 2.0e-13,
          message);
}

PatchComplex make_contact_constraint_fixture(
    std::span<const SpherePatchConstraint> constraints) {
  PatchComplex complex;
  complex.probe_radius = 1.0;
  complex.status = "pass";
  complex.atoms.push_back({0.0, 0.0, 0.0, 0.0, 1.0});
  complex.atom_neighbors.emplace_back();
  constexpr double neighbor_distance = 3.0;
  constexpr double expanded_radius = 2.0;
  for (const SpherePatchConstraint &constraint : constraints) {
    const double neighbor_expanded = std::sqrt(
        expanded_radius * expanded_radius +
        neighbor_distance * neighbor_distance +
        2.0 * expanded_radius * neighbor_distance * constraint.minimum_dot);
    const std::size_t neighbor = complex.atoms.size();
    complex.atoms.push_back({-neighbor_distance * constraint.normal.x,
                             -neighbor_distance * constraint.normal.y,
                             -neighbor_distance * constraint.normal.z, 0.0,
                             neighbor_expanded - complex.probe_radius});
    complex.atom_neighbors[0].push_back(neighbor);
    complex.atom_neighbors.push_back({0U});
  }
  complex.contacts.push_back({0U, {0.0, 0.0, 0.0}, 1.0, expanded_radius});
  return complex;
}

void check_two_active_cache_against_exhaustive_oracle() {
  const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
  const std::array<SpherePatchConstraint, 4> contact_constraints = {
      SpherePatchConstraint{{1.0, 0.0, 0.0}, 0.5},
      SpherePatchConstraint{{0.0, 1.0, 0.0}, 0.5},
      SpherePatchConstraint{{0.0, 0.0, 1.0}, -0.95},
      SpherePatchConstraint{{inverse_sqrt_two, inverse_sqrt_two, 0.0}, 0.60}};
  const PatchComplex contact_complex =
      make_contact_constraint_fixture(contact_constraints);
  const LocalPatchProjectionIndex contact_index =
      build_local_patch_projection_index(contact_complex);
  require(contact_index.bounded_domains_certified &&
              contact_index.contact_two_active_offsets.size() == 2U &&
              contact_index.contact_two_active_offsets[1] >
                  contact_index.contact_two_active_offsets[0],
          "contact two-active cache fixture was not certified");

  const Vec3 corner_query{0.0, 0.0, -1.2};
  const LocalProjectionHit cached_corner =
      project_patch_bounded_closest(contact_complex, contact_index,
                                    PatchKind::contact, 0U, corner_query, 3.0);
  const LocalProjectionHit exhaustive_corner =
      project_patch_bounded_closest_exhaustive_for_audit(
          contact_complex, contact_index, PatchKind::contact, 0U, corner_query,
          3.0);
  const Vec3 expected_corner{0.5, 0.5, -std::sqrt(0.5)};
  require_cached_matches_exhaustive(
      cached_corner, exhaustive_corner,
      "contact cached corner differs from exhaustive active pairs");
  require(cached_corner.found &&
              vnorm(vsub(cached_corner.point, expected_corner)) < 2.0e-12,
          "known two-active contact corner was not selected");

  std::uint64_t state = 0x6a09e667f3bcc909ULL;
  const auto next_direction = [&]() {
    const auto next_component = [&]() {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      const double unit = static_cast<double>(state >> 11U) /
                          static_cast<double>(std::uint64_t{1} << 53U);
      return 2.0 * unit - 1.0;
    };
    Vec3 direction_{next_component(), next_component(), next_component()};
    if (vnorm(direction_) < 1.0e-12)
      direction_ = {1.0, 0.0, 0.0};
    return vnormalized(direction_);
  };
  for (std::size_t query_index = 0; query_index < 64U; ++query_index) {
    const Vec3 query =
        vscale(1.1 + 0.01 * static_cast<double>(query_index), next_direction());
    require_cached_matches_exhaustive(
        project_patch_bounded_closest(contact_complex, contact_index,
                                      PatchKind::contact, 0U, query, 3.0),
        project_patch_bounded_closest_exhaustive_for_audit(
            contact_complex, contact_index, PatchKind::contact, 0U, query, 3.0),
        "contact cache disagrees with exhaustive fixed-seed oracle");
  }

  SeamFixture reentrant_fixture = make_torus_reentrant_seam_fixture();
  const LocalPatchProjectionIndex reentrant_index =
      build_local_patch_projection_index(reentrant_fixture.complex);
  require(reentrant_index.bounded_domains_certified &&
              reentrant_index.reentrant_two_active_offsets.size() == 2U &&
              reentrant_index.reentrant_two_active_offsets[1] >
                  reentrant_index.reentrant_two_active_offsets[0],
          "reentrant triangle cache fixture was not certified");
  for (std::size_t query_index = 0; query_index < 64U; ++query_index) {
    const Vec3 query =
        vadd(reentrant_fixture.complex.reentrant[0].probe_center,
             vscale(1.2 + 0.005 * static_cast<double>(query_index),
                    next_direction()));
    require_cached_matches_exhaustive(
        project_patch_bounded_closest(reentrant_fixture.complex,
                                      reentrant_index, PatchKind::reentrant, 0U,
                                      query, 3.0),
        project_patch_bounded_closest_exhaustive_for_audit(
            reentrant_fixture.complex, reentrant_index, PatchKind::reentrant,
            0U, query, 3.0),
        "reentrant cache disagrees with exhaustive fixed-seed oracle");
  }
}

void check_empty_contact_constraint_fails_closed() {
  PatchComplex complex;
  complex.probe_radius = 1.0;
  complex.status = "pass";
  complex.atoms = {{0.0, 0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 0.0, 0.0, 2.0}};
  complex.atom_neighbors = {{1U}, {0U}};
  complex.contacts.push_back({0U, {0.0, 0.0, 0.0}, 1.0, 2.0});
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  require(!index.contact_constraints.empty() &&
              index.contact_constraints[0].minimum_dot == 2.0 &&
              index.bounded_domains_certified &&
              index.contact_domain_empty.size() == 1U &&
              index.contact_domain_empty[0] == 1U,
          "empty contact sentinel was not certified explicitly");
  const LocalProjectionHit hit =
      project_local_patch_closest(complex, index, 0U, {1.05, 0.0, 0.0}, 0.10);
  require(!hit.found, "certified-empty contact domain returned a projection");
}

void check_near_parallel_constraints_fail_closed() {
  for (const double angle : {5.0e-8, 1.0e-13}) {
    PatchComplex complex;
    complex.probe_radius = 1.0;
    complex.status = "pass";
    complex.atoms = {
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {3.0, 0.0, 0.0, 0.0, 1.0},
        {3.0 * std::cos(angle), 3.0 * std::sin(angle), 0.0, 0.0, 1.0}};
    complex.atom_neighbors = {{1U, 2U}, {0U}, {0U}};
    complex.contacts.push_back({0U, atom_position(complex.atoms[0]), 1.0, 2.0});
    const LocalPatchProjectionIndex index =
        build_local_patch_projection_index(complex);
    require(index.bvh_invariants_certified && !index.bounded_domains_certified,
            "near-parallel active constraints did not fail closed");
    const LocalProjectionHit hit = project_local_patch_closest(
        complex, index, 0U, {-1.05, 0.0, 0.0}, 0.10);
    require(hit.found && !hit.bounded_closest_certified &&
                !hit.traversal_invariant_failed,
            "near-parallel query did not preserve an uncertified result");
  }
}

void check_strict_nearest_is_owner_order_independent() {
  constexpr double nearest_distance = 0.05;
  constexpr double former_tie_tolerance = 2.0e-10;
  const std::array<double, 3> distances = {
      nearest_distance + 1.8 * former_tie_tolerance,
      nearest_distance + 0.9 * former_tie_tolerance, nearest_distance};
  PatchComplex complex;
  complex.probe_radius = 1.0;
  complex.status = "pass";
  for (std::size_t patch = 0; patch < distances.size(); ++patch) {
    complex.atoms.push_back({1.0 + distances[patch], 0.0, 0.0, 0.0, 1.0});
    complex.atom_neighbors.emplace_back();
    complex.contacts.push_back(
        {patch, atom_position(complex.atoms.back()), 1.0, 2.0});
  }
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  require(index.bounded_domains_certified && index.spatial_nodes.size() == 1U,
          "strict-nearest leaf fixture was not certified");

  std::array<LocalProjectionHit, 3> hits;
  for (std::size_t owner = 0; owner < hits.size(); ++owner)
    hits[owner] = project_local_patch_closest(complex, index, owner,
                                              {0.0, 0.0, 0.0}, 0.10);
  for (const LocalProjectionHit &hit : hits) {
    require(hit.found && hit.bounded_closest_certified &&
                hit.kind == PatchKind::contact && hit.patch_index == 2U,
            "epsilon-chain ordering did not select the strict nearest patch");
    require(
        hit.distance == hits[0].distance && hit.point.x == hits[0].point.x &&
            hit.point.y == hits[0].point.y && hit.point.z == hits[0].point.z,
        "owner hint changed strict nearest distance or point");
  }
}

void check_bvh_proof_data_fail_closed() {
  PatchComplex malformed;
  malformed.probe_radius = 1.0;
  malformed.status = "pass";
  malformed.atoms = {{0.0, 0.0, 0.0, 0.0, 1.0}};
  malformed.atom_neighbors.resize(1U);
  malformed.contacts.push_back({0U, atom_position(malformed.atoms[0]),
                                std::numeric_limits<double>::quiet_NaN(), 2.0});
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(malformed);
  require(!index.bvh_invariants_certified && !index.bounded_domains_certified &&
              index.spatial_nodes.empty(),
          "malformed BVH proof bounds were certified");
  const LocalProjectionHit hit =
      project_local_patch_closest(malformed, index, 0U, {1.0, 0.0, 0.0}, 0.10);
  require(!hit.bounded_closest_certified && hit.traversal_invariant_failed,
          "malformed BVH proof data did not fail closed at query");
}

void check_near_axis_torus_bounds_and_bvh() {
  constexpr double axis_minor = 1.0e-8;
  constexpr double circle_radius = 1.25e8;
  constexpr double axis_offset = 1.0e4;
  constexpr double probe_radius = 1.0;
  const Vec3 axis{1.0, axis_minor, 0.0};
  const Vec3 basis_u{-axis_minor, 1.0, 0.0};
  const Vec3 basis_v{0.0, 0.0, 1.0};
  const double expanded_radius = std::hypot(circle_radius, axis_offset);

  PatchComplex complex;
  complex.probe_radius = probe_radius;
  complex.status = "pass";
  complex.atoms = {{-axis_offset * axis.x, -axis_offset * axis.y, 0.0, 0.0,
                    expanded_radius - probe_radius},
                   {axis_offset * axis.x, axis_offset * axis.y, 0.0, 0.0,
                    expanded_radius - probe_radius}};
  complex.atom_neighbors = {{1U}, {0U}};
  complex.tori.push_back(
      {0U,
       1U,
       {0.0, 0.0, 0.0},
       axis,
       basis_u,
       basis_v,
       circle_radius,
       axis_offset,
       expanded_radius,
       expanded_radius,
       {{std::numbers::pi - 0.10, std::numbers::pi + 0.10}}});

  // Separated contacts make atom ownership stale and force the BVH path.
  // Their x upper bounds remain far below the near-axis torus maximum.
  constexpr std::size_t far_count = 16U;
  for (std::size_t i = 0; i < far_count; ++i) {
    const std::size_t atom = complex.atoms.size();
    complex.atoms.push_back(
        {-1000.0, 0.0, 1.0e6 * static_cast<double>(i + 1U), 0.0, 0.25});
    complex.atom_neighbors.emplace_back();
    complex.contacts.push_back(
        {atom, atom_position(complex.atoms.back()), 0.25, 1.25});
  }

  const Vec3 radial = vscale(-1.0, basis_u);
  const Vec3 torus_point = vscale(circle_radius - probe_radius, radial);
  const PatchBounds bounds = patch_bounds(complex, PatchKind::toroidal, 0U);
  require(
      torus_point.x > complex.tori[0].circle_center.x + probe_radius + 0.10 &&
          torus_point.x <= bounds.upper.x && torus_point.y >= bounds.lower.y &&
          torus_point.y <= bounds.upper.y,
      "near-axis torus AABB does not contain its analytic surface");

  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  require(index.bvh_invariants_certified && index.bounded_domains_certified,
          "near-axis torus BVH was not certified");
  const LocalProjectionHit hit = project_local_patch_closest(
      complex, index, complex.atoms.size() - 1U, torus_point, 1.0e-5);
  require(hit.found && hit.kind == PatchKind::toroidal &&
              hit.patch_index == 0U && hit.distance < 1.0e-7 &&
              hit.spatial_fallback && !hit.owner_hint_hit &&
              hit.bounded_closest_certified,
          "near-axis torus was incorrectly pruned by branch-and-bound");
}

void check_spindle_cusp_is_not_normal_certified() {
  constexpr double circle_radius = 0.50;
  const double atom_distance =
      2.0 * std::sqrt(4.0 - circle_radius * circle_radius);
  const std::vector<Atom> atoms = {{-0.5 * atom_distance, 0.0, 0.0, 0.0, 1.0},
                                   {0.5 * atom_distance, 0.0, 0.0, 0.0, 1.0}};
  const PatchComplex complex = build_patch_complex(atoms, 1.0, nullptr, false);
  require(complex.status == "pass" && complex.tori.size() == 1U &&
              complex.tori[0].circle_radius < complex.probe_radius,
          "spindle cusp fixture failed");
  const ToroidalPatch &torus = complex.tori[0];
  const double beta = std::acos(-torus.circle_radius / complex.probe_radius);
  const Vec3 cusp =
      vadd(torus.circle_center,
           vscale(complex.probe_radius * std::sin(beta), torus.axis));
  const LocalPatchProjectionIndex index =
      build_local_patch_projection_index(complex);
  const LocalProjectionHit hit =
      project_local_patch_closest(complex, index, 0U, cusp, 0.10);
  const LocalProjectionHit repeated =
      project_local_patch_closest(complex, index, 0U, cusp, 0.10);
  require(index.bounded_domains_certified && hit.found &&
              hit.kind == PatchKind::toroidal && hit.distance < 2.0e-12 &&
              !hit.normal_certified && !hit.bounded_closest_certified,
          "spindle cusp incorrectly received a certified normal");
  require(repeated.found && repeated.point.x == hit.point.x &&
              repeated.point.y == hit.point.y &&
              repeated.point.z == hit.point.z &&
              repeated.normal.x == hit.normal.x &&
              repeated.normal.y == hit.normal.y &&
              repeated.normal.z == hit.normal.z,
          "spindle cusp handling is not deterministic");
}

void check_shared_seam_certificate() {
  const std::vector<Atom> pair = {{-1.5, 0.0, 0.0, 0.0, 1.0},
                                  {1.5, 0.0, 0.0, 0.0, 1.0}};
  const PatchComplex pair_complex =
      build_patch_complex(pair, 1.0, nullptr, false);
  const AnalyticPatchSeamCertificate pair_certificate =
      build_analytic_patch_seam_certificate(pair_complex);
  require(pair_certificate.source_complex == &pair_complex &&
              pair_certificate.structure_certified &&
              pair_certificate.all_seams_certified &&
              pair_certificate.status == "pass" &&
              pair_certificate.torus_contact_curve_certified.size() == 1U &&
              pair_certificate.torus_contact_curve_certified[0][0] == 1U &&
              pair_certificate.torus_contact_curve_certified[0][1] == 1U &&
              pair_certificate.reentrant_edges.empty(),
          "shared pair beta/contact seam certificate failed");

  PatchComplex one_bad_side = pair_complex;
  one_bad_side.atoms[0].radius += 0.01;
  for (ContactPatch &contact : one_bad_side.contacts) {
    if (contact.atom == 0U) {
      contact.atom_radius += 0.01;
      contact.expanded_radius += 0.01;
    }
  }
  const AnalyticPatchSeamCertificate side_certificate =
      build_analytic_patch_seam_certificate(one_bad_side);
  require(side_certificate.structure_certified &&
              !side_certificate.all_seams_certified &&
              side_certificate.torus_contact_curve_certified.size() == 1U &&
              side_certificate.torus_contact_curve_certified[0][0] == 0U &&
              side_certificate.torus_contact_curve_certified[0][1] == 1U,
          "one invalid beta/contact side contaminated the other side");

  SeamFixture welded = make_torus_reentrant_seam_fixture();
  const AnalyticPatchSeamCertificate welded_certificate =
      build_analytic_patch_seam_certificate(welded.complex);
  require(welded_certificate.structure_certified &&
              welded_certificate.all_seams_certified &&
              welded_certificate.reentrant_edges.size() == 3U &&
              std::all_of(welded_certificate.reentrant_edges.begin(),
                          welded_certificate.reentrant_edges.end(),
                          [](const ReentrantEdgeSeamCertificate &edge) {
                            return edge.reentrant.valid() &&
                                   edge.torus.valid() &&
                                   edge.relation ==
                                       PatchRelation::curve_adjacent;
                          }),
          "shared torus/reentrant endpoint certificate failed");

  PatchComplex missing_endpoint = welded.complex;
  missing_endpoint.tori[0].exposed_intervals[0].begin = 0.10;
  const AnalyticPatchSeamCertificate missing_certificate =
      build_analytic_patch_seam_certificate(missing_endpoint);
  require(missing_certificate.structure_certified &&
              !missing_certificate.all_seams_certified &&
              std::any_of(missing_certificate.reentrant_edges.begin(),
                          missing_certificate.reentrant_edges.end(),
                          [](const ReentrantEdgeSeamCertificate &edge) {
                            return edge.relation == PatchRelation::uncertified;
                          }),
          "missing exposed endpoint did not fail the shared seam certificate");

  const LocalPatchProjectionIndex projection =
      build_local_patch_projection_index(welded.complex);
  const AnalyticPatchSeamCertificate summary =
      build_analytic_patch_seam_certificate(
          welded.complex, SeamCertificateStorage::summary_only,
          projection.atom_offsets, projection.patches);
  require(summary.source_complex == &welded.complex &&
              summary.all_seams_certified &&
              summary.torus_contact_curve_certified.empty() &&
              summary.reentrant_edges.empty() &&
              summary.retained_payload_bytes() == 0U,
          "summary-only seam certificate retained provenance payload");

  PatchComplex nan_probe = pair_complex;
  nan_probe.probe_radius = std::numeric_limits<double>::quiet_NaN();
  const AnalyticPatchSeamCertificate nan_probe_certificate =
      build_analytic_patch_seam_certificate(nan_probe);
  require(!nan_probe_certificate.structure_certified &&
              !nan_probe_certificate.all_seams_certified &&
              nan_probe_certificate.status != "pass",
          "NaN probe radius received a seam certificate");

  PatchComplex infinite_radius = pair_complex;
  infinite_radius.atoms[0].radius = std::numeric_limits<double>::infinity();
  const AnalyticPatchSeamCertificate infinite_radius_certificate =
      build_analytic_patch_seam_certificate(infinite_radius);
  require(!infinite_radius_certificate.structure_certified &&
              !infinite_radius_certificate.all_seams_certified &&
              infinite_radius_certificate.status != "pass",
          "infinite atom radius received a seam certificate");

  PatchComplex shifted_contact = pair_complex;
  shifted_contact.contacts[0].center.x += 1.0e-6;
  const AnalyticPatchSeamCertificate shifted_contact_certificate =
      build_analytic_patch_seam_certificate(shifted_contact);
  require(!shifted_contact_certificate.structure_certified &&
              !shifted_contact_certificate.all_seams_certified &&
              shifted_contact_certificate.status != "pass",
          "inconsistent contact centre borrowed a torus seam certificate");

  PatchComplex wrong_contact_radius = pair_complex;
  wrong_contact_radius.contacts[0].expanded_radius += 1.0e-6;
  const AnalyticPatchSeamCertificate wrong_contact_certificate =
      build_analytic_patch_seam_certificate(wrong_contact_radius);
  require(!wrong_contact_certificate.structure_certified &&
              !wrong_contact_certificate.all_seams_certified &&
              wrong_contact_certificate.status != "pass",
          "inconsistent contact radius borrowed a torus seam certificate");
}

void check_patch_adjacency_oracle_truth_table() {
  const std::vector<Atom> pair = {{-1.5, 0.0, 0.0, 0.0, 1.0},
                                  {1.5, 0.0, 0.0, 0.0, 1.0}};
  PatchComplex pair_complex = build_patch_complex(pair, 1.0, nullptr, false);
  require(pair_complex.status == "pass" && pair_complex.contacts.size() == 2U &&
              pair_complex.tori.size() == 1U,
          "adjacency pair fixture failed");
  const LocalPatchProjectionIndex pair_index =
      build_local_patch_projection_index(
          pair_complex, LocalPatchProjectionIndexBuildOptions{true});
  const AnalyticPatchAdjacencyOracle pair_oracle =
      build_analytic_patch_adjacency_oracle(pair_complex, pair_index);
  require(pair_oracle.certified() && pair_oracle.status() == "pass" &&
              pair_oracle.build_ms() >= 0.0 &&
              pair_oracle.retained_payload_bytes() > 0U,
          "adjacency pair oracle was not certified");
  const PatchRefId contact_zero = make_patch_ref(PatchKind::contact, 0U);
  const PatchRefId contact_one = make_patch_ref(PatchKind::contact, 1U);
  const PatchRefId torus_zero = make_patch_ref(PatchKind::toroidal, 0U);
  require(pair_oracle.relation(contact_zero, contact_zero) ==
                  PatchRelation::same_surface &&
              pair_oracle.relation(contact_zero, torus_zero) ==
                  PatchRelation::curve_adjacent &&
              pair_oracle.relation(torus_zero, contact_zero) ==
                  PatchRelation::curve_adjacent &&
              pair_oracle.relation(contact_zero, contact_one) ==
                  PatchRelation::separated,
          "pair relation truth table failed");
  const PatchRefId invalid;
  require(pair_oracle.relation(invalid, contact_zero) ==
              PatchRelation::uncertified,
          "invalid packed ref did not fail closed");

  const auto require_symmetric = [](const AnalyticPatchAdjacencyOracle &oracle,
                                    const std::vector<PatchRefId> &refs) {
    for (const PatchRefId first : refs) {
      for (const PatchRefId second : refs) {
        require(oracle.relation(first, second) ==
                    oracle.relation(second, first),
                "patch relation is not symmetric");
      }
    }
  };
  require_symmetric(pair_oracle, {contact_zero, contact_one, torus_zero});

  const PatchComplex copied_pair = pair_complex;
  const AnalyticPatchAdjacencyOracle stale_identity =
      build_analytic_patch_adjacency_oracle(copied_pair, pair_index);
  require(!stale_identity.certified() &&
              stale_identity.relation(contact_zero, contact_zero) ==
                  PatchRelation::same_surface &&
              stale_identity.relation(contact_zero, torus_zero) ==
                  PatchRelation::uncertified,
          "copied-complex index identity did not fail closed");

  PatchComplex mutated_pair = pair_complex;
  const LocalPatchProjectionIndex pre_mutation_index =
      build_local_patch_projection_index(
          mutated_pair, LocalPatchProjectionIndexBuildOptions{true});
  mutated_pair.contacts[0].center.x += 1.0e-6;
  const AnalyticPatchAdjacencyOracle stale_geometry =
      build_analytic_patch_adjacency_oracle(mutated_pair, pre_mutation_index);
  require(!stale_geometry.certified(),
          "mutated complex retained a stale spatial certificate");

  LocalPatchProjectionIndex duplicate_spatial = pair_index;
  require(duplicate_spatial.spatial_entries.size() >= 2U,
          "spatial corruption fixture is too small");
  duplicate_spatial.spatial_entries[1] = duplicate_spatial.spatial_entries[0];
  require(
      !build_analytic_patch_adjacency_oracle(pair_complex, duplicate_spatial)
           .certified(),
      "duplicate/missing spatial canonical id was accepted");

  LocalPatchProjectionIndex partial_root = pair_index;
  partial_root.spatial_nodes[0].begin = 1U;
  require(!build_analytic_patch_adjacency_oracle(pair_complex, partial_root)
               .certified(),
          "partial BVH root coverage was accepted");

  LocalPatchProjectionIndex malformed_constraints = pair_index;
  ++malformed_constraints.contact_constraint_offsets.back();
  require(!build_analytic_patch_adjacency_oracle(pair_complex,
                                                 malformed_constraints)
               .certified(),
          "malformed constraint CSR was accepted");

  PatchComplex duplicate_contact = pair_complex;
  duplicate_contact.contacts.push_back(duplicate_contact.contacts.front());
  const LocalPatchProjectionIndex duplicate_contact_index =
      build_local_patch_projection_index(
          duplicate_contact, LocalPatchProjectionIndexBuildOptions{true});
  require(!build_analytic_patch_adjacency_oracle(duplicate_contact,
                                                 duplicate_contact_index)
               .certified(),
          "duplicate atom contact patch was accepted");

  PatchComplex missing_contact = pair_complex;
  missing_contact.contacts.pop_back();
  const LocalPatchProjectionIndex missing_contact_index =
      build_local_patch_projection_index(
          missing_contact, LocalPatchProjectionIndexBuildOptions{true});
  require(!build_analytic_patch_adjacency_oracle(missing_contact,
                                                 missing_contact_index)
               .certified(),
          "missing atom contact patch was accepted");

  PatchComplex invalid_atom_ref = pair_complex;
  invalid_atom_ref.contacts[0].atom = invalid_atom_ref.atoms.size() + 7U;
  const LocalPatchProjectionIndex invalid_atom_ref_index =
      build_local_patch_projection_index(
          invalid_atom_ref, LocalPatchProjectionIndexBuildOptions{true});
  require(!invalid_atom_ref_index.bounded_domains_certified &&
              !invalid_atom_ref_index.bvh_invariants_certified &&
              invalid_atom_ref_index.atom_offsets.empty() &&
              !build_analytic_patch_adjacency_oracle(invalid_atom_ref,
                                                     invalid_atom_ref_index)
                   .certified(),
          "out-of-range patch atom ref did not fail closed");

  PatchComplex bad_side = pair_complex;
  bad_side.atoms[0].radius += 0.01;
  bad_side.contacts[0].atom_radius += 0.01;
  bad_side.contacts[0].expanded_radius += 0.01;
  const LocalPatchProjectionIndex bad_side_index =
      build_local_patch_projection_index(
          bad_side, LocalPatchProjectionIndexBuildOptions{true});
  const AnalyticPatchAdjacencyOracle bad_side_oracle =
      build_analytic_patch_adjacency_oracle(bad_side, bad_side_index);
  require(!bad_side_oracle.certified() &&
              bad_side_oracle.relation(contact_zero, contact_zero) ==
                  PatchRelation::same_surface &&
              bad_side_oracle.relation(contact_zero, torus_zero) ==
                  PatchRelation::uncertified,
          "uncertified seam leaked a cross-patch relation");

  const double root3 = std::sqrt(3.0);
  const std::vector<Atom> triple_atoms = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                          {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                          {0.0, root3, 0.0, 0.0, 1.0}};
  PatchComplex triple = build_patch_complex(triple_atoms, 1.0, nullptr, false);
  const LocalPatchProjectionIndex triple_index =
      build_local_patch_projection_index(
          triple, LocalPatchProjectionIndexBuildOptions{true});
  const AnalyticPatchAdjacencyOracle triple_oracle =
      build_analytic_patch_adjacency_oracle(triple, triple_index);
  require(triple.status == "pass" && triple.reentrant.size() >= 2U &&
              triple_oracle.certified(),
          "production triple adjacency oracle was not certified");
  const std::size_t reentrant_index = 0U;
  const PatchRefId reentrant_ref =
      make_patch_ref(PatchKind::reentrant, reentrant_index);
  const std::uint32_t atom =
      static_cast<std::uint32_t>(triple.reentrant[reentrant_index].atoms[0]);
  const PatchRefId contact_ref = make_patch_ref(PatchKind::contact, atom);
  const auto &seams = triple_oracle.seam_certificate();
  const PatchRefId first_torus =
      seams.reentrant_edges[reentrant_index * 3U].torus;
  const PatchRefId second_torus =
      seams.reentrant_edges[reentrant_index * 3U + 2U].torus;
  require(triple_oracle.relation(first_torus, reentrant_ref) ==
                  PatchRelation::curve_adjacent &&
              triple_oracle.relation(reentrant_ref, first_torus) ==
                  PatchRelation::curve_adjacent &&
              triple_oracle.relation(contact_ref, reentrant_ref) ==
                  PatchRelation::point_only &&
              triple_oracle.relation(first_torus, second_torus) ==
                  PatchRelation::uncertified,
          "triple seam/junction relation truth table failed");
  const PatchRefId opposite_reentrant =
      make_patch_ref(PatchKind::reentrant, 1U);
  require(triple_oracle.relation(reentrant_ref, opposite_reentrant) ==
              PatchRelation::uncertified,
          "tangent reentrant domains were certified without a bounded "
          "tangent-point proof");

  std::vector<PatchRefId> triple_refs;
  for (std::size_t i = 0U; i < triple.contacts.size(); ++i)
    triple_refs.push_back(make_patch_ref(PatchKind::contact, i));
  for (std::size_t i = 0U; i < triple.tori.size(); ++i)
    triple_refs.push_back(make_patch_ref(PatchKind::toroidal, i));
  for (std::size_t i = 0U; i < triple.reentrant.size(); ++i)
    triple_refs.push_back(make_patch_ref(PatchKind::reentrant, i));
  require_symmetric(triple_oracle, triple_refs);

  PatchComplex duplicate_reentrant = triple;
  duplicate_reentrant.reentrant.push_back(
      duplicate_reentrant.reentrant.front());
  const LocalPatchProjectionIndex duplicate_reentrant_index =
      build_local_patch_projection_index(
          duplicate_reentrant, LocalPatchProjectionIndexBuildOptions{true});
  const AnalyticPatchAdjacencyOracle duplicate_reentrant_oracle =
      build_analytic_patch_adjacency_oracle(duplicate_reentrant,
                                            duplicate_reentrant_index);
  require(duplicate_reentrant_oracle.certified() &&
              duplicate_reentrant_oracle.relation(
                  make_patch_ref(PatchKind::reentrant, 0U),
                  make_patch_ref(PatchKind::reentrant,
                                 duplicate_reentrant.reentrant.size() - 1U)) ==
                  PatchRelation::same_surface,
          "compatible coincident reentrant domains were not recognized");

  std::vector<Atom> separated_triples = triple_atoms;
  for (Atom atom_ : triple_atoms) {
    atom_.x += 20.0;
    separated_triples.push_back(atom_);
  }
  PatchComplex separated_complex =
      build_patch_complex(separated_triples, 1.0, nullptr, false);
  const LocalPatchProjectionIndex separated_index =
      build_local_patch_projection_index(
          separated_complex, LocalPatchProjectionIndexBuildOptions{true});
  const AnalyticPatchAdjacencyOracle separated_oracle =
      build_analytic_patch_adjacency_oracle(separated_complex, separated_index);
  require(separated_oracle.certified() &&
              separated_complex.reentrant.size() >= 4U,
          "separated reentrant fixture was not certified");
  std::size_t left_reentrant = 0U;
  std::size_t right_reentrant = 0U;
  for (std::size_t patch = 1U; patch < separated_complex.reentrant.size();
       ++patch) {
    if (separated_complex.reentrant[patch].probe_center.x <
        separated_complex.reentrant[left_reentrant].probe_center.x)
      left_reentrant = patch;
    if (separated_complex.reentrant[patch].probe_center.x >
        separated_complex.reentrant[right_reentrant].probe_center.x)
      right_reentrant = patch;
  }
  const PatchRefId left_ref =
      make_patch_ref(PatchKind::reentrant, left_reentrant);
  const PatchRefId right_ref =
      make_patch_ref(PatchKind::reentrant, right_reentrant);
  require(separated_oracle.relation(left_ref, right_ref) ==
                  PatchRelation::separated &&
              separated_oracle.relation(right_ref, left_ref) ==
                  PatchRelation::separated,
          "geometrically disjoint reentrant spheres were not separated");
}

} // namespace

int main() {
  try {
    const std::vector<Atom> pair = {{-1.5, 0.0, 0.0, 0.0, 1.0},
                                    {1.5, 0.0, 0.0, 0.0, 1.0}};
    const PatchComplex base = build_patch_complex(pair, 1.0, nullptr, false);
    require(base.status == "pass" && base.tori.size() == 1,
            "pair patch build failed");
    require(base.ray_buckets.empty(),
            "local projection fixture unexpectedly built ray buckets");
    compare_shifted_roots(
        base,
        build_patch_complex(translate(pair, 1.0e3, -17.0, 29.0), 1.0, nullptr,
                            false),
        1.0e3);
    compare_shifted_roots(
        base,
        build_patch_complex(translate(pair, 1.0e6, -17.0, 29.0), 1.0, nullptr,
                            false),
        1.0e6);

    const ReentrantProbeIndex index = build_reentrant_probe_index(base);
    const auto x = intersect_patch_axis(
        base, index, PatchKind::toroidal, 0, RayAxis::x,
        base.tori[0].circle_center.y + 0.37, base.tori[0].circle_center.z);
    const auto y = intersect_patch_axis(
        base, index, PatchKind::toroidal, 0, RayAxis::y,
        base.tori[0].circle_center.x + 0.37, base.tori[0].circle_center.z);
    const auto z = intersect_patch_axis(
        base, index, PatchKind::toroidal, 0, RayAxis::z,
        base.tori[0].circle_center.x + 0.37, base.tori[0].circle_center.y);
    require(x.certified && y.certified && z.certified,
            "three-axis torus roots were not certified");
    require(x.count > 0, "torus fixture has no x projection root");
    QuarticDiagnostics quartic_diagnostics;
    TorusPatchAxisDiagnostics filter_diagnostics;
    const PatchAxisRoots diagnostic_x =
        intersect_patch_axis_with_quartic_diagnostics(
            base, index, PatchKind::toroidal, 0U, RayAxis::x,
            base.tori[0].circle_center.y + 0.37, base.tori[0].circle_center.z,
            &quartic_diagnostics, &filter_diagnostics);
    require(diagnostic_x.count == x.count &&
                diagnostic_x.certified == x.certified &&
                filter_diagnostics.root_count ==
                    quartic_diagnostics.output_root_count,
            "torus filter diagnostics changed the intersection result");
    for (std::size_t root = 0U; root < x.count; ++root) {
      require(std::bit_cast<std::uint64_t>(diagnostic_x.coordinates[root]) ==
                      std::bit_cast<std::uint64_t>(x.coordinates[root]) &&
                  diagnostic_x.primitive_simplicity_certified[root] ==
                      x.primitive_simplicity_certified[root],
              "torus filter diagnostics changed retained root bits");
    }
    std::size_t audited_filtered_roots = 0U;
    for (int ui = -8; ui <= 8; ++ui) {
      for (int vi = -8; vi <= 8; ++vi) {
        QuarticDiagnostics scan_quartic;
        TorusPatchAxisDiagnostics scan_filters;
        (void)intersect_patch_axis_with_quartic_diagnostics(
            base, index, PatchKind::toroidal, 0U, RayAxis::x,
            base.tori[0].circle_center.y + 0.31 * static_cast<double>(ui),
            base.tori[0].circle_center.z + 0.31 * static_cast<double>(vi),
            &scan_quartic, &scan_filters);
        require(scan_filters.root_count == scan_quartic.output_root_count,
                "torus filter diagnostics lost a primitive root");
        for (std::size_t root = 0U; root < scan_filters.root_count; ++root) {
          if (scan_filters.roots[root].reason !=
              TorusRootFilterReason::accepted)
            ++audited_filtered_roots;
        }
      }
    }
    require(audited_filtered_roots > 0U,
            "torus filter diagnostic fixture exercised no rejected roots");

    constexpr std::size_t kNullPathCalls = 20000U;
    const auto time_axis_calls = [&](bool explicit_null_audit) {
      std::uint64_t checksum = 0U;
      const auto started = std::chrono::steady_clock::now();
      for (std::size_t call = 0U; call < kNullPathCalls; ++call) {
        const PatchAxisRoots roots =
            explicit_null_audit
                ? intersect_patch_axis_with_quartic_diagnostics(
                      base, index, PatchKind::toroidal, 0U, RayAxis::x,
                      base.tori[0].circle_center.y + 0.37,
                      base.tori[0].circle_center.z, nullptr, nullptr)
                : intersect_patch_axis(base, index, PatchKind::toroidal, 0U,
                                       RayAxis::x,
                                       base.tori[0].circle_center.y + 0.37,
                                       base.tori[0].circle_center.z);
        checksum ^= static_cast<std::uint64_t>(roots.count) + call;
        for (const double coordinate : roots.coordinates)
          checksum ^= std::bit_cast<std::uint64_t>(coordinate) +
                      0x9e3779b97f4a7c15ULL + (checksum << 6U) +
                      (checksum >> 2U);
      }
      return std::pair{std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - started)
                           .count(),
                       checksum};
    };
    (void)time_axis_calls(false);
    (void)time_axis_calls(true);
    std::array<double, 7> axis_normal_samples{};
    std::array<double, 7> axis_null_samples{};
    for (std::size_t sample = 0U; sample < axis_normal_samples.size();
         ++sample) {
      const auto [normal_ms, normal_checksum] = time_axis_calls(false);
      const auto [null_ms, null_checksum] = time_axis_calls(true);
      require(normal_checksum == null_checksum,
              "null patch-axis diagnostics changed result checksum");
      axis_normal_samples[sample] = normal_ms;
      axis_null_samples[sample] = null_ms;
    }
    std::sort(axis_normal_samples.begin(), axis_normal_samples.end());
    std::sort(axis_null_samples.begin(), axis_null_samples.end());
    const double axis_normal_ms =
        axis_normal_samples[axis_normal_samples.size() / 2U];
    const double axis_null_ms =
        axis_null_samples[axis_null_samples.size() / 2U];
    const LocalPatchProjectionIndex local =
        build_local_patch_projection_index(base);
    require(local.retained_bytes() > 0,
            "local projection index retained no compact data");
    const Vec3 query{x.coordinates[0] + 0.05,
                     base.tori[0].circle_center.y + 0.37,
                     base.tori[0].circle_center.z};
    const LocalProjectionHit projected =
        project_local_patch_x(base, local, 0, query, 0.10);
    const LocalProjectionHit repeated =
        project_local_patch_x(base, local, 0, query, 0.10);
    const LocalProjectionHit closest =
        project_local_patch_closest(base, local, 0, query, 0.10);
    require(projected.found &&
                std::abs(projected.point.x - x.coordinates[0]) < 1.0e-11,
            "atom-local projection missed the analytic torus root");
    require(repeated.found && repeated.kind == projected.kind &&
                repeated.patch_index == projected.patch_index &&
                repeated.point.x == projected.point.x,
            "atom-local projection is not deterministic");
    require(closest.found && closest.distance <= projected.axis_distance &&
                std::hypot(closest.normal.x, closest.normal.y,
                           closest.normal.z) > 0.9,
            "closest-point projection missed the local analytic patch");
    require(closest.bounded_closest_certified && local.bvh_invariants_certified,
            "base closest projection was not certified");
    require(closest.owner_hint_candidates == 2U,
            "production owner-only patch hints were not deduplicated");
    const PatchComplex shifted_closest_complex = build_patch_complex(
        translate(pair, 1.0e6, -17.0, 29.0), 1.0, nullptr, false);
    const LocalPatchProjectionIndex shifted_closest_index =
        build_local_patch_projection_index(shifted_closest_complex);
    const Vec3 shifted_query{query.x + 1.0e6, query.y - 17.0, query.z + 29.0};
    const LocalProjectionHit shifted_closest = project_local_patch_closest(
        shifted_closest_complex, shifted_closest_index, 0U, shifted_query,
        0.10);
    require(shifted_closest.found &&
                shifted_closest.bounded_closest_certified &&
                shifted_closest.kind == closest.kind &&
                shifted_closest.patch_index == closest.patch_index &&
                std::abs((shifted_closest.point.x - 1.0e6) - closest.point.x) <
                    2.0e-9 &&
                std::abs((shifted_closest.point.y + 17.0) - closest.point.y) <
                    2.0e-9 &&
                std::abs((shifted_closest.point.z - 29.0) - closest.point.z) <
                    2.0e-9 &&
                std::abs(shifted_closest.distance - closest.distance) < 2.0e-9,
            "outward-padded BVH changed closest projection at 1e6 shift");
    const LocalProjectionHit nonfinite = project_local_patch_closest(
        base, local, 0U,
        {std::numeric_limits<double>::quiet_NaN(), query.y, query.z}, 0.10);
    require(!nonfinite.found && nonfinite.candidate_patches == 0U,
            "non-finite closest query did not fail closed");
    const PatchComplex bucket_reference =
        build_patch_complex(pair, 1.0, nullptr, true);
    const std::vector<RayHit> reference_hits =
        intersect_ray_x(bucket_reference, query.y, query.z);
    require(!reference_hits.empty(),
            "bucket-reference projection has no analytic hit");
    const auto closest_reference = std::min_element(
        reference_hits.begin(), reference_hits.end(),
        [&](const RayHit &a, const RayHit &b) {
          return std::abs(a.x - query.x) < std::abs(b.x - query.x);
        });
    require(std::abs(closest_reference->x - projected.point.x) < 1.0e-9,
            "compact local and bucket-reference projections disagree");

    const LocalProjectionHit stale_owner =
        check_spatial_projection_outside_owner_hint();
    const ConstrainedSeamHits seam_hits = check_constrained_seam_minima();
    check_two_active_cache_against_exhaustive_oracle();
    check_empty_contact_constraint_fails_closed();
    check_near_parallel_constraints_fail_closed();
    check_strict_nearest_is_owner_order_independent();
    check_bvh_proof_data_fail_closed();
    check_near_axis_torus_bounds_and_bvh();
    check_spindle_cusp_is_not_normal_certified();
    check_shared_seam_certificate();
    check_patch_adjacency_oracle_truth_table();
    const LargeIndexTiming large_index = measure_large_index_setup();

    std::cout << "tabi_ses_patch_intersections_smoke=pass x=" << x.count
              << " y=" << y.count << " z=" << z.count
              << " projection_bytes=" << local.retained_bytes()
              << " owner_hint_candidates=" << closest.owner_hint_candidates
              << " stale_owner_candidates=" << stale_owner.candidate_patches
              << " stale_owner_bound_pruned="
              << stale_owner.bound_pruned_patches
              << " stale_owner_nodes=" << stale_owner.spatial_nodes_visited
              << " seam_candidates="
              << seam_hits.torus_reentrant.candidate_patches
              << " contact_seam_distance=" << seam_hits.contact_torus.distance
              << " torus_reentrant_audit_pairs="
              << seam_hits.torus_reentrant_pairs_audited
              << " torus_reentrant_audit_ms="
              << seam_hits.torus_reentrant_audit_ms << " reentrant_fail_closed="
              << (!seam_hits.reentrant_only.bounded_closest_certified)
              << " large_patches=4096"
              << " large_patch_build_ms=" << large_index.patch_build_ms
              << " large_index_build_ms=" << large_index.index_build_ms
              << " large_index_bytes=" << large_index.retained_bytes
              << " audited_filtered_roots=" << audited_filtered_roots
              << " patch_axis_normal_ms_median=" << axis_normal_ms
              << " patch_axis_null_ms_median=" << axis_null_ms
              << " patch_axis_null_ratio="
              << (axis_normal_ms > 0.0 ? axis_null_ms / axis_normal_ms : 0.0)
              << " stamp_bytes_per_worker="
              << large_index.stamp_bytes_per_worker << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_patch_intersections_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
