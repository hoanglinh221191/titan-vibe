// Checks of the curved PN lift and the near-pair quadrature on small
// synthetic cases (RSOS revision; the counterexamples of the code review):
//   fold      random PN triangles accepted by the lift must keep their
//             orientation everywhere (dense lattice), and a known patch that
//             folds between the 6x6 lattice nodes must be refused;
//   crossing  two patches without a shared vertex that bend into each other
//             (flat triangles 0.1 A apart) must be repaired apart by their
//             vertex normals;
//   self      a vertex whose own analytic normal makes its patches steep:
//             the repair replaces it, every patch is then tangent at the
//             vertex to one plane, and the self integral with the repaired
//             collocation normal converges with the order;
//   fallback  a fan with a triangle turned against its vertices: the lift
//             must refuse it;
//   quadrature the double-layer coefficient of a flat dual cell against the
//             solid-angle value, for targets 0.1 A down to 1e-14 A above it:
//             resolved within 1e-8 down to 1e-10 A, reported unresolved below.
// Built against the translation unit itself:
//   clang++ -std=c++23 -O2 -march=native -Isrc test/tabi_pn_certificate_probe.cpp
//     build/native/src/titan_simd_dispatch.o ... -o pn_certificate_probe.exe
#include "../src/tabi/tabi_node_patch_curved.cpp"
#include <chrono>
#include <iomanip>
#include <random>
using namespace TitanPBE::Tabi;
namespace {
int failures = 0;
void check(bool ok, const std::string &what) {
  if (!ok) {
    ++failures;
    std::cout << "FAIL " << what << "\n";
  }
}
void fill_incident(NodePatchGeometry &g) {
  g.incident_offsets = {0};
  g.incident_triangles.clear();
  for (std::size_t v = 0; v < g.vertex_xyz.size(); ++v) {
    for (std::size_t t = 0; t < g.triangles.size(); ++t)
      for (auto a : g.triangles[t].vertices)
        if (a == v)
          g.incident_triangles.push_back(t);
    g.incident_offsets.push_back(g.incident_triangles.size());
  }
}
double min_jacobian(const PnPatch &p, int lattice) {
  double lowest = std::numeric_limits<double>::infinity();
  for (int i = 0; i <= lattice; ++i)
    for (int j = 0; i + j <= lattice; ++j) {
      const double a = double(i) / lattice, b = double(j) / lattice;
      V x, d1, d2, d3;
      pn_evaluate(p, a, b, 1 - a - b, x, d1, d2, d3);
      lowest = std::min(lowest, cross(sub(d1, d3), sub(d2, d3))[2]);
    }
  return lowest;
}
void fold() {
  std::mt19937_64 rng(991710);
  std::uniform_real_distribution<double> u(0, 1);
  int accepted = 0, folded_accepted = 0, lattice_only_accepted = 0;
  for (int it = 0; it < 200000; ++it) {
    NodePatchGeometry g;
    g.vertex_xyz = {{0, 0, 0}, {1, 0, 0},
                    {2 * u(rng) - .5, std::pow(10., -2 + 3 * u(rng)), 0}};
    NodePatchTriangle t{};
    t.vertices = {0, 1, 2};
    t.nz = 1;
    g.triangles = {t};
    std::vector<std::array<double, 3>> n;
    for (int j = 0; j < 3; ++j) {
      const double z = .3001 + .6999 * u(rng), a = 2 * pi * u(rng),
                   r = std::sqrt(1 - z * z);
      n.push_back({r * std::cos(a), r * std::sin(a), z});
    }
    const PnPatch p = pn_triangle_patch(n, g, 0, 0);
    if (min_jacobian(p, 6) > 0)
      ++lattice_only_accepted;
    if (pn_triangle_steep(n, g, 0, 0))
      continue;
    ++accepted;
    if (min_jacobian(p, 120) < 0)
      ++folded_accepted;
  }
  std::cout << "FOLD trials=200000 accepted=" << accepted
            << " positive_on_6x6_lattice=" << lattice_only_accepted
            << " accepted_but_folded_on_120x120=" << folded_accepted << "\n";
  check(folded_accepted == 0, "an accepted patch folds");
  // The review's patch: accepted by the 6x6 lattice, folded at (1/30, 29/30, 0).
  NodePatchGeometry g;
  g.vertex_xyz = {{0, 0, 0}, {1, 0, 0}, {0.15264677138308202, 0.014783066384448002, 0}};
  NodePatchTriangle t{};
  t.vertices = {0, 1, 2};
  t.nz = 1;
  g.triangles = {t};
  const std::vector<std::array<double, 3>> n = {
      {-0.68671764920558143, -0.36444861301252346, 0.62896429051481761},
      {0.92191028271830677, -0.20156614736789535, 0.33083609061514979},
      {0.4000067371771997, 0.036691439037907163, 0.91577745577961145}};
  const bool steep = pn_triangle_steep(n, g, 0, 0);
  std::cout << "FOLD review_patch steep=" << steep
            << " min_jacobian_120=" << min_jacobian(pn_triangle_patch(n, g, 0, 0), 120)
            << "\n";
  check(steep, "the review's folded patch is accepted");
}
void crossing() {
  NodePatchGeometry g;
  std::vector<std::array<double, 3>> n, star;
  for (int top = 0; top < 2; ++top)
    for (int k = 0; k < 3; ++k) {
      const double a = 2 * pi * k / 3;
      g.vertex_xyz.push_back({std::cos(a), std::sin(a), top * .1});
      n.push_back({.6 * std::cos(a), .6 * std::sin(a), top ? -.8 : .8});
      star.push_back({0, 0, top ? -1.0 : 1.0});
    }
  NodePatchTriangle a{}, b{};
  a.vertices = {0, 1, 2};
  b.vertices = {3, 5, 4};
  a.nz = 1;
  b.nz = -1;
  g.triangles = {a, b};
  fill_incident(g);
  std::size_t before = 0, after = 0;
  const std::vector<std::uint8_t> zero(2, 0U);
  const auto crossing_before = pn_unproven_pairs(n, zero, zero, g, before);
  PnRepairReport report;
  pn_repair_normals(n, star, g, report);
  const auto unproven = pn_unproven_pairs(n, zero, zero, g, after);
  std::cout << "CROSSING unproven_before=" << crossing_before.size()
            << " star_normals=" << report.star << " turned=" << report.turned
            << " rounds=" << report.rounds << " unproven_after=" << unproven.size()
            << " flat_contacts=" << after << "\n";
  check(!crossing_before.empty() && unproven.empty(), "crossing patches kept");
}
void self_term() {
  const double R = 1.6, cap = 0.55;
  NodePatchGeometry g;
  std::vector<std::array<double, 3>> n;
  g.vertex_xyz.push_back({0, 0, R});
  n.push_back({0, 0, 1});
  for (int k = 0; k < 6; ++k) {
    const double phi = 2 * pi * k / 6 + 0.15 * std::sin(3.0 * k),
                 th = cap * (1 + 0.12 * std::cos(2.0 * k));
    const V p{R * std::sin(th) * std::cos(phi), R * std::sin(th) * std::sin(phi),
              R * std::cos(th)};
    g.vertex_xyz.push_back(p);
    n.push_back(mul(p, 1 / R));
  }
  for (int k = 0; k < 6; ++k) {
    NodePatchTriangle t{};
    t.vertices = {0, std::size_t(1 + k), std::size_t(1 + (k + 1) % 6)};
    V f = cross(sub(g.vertex_xyz[1 + k], g.vertex_xyz[0]),
                sub(g.vertex_xyz[1 + (k + 1) % 6], g.vertex_xyz[0]));
    f = mul(f, 1 / norm(f));
    t.nx = f[0];
    t.ny = f[1];
    t.nz = f[2];
    g.triangles.push_back(t);
  }
  fill_incident(g);
  std::vector<std::array<double, 3>> star(7, V{});
  for (const auto &t : g.triangles) {
    const double area = norm(cross(sub(g.vertex_xyz[t.vertices[1]], g.vertex_xyz[t.vertices[0]]),
                                   sub(g.vertex_xyz[t.vertices[2]], g.vertex_xyz[t.vertices[0]])));
    for (auto v : t.vertices)
      star[v] = add(star[v], mul(V{t.nx, t.ny, t.nz}, area));
  }
  for (auto &s : star)
    s = mul(s, 1 / norm(s));
  // The vertex itself carries an analytic normal 75 deg off, so the repair
  // must replace the collocation normal the self term is taken with.
  n[0] = {std::sin(75 * pi / 180), 0, std::cos(75 * pi / 180)};
  const V analytic = n[0];
  int steep_before = 0;
  for (std::size_t ti = 0; ti < 6; ++ti)
    steep_before += pn_triangle_steep(n, g, ti, 0U);
  PnRepairReport report;
  pn_repair_normals(n, star, g, report);
  const bool replaced = norm(sub(n[0], analytic)) > 1e-3;
  double tilt = 0.0;
  for (std::size_t ti = 0; ti < 6; ++ti) {
    const PnPatch p = pn_triangle_patch(n, g, ti, 0U);
    V m = cross(sub(p.b[3], p.b[0]), sub(p.b[8], p.b[0]));
    m = mul(m, 1 / norm(m));
    tilt = std::max(tilt, std::acos(std::min(1.0, std::abs(dot(m, n[0])))) * 180 / pi);
  }
  // The production target: the vertex with the normal the repair left there.
  Panel target;
  target.z = R;
  target.nx = n[0][0];
  target.ny = n[0][1];
  target.nz = n[0][2];
  double previous = 0.0, change = 0.0;
  for (int order : {16, 32, 64, 128}) {
    const B v = integrate(target, pn_samples(n, {}, g, 0, order), 80.0, 0.1);
    if (order > 16)
      change = std::max(change, std::abs(v[2] - previous) / std::abs(v[2]));
    previous = v[2];
  }
  std::cout << "SELF steep_before=" << steep_before << " star_normals=" << report.star
            << " vertex_normal_replaced=" << replaced
            << " largest_tangent_tilt_deg=" << tilt
            << " K_tgt_change_16_to_128=" << change << "\n";
  check(steep_before > 0 && replaced && tilt < 1e-9 && change < 1e-8,
        "self term at a repaired vertex");
}
// A fan with one triangle turned against its vertices (a fold of the flat
// mesh): no normal per vertex mends it.  It must come back flat, its vertex
// a corner, and the corner cell must collocate at a smooth point inside its
// host patch, where its self integral converges with the order; at the
// vertex itself the same integral drifts with the order.
void fallback() {
  const double R = 1.6;
  NodePatchGeometry g;
  std::vector<std::array<double, 3>> n;
  g.vertex_xyz.push_back({0, 0, R});
  n.push_back({0, 0, 1});
  for (int k = 0; k < 6; ++k) {
    const double phi = 2 * pi * k / 6, th = 0.5;
    const V p{R * std::sin(th) * std::cos(phi), R * std::sin(th) * std::sin(phi),
              R * std::cos(th)};
    g.vertex_xyz.push_back(p);
    n.push_back(mul(p, 1 / R));
  }
  for (int k = 0; k < 6; ++k) {
    NodePatchTriangle t{};
    t.vertices = {0, std::size_t(1 + k), std::size_t(1 + (k + 1) % 6)};
    V f = cross(sub(g.vertex_xyz[1 + k], g.vertex_xyz[0]),
                sub(g.vertex_xyz[1 + (k + 1) % 6], g.vertex_xyz[0]));
    f = mul(f, (k == 2 ? -1.0 : 1.0) / norm(f));
    t.nx = f[0];
    t.ny = f[1];
    t.nz = f[2];
    g.triangles.push_back(t);
  }
  fill_incident(g);
  std::vector<std::array<double, 3>> star(7, V{});
  for (const auto &t : g.triangles)
    for (auto v : t.vertices)
      star[v] = add(star[v], V{t.nx, t.ny, t.nz});
  for (auto &s : star)
    s = mul(s, 1 / norm(s));
  PnRepairReport report;
  const auto residual = pn_repair_normals(n, star, g, report);
  PnLiftReport lift;
  const auto masks = residual.empty() ? std::vector<std::uint8_t>{}
                                      : pn_flatten_residual(n, g, residual, lift);
  // The vertex is a corner when an edge through it is straight.
  bool corner = false;
  std::size_t host = 0;
  double best = -1.0;
  for (auto pos = g.incident_offsets[0]; pos < g.incident_offsets[1]; ++pos) {
    const auto ti = g.incident_triangles[pos];
    corner = corner || (!masks.empty() && ((masks[ti] & 1U) || (masks[ti] & 4U)));
    const auto &t = g.triangles[ti];
    const double area = norm(cross(sub(g.vertex_xyz[t.vertices[1]], g.vertex_xyz[t.vertices[0]]),
                                   sub(g.vertex_xyz[t.vertices[2]], g.vertex_xyz[t.vertices[0]])));
    if (area > best) {
      best = area;
      host = pos;
    }
  }
  CurvedNodePatchGeometry c;
  c.pn_node_normals = n;
  c.pn_straight_edges = masks;
  c.cells.resize(7);
  c.cells[0].offset = true;
  c.cells[0].host_pos = host;
  const PnQuad q = pn_quad(n, masks, g, 0, host, 0);
  V P, normal;
  double area;
  pn_quad_point(q, 0.5, 1.0, P, normal, area);
  Panel off;
  off.x = P[0]; off.y = P[1]; off.z = P[2];
  off.nx = normal[0]; off.ny = normal[1]; off.nz = normal[2];
  Panel vertex;
  vertex.z = R;
  vertex.nx = n[0][0]; vertex.ny = n[0][1]; vertex.nz = n[0][2];
  double previous = 0.0, change = 0.0, at_vertex_prev = 0.0, at_vertex_change = 0.0;
  for (int order : {16, 32, 64}) {
    const B v = offset_self_block(off, g, c, 0, order, 80.0, 0.1);
    const B w = integrate(vertex, pn_samples(n, masks, g, 0, order), 80.0, 0.1);
    if (order > 16) {
      change = std::max(change, std::abs(v[2] - previous) / std::abs(v[2]));
      at_vertex_change = std::max(at_vertex_change, std::abs(w[2] - at_vertex_prev) / std::abs(w[2]));
    }
    previous = v[2];
    at_vertex_prev = w[2];
  }
  std::cout << "FALLBACK folded_fan residual=" << std::count(residual.begin(), residual.end(), 1U)
            << " flat=" << lift.flattened << " corner=" << corner
            << " offset_A=" << norm(sub(P, g.vertex_xyz[0]))
            << " K_tgt_change_16_to_64 off_vertex=" << change
            << " at_vertex=" << at_vertex_change << "\n";
  check(!residual.empty() && corner && change < 1e-6 && at_vertex_change > 1e-3,
        "corner cell off-vertex self term");
  // The production order sequence: off the vertex the block is accepted; at
  // the vertex it is still moving at 128, which the lift refuses rather than
  // keeping the order-128 value.
  const SelfBlock accepted = converged_self_block(off, g, c, 0, false, 16, 80.0, 0.1);
  c.cells[0].offset = false;
  const SelfBlock refused = converged_self_block(vertex, g, c, 0, false, 16, 80.0, 0.1);
  std::cout << "SELF_SEQUENCE off_vertex converged=" << accepted.converged
            << " order=" << accepted.order << " change=" << accepted.change
            << " at_vertex converged=" << refused.converged
            << " order=" << refused.order << " change=" << refused.change << "\n";
  check(accepted.converged && !refused.converged && refused.order == 128,
        "self order sequence refuses the corner vertex");
}
// A large kernel must not hide the quadrature error of a small one: K1 = 1
// fixed, K3 = 1e-4 moving by 1e-8 (1e-4 of itself) is not converged, though
// it moved by only 1e-8 of the block's largest coefficient, and however
// large the other kernels' rounding sizes.  A coefficient that cancels to
// rounding - terms of size 1 summing to 1e-15 - is judged against its own
// rounding size and accepted.
void self_criterion() {
  const auto value = [](B v, B size) { return SelfValue{v, size}; };
  const B plain{1.0, 0.5, 1e-4, 0.25};
  double change = 0.0;
  const bool masked =
      self_orders_agree(value({1.0, 0.5, 1e-4 + 1e-8, 0.25}, {1.0, 1e6, 1e-4, 1e6}),
                        value({1.0, 0.5, 1e-4, 0.25}, {1.0, 1e6, 1e-4, 1e6}), change);
  const double masked_change = change;
  const bool rounding =
      self_orders_agree(value({1.0, 0.5, 2e-15, 0.25}, {1.0, 0.5, 1.0, 0.25}),
                        value({1.0, 0.5, -1e-15, 0.25}, {1.0, 0.5, 1.0, 0.25}), change);
  const bool moved = self_orders_agree(value({1.0 + 1e-7, 0.5, 1e-4, 0.25}, plain),
                                       value(plain, plain), change);
  std::cout << "SELF_CRITERION small_kernel_moving_accepted=" << masked
            << " change=" << masked_change << " rounding_accepted=" << rounding
            << " large_kernel_1e-7_accepted=" << moved << "\n";
  check(!masked && masked_change > 1e-5 && rounding && moved,
        "self block criterion per kernel");
}
double solid_triangle(V a, V b, V c, V p) {
  a = sub(a, p);
  b = sub(b, p);
  c = sub(c, p);
  return 2 * std::atan2(dot(a, cross(b, c)),
                        norm(a) * norm(b) * norm(c) + dot(a, b) * norm(c) +
                            dot(b, c) * norm(a) + dot(c, a) * norm(b));
}
void quadrature() {
  NodePatchGeometry g;
  g.vertex_xyz = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  NodePatchTriangle t{};
  t.vertices = {0, 1, 2};
  t.nz = 1;
  g.triangles = {t};
  fill_incident(g);
  CurvedNodePatchGeometry c;
  c.pn_node_normals = {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
  c.pn_straight_edges = {0};
  c.cells.resize(3);
  c.one_point_beyond_radii = 6;
  c.cells[0].radius = .5;
  c.cells[0].samples = pn_samples(c.pn_node_normals, c.pn_straight_edges, g, 0, 2);
  for (double gap : {.1, .01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14}) {
    Panel p;
    p.x = .2;
    p.y = .1;
    p.z = gap;
    p.nz = 1;
    const V pos{p.x, p.y, p.z};
    const double omega =
        solid_triangle({0, 0, 0}, {.5, 0, 0}, {1. / 3, 1. / 3, 0}, pos) +
        solid_triangle({0, 0, 0}, {1. / 3, 1. / 3, 0}, {0, .5, 0}, pos);
    const double reference = -(80. - 1) * omega / (4 * pi);
    const std::size_t unresolved_before = c.near_unresolved_evaluations.load();
    const auto start = std::chrono::steady_clock::now();
    B value{};
    bool rejected = false;
    try {
      value = geometric_pair(p, g, c, 0, 80, 0);
    } catch (const std::runtime_error &e) {
      rejected = std::string(e.what()) == "curved_near_quadrature_unresolved";
      if (!rejected) throw;
    }
    const double ms = std::chrono::duration<double, std::milli>(
                          std::chrono::steady_clock::now() - start)
                          .count();
    const bool unresolved = c.near_unresolved_evaluations.load() > unresolved_before;
    const double error = std::abs(value[1] / reference - 1);
    std::cout << "QUAD gap_A=" << gap << " reference=" << reference
              << " value=" << value[1] << " relative_error=" << error
              << " unresolved=" << unresolved << " rejected=" << rejected
              << " ms=" << ms << "\n";
    if (gap >= 1e-10)
      check(!rejected && !unresolved && error < 1e-8,
            "quadrature at gap " + std::to_string(gap));
    else
      check((rejected && unresolved) || error < 1e-8,
            "unresolved integral escaped at gap " + std::to_string(gap));
  }
}
} // namespace
int main() {
  std::cout << std::setprecision(10);
  try {
    fold();
    crossing();
    self_term();
    fallback();
    self_criterion();
    quadrature();
  } catch (const std::exception &e) {
    std::cout << "EXCEPTION " << e.what() << "\n";
    return 2;
  }
  std::cout << (failures ? "FAILED " : "PASSED ") << failures << "\n";
  return failures ? 1 : 0;
}
