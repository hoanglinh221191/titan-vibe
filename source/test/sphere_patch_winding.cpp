// Constant-density Laplace double layer at charge locations. Test-only:
// production centroid projection versus a converged curved-boundary rule.
#include "tabi/tabi_types.h"
#include "titan_executor.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>
using namespace TitanPBE::Tabi;
namespace {
using V = std::array<double, 3>;
using Tri = std::array<V, 3>;
constexpr double pi = std::numbers::pi;
V add(V a, V b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V sub(V a, V b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V mul(V a, double s) { return {s * a[0], s * a[1], s * a[2]}; }
double dot(V a, V b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V cross(V a, V b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
double norm(V a) { return std::sqrt(dot(a, a)); }
V unit(V a) { return mul(a, 1 / norm(a)); }
V xyz(const Atom &a) { return {a.x, a.y, a.z}; }
using Rule = std::vector<std::array<double, 2>>;
Rule gauss(int order) {
  Rule r;
  for (int i = 0; i < order; ++i) {
    double z = std::cos(pi * (i + .75) / (order + .5)), derivative = 0;
    for (int it = 0; it < 64; ++it) {
      double p = 1, previous = 0;
      for (int k = 1; k <= order; ++k) {
        double next = ((2 * k - 1) * z * p - (k - 1) * previous) / k;
        previous = p;
        p = next;
      }
      derivative = order * (z * p - previous) / (z * z - 1);
      double step = p / derivative;
      z -= step;
      if (std::abs(step) < 2e-15)
        break;
    }
    r.push_back({(1 + z) / 2, 1 / ((1 - z * z) * derivative * derivative)});
  }
  return r;
}
struct Shape {
  Tri v;
  V center;
  double radius;
};
Shape shape(Tri v) {
  if (dot(cross(sub(v[1], v[0]), sub(v[2], v[0])), add(add(v[0], v[1]), v[2])) <
      0)
    std::swap(v[1], v[2]);
  const V center = unit(add(add(v[0], v[1]), v[2]));
  return {v, center,
          std::max({norm(sub(v[0], center)), norm(sub(v[1], center)),
                    norm(sub(v[2], center))})};
}
Shape shape(const Panel &p, const Atom &a) {
  const V c = xyz(a);
  return shape({unit(sub({p.v0x, p.v0y, p.v0z}, c)),
                unit(sub({p.v1x, p.v1y, p.v1z}, c)),
                unit(sub({p.v2x, p.v2y, p.v2z}, c))});
}
// Solid angle is the circulation of a monopole vector potential. The gauge
// pole is local to the patch. A conservative spherical-cap bound guarantees
// the entire patch is in its nonsingular hemisphere. Otherwise partition the
// integration domain (not the production mesh) and repeat.
double curved(const Shape &s, V x, const Rule &r, int depth = 0) {
  const V center_delta = sub(s.center, x);
  const double d = norm(center_delta);
  if (d <= 1.001 * s.radius) {
    if (depth == 12)
      throw std::runtime_error("winding observation too close to surface");
    const auto &v = s.v;
    const V ab = unit(add(v[0], v[1])), bc = unit(add(v[1], v[2])),
            ca = unit(add(v[2], v[0]));
    return curved(shape({v[0], ab, ca}), x, r, depth + 1) +
           curved(shape({v[1], bc, ab}), x, r, depth + 1) +
           curved(shape({v[2], ca, bc}), x, r, depth + 1) +
           curved(shape({ab, bc, ca}), x, r, depth + 1);
  }
  const V pole = mul(center_delta, 1 / d);
  long double omega = 0;
  for (int e = 0; e < 3; ++e) {
    const V a = s.v[e], edge = sub(s.v[(e + 1) % 3], a);
    for (auto q : r) {
      const V p = add(a, mul(edge, q[0]));
      const double length = norm(p);
      const V y = mul(p, 1 / length),
              dy = mul(sub(edge, mul(y, dot(y, edge))), 1 / length);
      const V delta = sub(y, x);
      const double distance = norm(delta);
      const double denominator = distance * (distance + dot(pole, delta));
      if (!(denominator > 0))
        throw std::runtime_error("winding gauge singularity");
      omega += q[1] * dot(pole, cross(delta, dy)) / denominator;
    }
  }
  return -static_cast<double>(omega) / (4 * pi);
}
// Independent 2D radial-map area quadrature, used only to validate the 1D
// rule on selected regular/thin patches and arbitrary observation points.
double area_rule(const Shape &s, V x, const Rule &r) {
  const V center = mul(add(add(s.v[0], s.v[1]), s.v[2]), 1. / 3.);
  long double sum = 0;
  for (int e = 0; e < 3; ++e)
    for (auto u : r)
      for (auto v : r) {
        const V a = s.v[e], b = s.v[(e + 1) % 3];
        const V radial = sub(add(mul(a, 1 - v[0]), mul(b, v[0])), center);
        const V p = add(center, mul(radial, u[0]));
        const double len = norm(p);
        const V y = mul(p, 1 / len), delta = sub(x, y);
        const double distance = norm(delta);
        const double jac = u[0] * std::abs(dot(p, cross(radial, sub(b, a)))) /
                           (len * len * len);
        sum += u[1] * v[1] * jac * dot(y, delta) /
               (4 * pi * distance * distance * distance);
      }
  return static_cast<double>(sum);
}
struct Sum {
  long double centroid = 0, q8 = 0, q16 = 0, q32 = 0;
};
Sum evaluate(const std::vector<Atom> &atoms, const std::vector<Panel> &panels,
             V point, const std::array<Rule, 3> &rules) {
  Sum sum;
  for (const auto &p : panels) {
    const auto &a = atoms[p.atom_index];
    const V x = mul(sub(point, xyz(a)), 1 / a.radius);
    const Shape s = shape(p, a);
    const double radial_orientation = dot({p.nx, p.ny, p.nz}, s.center);
    if (std::abs(std::abs(radial_orientation) - 1) > 1e-10)
      throw std::runtime_error(
          "winding audit requires radial sphere-patch normals");
    const double orientation = radial_orientation >= 0 ? 1 : -1;
    const V delta = sub(point, {p.x, p.y, p.z});
    const double d = norm(delta);
    sum.centroid +=
        p.area * dot({p.nx, p.ny, p.nz}, delta) / (4 * pi * d * d * d);
    sum.q8 += orientation * curved(s, x, rules[0]);
    sum.q16 += orientation * curved(s, x, rules[1]);
    sum.q32 += orientation * curved(s, x, rules[2]);
  }
  return sum;
}
} // namespace
int main(int argc, char **argv) {
  try {
    std::cout << std::setprecision(17) << std::unitbuf;
    const std::array<Rule, 3> rules{gauss(8), gauss(16), gauss(32)};
    TitanExecutor executor(8);
    Options o;
    o.surface_only = true;
    o.surface_mode = "vdw";
    o.highres_surface = true;
    if (argc == 1) {
      int checks = 0;
      auto require = [&](bool v, const char *m) {
        ++checks;
        if (!v)
          throw std::runtime_error(m);
      };
      for (auto mode : {"icosphere", "latlong"}) {
        o.sphere_mesh_mode = mode;
        o.surface_point_count = 96;
        std::vector<Atom> a{{0, 0, 0, 0, 1.7}};
        auto p = build_surface(a, o, nullptr, &executor);
        for (V x : std::vector<V>{
                 {0, 0, 0}, {1.36, 0, 0}, {3.4, 0, 0}, {100, 200, 300}}) {
          auto s = evaluate(a, p, x, rules);
          const double expected = norm(x) < 1.7 ? -1 : 0;
          require(std::abs(s.q32 - expected) < 2e-11,
                  "closed sphere inside/outside identity");
          require(std::abs(s.q32 - s.q16) < 2e-10,
                  "sphere quadrature convergence");
        }
        const auto s = shape(p[0], a[0]);
        for (V x : std::vector<V>{{0, 0, 0}, {.4, .2, -.1}, {2, .3, .2}})
          require(std::abs(curved(s, x, rules[2]) - area_rule(s, x, rules[2])) <
                      1e-12,
                  "boundary/area integration agree");
        auto removed = p;
        removed.erase(removed.begin());
        const auto cut = evaluate(a, removed, {0, 0, 0}, rules);
        require(std::abs(cut.q32 + 1 - p[0].area / (4 * pi * 1.7 * 1.7)) <
                    1e-12,
                "known removed patch solid-angle defect");
        auto flipped = p;
        flipped[0].nx = -flipped[0].nx;
        flipped[0].ny = -flipped[0].ny;
        flipped[0].nz = -flipped[0].nz;
        const auto reverse = evaluate(a, flipped, {0, 0, 0}, rules);
        require(std::abs(reverse.q32 + 1 -
                         2 * p[0].area / (4 * pi * 1.7 * 1.7)) < 1e-12,
                "reversed production normal changes the signed winding");
      }
      const auto thin = shape(
          {unit({-.1, -1e-5, 1}), unit({.1, -1e-5, 1}), unit({0, 2e-5, 1})});
      require(std::abs(curved(thin, {0, 0, 0}, rules[2]) -
                       area_rule(thin, {0, 0, 0}, rules[2])) < 1e-12,
              "thin patch independent integration");
      std::cout << "PASS " << checks << " winding controls\n";
      return 0;
    }
    if (argc != 4)
      return 2;
    std::ifstream in(argv[1]);
    std::size_t n{};
    if (!(in >> n))
      return 2;
    std::vector<Atom> a(n);
    for (auto &x : a)
      if (!(in >> x.x >> x.y >> x.z >> x.radius >> x.q))
        return 2;
    o.surface_point_count = std::stoi(argv[2]);
    o.sphere_mesh_mode = argv[3];
    auto p = build_surface(a, o, nullptr, &executor);
    std::vector<Sum> sums(n);
    std::vector<double> owned_area(n), max_edge(n);
    std::vector<std::size_t> owned_count(n);
    for (const auto &s : p) {
      owned_area[s.atom_index] += s.area;
      ++owned_count[s.atom_index];
      const V v0{s.v0x, s.v0y, s.v0z}, v1{s.v1x, s.v1y, s.v1z},
          v2{s.v2x, s.v2y, s.v2z};
      max_edge[s.atom_index] =
          std::max({max_edge[s.atom_index], norm(sub(v0, v1)),
                    norm(sub(v1, v2)), norm(sub(v2, v0))});
    }
    const auto start = std::chrono::steady_clock::now();
    executor.parallel_for(
        n, 8, 1, [&](std::size_t, std::size_t first, std::size_t last) {
          for (std::size_t i = first; i < last; ++i)
            sums[i] = evaluate(a, p, xyz(a[i]), rules);
        });
    std::cout << "atom_index0,charge,radius,centroid,q8,q16,q32,centroid_"
                 "defect,curved_defect,q16_q32_difference,owned_panels,owned_"
                 "area,max_owned_edge\n";
    for (std::size_t i = 0; i < n; ++i) {
      const auto &s = sums[i];
      std::cout << i << ',' << a[i].q << ',' << a[i].radius << ',' << s.centroid
                << ',' << s.q8 << ',' << s.q16 << ',' << s.q32 << ','
                << s.centroid + 1 << ',' << s.q32 + 1 << ','
                << std::abs(s.q32 - s.q16) << ',' << owned_count[i] << ','
                << owned_area[i] << ',' << max_edge[i] << '\n';
    }
    std::cerr << "panels=" << p.size() << " seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count()
              << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
