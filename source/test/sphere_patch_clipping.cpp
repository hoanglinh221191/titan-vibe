#include "tabi/tabi_charge_resolution.h"
#include "tabi/tabi_near.h"
#include "tabi/tabi_sphere_clip.h"
#include "titan_executor.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
using namespace TitanPBE::Tabi;
int main(int argc, char **argv) {
  try {
    std::cout << std::unitbuf << std::setprecision(17);
    int checks = 0;
    auto require = [&](bool pass, const char *m) {
      ++checks;
      if (!pass)
        throw std::runtime_error(m);
    };
    TitanExecutor executor(4);
    Options o;
    o.surface_mode = "vdw";
    o.surface_only = true;
    o.surface_point_count = 96;
    o.highres_surface = true;
    auto audit = [&](const std::vector<Atom> &a, int points) {
      o.surface_point_count = points;
      auto p = build_surface(a, o, nullptr, &executor);
      double area = 0, depth = 0;
      for (const auto &s : p) {
        require(s.area > 0 && std::isfinite(s.area),
                "finite positive panel area");
        area += s.area;
        for (const auto &b : a)
          depth = std::max(depth, b.radius - charge_to_spherical_panel_distance(
                                                 b, a[s.atom_index], s));
      }
      std::cout << "GEOMETRY " << a.size() << ' ' << points << ' ' << p.size()
                << ' ' << area << ' ' << depth << '\n';
      require(depth < 2e-10,
              "no curved patch enters any atom, including edges/interiors");
      return area;
    };
    if (argc > 1) {
      if (argc > 2 && std::string(argv[2]) == "prepare") {
        o.surface_only = false;
        o.salt_m = 0;
        o.near_integration = "auto";
      }
      std::ifstream in(argv[1]);
      std::size_t n{};
      if (!(in >> n))
        return 2;
      std::vector<Atom> a(n);
      for (auto &b : a)
        if (!(in >> b.x >> b.y >> b.z >> b.radius >> b.q))
          return 2;
      for (int points : {96, 384, 1536})
        audit(a, points);
      return 0;
    }
    const double pi = std::numbers::pi, R = 1.7;
    const std::vector<Atom> two = {{0, 0, 0, 1, R}, {2.5, 0, 0, -1, R}};
    const double exact = 4 * pi * R * R + 2 * pi * R * 2.5;
    double previous = 1;
    for (int points : {96, 384, 1536, 6144}) {
      const double area = audit(two, points), error = (exact - area) / exact;
      require(error >= -1e-12 && error < previous,
              "two-sphere area approaches exact union from below");
      previous = error;
    }
    require(previous < 2e-5, "fine two-sphere union area");
    for (auto atoms : std::vector<std::vector<Atom>>{
             {{0, 0, 0, 1, R}, {0, 0, 0, -1, R}},
             {{0, 0, 0, 1, R}, {.5, 0, 0, -1, .3}},
             {{0, 0, 0, 1, R}, {3.4, 0, 0, -1, R}},
             {{0, 0, 0, 1, R}, {1.4, 0, 0, -1, .3}},
             {{0, 0, 0, 1, R}, {1.6, 0, 0, -1, .3}},
             {{0, 0, 0, 1, R}, {1e-7, 0, 0, -1, R}},
             {{0, 0, 0, 1, R}, {2.5, 0, 0, -1, R}, {1.25, 2, 0, .5, R}}})
      audit(atoms, 96);
    auto shifted = two;
    for (auto &a : shifted) {
      a.x += 100;
      a.y -= 50;
      a.z += 25;
    }
    require(std::abs(audit(two, 384) - audit(shifted, 384)) < 1e-10,
            "translated clipping area invariant");
    auto reversed = two;
    std::reverse(reversed.begin(), reversed.end());
    require(std::abs(audit(two, 384) - audit(reversed, 384)) < 1e-10,
            "two-sphere atom permutation area");
    // Hidden cap / exposed island: all three base vertices have the same
    // classification although the interior has the opposite classification.
    {
      using namespace SphereClip;
      const SphereClip::Triangle t{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
      const V n = unit({1, 1, 1});
      for (bool island : {false, true}) {
        const V axis = island ? mul(n, -1) : n;
        const double c = (island ? -1 : 1) * std::cos(.01);
        std::vector<SphereClip::Triangle> pieces;
        SphereClip::clip(t, {{axis, c}}, pieces);
        double kept = 0;
        for (const auto &piece : pieces) {
          Panel panel;
          panel.has_triangle_geometry = true;
          panel.v0x = piece[0][0];
          panel.v0y = piece[0][1];
          panel.v0z = piece[0][2];
          panel.v1x = piece[1][0];
          panel.v1y = piece[1][1];
          panel.v1z = piece[1][2];
          panel.v2x = piece[2][0];
          panel.v2y = piece[2][1];
          panel.v2z = piece[2][2];
          Atom owner{0, 0, 0, 0, 1},
              cover{axis[0], axis[1], axis[2], 0, std::sqrt(2 - 2 * c)};
          require(charge_to_spherical_panel_distance(cover, owner, panel) >=
                      cover.radius - 1e-12,
                  "hidden-component clipping excludes the covering ball");
          kept += area(piece);
        }
        const double cap = 2 * pi * (1 - std::cos(.01));
        const double reference = island ? cap : pi / 2 - cap;
        require(kept > 0 && kept <= reference + 1e-12 &&
                    reference - kept < .1 * cap,
                "hidden cap and island retain their exposed component");
      }
    }
    // Independent planar weak-singular formula in a shrinking tangent patch.
    // Curvature corrections vanish as h^2; includes severe clipping slivers.
    {
      using namespace SphereClip;
      bool rejected = false;
      try {
        const double h = 1e-12;
        prepare_sphere_patch_quadrature(unit({-h, 0, 1}), unit({h, 0, 1}),
                                        unit({0, h, 1}), nullptr, true);
      } catch (const std::runtime_error &e) {
        rejected =
            std::string(e.what()) == "sphere_patch_quadrature_nodes_invalid";
      }
      require(rejected, "a truncated node rule must not be accepted by "
                        "boundary self preparation");
    }
    {
      using namespace SphereClip;
      const SphereClip::Triangle t{unit({-.2, -.1, 1}), unit({.2, -.1, 1}),
                                   unit({0, .2, 1})};
      const auto duffy =
          prepare_sphere_patch_quadrature(t[0], t[1], t[2], nullptr, false);
      const auto boundary =
          prepare_sphere_patch_quadrature(t[0], t[1], t[2], nullptr, true);
      require(std::abs(boundary->self_double_layer / duffy->self_double_layer -
                       1) < 2e-8,
              "boundary self agrees with independent radial Duffy rule");
    }
    for (double aspect : {5., 50., 500., 5000.}) {
      using namespace SphereClip;
      const double h = 1e-3, w = h / aspect;
      const std::array<V, 3> flat{{{-h, -w, 0}, {h, -w, 0}, {0, 2 * w, 0}}};
      double planar = 0;
      for (int e = 0; e < 3; ++e) {
        const V a = flat[e], b = flat[(e + 1) % 3], d = add(b, mul(a, -1));
        const double length = norm(d), ra = norm(a), rb = norm(b),
                     height = norm(cross(a, b)) / length;
        planar += height * std::log((ra + rb + length) / (ra + rb - length));
      }
      SphereClip::Triangle t;
      for (int i = 0; i < 3; ++i)
        t[i] = unit({flat[i][0], flat[i][1], 1});
      auto q = prepare_sphere_patch_quadrature(t[0], t[1], t[2], nullptr, true);
      const double ref = -planar / (8 * pi),
                   error = std::abs(q->self_double_layer / ref - 1);
      std::cout << "SELF " << aspect << ' ' << q->self_double_layer << ' '
                << ref << ' ' << error << '\n';
      require(error < 2e-6,
              "curved sliver self approaches independent planar formula");
    }
    std::cout << "PASS " << checks << " checks\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
