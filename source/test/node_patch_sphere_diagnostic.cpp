// Frozen baseline diagnosis: the pre-cache node conversion versus exact duals.
// The historical "production" columns intentionally retain that baseline.
// Use node_patch_cached_rows.cpp to audit the current prepared production path.
// No fitted coefficients; no production geometry or quadrature changes.
#include "tabi/tabi_near.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

using namespace TitanPBE::Tabi;
constexpr double pi = std::numbers::pi;
constexpr double coulomb = 332.063713299;
struct V {
  double x, y, z;
};
V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
V operator*(V a, double b) { return {a.x * b, a.y * b, a.z * b}; }
double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
V cross(V a, V b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(V a) { return std::sqrt(dot(a, a)); }

// Kirkwood boundary trace in the solver's 1/(4*pi) normalization.
double exact_trace(double cosine, double offset, double radius, double q,
                   double ei, double eo) {
  long double sum = 0, pprev = 1, p = 1, power = 1;
  for (int l = 0; l < 400; ++l) {
    if (l == 1)
      p = cosine;
    else if (l > 1) {
      long double next = ((2 * l - 1) * cosine * p - (l - 1) * pprev) / l;
      pprev = p;
      p = next;
    }
    sum += (2 * l + 1) * power * p / (ei * l + eo * (l + 1));
    power *= offset;
  }
  return q * static_cast<double>(sum) / (4 * pi * radius);
}
double exact_energy(double offset, double radius, double q, double ei,
                    double eo) {
  long double sum = 0, power = 1;
  for (int l = 0; l < 400; ++l) {
    sum += (ei - eo) * (l + 1) * power / (ei * l + eo * (l + 1));
    power *= offset * offset;
  }
  return 0.5 * coulomb * q * q * static_cast<double>(sum) / (ei * radius);
}

// Screened (linearised PB) Kirkwood energy, for salt > 0.
// Matching phi and eps*dphi/dr on the sphere gives the exterior response
// gamma_l = kappa*a*k_l'(kappa*a)/k_l(kappa*a), with k_l the modified spherical
// Bessel function of the second kind.  Since k_l'(z) = -k_{l+1}(z)+(l/z)k_l(z),
// gamma_l = l - z*r_l with r_l = k_{l+1}(z)/k_l(z).  That ratio obeys
// r_{l+1} = (2l+3)/z + 1/r_l, r_0 = 1 + 1/z, which stays bounded where the
// functions themselves overflow.  As z -> 0, gamma_l -> -(l+1) and the series
// reduces term by term to the unscreened exact_energy above.
double exact_energy_screened(double offset, double radius, double q, double ei,
                             double eo, double kappa) {
  if (kappa <= 0)
    return exact_energy(offset, radius, q, ei, eo);
  const long double z = kappa * radius;
  long double sum = 0, power = 1, r = 1 + 1 / z;
  for (int l = 0; l < 400; ++l) {
    const long double gamma = l - z * r;
    sum += (eo * gamma + ei * (l + 1)) * power / (ei * l - eo * gamma);
    power *= offset * offset;
    r = (2.0L * l + 3.0L) / z + 1 / r;
  }
  return 0.5 * coulomb * q * q * static_cast<double>(sum) / (ei * radius);
}

struct Rule {
  std::vector<double> x, w;
};
Rule gauss(int n) {
  Rule r;
  for (int i = 0; i < n; ++i) {
    double z = std::cos(pi * (i + .75) / (n + .5)), d = 0;
    for (int it = 0; it < 30; ++it) {
      double p = 1, prev = 0;
      for (int k = 1; k <= n; ++k) {
        double old = p;
        p = ((2 * k - 1) * z * p - (k - 1) * prev) / k;
        prev = old;
      }
      d = n * (z * p - prev) / (z * z - 1);
      double next = z - p / d;
      if (std::abs(next - z) < 1e-15) {
        z = next;
        break;
      }
      z = next;
    }
    r.x.push_back((z + 1) / 2);
    r.w.push_back(1 / ((1 - z * z) * d * d));
  }
  return r;
}

#include "tabi/tabi_node_patch.h"
#include <Eigen/Dense>
#include <array>
#include <fstream>
#include <set>

struct Trace {
  std::vector<long double> f, df;
  Trace(double offset, double radius, double ei, double eo, double kappa) {
    long double power = 1, z = kappa * radius, ratio = z > 0 ? 1 + 1 / z : 0;
    for (int l = 0; l < 200; ++l) {
      long double gamma = z > 0 ? l - z * ratio : -(l + 1);
      long double a =
          (2 * l + 1) * power / (4 * pi * radius * (ei * l - eo * gamma));
      f.push_back(a);
      df.push_back(a * eo / ei * gamma / radius);
      if (l > 2 && std::abs(power) < 1e-20)
        break;
      power *= offset;
      if (z > 0)
        ratio = (2 * l + 3) / z + 1 / ratio;
    }
  }
  std::array<double, 2> operator()(double cosine) const {
    long double p = 1, prev = 1, s = 0, d = 0;
    for (std::size_t l = 0; l < f.size(); ++l) {
      if (l == 1)
        p = cosine;
      else if (l > 1) {
        auto next = ((2 * l - 1) * cosine * p - (l - 1) * prev) / l;
        prev = p;
        p = next;
      }
      s += f[l] * p;
      d += df[l] * p;
    }
    return {static_cast<double>(s), static_cast<double>(d)};
  }
};
struct QNode {
  V y;
  double w;
  std::array<double, 2> trace;
};
using QCell = std::vector<QNode>;
V xyz(const std::array<double, 3> &a) { return {a[0], a[1], a[2]}; }
V radial(V p, double R) { return p * (R / norm(p)); }
QCell dual_rule(const NodePatchGeometry &g, std::size_t node, double R, V axis,
                const Trace &trace, int order) {
  QCell out;
  const Rule rule = gauss(order);
  for (auto pos = g.incident_offsets[node]; pos < g.incident_offsets[node + 1];
       ++pos) {
    const auto &t = g.triangles[g.incident_triangles[pos]];
    int k = 0;
    while (t.vertices[k] != node)
      ++k;
    const V a = xyz(g.vertex_xyz[node]),
            b = xyz(g.vertex_xyz[t.vertices[(k + 1) % 3]]),
            c = xyz(g.vertex_xyz[t.vertices[(k + 2) % 3]]),
            mid = (a + b + c) * (1. / 3);
    for (const auto &edge : std::array<std::array<V, 2>, 2>{
             {{{(a + b) * .5, mid}}, {{mid, (a + c) * .5}}}}) {
      for (std::size_t i = 0; i < rule.x.size(); ++i)
        for (std::size_t j = 0; j < rule.x.size(); ++j) {
          const double u = rule.x[i], v = rule.x[j];
          const V du = edge[0] * (1 - v) + edge[1] * v - a,
                  dv = (edge[1] - edge[0]) * u, p = a + du * u;
          const double len = norm(p);
          const V y = p * (R / len);
          const double w = rule.w[i] * rule.w[j] * R * R *
                           std::abs(dot(p, cross(du, dv))) / (len * len * len);
          out.push_back({y, w, trace(dot(y, axis) / R)});
        }
    }
  }
  return out;
}
using Blocks = std::array<double, 4>;
Blocks kernel_sphere(V x, V y, double R, double ratio, double kappa) {
  const long double r = norm(x - y), s = kappa * r, e = std::exp(-s);
  long double E1;
  if (s < .5L) {
    E1 = 0;
    long double power = s * s, factorial = 1;
    for (int m = 0; m < 24; ++m) {
      if (m)
        factorial *= m;
      E1 += (m % 2 ? -1 : 1) * power / (factorial * (m + 2));
      power *= s;
    }
  } else
    E1 = 1 - (1 + s) * e;
  const long double c = 1 / (4 * pi), r2 = r * r, rr = r2 / (4 * R * R);
  return {static_cast<double>(-std::expm1(-s) * c / r),
          static_cast<double>(-(ratio * (1 - E1) - 1) * c / (2 * R * r)),
          static_cast<double>(-(1 - (1 - E1) / ratio) * c / (2 * R * r)),
          static_cast<double>((s * s * e * rr - E1 * (1 + rr)) * c / (r * r2))};
}
Blocks weight_trace(const Blocks &k, const std::array<double, 2> &u) {
  return {k[0] * u[1], k[1] * u[0], k[2] * u[1], k[3] * u[0]};
}
void add(Blocks &a, const Blocks &b, double w = 1) {
  for (int k = 0; k < 4; ++k)
    a[k] += w * b[k];
}
double project(V charge, V y, V normal, double ratio, double kappa,
               std::array<double, 2> u) {
  const V delta = charge - y;
  const double r = norm(delta), s = kappa * r, e = std::exp(-s);
  const double k1 = -std::expm1(-s) / (4 * pi * r),
               k2 = dot(normal, delta) / (4 * pi * r * r * r) *
                    (ratio * e * (1 + s) - 1);
  return k1 * u[1] + k2 * u[0];
}
void solve_variants(const std::vector<Panel> &panels,
                    const std::vector<QCell> &far, const Trace &trace, V axis,
                    V charge, double R, double ratio, double kappa,
                    int near_order, TitanExecutor &ex) {
  const auto &g = *panels.front().node_patch_geometry;
  const std::size_t n = panels.size();
  const bool coupled = kappa > 1e-14;
  const std::size_t dim = coupled ? 2 * n : n;
  Eigen::MatrixXd production = Eigen::MatrixXd::Zero(dim, dim),
                  curved = production;
  Eigen::VectorXd rhs(dim), radial_rhs(dim);
  std::vector<int> groups(n * n);
  ex.parallel_for(
      n, 8, 1, [&](std::size_t, std::size_t begin, std::size_t end) {
        for (auto i = begin; i < end; ++i) {
          const auto &p = panels[i];
          const V x = xyz(g.vertex_xyz[i]), delta = x - charge;
          const double d = norm(delta);
          rhs[i] = radial_rhs[i] = 1 / (4 * pi * d);
          if (coupled) {
            rhs[n + i] = -dot({p.nx, p.ny, p.nz}, delta) / (4 * pi * d * d * d);
            radial_rhs[n + i] = -dot(x * (1 / R), delta) / (4 * pi * d * d * d);
          }
          production(i, i) = curved(i, i) = .5 * (1 + ratio);
          if (coupled)
            production(n + i, n + i) = curved(n + i, n + i) =
                .5 * (1 + 1 / ratio);
          for (std::size_t j = 0; j < n; ++j) {
            const auto k = evaluate_jcp_near_pair_kernel(
                p, i, panels[j], j, ratio, kappa, coupled, "analytic_singular",
                "off");
            if (!k.valid && !is_intentional_zero_self_kernel(
                                k, coupled, "analytic_singular"))
              throw std::runtime_error("invalid matrix kernel");
            const int group = i == j ? 0 : k.one_ring ? 1 : 2;
            groups[i * n + j] = group;
            const auto fine = group < 2
                                  ? dual_rule(g, j, R, axis, trace, near_order)
                                  : QCell{};
            Blocks b{};
            for (const auto &q : group < 2 ? fine : far[j])
              add(b, kernel_sphere(x, q.y, R, ratio, kappa), q.w);
            production(i, j) -= k.k2;
            curved(i, j) -= b[1];
            if (coupled) {
              production(i, n + j) -= k.k1;
              production(n + i, j) -= k.k4;
              production(n + i, n + j) -= k.k3;
              curved(i, n + j) -= b[0];
              curved(n + i, j) -= b[3];
              curved(n + i, n + j) -= b[2];
            }
          }
        }
      });
  for (int mode = 0; mode < 4; ++mode) {
    Eigen::MatrixXd a = production;
    Eigen::VectorXd b = rhs;
    if (mode == 3) {
      a = curved;
      b = radial_rhs;
    } else if (mode)
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
          if (groups[i * n + j] < mode) {
            a(i, j) = curved(i, j);
            if (coupled) {
              a(i, n + j) = curved(i, n + j);
              a(n + i, j) = curved(n + i, j);
              a(n + i, n + j) = curved(n + i, n + j);
            }
          }
    Eigen::VectorXd u = a.partialPivLu().solve(b);
    const double residual = (a * u - b).norm() / b.norm();
    double ep = 0, ec = 0;
    for (std::size_t j = 0; j < n; ++j) {
      const auto &p = panels[j];
      std::array<double, 2> t{u[j], coupled ? u[n + j] : 0};
      ep += p.area * project(charge, {p.source_x, p.source_y, p.source_z},
                             {p.nx, p.ny, p.nz}, ratio, kappa, t);
      for (const auto &q : far[j])
        ec += q.w * project(charge, q.y, q.y * (1 / R), ratio, kappa, t);
    }
    const double factor = .5 * coulomb * 4 * pi;
    std::cout << "SOLVE {\"mode\":" << mode << ",\"residual\":" << residual
              << ",\"energy_production_projection\":" << ep * factor
              << ",\"energy_curved_projection\":" << ec * factor << "}\n";
    if (!u.allFinite() || residual > 1e-10)
      throw std::runtime_error("dense reference solve failed");
  }
}
int main(int argc, char **argv) {
  try {
    if (argc == 1) {
      int checks = 0;
      auto require = [&](bool pass) {
        ++checks;
        if (!pass)
          throw std::runtime_error("diagnostic control failed");
      };
      for (double R : {.7, 1.7, 3.4})
        for (double ratio : {1., 3., 78.5})
          for (double kappa : {0., .1, 1.}) {
            V x = radial({.2, .3, 1}, R), y = radial({-.4, .2, .8}, R);
            Panel target, source;
            target.x = x.x;
            target.y = x.y;
            target.z = x.z;
            target.nx = x.x / R;
            target.ny = x.y / R;
            target.nz = x.z / R;
            source.x = y.x;
            source.y = y.y;
            source.z = y.z;
            source.nx = y.x / R;
            source.ny = y.y / R;
            source.nz = y.z / R;
            source.area = 1;
            const auto p = evaluate_jcp_centroid_monopole_kernel(target, source,
                                                                 ratio, kappa);
            const Blocks generic{p.k1, p.k2, p.k3, p.k4},
                spherical = kernel_sphere(x, y, R, ratio, kappa);
            for (int b = 0; b < 4; ++b)
              require(std::abs(generic[b] - spherical[b]) <
                      1e-11 * std::max(1., std::abs(generic[b])));
            Trace center(0, R, 1, ratio, kappa);
            const auto u = center(.3);
            require(std::abs(u[0] - 1 / (4 * pi * ratio * R *
                                         (1 + kappa * R))) < 1e-13);
            require(std::abs(u[1] + 1 / (4 * pi * R * R)) < 1e-13);
          }
      std::cout << "PASS " << checks << " kernel and boundary trace controls\n";
      return 0;
    }
    if (argc != 7 && argc != 8)
      throw std::runtime_error(
          "points offset salt direction far_order near_order [solve]");
    const int points = std::stoi(argv[1]), direction = std::stoi(argv[4]),
              fo = std::stoi(argv[5]), no = std::stoi(argv[6]);
    const double off = std::stod(argv[2]), salt = std::stod(argv[3]), R = 1.7,
                 ei = 1, eo = 78.5, ratio = eo / ei;
    const double kappa = .3291 * std::sqrt(salt);
    const bool coupled = kappa > 1e-14;
    if (off < 0 || off > .8 || salt < 0 || fo < 2 || no < fo)
      throw std::runtime_error("unsupported domain");
    const V axis = direction ? V{1 / std::sqrt(3.), 1 / std::sqrt(3.),
                                 1 / std::sqrt(3.)}
                             : V{1, 0, 0},
            charge = axis * (off * R);
    Options o;
    o.primary_surface = "sphere-patch";
    o.surface_mode = "vdw";
    o.sphere_mesh_mode = "icosphere";
    o.surface_point_count = points;
    o.salt_m = salt;
    o.surface_only = true;
    TitanExecutor ex(8);
    auto triangles = build_surface({{0, 0, 0, 0, R}}, o, nullptr, &ex);
    auto panels = build_node_patch_panels(triangles, "solid-angle");
    const std::size_t n = panels.size();
    const auto &g = *panels.front().node_patch_geometry;
    Trace trace(off, R, ei, eo, kappa);
    std::vector<std::array<double, 2>> u(n);
    std::vector<QCell> far(n);
    ex.parallel_for(n, 8, 1,
                    [&](std::size_t, std::size_t begin, std::size_t end) {
                      for (auto j = begin; j < end; ++j) {
                        u[j] = trace(dot(xyz(g.vertex_xyz[j]), axis) / R);
                        far[j] = dual_rule(g, j, R, axis, trace, fo);
                      }
                    });
    double area = 0, area_difference = 0, normal_error = 0, projprod = 0,
           projconstant = 0, projexact = 0;
    for (std::size_t j = 0; j < n; ++j) {
      const auto &p = panels[j];
      double a = 0;
      for (const auto &q : far[j]) {
        a += q.w;
        projconstant +=
            q.w * project(charge, q.y, q.y * (1 / R), ratio, kappa, u[j]);
        projexact +=
            q.w * project(charge, q.y, q.y * (1 / R), ratio, kappa, q.trace);
      }
      area += a;
      area_difference = std::max(area_difference, std::abs(a / p.area - 1));
      normal_error =
          std::max(normal_error,
                   norm(V{p.nx, p.ny, p.nz} - xyz(g.vertex_xyz[j]) * (1 / R)));
      projprod += p.area * project(charge, {p.source_x, p.source_y, p.source_z},
                                   {p.nx, p.ny, p.nz}, ratio, kappa, u[j]);
    }
    const double factor = .5 * coulomb * 4 * pi;
    std::cout << std::setprecision(17) << std::unitbuf;
    std::cout << "META {\"triangles\":" << triangles.size()
              << ",\"unknown_nodes\":" << n
              << ",\"scalar_unknowns\":" << (coupled ? 2 * n : n)
              << ",\"area_relative_error\":" << area / (4 * pi * R * R) - 1
              << ",\"max_cell_area_relative_difference\":" << area_difference
              << ",\"max_normal_vector_error\":" << normal_error
              << ",\"projection_production\":" << projprod * factor
              << ",\"projection_curved_constant\":" << projconstant * factor
              << ",\"projection_curved_exact\":" << projexact * factor
              << ",\"energy_exact\":"
              << exact_energy_screened(off, R, 1, ei, eo, kappa) << "}\n";
    std::set<std::size_t> targets;
    for (std::size_t i = 0; i < std::min<std::size_t>(n, 16); ++i)
      targets.insert(i * n / std::min<std::size_t>(n, 16));
    std::size_t closest = 0;
    double best = 1e300;
    for (std::size_t i = 0; i < n; ++i) {
      double d = norm(xyz(g.vertex_xyz[i]) - charge);
      if (d < best) {
        best = d;
        closest = i;
      }
    }
    targets.insert(closest);
    std::cout << "ROW_HEADER "
                 "target,group,block,production,curved_constant,curved_exact\n";
    for (const auto i : targets) {
      std::array<Blocks, 3> prod{}, constant{}, exact{};
      double winding = 0;
      const auto &p = panels[i];
      const V x = xyz(g.vertex_xyz[i]);
      for (std::size_t j = 0; j < n; ++j) {
        const auto k =
            evaluate_jcp_near_pair_kernel(p, i, panels[j], j, ratio, kappa,
                                          coupled, "analytic_singular", "off");
        if (!k.valid &&
            !is_intentional_zero_self_kernel(k, coupled, "analytic_singular"))
          throw std::runtime_error("invalid production kernel");
        const int group = i == j ? 0 : k.one_ring ? 1 : 2;
        add(prod[group], weight_trace({k.k1, k.k2, k.k3, k.k4}, u[j]));
        const QCell fine =
            group < 2 ? dual_rule(g, j, R, axis, trace, no) : QCell{};
        for (const auto &q : group < 2 ? fine : far[j]) {
          auto v = kernel_sphere(x, q.y, R, ratio, kappa);
          if (!coupled) {
            v[0] = 0;
            v[2] = 0;
            v[3] = 0;
          }
          add(constant[group], weight_trace(v, u[j]), q.w);
          add(exact[group], weight_trace(v, q.trace), q.w);
          winding -= q.w / (8 * pi * R * norm(x - q.y));
        }
      }
      for (int group = 0; group < 3; ++group)
        for (int block = 0; block < 4; ++block)
          std::cout << "ROW " << i << ',' << group << ',' << block << ','
                    << prod[group][block] << ',' << constant[group][block]
                    << ',' << exact[group][block] << '\n';
      const V delta = x - charge;
      const double d = norm(delta), rhs1 = 1 / (4 * pi * ei * d),
                   rhs2 = -dot({p.nx, p.ny, p.nz}, delta) /
                          (4 * pi * ei * d * d * d),
                   rhs2radial =
                       -dot(x * (1 / R), delta) / (4 * pi * ei * d * d * d);
      std::cout << "RES " << i << ',' << rhs1 << ',' << rhs2;
      int representation = 0;
      for (const auto &v : {prod, constant, exact}) {
        double s1 = 0, s2 = 0;
        for (const auto &b : v) {
          s1 += b[0] + b[1];
          s2 += b[2] + b[3];
        }
        std::cout << ',' << (.5 * (1 + ratio) * u[i][0] - s1 - rhs1) << ','
                  << (coupled ? .5 * (1 + 1 / ratio) * u[i][1] - s2 -
                                    (representation ? rhs2radial : rhs2)
                              : 0);
        ++representation;
      }
      std::cout << ',' << winding + .5 << ',' << rhs2radial << '\n';
    }
    if (argc == 8) {
      if (std::string(argv[7]) != "solve")
        throw std::runtime_error("unknown extra option");
      solve_variants(panels, far, trace, axis, charge, R, ratio, kappa, no, ex);
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
