// Geometry audit for TITAN's production-default TABI surface, sphere-patch.
//
// The 2026-09-14 SES audit covered grid-ses, analytic-ses/current, patch-ray
// and the PBE grid surface. It did not cover sphere-patch, which is the
// tabi_primary_surface default and therefore the surface under the FAGBI
// manuscript evidence. "Not in the failing set" is not evidence of being
// correct, so this driver applies positive geometric checks to it.
//
// sphere-patch is not a stitched watertight triangulation; it is a union of
// per-atom spherical patches with buried patches discarded. The SES
// self-intersection auditor, which walks indexed vertex adjacency, is therefore
// the wrong instrument. The equivalent and stricter question for a patch union
// is whether any retained panel lies inside another atom's sphere: such a panel
// is interior geometry that the boundary integral must never see.
//
// Checks, all against closed-form references where one exists:
//   1. single sphere   - panel corner radii, total area vs 4*pi*R^2, outward
//                        normals, and convergence with the point count
//   2. two spheres     - total area vs the analytic union area, which is the
//                        pair of buried spherical caps removed
//   3. any input       - burial violations: panel corners and centroids that
//                        sit strictly inside another atom's sphere, with the
//                        deepest penetration reported
//   4. any input       - panel area and shape quality floors
//
// Usage: sphere_patch_surface_audit <atoms file> [point_count] [mesh mode]
//                                   [surface mode] [probe] [threads]
// Atom file: first line is the count, then "x y z radius charge" per atom,
// matching tmp/ses_audit_20260914/*.atoms.

#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Panel;

struct Point {
  double x, y, z;
};

double sphere_area(double r) { return 4.0 * std::numbers::pi * r * r; }

// Exposed area of the union of two spheres: each sphere loses the cap that the
// other one buries. Returns a negative value when the configuration is not a
// proper partial overlap, so the caller can skip the comparison.
double two_sphere_union_area(double r1, double r2, double d) {
  if (!(d > 0.0) || d >= r1 + r2 || d <= std::abs(r1 - r2))
    return -1.0;
  const double x1 = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
  const double x2 = (d * d + r2 * r2 - r1 * r1) / (2.0 * d);
  const double cap1 = 2.0 * std::numbers::pi * r1 * (r1 - x1);
  const double cap2 = 2.0 * std::numbers::pi * r2 * (r2 - x2);
  return sphere_area(r1) - cap1 + sphere_area(r2) - cap2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2 || argc > 7) {
      std::cerr << "usage: sphere_patch_surface_audit <atoms> [points] "
                   "[mesh mode] [surface mode] [probe] [threads]\n";
      return 2;
    }
    std::ifstream in(argv[1]);
    std::size_t n = 0;
    if (!(in >> n) || n == 0)
      return 3;
    std::vector<Atom> atoms(n);
    for (Atom &a : atoms)
      if (!(in >> a.x >> a.y >> a.z >> a.radius >> a.q))
        return 3;

    const int point_count = argc > 2 ? std::stoi(argv[2]) : 96;
    const std::string mesh_mode = argc > 3 ? argv[3] : "auto";
    const std::string surface_mode = argc > 4 ? argv[4] : "vdw";
    const double probe = argc > 5 ? std::stod(argv[5]) : 1.4;
    const std::size_t threads = argc > 6 ? std::stoull(argv[6]) : 1;

    TitanPBE::Tabi::Options options;
    options.primary_surface = "sphere-patch";
    options.surface_mode = surface_mode;
    options.probe_radius = probe;
    options.surface_point_count = point_count;
    options.sphere_mesh_mode = mesh_mode;
    options.surface_only = true; // geometry only, no solid-angle diagonals
    // build_triangle_mesh clamps the point count to 4096 unless the high
    // resolution path is enabled, so a refinement study silently plateaus
    // there. Enable it automatically once the request exceeds the clamp.
    options.highres_surface = point_count > 4096;

    TitanExecutor executor(threads);
    const auto begin = std::chrono::steady_clock::now();
    const std::vector<Panel> panels =
        TitanPBE::Tabi::build_surface(atoms, options, nullptr, &executor);
    const double wall_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin)
            .count();

    // The effective sphere radius is what the surface is actually built on:
    // "sas" inflates every atom by the probe, "vdw" does not.
    const double inflate = surface_mode == "sas" ? std::max(0.0, probe) : 0.0;

    double area = 0.0, min_area = std::numeric_limits<double>::infinity();
    double min_quality = std::numeric_limits<double>::infinity();
    double volume = 0.0;
    double centre_x = 0.0, centre_y = 0.0, centre_z = 0.0;
    for (const Atom &a : atoms) {
      centre_x += a.x / static_cast<double>(n);
      centre_y += a.y / static_cast<double>(n);
      centre_z += a.z / static_cast<double>(n);
    }

    double rad_min = std::numeric_limits<double>::infinity();
    double rad_max = 0.0;
    double normal_cos_min = 1.0;
    std::size_t burial_points = 0;
    std::size_t burial_panels = 0;
    double burial_depth_max = 0.0;
    // A point is only flagged when it is inside another sphere by more than a
    // tolerance, so shared boundary circles between touching spheres, where
    // points legitimately sit on two surfaces at once, are not counted.
    constexpr double kBurialTolerance = 1.0e-9;

    for (const Panel &p : panels) {
      area += p.area;
      min_area = std::min(min_area, p.area);
      const double ab2 = (p.v1x - p.v0x) * (p.v1x - p.v0x) +
                         (p.v1y - p.v0y) * (p.v1y - p.v0y) +
                         (p.v1z - p.v0z) * (p.v1z - p.v0z);
      const double ac2 = (p.v2x - p.v0x) * (p.v2x - p.v0x) +
                         (p.v2y - p.v0y) * (p.v2y - p.v0y) +
                         (p.v2z - p.v0z) * (p.v2z - p.v0z);
      const double bc2 = (p.v2x - p.v1x) * (p.v2x - p.v1x) +
                         (p.v2y - p.v1y) * (p.v2y - p.v1y) +
                         (p.v2z - p.v1z) * (p.v2z - p.v1z);
      if (ab2 + ac2 + bc2 > 0.0)
        min_quality = std::min(
            min_quality, 4.0 * std::sqrt(3.0) * p.area / (ab2 + ac2 + bc2));
      volume += p.area * ((p.x - centre_x) * p.nx + (p.y - centre_y) * p.ny +
                          (p.z - centre_z) * p.nz) / 3.0;

      const Point probe_points[4] = {{p.v0x, p.v0y, p.v0z},
                                     {p.v1x, p.v1y, p.v1z},
                                     {p.v2x, p.v2y, p.v2z},
                                     {p.x, p.y, p.z}};
      bool panel_flagged = false;
      for (const Point &q : probe_points) {
        if (n == 1) {
          const double d = std::hypot(q.x - atoms[0].x, q.y - atoms[0].y,
                                      q.z - atoms[0].z);
          rad_min = std::min(rad_min, d);
          rad_max = std::max(rad_max, d);
        }
        for (std::size_t k = 0; k < n; ++k) {
          if (k == p.atom_index)
            continue;
          const double rk = atoms[k].radius + inflate;
          if (!(rk > 0.0))
            continue;
          const double d = std::hypot(q.x - atoms[k].x, q.y - atoms[k].y,
                                      q.z - atoms[k].z);
          const double depth = rk - d;
          if (depth > kBurialTolerance) {
            ++burial_points;
            burial_depth_max = std::max(burial_depth_max, depth);
            panel_flagged = true;
          }
        }
      }
      if (panel_flagged)
        ++burial_panels;

      if (n == 1) {
        const double dx = p.x - atoms[0].x, dy = p.y - atoms[0].y,
                     dz = p.z - atoms[0].z;
        const double len = std::hypot(dx, dy, dz);
        if (len > 0.0)
          normal_cos_min = std::min(
              normal_cos_min, (dx * p.nx + dy * p.ny + dz * p.nz) / len);
      }
    }

    double area_reference = -1.0;
    if (n == 1) {
      area_reference = sphere_area(atoms[0].radius + inflate);
    } else if (n == 2) {
      const double d = std::hypot(atoms[0].x - atoms[1].x,
                                  atoms[0].y - atoms[1].y,
                                  atoms[0].z - atoms[1].z);
      area_reference = two_sphere_union_area(atoms[0].radius + inflate,
                                             atoms[1].radius + inflate, d);
    }
    const double area_rel_error =
        area_reference > 0.0 ? (area - area_reference) / area_reference : 0.0;

    std::cout << std::setprecision(17);
    std::cout << "mode,points,mesh_mode,surface_mode,atoms,panels,wall_ms,"
                 "area,area_reference,area_rel_error,volume,min_area,"
                 "min_quality,sphere_rad_min,sphere_rad_max,"
                 "sphere_normal_cos_min,burial_points,burial_panels,"
                 "burial_depth_max\n";
    std::cout << "sphere-patch," << point_count << ',' << mesh_mode << ','
              << surface_mode << ',' << n << ',' << panels.size() << ','
              << wall_ms << ',' << area << ',' << area_reference << ','
              << area_rel_error << ',' << volume << ',' << min_area << ','
              << min_quality << ','
              << (n == 1 ? rad_min : 0.0) << ','
              << (n == 1 ? rad_max : 0.0) << ','
              << (n == 1 ? normal_cos_min : 0.0) << ',' << burial_points << ','
              << burial_panels << ',' << burial_depth_max << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
