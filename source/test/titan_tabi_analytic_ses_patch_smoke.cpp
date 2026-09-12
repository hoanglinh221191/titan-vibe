#include "tabi/tabi_analytic_ses.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::analytic_ses::build_patch_complex;
using TitanPBE::Tabi::analytic_ses::intersect_ray_x;
using TitanPBE::Tabi::analytic_ses::PatchComplex;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

double exposed_angle(const PatchComplex &complex) {
  double result = 0.0;
  for (const auto &torus : complex.tori) {
    for (const auto &interval : torus.exposed_intervals)
      result += interval.end - interval.begin;
  }
  return result;
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

std::vector<Atom> rotate_z(std::vector<Atom> atoms, double angle) {
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  for (Atom &atom : atoms) {
    const double x = atom.x;
    const double y = atom.y;
    atom.x = cosine * x - sine * y;
    atom.y = sine * x + cosine * y;
  }
  return atoms;
}

} // namespace

int main() {
  try {
    const PatchComplex sphere =
        build_patch_complex(std::vector<Atom>{{0.0, 0.0, 0.0, 0.0, 1.5}}, 1.4);
    require(sphere.status == "pass", "single sphere build failed");
    require(sphere.contacts.size() == 1, "single sphere contact count");
    require(sphere.tori.empty(), "single sphere has no torus");
    require(sphere.reentrant.empty(), "single sphere has no reentrant probe");
    const auto sphere_hits = intersect_ray_x(sphere, 0.0, 0.0);
    require(sphere_hits.size() == 2, "central ray must hit sphere twice");
    require(std::abs(sphere_hits[0].x + 1.5) < 1.0e-10 &&
                std::abs(sphere_hits[1].x - 1.5) < 1.0e-10,
            "contact ray intersections must be analytic");
    require(sphere_hits[0].normal.x < 0.0 && sphere_hits[1].normal.x > 0.0,
            "contact ray normals must point outwards");

    const std::vector<Atom> pair = {{-1.5, 0.0, 0.0, 0.0, 1.0},
                                    {1.5, 0.0, 0.0, 0.0, 1.0}};
    const PatchComplex two = build_patch_complex(pair, 1.0);
    require(two.status == "pass", "two-atom build failed");
    require(two.tori.size() == 1, "two atoms require one toroidal patch");
    require(two.tori.front().exposed_intervals.size() == 1,
            "isolated pair circle must be fully exposed");
    require(std::abs(exposed_angle(two) - 2.0 * std::acos(-1.0)) < 1.0e-10,
            "isolated torus must cover a full probe-centre circle");
    const double torus_major = two.tori.front().circle_radius;
    const auto torus_hits = intersect_ray_x(two, torus_major - 0.5, 0.0);
    require(torus_hits.size() >= 2,
            "ray through the saddle band must hit analytic torus");

    const double root3 = std::sqrt(3.0);
    const std::vector<Atom> triple = {{-1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                      {1.5, -0.5 * root3, 0.0, 0.0, 1.0},
                                      {0.0, root3, 0.0, 0.0, 1.0}};
    const PatchComplex three = build_patch_complex(triple, 1.0);
    require(three.status == "pass", "three-atom build failed");
    require(three.candidate_triples == 1, "one geometric triple expected");
    require(three.reentrant.size() == 2,
            "equilateral triple requires two empty probe centres");

    const PatchComplex shifted =
        build_patch_complex(translate(triple, 31.25, -17.5, 9.75), 1.0);
    require(shifted.tori.size() == three.tori.size(),
            "translation changed torus count");
    require(shifted.reentrant.size() == three.reentrant.size(),
            "translation changed reentrant count");
    require(std::abs(exposed_angle(shifted) - exposed_angle(three)) < 1.0e-10,
            "translation changed exposed alpha-edge measure");
    const PatchComplex rotated =
        build_patch_complex(rotate_z(triple, 0.371), 1.0);
    require(rotated.tori.size() == three.tori.size() &&
                rotated.reentrant.size() == three.reentrant.size(),
            "rotation changed analytic patch topology");
    require(std::abs(exposed_angle(rotated) - exposed_angle(three)) < 1.0e-10,
            "rotation changed exposed alpha-edge measure");

    // Four expanded atom spheres meet at the origin.  Each of the four atom
    // triples owns a different reentrant face there; centre-only deduplication
    // would silently retain just one quarter of this degenerate polygon.
    const double tetrahedron_scale = 2.0 / std::sqrt(3.0);
    const std::vector<Atom> four_contact = {
        {tetrahedron_scale, tetrahedron_scale, tetrahedron_scale, 0.0, 1.0},
        {tetrahedron_scale, -tetrahedron_scale, -tetrahedron_scale, 0.0, 1.0},
        {-tetrahedron_scale, tetrahedron_scale, -tetrahedron_scale, 0.0, 1.0},
        {-tetrahedron_scale, -tetrahedron_scale, tetrahedron_scale, 0.0, 1.0}};
    const PatchComplex four = build_patch_complex(four_contact, 1.0);
    std::size_t shared_centre_faces = 0U;
    for (const auto &patch : four.reentrant) {
      if (std::hypot(patch.probe_center.x, patch.probe_center.y,
                     patch.probe_center.z) < 1.0e-10)
        ++shared_centre_faces;
    }
    require(four.status == "pass" && four.candidate_triples == 4U &&
                shared_centre_faces == 4U,
            "four-atom co-contact lost reentrant spherical-polygon faces");

    TitanExecutor executor(8);
    const PatchComplex parallel =
        build_patch_complex(triple, 1.0, &executor, false);
    require(parallel.status == "pass", "parallel patch build failed");
    require(parallel.ray_buckets.empty(),
            "patch-ray build retained unused legacy ray buckets");
    require(parallel.candidate_pairs == three.candidate_pairs &&
                parallel.candidate_triples == three.candidate_triples &&
                parallel.tori.size() == three.tori.size() &&
                parallel.reentrant.size() == three.reentrant.size(),
            "parallel slab/halo patch counts differ from serial build");
    for (std::size_t index = 0; index < three.reentrant.size(); ++index) {
      const auto &a = three.reentrant[index].probe_center;
      const auto &b = parallel.reentrant[index].probe_center;
      require(std::hypot(a.x - b.x, a.y - b.y, a.z - b.z) < 1.0e-12,
              "parallel slab/halo patch order is not deterministic");
    }

    std::cout << "tabi_analytic_ses_patch_smoke=pass contacts="
              << three.contacts.size() << " tori=" << three.tori.size()
              << " reentrant=" << three.reentrant.size()
              << " four_contact_faces=" << shared_centre_faces
              << " parallel_tori=" << parallel.tori.size() << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_analytic_ses_patch_smoke=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
