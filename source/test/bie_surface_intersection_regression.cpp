// SES self-intersection audit on triangle pairs that share a vertex.
//
// 1. A true overlap (pair 157908/157910 of the 13,541-atom frame-0 mesh;
//    exact rational plane slicing finds 0.27774364742302693 A of overlap
//    beyond the shared vertex) must still be counted.
// 2. Two pairs that meet only at their shared vertex used to be counted when
//    an ill-conditioned line reconstruction returned an invalid interval.
//    Exact rational axis separation proves the contact is the vertex alone;
//    the interval certificate (tabi_shared_vertex_separation.h) must clear
//    them.
// Each case is repeated under rotation, scaling and translation.
//
// Build: link with the engine objects (test/pbe_p0_oracles.py --case ses).
#include "tabi/tabi_calderon_aux_surface.h"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace TitanPBE::Tabi;

namespace {

// Rotate about z, scale and shift, which must not change any classification.
void transform(SurfaceMesh &mesh, double scale, double shift) {
  for (auto &v : mesh.vertices) {
    const double x = v.x, y = v.y, z = v.z;
    v.x = shift + scale * (0.6 * x - 0.8 * y);
    v.y = -shift + scale * (0.8 * x + 0.6 * y);
    v.z = shift + scale * z;
  }
}

} // namespace

int main() {
  try {
    int checks = 0;

    SurfaceMesh overlap;
    overlap.vertices.push_back({132.4748905848725, 116.83094530036095,
                                99.47755698628418, 0.09312726855093977,
                                0.6948641996030784, 0.7130855881044115, 0});
    overlap.vertices.push_back({133.34437687212417, 116.80163679883069,
                                99.54127516661086, -0.09264683512887889,
                                0.27647491578767797, 0.9565449204720075, 0});
    overlap.vertices.push_back({133.4023951607255, 116.91161190376646,
                                99.52103340300323, -0.13408846984410633,
                                0.19792126940498528, 0.9710033230488891, 0});
    overlap.vertices.push_back({132.66793402724315, 116.1499677469271,
                                99.0503050378997, 0.7051925041821494,
                                -0.6639527003952129, -0.2487374995516636, 0});
    overlap.vertices.push_back({132.69495367929557, 116.7173717968462,
                                99.59739643954079, -0.06406065603696559,
                                0.7759881306850481, 0.62748597863541, 0});
    overlap.triangles.push_back({0, 1, 2, 0});
    overlap.triangles.push_back({3, 4, 2, 0});
    for (const double scale : {0.001, 1.0, 1000.0})
      for (const double shift : {0.0, 10000.0}) {
        SurfaceMesh mesh = overlap;
        transform(mesh, scale, shift);
        if (calderon::audit_surface_self_intersections(mesh, 0.75 * scale)
                .count != 1)
          throw std::runtime_error("real intersection was missed");
        ++checks;
      }

    const std::array<std::array<Vertex, 5>, 2> separated{{
        {{{156.06712362033372, 102.36998970359356, 90.23659981699544},
          {155.95144221103467, 102.47165170834094, 90.38705372618462},
          {155.95144221103467, 102.4716621284233, 90.38708307037803},
          {155.95144221103467, 102.39960681876636, 90.15884007062257},
          {156.24131894620933, 102.24341912866785, 90.15884007062257}}},
        {{{1.6827148437500001, -1.0376998771057635, -0.5777389396057632},
          {1.6825123398858237, -1.0377025038641767, -0.5777415663641764},
          {1.6827148437500001, -1.037812906054384, -0.5775390624999996},
          {1.6832542808301802, -1.0376929847948198, -0.5777320472948199},
          {1.6830619424065199, -1.037695424770718, -0.5777344872707167}}},
    }};
    for (const auto &points : separated)
      for (const double scale : {1.0, 1000.0})
        for (const double shift : {0.0, 10000.0}) {
          SurfaceMesh mesh;
          mesh.vertices.assign(points.begin(), points.end());
          transform(mesh, scale, shift);
          mesh.triangles.push_back({0, 1, 2, 0});
          mesh.triangles.push_back({0, 3, 4, 0});
          if (calderon::audit_surface_self_intersections(mesh, 0.75 * scale)
                  .count != 0)
            throw std::runtime_error(
                "shared-vertex-only contact misclassified");
          ++checks;
        }

    std::cout << "PASS " << checks
              << " true-intersection and certified-separation controls\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
  }
}
