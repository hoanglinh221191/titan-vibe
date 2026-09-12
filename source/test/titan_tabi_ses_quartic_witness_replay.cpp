#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_intersections.h"
#include "tabi/tabi_ses_quartic.h"
#include "titan_executor.h"

#include <bit>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

constexpr std::size_t kAtomCap = 8192U;
constexpr std::size_t kTorusIndex = 6693U;
constexpr double kProbeRadius = 1.4;
constexpr double kRayU = -129.77699999999999;
constexpr double kRayV = -128.863;

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

std::string trimmed(std::string value) {
  const auto first = value.find_first_not_of(' ');
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(' ');
  return value.substr(first, last - first + 1U);
}

double pdb_radius(std::string element) {
  element = trimmed(std::move(element));
  if (element == "H")
    return 1.20;
  if (element == "N")
    return 1.55;
  if (element == "O")
    return 1.52;
  if (element == "S" || element == "P")
    return 1.80;
  return 1.70;
}

std::vector<Atom> read_atoms(const std::string &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open input");
  std::vector<Atom> atoms;
  atoms.reserve(kAtomCap);
  std::string line;
  while (std::getline(input, line) && atoms.size() < kAtomCap) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM")))
      continue;
    if (line.size() < 54U)
      throw std::runtime_error("malformed selected PDB atom");
    Atom atom;
    atom.x = std::stod(line.substr(30U, 8U));
    atom.y = std::stod(line.substr(38U, 8U));
    atom.z = std::stod(line.substr(46U, 8U));
    atom.q = 0.0;
    atom.radius =
        pdb_radius(line.size() >= 78U ? line.substr(76U, 2U) : std::string{});
    if (!std::isfinite(atom.x) || !std::isfinite(atom.y) ||
        !std::isfinite(atom.z) || !std::isfinite(atom.radius) ||
        atom.radius <= 0.0)
      throw std::runtime_error("invalid selected PDB atom");
    atoms.push_back(atom);
  }
  require(atoms.size() == kAtomCap, "PDB prefix atom count mismatch");
  return atoms;
}

const char *filter_name(TorusRootFilterReason reason) noexcept {
  switch (reason) {
  case TorusRootFilterReason::accepted:
    return "accepted";
  case TorusRootFilterReason::radial_axis_singular:
    return "radial_axis_singular";
  case TorusRootFilterReason::theta_interval:
    return "theta_interval";
  case TorusRootFilterReason::tube_consistency:
    return "tube_consistency";
  case TorusRootFilterReason::beta_minor_arc:
    return "beta_minor_arc";
  case TorusRootFilterReason::reentrant_coverage:
    return "reentrant_coverage";
  case TorusRootFilterReason::tangential:
    return "tangential";
  }
  return "invalid";
}

bool same_roots(const PatchAxisRoots &left,
                const PatchAxisRoots &right) noexcept {
  if (left.count != right.count || left.certified != right.certified)
    return false;
  for (std::size_t index = 0U; index < left.coordinates.size(); ++index) {
    if (std::bit_cast<std::uint64_t>(left.coordinates[index]) !=
            std::bit_cast<std::uint64_t>(right.coordinates[index]) ||
        left.primitive_simplicity_certified[index] !=
            right.primitive_simplicity_certified[index])
      return false;
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  try {
    require(argc == 2, "usage: witness_replay <huge.pdb>");
    const std::vector<Atom> atoms = read_atoms(argv[1]);
    TitanExecutor executor(8U);
    const PatchComplex complex =
        build_patch_complex(atoms, kProbeRadius, &executor, false);
    require(complex.status == "pass" && complex.tori.size() > kTorusIndex,
            "prefix patch complex did not reproduce torus index");
    const ReentrantProbeIndex probes = build_reentrant_probe_index(complex);
    const PatchAxisRoots normal =
        intersect_patch_axis(complex, probes, PatchKind::toroidal, kTorusIndex,
                             RayAxis::y, kRayU, kRayV);
    QuarticDiagnostics quartic;
    TorusPatchAxisDiagnostics filters;
    const PatchAxisRoots audited =
        intersect_patch_axis_with_quartic_diagnostics(
            complex, probes, PatchKind::toroidal, kTorusIndex, RayAxis::y,
            kRayU, kRayV, &quartic, &filters);
    std::cout << std::setprecision(17) << "replay_quartic_coefficients";
    for (std::size_t index = 0U; index < quartic.input_coefficients.size();
         ++index) {
      std::cout << " c" << index << '=' << quartic.input_coefficients[index]
                << ":0x" << std::hex
                << std::bit_cast<std::uint64_t>(
                       quartic.input_coefficients[index])
                << std::dec;
    }
    std::cout << " roots";
    for (std::size_t root = 0U; root < quartic.output_root_count; ++root) {
      std::cout << " r" << root << '=' << quartic.output_roots[root] << ":0x"
                << std::hex
                << std::bit_cast<std::uint64_t>(quartic.output_roots[root])
                << std::dec;
    }
    const ToroidalPatch &torus = complex.tori[kTorusIndex];
    std::cout << " contacts=" << complex.contacts.size()
              << " tori=" << complex.tori.size()
              << " reentrant=" << complex.reentrant.size()
              << " atom_i=" << torus.atom_i << " atom_j=" << torus.atom_j
              << " exposed_intervals=" << torus.exposed_intervals.size()
              << '\n';
    require(same_roots(normal, audited),
            "audit overload changed patch-axis result bits");
    require(quartic.output_root_count == 2U && quartic.failure_mask == 0U,
            "target quartic was not certified after refinement fix");
    const std::array<std::uint64_t, 5> expected_coefficient_bits = {
        0xbf00555c1486006bULL, 0xc04508da863a9e0eULL, 0x4033d130970fc129ULL,
        0x0000000000000000ULL, 0x3ff0000000000000ULL};
    for (std::size_t index = 0U; index < expected_coefficient_bits.size();
         ++index) {
      require(std::bit_cast<std::uint64_t>(quartic.input_coefficients[index]) ==
                  expected_coefficient_bits[index],
              "target quartic coefficient snapshot changed");
    }
    require(std::bit_cast<std::uint64_t>(quartic.output_roots[0]) ==
                    0xbea8d91e15703802ULL &&
                std::bit_cast<std::uint64_t>(quartic.output_roots[1]) ==
                    0x3ffd1b0f13aa87d8ULL,
            "target quartic output roots changed");
    require(normal.certified && normal.count == 0U,
            "target patch-axis filter result changed");
    require(filters.root_count == quartic.output_root_count,
            "target filter diagnostics lost primitive roots");

    std::cout << std::setprecision(17)
              << "tabi_ses_quartic_witness_replay=pass atoms=" << atoms.size()
              << " contacts=" << complex.contacts.size()
              << " tori=" << complex.tori.size()
              << " reentrant=" << complex.reentrant.size()
              << " patch_axis_count=" << normal.count
              << " patch_axis_certified=" << (normal.certified ? 1 : 0)
              << " quartic_count=" << quartic.output_root_count
              << " quartic_max_residual="
              << quartic.output_max_normalized_residual << '\n';
    for (std::size_t root = 0U; root < filters.root_count; ++root) {
      const TorusRootFilterDiagnostic &filter = filters.roots[root];
      std::cout << "root=" << root << " value=" << quartic.output_roots[root]
                << " value_bits=0x" << std::hex
                << std::bit_cast<std::uint64_t>(quartic.output_roots[root])
                << std::dec << " pre_residual="
                << quartic.pre_cast_normalized_residuals[root]
                << " post_residual="
                << quartic.post_cast_normalized_residuals[root]
                << " filter=" << filter_name(filter.reason)
                << " coordinate=" << filter.coordinate
                << " theta=" << filter.theta << " beta=" << filter.beta
                << " tube_residual=" << filter.tube_residual
                << " directional_normal=" << filter.directional_normal << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "tabi_ses_quartic_witness_replay=fail reason=" << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
