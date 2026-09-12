#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_transfer.h"
#include "tabi/tabi_types.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Options;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::Result;
using TitanPBE::Tabi::build_surface;
using TitanPBE::Tabi::calderon::GeometryNormalMatching;
using TitanPBE::Tabi::calderon::PrimaryAuxiliaryTransferConfig;
using TitanPBE::Tabi::calderon::PrimaryAuxiliaryTransferResult;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;
using TitanPBE::Tabi::calderon::build_primary_auxiliary_trace_transfer;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::vector<Atom> read_pqr(const std::string &path) {
  std::ifstream input(path);
  if (!input.is_open()) {
    throw std::runtime_error("failed to open PQR fixture");
  }
  std::vector<Atom> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM"))) {
      continue;
    }
    std::istringstream parser(line);
    std::vector<std::string> token;
    std::string value;
    while (parser >> value) {
      token.push_back(value);
    }
    std::vector<double> numeric;
    for (const std::string &field : token) {
      char *end = nullptr;
      const double parsed = std::strtod(field.c_str(), &end);
      if (end != field.c_str() && *end == '\0') {
        numeric.push_back(parsed);
      }
    }
    if (numeric.size() < 7) {
      throw std::runtime_error("invalid PQR coordinate/charge/radius row");
    }
    const std::size_t offset = numeric.size() - 5;
    atoms.push_back({numeric[offset], numeric[offset + 1], numeric[offset + 2],
                     numeric[offset + 3], numeric[offset + 4]});
  }
  require(!atoms.empty(), "PQR fixture contains no atoms");
  return atoms;
}

void write_mismatches(const std::string &path,
                      const PrimaryAuxiliaryTransferResult &transfer) {
  if (path.empty()) {
    return;
  }
  std::ofstream output(path);
  require(output.is_open(), "failed to create mismatch CSV");
  output << "direction,source_index,source_atom_index,matched_index,"
            "matched_atom_index,distance,normal_dot,same_owner,"
            "owner_ambiguous,positive_alternative\n";
  output << std::setprecision(17);
  for (const auto &sample : transfer.normal_mismatch_samples) {
    output << sample.direction << ',' << sample.source_index << ','
           << sample.source_atom_index << ',' << sample.matched_index << ','
           << sample.matched_atom_index << ',' << sample.distance << ','
           << sample.normal_dot << ',' << (sample.same_owner ? "true" : "false")
           << ',' << (sample.owner_ambiguous ? "true" : "false") << ','
           << (sample.positive_alternative ? "true" : "false") << '\n';
  }
}

void print_result(const PrimaryAuxiliaryTransferResult &transfer,
                  std::size_t atom_count, double slack) {
  std::cout << std::setprecision(10)
            << "calderon_c5c7_normal_diagnostic"
            << " status=" << transfer.status
            << " policy=" << transfer.normal_matching_policy
            << " slack=" << slack << " input_atoms=" << atom_count
            << " primary_panels=" << transfer.primary_panels
            << " auxiliary_vertices=" << transfer.auxiliary_vertices
            << " area_rel=" << transfer.relative_area_mismatch
            << " p2a_rms=" << transfer.primary_to_auxiliary_rms_distance
            << " p2a_max=" << transfer.primary_to_auxiliary_max_distance
            << " a2p_rms=" << transfer.auxiliary_to_primary_rms_distance
            << " a2p_max=" << transfer.auxiliary_to_primary_max_distance
            << " p2a_coverage=" << transfer.primary_coverage_fraction
            << " a2p_coverage=" << transfer.auxiliary_coverage_fraction
            << " p2a_normal_fraction="
            << transfer.primary_normal_alignment_fraction
            << " a2p_normal_fraction="
            << transfer.auxiliary_normal_alignment_fraction
            << " p2a_weighted_fraction="
            << transfer.primary_area_weighted_normal_alignment_fraction
            << " a2p_weighted_fraction="
            << transfer.auxiliary_area_weighted_normal_alignment_fraction
            << " p2a_same_owner=" << transfer.primary_same_owner_fraction
            << " a2p_same_owner=" << transfer.auxiliary_same_owner_fraction
            << " p2a_ambiguity=" << transfer.primary_owner_ambiguity_fraction
            << " a2p_ambiguity=" << transfer.auxiliary_owner_ambiguity_fraction
            << " p2a_negative=" << transfer.primary_negative_normal_samples
            << " a2p_negative=" << transfer.auxiliary_negative_normal_samples
            << " p2a_positive_alternative="
            << transfer.primary_positive_alternative_samples
            << " a2p_positive_alternative="
            << transfer.auxiliary_positive_alternative_samples
            << " mismatch_records=" << transfer.normal_mismatch_samples.size()
            << " forward_constant=" << transfer.forward_constant_relative
            << " adjoint=" << transfer.adjoint_relative
            << " build_ms=" << transfer.build_ms << '\n';
}

PrimaryAuxiliaryTransferResult run_transfer(
    const std::vector<Panel> &primary,
    const TitanPBE::Tabi::calderon::AuxiliarySurfaceResult &auxiliary,
    GeometryNormalMatching policy, double slack, bool fixture_mode) {
  const double primary_area = std::accumulate(
      primary.begin(), primary.end(), 0.0,
      [](double sum, const Panel &panel) { return sum + panel.area; });
  const double auxiliary_area = auxiliary.trace_spaces.quality.total_area;
  const double area_mismatch = std::abs(primary_area - auxiliary_area) /
                               std::max(1.0e-300, auxiliary_area);
  PrimaryAuxiliaryTransferConfig config;
  config.nearest_primary_panels =
      !fixture_mode ? 8 : (area_mismatch > 0.15 ? 64 : 32);
  config.geometry_sample_cap = fixture_mode ? 512 : 0;
  config.max_primary_to_auxiliary_distance = fixture_mode ? 2.50 : 0.75;
  config.max_auxiliary_to_primary_distance = fixture_mode ? 2.50 : 0.75;
  config.min_normal_dot = 0.0;
  config.min_normal_alignment_fraction = fixture_mode ? 0.90 : 1.0;
  config.max_relative_area_mismatch = fixture_mode ? 0.75 : 0.25;
  config.normal_matching = policy;
  config.normal_candidate_distance_slack = slack;
  config.require_area_weighted_normal_alignment = true;
  config.mismatch_sample_cap = 2048;
  return build_primary_auxiliary_trace_transfer(primary, auxiliary, config);
}

} // namespace

int main(int argc, char **argv) {
  try {
    const bool fixture_mode = argc >= 2;
    const std::vector<Atom> atoms = fixture_mode
                                        ? read_pqr(argv[1])
                                        : std::vector<Atom>{{-0.9, 0.0, 0.0, 0.0, 1.5},
                                                            {0.9, 0.0, 0.0, 0.0, 1.5}};
    const std::string policy_name = argc >= 3 ? argv[2] : "compare";
    const double slack = argc >= 4 ? std::stod(argv[3]) : 0.25;
    const std::string mismatch_path = argc >= 5 ? argv[4] : "";
    require(slack >= 0.0 && slack <= 1.0, "invalid matching distance slack");

    Options options;
    options.surface_geometry = "triangle";
    options.surface_mode = "vdw";
    options.surface_point_count = fixture_mode ? 8 : 32;
    options.sphere_mesh_mode = "icosphere";
    options.diagonal_mode = "solid-angle";
    Result surface_result;
    const std::vector<Panel> primary =
        build_surface(atoms, options, &surface_result, nullptr);
    const auto auxiliary =
        build_implicit_vdw_auxiliary_surface(atoms, fixture_mode ? 1.25 : 0.30);
    require(!primary.empty(), "production primary surface must contain panels");
    require(auxiliary.status == "pass", "auxiliary surface must be watertight");

    if (policy_name == "compare") {
      const auto legacy = run_transfer(primary, auxiliary,
                                       GeometryNormalMatching::FacetNearest,
                                       slack, fixture_mode);
      const auto corrected = run_transfer(
          primary, auxiliary, GeometryNormalMatching::OwnerSdfAware, slack,
          fixture_mode);
      const auto repeated = run_transfer(
          primary, auxiliary, GeometryNormalMatching::OwnerSdfAware, slack,
          fixture_mode);
      print_result(legacy, atoms.size(), slack);
      print_result(corrected, atoms.size(), slack);
      require(corrected.primary_normal_alignment_fraction + 1.0e-15 >=
                  legacy.primary_normal_alignment_fraction,
              "owner/SDF matching must not reduce primary alignment");
      require(corrected.auxiliary_normal_alignment_fraction + 1.0e-15 >=
                  legacy.auxiliary_normal_alignment_fraction,
              "owner/SDF matching must not reduce auxiliary alignment");
      require(corrected.forward_constant_relative <= 1.0e-13 &&
                  corrected.adjoint_relative <= 1.0e-13,
              "geometry matching must preserve transfer algebra");
      require(corrected.status == repeated.status &&
                  corrected.primary_normal_alignment_fraction ==
                      repeated.primary_normal_alignment_fraction &&
                  corrected.auxiliary_normal_alignment_fraction ==
                      repeated.auxiliary_normal_alignment_fraction &&
                  corrected.primary_to_auxiliary_max_distance ==
                      repeated.primary_to_auxiliary_max_distance &&
                  corrected.auxiliary_to_primary_max_distance ==
                      repeated.auxiliary_to_primary_max_distance &&
                  corrected.normal_mismatch_samples.size() ==
                      repeated.normal_mismatch_samples.size(),
              "owner/SDF geometry attribution must be deterministic");
    } else {
      const GeometryNormalMatching policy =
          policy_name == "owner-sdf-aware"
              ? GeometryNormalMatching::OwnerSdfAware
              : GeometryNormalMatching::FacetNearest;
      require(policy_name == "owner-sdf-aware" || policy_name == "facet-nearest",
              "unsupported normal matching policy");
      const auto transfer =
          run_transfer(primary, auxiliary, policy, slack, fixture_mode);
      print_result(transfer, atoms.size(), slack);
      write_mismatches(mismatch_path, transfer);
      require(transfer.forward_constant_relative <= 1.0e-13 &&
                  transfer.adjoint_relative <= 1.0e-13,
              "geometry matching must preserve transfer algebra");
    }
    std::cout << "calderon_c5c7_normal_attribution_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c7_normal_attribution_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
