#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_transfer.h"
#include "tabi/tabi_types.h"

#include <cstdlib>
#include <fstream>
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
using TitanPBE::Tabi::calderon::MassCompatibleTransferConfig;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;
using TitanPBE::Tabi::calderon::build_mass_compatible_primary_auxiliary_trace_transfer;

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
    if (token.size() < 10) {
      continue;
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
  if (atoms.empty()) {
    throw std::runtime_error("PQR fixture contains no atoms");
  }
  return atoms;
}

} // namespace

int main(int argc, char **argv) {
  try {
    const bool fixture_mode = argc >= 2;
    const std::vector<Atom> atoms = fixture_mode
                                        ? read_pqr(argv[1])
                                        : std::vector<Atom>{{-0.9, 0.0, 0.0, 0.0, 1.5},
                                                            {0.9, 0.0, 0.0, 0.0, 1.5}};
    const double grid_spacing = fixture_mode ? 1.25 : 0.30;
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
        build_implicit_vdw_auxiliary_surface(atoms, grid_spacing);
    require(!primary.empty(), "production primary surface must contain panels");
    require(auxiliary.status == "pass", "auxiliary surface must be watertight");

    MassCompatibleTransferConfig config;
    const double primary_area = std::accumulate(
        primary.begin(), primary.end(), 0.0,
        [](double sum, const Panel &panel) { return sum + panel.area; });
    const double auxiliary_area = auxiliary.trace_spaces.quality.total_area;
    const double pre_transfer_area_mismatch =
        std::abs(primary_area - auxiliary_area) /
        std::max(1.0e-300, auxiliary_area);
    const std::size_t support_k =
        !fixture_mode ? 8 : (pre_transfer_area_mismatch > 0.15 ? 64 : 32);
    config.geometry.nearest_primary_panels = support_k;
    config.geometry.geometry_sample_cap = fixture_mode ? 512 : 0;
    config.geometry.max_primary_to_auxiliary_distance =
        fixture_mode ? 2.50 : 0.75;
    config.geometry.max_auxiliary_to_primary_distance =
        fixture_mode ? 2.50 : 0.75;
    config.geometry.min_normal_dot = 0.0;
    config.geometry.min_normal_alignment_fraction = fixture_mode ? 0.90 : 1.0;
    config.geometry.max_relative_area_mismatch = fixture_mode ? 0.75 : 0.25;
    config.max_balance_iterations = 8000;
    config.balance_relative_tolerance = 1.0e-9;
    const auto transfer = build_mass_compatible_primary_auxiliary_trace_transfer(
        primary, auxiliary, config);
    std::cout << "calderon_c5c6_real_primary_diagnostic status="
              << transfer.status << " geometry_status=" << transfer.geometry_status
              << " input_atoms=" << atoms.size()
              << " primary_panels=" << transfer.primary_panels
              << " auxiliary_vertices=" << transfer.auxiliary_vertices
              << " support_edges=" << transfer.support_edges
              << " support_k=" << support_k
              << " support_components=" << transfer.support_components
              << " repair_edges=" << transfer.support_repair_edges
              << " balance_iterations=" << transfer.balance_iterations
              << " balance_relative=" << transfer.balance_relative
              << " mass_scale=" << transfer.primary_mass_global_scale
              << " component_scale_min=" << transfer.component_mass_scale_min
              << " component_scale_max=" << transfer.component_mass_scale_max
              << " forward_constant=" << transfer.forward_constant_relative
              << " reverse_constant=" << transfer.reverse_constant_relative
              << " adjoint=" << transfer.adjoint_relative
              << " area_rel=" << transfer.geometry.relative_area_mismatch
              << " p2a_rms="
              << transfer.geometry.primary_to_auxiliary_rms_distance
              << " p2a_max="
              << transfer.geometry.primary_to_auxiliary_max_distance
              << " a2p_rms="
              << transfer.geometry.auxiliary_to_primary_rms_distance
              << " a2p_max="
              << transfer.geometry.auxiliary_to_primary_max_distance
              << " p2a_coverage=" << transfer.geometry.primary_coverage_fraction
              << " a2p_coverage=" << transfer.geometry.auxiliary_coverage_fraction
              << " p2a_normal_min="
              << transfer.geometry.primary_to_auxiliary_min_normal_dot
              << " a2p_normal_min="
              << transfer.geometry.auxiliary_to_primary_min_normal_dot
              << " p2a_normal_fraction="
              << transfer.geometry.primary_normal_alignment_fraction
              << " a2p_normal_fraction="
              << transfer.geometry.auxiliary_normal_alignment_fraction
              << " build_ms=" << transfer.build_ms << '\n';
    require(transfer.balance_relative <= 1.0e-9 &&
                transfer.forward_constant_relative <= 1.0e-8 &&
                transfer.reverse_constant_relative <= 1.0e-8 &&
                transfer.adjoint_relative <= 1.0e-8,
            "real-primary mass-compatible transfer algebra must pass");
    if (!fixture_mode) {
      require(transfer.status == "pass",
              "synthetic production-primary geometry must pass");
    }
    std::cout << "calderon_c5c6_real_primary_transfer_smoke=pass status="
              << transfer.status << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5c6_real_primary_transfer_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
