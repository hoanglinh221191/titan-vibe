#include "tabi/tabi_calderon_aux_surface.h"
#include "tabi/tabi_calderon_prepared.h"
#include "tabi/tabi_calderon_weak_form.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
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
using TitanPBE::Tabi::calderon::CompleteMixedJcpConfig;
using TitanPBE::Tabi::calderon::PreparedCalderonConfig;
using TitanPBE::Tabi::calderon::apply_complete_mixed_jcp_dense_reference;
using TitanPBE::Tabi::calderon::build_implicit_vdw_auxiliary_surface;
using TitanPBE::Tabi::calderon::prepare_fmm_native_calderon_operator;

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::vector<Atom> read_pqr(const std::string &path) {
  std::ifstream input(path);
  require(input.is_open(), "failed to open PQR fixture");
  std::vector<Atom> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM"))) {
      continue;
    }
    std::istringstream parser(line);
    std::vector<double> numeric;
    std::string field;
    while (parser >> field) {
      char *end = nullptr;
      const double value = std::strtod(field.c_str(), &end);
      if (end != field.c_str() && *end == '\0') {
        numeric.push_back(value);
      }
    }
    require(numeric.size() >= 7, "invalid PQR row");
    const std::size_t offset = numeric.size() - 5;
    atoms.push_back({numeric[offset], numeric[offset + 1], numeric[offset + 2],
                     numeric[offset + 3], numeric[offset + 4]});
  }
  require(!atoms.empty(), "PQR fixture contains no atoms");
  return atoms;
}

double relative_l2(const std::vector<double> &left,
                   const std::vector<double> &right) {
  require(left.size() == right.size(), "relative L2 shape mismatch");
  double numerator = 0.0;
  double denominator = 0.0;
  for (std::size_t i = 0; i < left.size(); ++i) {
    const double difference = left[i] - right[i];
    numerator += difference * difference;
    denominator += right[i] * right[i];
  }
  return std::sqrt(numerator / std::max(1.0e-300, denominator));
}

std::vector<double> deterministic_probe(std::size_t count, int seed) {
  std::vector<double> values(count, 0.0);
  for (std::size_t i = 0; i < count; ++i) {
    const double phase = static_cast<double>((i + 1) * (seed + 3));
    values[i] = std::sin(0.173 * phase) + 0.35 * std::cos(0.071 * phase);
  }
  return values;
}

PreparedCalderonConfig make_config(std::size_t primary_count,
                                   double area_mismatch) {
  PreparedCalderonConfig config;
  config.eps_ratio = 80.0;
  config.kappa = 0.125;
  config.leaf_size = 64;
  config.opening_theta = 0.6;
  config.expansion_order = 4;
  config.transfer.geometry.nearest_primary_panels =
      std::min(primary_count, area_mismatch > 0.15 ? std::size_t{64}
                                                   : std::size_t{32});
  config.transfer.geometry.geometry_sample_cap = 512;
  config.transfer.geometry.max_primary_to_auxiliary_distance = 2.5;
  config.transfer.geometry.max_auxiliary_to_primary_distance = 2.5;
  config.transfer.geometry.min_normal_dot = 0.0;
  config.transfer.geometry.min_normal_alignment_fraction = 0.90;
  config.transfer.geometry.max_relative_area_mismatch = 0.75;
  config.transfer.geometry.normal_candidate_distance_slack = 0.25;
  config.transfer.max_balance_iterations = 8000;
  config.transfer.balance_relative_tolerance = 1.0e-9;
  return config;
}

} // namespace

int main(int argc, char **argv) {
  try {
    const bool protein = argc >= 2;
    const std::vector<Atom> atoms =
        protein ? read_pqr(argv[1])
                : std::vector<Atom>{{0.0, 0.0, 0.0, 0.0, 1.5}};
    Options options;
    options.surface_geometry = "triangle";
    options.surface_mode = "vdw";
    options.surface_point_count = protein ? 8 : 32;
    options.sphere_mesh_mode = "icosphere";
    options.diagonal_mode = "solid-angle";
    Result surface_result;
    const std::vector<Panel> primary =
        build_surface(atoms, options, &surface_result, nullptr);
    const auto auxiliary =
        build_implicit_vdw_auxiliary_surface(atoms, protein ? 1.25 : 1.0);
    require(!primary.empty(), "primary surface must not be empty");
    require(auxiliary.status == "pass", "auxiliary surface must pass");
    const double primary_area = std::accumulate(
        primary.begin(), primary.end(), 0.0,
        [](double sum, const Panel &panel) { return sum + panel.area; });
    const double area_mismatch =
        std::abs(primary_area - auxiliary.trace_spaces.quality.total_area) /
        auxiliary.trace_spaces.quality.total_area;
    const PreparedCalderonConfig config =
        make_config(primary.size(), area_mismatch);

    TitanExecutor executor4(4);
    auto prepared = prepare_fmm_native_calderon_operator(
        primary, auxiliary, config, executor4);
    require(prepared.status() == "prepared", "prepared operator must build");
    require(prepared.geometry_status() == "pass",
            "prepared operator geometry must pass");
    require(prepared.prepare_build_count() == 1 &&
                prepared.context_build_count() == 1,
            "prepared lifecycle must build transfer/context once");

    const std::vector<double> primary_probe =
        deterministic_probe(2 * primary.size(), 11);
    const auto first = prepared.apply_primary(primary_probe, executor4);
    const auto second = prepared.apply_primary(primary_probe, executor4);
    require(first.status == "ok" && second.status == "ok" && first.finite &&
                second.finite,
            "prepared primary applies must be finite");
    require(relative_l2(first.values, second.values) <= 1.0e-14,
            "repeated prepared applies must be deterministic");
    require(prepared.prepare_build_count() == 1 &&
                prepared.context_build_count() == 1 &&
                prepared.apply_count() == 2,
            "repeated applies must not rebuild prepared state");

    double dense_relative = -1.0;
    double dense_phi_relative = -1.0;
    double dense_phin_relative = -1.0;
    std::string dense_status = "not_run_on_protein";
    if (!protein) {
      const std::size_t n = auxiliary.mesh.vertices.size();
      const std::vector<double> auxiliary_probe =
          deterministic_probe(2 * n, 23);
      const auto candidate = prepared.apply_auxiliary(auxiliary_probe, executor4);
      CompleteMixedJcpConfig dense_config;
      dense_config.eps_ratio = config.eps_ratio;
      dense_config.kappa = config.kappa;
      dense_config.adjacent_pair_subdivision_depth = 2;
      dense_config.direct_action_relative_tolerance = 0.025;
      const auto dense = apply_complete_mixed_jcp_dense_reference(
          auxiliary.mesh, auxiliary.trace_spaces, auxiliary_probe, dense_config);
      dense_status = dense.status;
      require(candidate.status == "ok" && dense.status == "pass",
              "candidate and dense auxiliary actions must run");
      dense_relative = relative_l2(candidate.values, dense.action);
      dense_phi_relative = relative_l2(
          std::vector<double>(candidate.values.begin(),
                              candidate.values.begin() +
                                  static_cast<std::ptrdiff_t>(n)),
          std::vector<double>(dense.action.begin(),
                              dense.action.begin() +
                                  static_cast<std::ptrdiff_t>(n)));
      dense_phin_relative = relative_l2(
          std::vector<double>(candidate.values.begin() +
                                  static_cast<std::ptrdiff_t>(n),
                              candidate.values.end()),
          std::vector<double>(dense.action.begin() +
                                  static_cast<std::ptrdiff_t>(n),
                              dense.action.end()));

      const std::vector<double> a = deterministic_probe(2 * n, 31);
      const std::vector<double> b = deterministic_probe(2 * n, 47);
      std::vector<double> combined(2 * n, 0.0);
      for (std::size_t i = 0; i < combined.size(); ++i) {
        combined[i] = 0.7 * a[i] - 1.2 * b[i];
      }
      const auto aa = prepared.apply_auxiliary(a, executor4);
      const auto bb = prepared.apply_auxiliary(b, executor4);
      const auto cc = prepared.apply_auxiliary(combined, executor4);
      std::vector<double> expected(cc.values.size(), 0.0);
      for (std::size_t i = 0; i < expected.size(); ++i) {
        expected[i] = 0.7 * aa.values[i] - 1.2 * bb.values[i];
      }
      require(relative_l2(cc.values, expected) <= 1.0e-12,
              "prepared auxiliary action must be linear");

      TitanExecutor executor1(1);
      auto serial = prepare_fmm_native_calderon_operator(
          primary, auxiliary, config, executor1);
      require(serial.status() == "prepared",
              "serial prepared operator must build");
      const auto serial_action = serial.apply_auxiliary(auxiliary_probe, executor1);
      require(serial_action.status == "ok" &&
                  relative_l2(serial_action.values, candidate.values) <= 1.0e-12,
              "prepared action must be deterministic for nt=1/4");
    }

    std::cout << "calderon_c5d1_prepared_diagnostic"
              << " status=" << prepared.status()
              << " geometry_status=" << prepared.geometry_status()
              << " input_atoms=" << atoms.size()
              << " primary_panels=" << prepared.primary_panels()
              << " auxiliary_vertices=" << prepared.auxiliary_vertices()
              << " fmm_quadrature_panels=" << prepared.fmm_quadrature_panels()
              << " prepare_build_count=" << prepared.prepare_build_count()
              << " context_build_count=" << prepared.context_build_count()
              << " apply_count=" << prepared.apply_count()
              << " prepare_seconds=" << prepared.prepare_seconds()
              << " transfer_seconds=" << prepared.transfer_seconds()
              << " fmm_setup_seconds=" << prepared.fmm_setup_seconds()
              << " first_apply_seconds=" << first.apply_seconds
              << " second_apply_seconds=" << second.apply_seconds
              << " pairing_seconds=" << second.pairing_seconds
              << " pairing_iterations=" << second.pairing_iterations_max
              << " near_pairs=" << second.near_pairs
              << " far_pairs=" << second.far_pairs
              << " estimated_bytes=" << prepared.estimated_bytes()
              << " dense_status=" << dense_status
              << " dense_relative=" << dense_relative
              << " dense_phi_relative=" << dense_phi_relative
              << " dense_phin_relative=" << dense_phin_relative << '\n';
    std::cout << "calderon_c5d1_prepared_apply_smoke=pass\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "calderon_c5d1_prepared_apply_smoke=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
