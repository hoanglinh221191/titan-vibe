#include "tabi/tabi_fmm_kifmm_lite.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmLiteConfig;
using titan::tabi::fmm_kifmm::KifmmLiteResult;
using titan::tabi::fmm_kifmm::KifmmLiteTargetPanel;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::run_kifmm_lite_matvec;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

double dot(Vec3 a, Vec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

KifmmPanelSource source(Vec3 p, int i) {
  return {p,
          unit({0.37 + 0.11 * i, -0.22 + 0.07 * i, 0.51 - 0.03 * i}),
          0.62 + 0.03 * (i % 5),
          ((i % 2) == 0 ? 0.35 : -0.27) + 0.015 * i,
          ((i % 3) == 0 ? -0.44 : 0.31) - 0.012 * i};
}

KifmmLiteTargetPanel target(Vec3 p, int i) {
  return {p, unit({-0.18 + 0.09 * i, 0.61 - 0.04 * i, 0.43 + 0.05 * i})};
}

std::vector<KifmmPanelSource> cluster_sources(Vec3 center, double spacing,
                                              int count, int offset) {
  std::vector<KifmmPanelSource> panels;
  panels.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const int ix = i % 3;
    const int iy = (i / 3) % 3;
    const int iz = i / 9;
    panels.push_back(source({center.x + spacing * (ix - 1),
                             center.y + spacing * (iy - 1),
                             center.z + spacing * (iz - 0.5)},
                            offset + i));
  }
  return panels;
}

std::vector<KifmmLiteTargetPanel> cluster_targets(Vec3 center, double spacing,
                                                  int count, int offset) {
  std::vector<KifmmLiteTargetPanel> panels;
  panels.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const int ix = i % 3;
    const int iy = (i / 3) % 3;
    const int iz = i / 9;
    panels.push_back(target({center.x + spacing * (ix - 1),
                             center.y + spacing * (iy - 1),
                             center.z + spacing * (iz - 0.5)},
                            offset + i));
  }
  return panels;
}

KifmmLiteConfig base_config() {
  KifmmLiteConfig config;
  config.order = 8;
  config.leaf_size = 8;
  config.opening_theta = 0.6;
  config.equiv_scale = 1.05;
  config.check_scale = 2.40;
  config.kappa = 0.12746;
  config.eps_ratio = 80.0;
  config.max_aggregate_relative_rms = 1.0e-5;
  config.max_channel_relative_rms = 1.0e-5;
  config.solver = "qr";
  return config;
}

void require_parity(const KifmmLiteResult &result, double gate,
                    const std::string &name,
                    const std::string &expected_status = "pass") {
  std::cout << std::setprecision(10) << "KIFMM_LITE_CASE name=" << name
            << " aggregate=" << result.aggregate_relative_rms
            << " k1=" << result.k1_relative_rms
            << " k2=" << result.k2_relative_rms
            << " k3=" << result.k3_relative_rms
            << " k4=" << result.k4_relative_rms
            << " far_pairs=" << result.accepted_far_cell_pairs
            << " fallback_pairs=" << result.direct_fallback_panel_pairs
            << " equivalent_system_builds="
            << result.equivalent_system_builds
            << " equivalent_rhs_solves=" << result.equivalent_rhs_solves
            << " status=" << result.status << '\n';
  require(result.finite, name + " finite");
  require(result.status == expected_status, name + " status");
  require(result.aggregate_relative_rms <= gate, name + " aggregate RMS gate");
  require(result.k1_relative_rms <= gate, name + " K1 RMS gate");
  require(result.k2_relative_rms <= gate, name + " K2 RMS gate");
  require(result.k3_relative_rms <= gate, name + " K3 RMS gate");
  require(result.k4_relative_rms <= gate, name + " K4 RMS gate");
}

} // namespace

int main() {
  const KifmmLiteConfig config = base_config();

  const auto near_sources = cluster_sources({0.0, 0.0, 0.0}, 0.03, 8, 0);
  const auto near_targets = cluster_targets({0.06, -0.03, 0.04}, 0.025, 8, 0);
  KifmmLiteConfig all_near_config = config;
  all_near_config.require_far_work = false;
  const KifmmLiteResult all_near =
      run_kifmm_lite_matvec(near_sources, near_targets, all_near_config);
  require_parity(all_near, 1.0e-12, "all-near",
                 "direct_only_reference");
  require(all_near.accepted_far_cell_pairs == 0, "all-near has no far pairs");
  require(all_near.direct_fallback_panel_pairs ==
              near_sources.size() * near_targets.size(),
          "all-near direct fallback coverage");
  std::vector<KifmmPanelSource> far_sources =
      cluster_sources({-8.0, 0.0, 0.0}, 0.05, 16, 100);
  const auto far_targets = cluster_targets({8.0, 0.2, -0.1}, 0.045, 14, 200);
  KifmmLiteConfig direct_only_config = config;
  direct_only_config.opening_theta = 0.0;
  const KifmmLiteResult direct_only =
      run_kifmm_lite_matvec(far_sources, far_targets, direct_only_config);
  require(!direct_only.finite, "direct-only default is not finite/pass");
  require(direct_only.status == "no_far_work",
          "direct-only default reports no_far_work");
  KifmmLiteConfig bad_eps_config = config;
  bad_eps_config.eps_ratio = -80.0;
  const KifmmLiteResult bad_eps =
      run_kifmm_lite_matvec(far_sources, far_targets, bad_eps_config);
  require(!bad_eps.finite && bad_eps.status == "invalid_config",
          "negative eps_ratio is rejected");

  const KifmmLiteResult separated =
      run_kifmm_lite_matvec(far_sources, far_targets, config);
  require_parity(separated, 1.0e-5, "separated");
  require(separated.accepted_far_cell_pairs > 0, "separated far pairs");
  require(separated.accepted_far_panel_pairs > 0, "separated far panel pairs");
  require(separated.equivalent_rhs_solves ==
              4 * separated.equivalent_system_builds,
          "separated equivalent RHS/build accounting");
  require(separated.equivalent_system_builds <
              2 * separated.accepted_far_cell_pairs,
          "separated source-node channel cache reused");

  auto mixed_sources = cluster_sources({-8.0, 0.0, 0.0}, 0.05, 12, 300);
  const auto near_tail = cluster_sources({3.0, 0.1, 0.0}, 0.04, 6, 400);
  mixed_sources.insert(mixed_sources.end(), near_tail.begin(), near_tail.end());
  const auto mixed_targets = cluster_targets({3.05, 0.0, 0.0}, 0.04, 12, 500);
  const KifmmLiteResult mixed =
      run_kifmm_lite_matvec(mixed_sources, mixed_targets, config);
  require_parity(mixed, 1.0e-6, "mixed");
  require(mixed.accepted_far_cell_pairs > 0, "mixed far pairs");
  require(mixed.accepted_far_panel_pairs > 0, "mixed far panel pairs");
  require(mixed.direct_fallback_panel_pairs > 0, "mixed fallback pairs");
  require(mixed.worklist_hash != 0, "mixed worklist hash is populated");
  const KifmmLiteResult mixed_repeat =
      run_kifmm_lite_matvec(mixed_sources, mixed_targets, config);
  require(mixed_repeat.worklist_hash == mixed.worklist_hash,
          "mixed worklist hash stable across repeated runs");
  require(mixed.equivalent_rhs_solves == 4 * mixed.equivalent_system_builds,
          "mixed equivalent RHS/build accounting");

  std::vector<KifmmPanelSource> invalid = near_sources;
  invalid.front().area = -1.0;
  const KifmmLiteResult bad =
      run_kifmm_lite_matvec(invalid, near_targets, config);
  require(!bad.finite, "invalid panel is rejected");
  require(bad.status != "pass", "invalid panel status is not pass");

  std::cout << std::setprecision(10)
            << "KIFMM_LITE all_near_rms=" << all_near.aggregate_relative_rms
            << " separated_rms=" << separated.aggregate_relative_rms
            << " mixed_rms=" << mixed.aggregate_relative_rms
            << " separated_far_pairs=" << separated.accepted_far_cell_pairs
            << " mixed_far_pairs=" << mixed.accepted_far_cell_pairs
            << " mixed_fallback_pairs=" << mixed.direct_fallback_panel_pairs
            << " mixed_hash=" << mixed.worklist_hash
            << " max_condition=" << mixed.max_condition_estimate << '\n';
  std::cout << "PASS titan_fmm_kifmm_lite_matvec_smoke\n";
  return EXIT_SUCCESS;
}
