#include "tabi/tabi_fmm_kifmm_hierarchy_jcp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmHierarchyJcpConfig;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::Vec3;

namespace {

double dot(Vec3 a, Vec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(Vec3 a) {
  return std::sqrt(dot(a, a));
}

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

KifmmPanelSource panel(Vec3 position, Vec3 normal, double area, double phi,
                       double phin) {
  return {position, unit(normal), area, phi, phin};
}

std::vector<KifmmPanelSource> compact_sources() {
  return {panel({-0.52, -0.44, -0.38}, {0.62, -0.37, 0.44}, 0.79, 0.31,
                -0.47),
          panel({-0.18, 0.12, 0.08}, {-0.35, 0.81, -0.27}, 0.96,
                -0.52, 0.58),
          panel({0.22, -0.21, 0.16}, {0.21, 0.44, 0.87}, 0.68, 0.49,
                -0.26),
          panel({0.58, 0.50, -0.29}, {-0.67, -0.28, 0.49}, 0.88,
                -0.34, 0.39)};
}

std::vector<Vec3> compact_targets() {
  return {{1.38, 0.40, -0.16}, {-1.31, -0.55, 0.24}, {0.16, 1.42, 0.46}};
}

std::vector<Vec3> compact_normals() {
  return {unit({0.96, 0.15, -0.12}), unit({-0.92, -0.08, 0.22}),
          unit({0.12, 0.98, 0.09})};
}

std::vector<KifmmPanelSource> separated_sources() {
  return {panel({-2.20, -0.34, -0.33}, {0.62, -0.37, 0.44}, 0.79,
                0.31, -0.47),
          panel({-2.05, -0.29, -0.31}, {-0.35, 0.81, -0.27}, 0.96,
                -0.52, 0.58),
          panel({-1.86, -0.36, -0.27}, {0.21, 0.44, 0.87}, 0.68,
                0.49, -0.26),
          panel({-1.70, -0.28, -0.37}, {-0.67, -0.28, 0.49}, 0.88,
                -0.34, 0.39),
          panel({-2.18, -0.24, -0.28}, {0.14, -0.93, 0.28}, 0.74,
                0.44, -0.22),
          panel({-1.78, -0.39, -0.35}, {-0.18, 0.21, -0.96}, 0.82,
                -0.28, 0.36)};
}

std::vector<Vec3> separated_targets() {
  return {{2.65, 0.28, 0.31},
          {2.78, 0.34, 0.27},
          {2.43, 0.30, 0.36},
          {2.94, 0.23, 0.25}};
}

std::vector<Vec3> separated_normals() {
  return {unit({0.96, 0.15, -0.12}), unit({-0.92, -0.08, 0.22}),
          unit({0.12, 0.98, 0.09}), unit({-0.14, -0.95, -0.18})};
}

std::vector<KifmmPanelSource> mixed_sources() {
  std::vector<KifmmPanelSource> sources = separated_sources();
  for (KifmmPanelSource &source : sources) {
    source.position.x -= 18.0;
  }
  const auto local = compact_sources();
  sources.insert(sources.end(), local.begin(), local.end());
  return sources;
}

std::vector<Vec3> mixed_targets() {
  return {{20.45, 0.28, 0.31}, {20.58, 0.34, 0.27}, {0.16, 1.42, 0.46},
          {-1.31, -0.55, 0.24}};
}

std::vector<Vec3> mixed_normals() {
  return {unit({0.96, 0.15, -0.12}), unit({-0.92, -0.08, 0.22}),
          unit({0.12, 0.98, 0.09}), unit({-0.14, -0.95, -0.18})};
}

std::vector<KifmmPanelSource> stress_sources() {
  std::vector<KifmmPanelSource> sources;
  sources.reserve(128);
  for (int ix = 0; ix < 8; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      for (int iz = 0; iz < 4; ++iz) {
        const double raw = ((ix + 2 * iy + 3 * iz) % 7 - 3) / 5.0;
        sources.push_back(panel(
            {-2.1 + 0.035 * ix, -0.38 + 0.055 * iy,
             -0.34 + 0.052 * iz},
            {0.3 + 0.07 * ix, -0.5 + 0.11 * iy, 0.4 + 0.09 * iz},
            0.55 + 0.01 * ((ix + iy + iz) % 5),
            raw == 0.0 ? 0.2 : raw,
            0.17 + 0.03 * ((2 * ix + iy + iz) % 9)));
      }
    }
  }
  return sources;
}

std::vector<Vec3> stress_targets() {
  std::vector<Vec3> targets;
  targets.reserve(128);
  for (int ix = 0; ix < 8; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      for (int iz = 0; iz < 4; ++iz) {
        targets.push_back({2.45 + 0.04 * ix, 0.22 + 0.057 * iy,
                           0.24 + 0.049 * iz});
      }
    }
  }
  return targets;
}

std::vector<Vec3> stress_normals() {
  std::vector<Vec3> normals;
  normals.reserve(128);
  for (int ix = 0; ix < 8; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      for (int iz = 0; iz < 4; ++iz) {
        normals.push_back(unit({0.8 + 0.02 * ix, 0.1 + 0.04 * iy,
                                -0.3 + 0.05 * iz}));
      }
    }
  }
  return normals;
}

int parse_int_arg(const std::string &arg, const char *prefix, int fallback) {
  const std::string p(prefix);
  if (arg.rfind(p, 0) == 0) {
    return std::stoi(arg.substr(p.size()));
  }
  return fallback;
}

double parse_double_arg(const std::string &arg, const char *prefix,
                        double fallback) {
  const std::string p(prefix);
  if (arg.rfind(p, 0) == 0) {
    return std::stod(arg.substr(p.size()));
  }
  return fallback;
}

} // namespace

int main(int argc, char **argv) {
  std::string case_id = "jcp_laplace_yukawa_separated_q6_check3_local3";
  KifmmHierarchyJcpConfig config;
  config.hierarchy.max_depth = 5;
  config.hierarchy.leaf_size = 2;
  config.hierarchy.order = 6;
  config.hierarchy.theta = 0.6;
  config.hierarchy.upward_check_scale = 3.0;
  config.hierarchy.local_equiv_scale = 3.0;
  config.hierarchy.local_check_scale = 1.05;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg.rfind("--case=", 0) == 0) {
      case_id = arg.substr(std::string("--case=").size());
    } else if (arg.rfind("--order=", 0) == 0) {
      config.hierarchy.order = parse_int_arg(arg, "--order=", 6);
    } else if (arg.rfind("--theta=", 0) == 0) {
      config.hierarchy.theta = parse_double_arg(arg, "--theta=", 0.6);
    } else if (arg.rfind("--upward-check-scale=", 0) == 0) {
      config.hierarchy.upward_check_scale =
          parse_double_arg(arg, "--upward-check-scale=", 3.0);
    } else if (arg.rfind("--equiv-scale=", 0) == 0) {
      config.hierarchy.equiv_scale =
          parse_double_arg(arg, "--equiv-scale=", 1.05);
    } else if (arg.rfind("--local-equiv-scale=", 0) == 0) {
      config.hierarchy.local_equiv_scale =
          parse_double_arg(arg, "--local-equiv-scale=", 3.0);
    } else if (arg.rfind("--local-check-scale=", 0) == 0) {
      config.hierarchy.local_check_scale =
          parse_double_arg(arg, "--local-check-scale=", 1.05);
    } else if (arg.rfind("--kappa=", 0) == 0) {
      config.kappa = parse_double_arg(arg, "--kappa=", 0.12746);
    }
  }

  std::vector<KifmmPanelSource> sources = separated_sources();
  std::vector<Vec3> targets = separated_targets();
  std::vector<Vec3> normals = separated_normals();
  if (case_id.find("all_near") != std::string::npos) {
    sources = compact_sources();
    targets = compact_targets();
    normals = compact_normals();
    config.hierarchy.theta = 1.0e-12;
  } else if (case_id.find("mixed") != std::string::npos) {
    sources = mixed_sources();
    targets = mixed_targets();
    normals = mixed_normals();
  } else if (case_id.find("stress") != std::string::npos) {
    sources = stress_sources();
    targets = stress_targets();
    normals = stress_normals();
    config.hierarchy.max_depth = 6;
    config.hierarchy.leaf_size = 8;
  }
  if (case_id.find("_q8") != std::string::npos) {
    config.hierarchy.order = 8;
  }
  if (case_id.find("check2p4") != std::string::npos) {
    config.hierarchy.upward_check_scale = 2.4;
    config.hierarchy.local_equiv_scale = 2.4;
  }

  const auto result = titan::tabi::fmm_kifmm::run_full_hierarchical_jcp_kifmm(
      sources, targets, normals, config);
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << "E253_RESULT"
            << " case=" << case_id << " status=" << result.status
            << " order=" << config.hierarchy.order
            << " theta=" << config.hierarchy.theta
            << " upward_check_scale=" << config.hierarchy.upward_check_scale
            << " local_equiv_scale=" << config.hierarchy.local_equiv_scale
            << " local_check_scale=" << config.hierarchy.local_check_scale
            << " condition_policy=" << config.hierarchy.condition_policy
            << " source_box_policy=" << config.hierarchy.source_box_policy
            << " aggregate_rms=" << result.aggregate_relative_rms
            << " k1_rms=" << result.k1_relative_rms
            << " k2_rms=" << result.k2_relative_rms
            << " k3_rms=" << result.k3_relative_rms
            << " k4_rms=" << result.k4_relative_rms
            << " max_channel_rms=" << result.max_channel_relative_rms
            << " max_channel_name=" << result.max_channel_name
            << " far_fraction=" << result.stats.far_fraction
            << " far_pairs=" << result.stats.far_panel_pairs
            << " direct_pairs=" << result.stats.direct_panel_pairs
            << " coverage=" << result.stats.exact_pair_coverage
            << " source_count=" << result.stats.source_count
            << " target_count=" << result.stats.target_count
            << " node_count=" << result.stats.node_count
            << " leaf_count=" << result.stats.leaf_count
            << " unique_half_widths=" << result.stats.unique_half_widths
            << " observed_level_count=" << result.stats.observed_level_count
            << " max_depth=" << result.stats.max_observed_depth
            << " far_pair_count=" << result.stats.far_pair_count
            << " direct_pair_count=" << result.stats.direct_pair_count
            << " p2m_leaf_solves=" << result.p2m_leaf_solves
            << " m2m_translations=" << result.m2m_translations
            << " m2m_parent_solves=" << result.m2m_parent_solves
            << " m2l_translations=" << result.m2l_translations
            << " local_node_solves=" << result.local_node_solves
            << " l2l_translations=" << result.l2l_translations
            << " l2p_target_count=" << result.l2p_target_count
            << " local_rhs_adds=" << result.local_rhs_adds
            << " upward_rhs_adds=" << result.upward_rhs_adds
            << " cache_entries=" << result.equivalent_factor_cache_entries
            << " cache_hits=" << result.equivalent_factor_cache_hits
            << " cache_misses=" << result.equivalent_factor_cache_misses
            << " p2m_seconds=" << result.p2m_seconds
            << " m2m_seconds=" << result.m2m_seconds
            << " m2l_seconds=" << result.m2l_seconds
            << " local_solve_seconds=" << result.local_solve_seconds
            << " l2l_seconds=" << result.l2l_seconds
            << " l2p_seconds=" << result.l2p_seconds
            << " direct_seconds=" << result.direct_fallback_seconds
            << " oracle_seconds=" << result.dense_oracle_seconds
            << " total_seconds=" << result.total_seconds
            << " tree_hash=" << result.tree_hash
            << " worklist_hash=" << result.worklist_hash
            << " policy_hash=" << result.policy_hash
            << " output_hash=" << result.output_hash << '\n';

  return result.status == "ok" ? 0 : 1;
}
