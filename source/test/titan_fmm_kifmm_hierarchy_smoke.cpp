#include "tabi/tabi_fmm_kifmm_hierarchy.h"
#include "tabi/tabi_fmm_kifmm_hierarchy_scalar.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmHierarchyPolicy;
using titan::tabi::fmm_kifmm::KifmmHierarchyScalarConfig;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmScalarSource;
using titan::tabi::fmm_kifmm::Vec3;

namespace {

void require(bool ok, const char *message) {
  if (!ok) {
    throw std::runtime_error(message);
  }
}

std::vector<KifmmScalarSource> separated_sources() {
  return {{{-2.20, -0.34, -0.33}, 1.0},
          {{-2.05, -0.29, -0.31}, -0.7},
          {{-1.86, -0.36, -0.27}, 0.5},
          {{-1.70, -0.28, -0.37}, -0.3},
          {{-2.18, -0.24, -0.28}, 0.4},
          {{-1.78, -0.39, -0.35}, -0.2}};
}

std::vector<Vec3> separated_targets() {
  return {{2.65, 0.28, 0.31},
          {2.78, 0.34, 0.27},
          {2.43, 0.30, 0.36},
          {2.94, 0.23, 0.25}};
}

std::vector<KifmmScalarSource> compact_sources() {
  return {{{-0.85, -0.75, -0.65}, 1.0},
          {{-0.40, 0.20, 0.10}, -0.7},
          {{0.30, -0.35, 0.25}, 0.5},
          {{0.95, 0.85, -0.45}, -0.3}};
}

std::vector<Vec3> compact_targets() {
  return {{1.8, 0.7, -0.2}, {-1.7, -0.9, 0.4}, {0.2, 1.9, 0.8}};
}

std::vector<KifmmScalarSource> stress_sources() {
  std::vector<KifmmScalarSource> sources;
  sources.reserve(128);
  for (int ix = 0; ix < 8; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      for (int iz = 0; iz < 4; ++iz) {
        const double strength =
            ((ix + 2 * iy + 3 * iz) % 7 - 3) / 5.0;
        sources.push_back(
            {{-2.1 + 0.035 * ix, -0.38 + 0.055 * iy,
              -0.34 + 0.052 * iz},
             strength == 0.0 ? 0.2 : strength});
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

KifmmHierarchyScalarConfig scalar_config(int order, double theta,
                                         KifmmKernel kernel, double kappa) {
  KifmmHierarchyScalarConfig config;
  config.hierarchy.max_depth = 5;
  config.hierarchy.leaf_size = 2;
  config.hierarchy.order = order;
  config.hierarchy.theta = theta;
  config.kernel = kernel;
  config.kappa = kappa;
  return config;
}

double parse_double_arg(const std::string &arg, const char *prefix,
                        double fallback) {
  const std::string p(prefix);
  if (arg.rfind(p, 0) == 0) {
    return std::stod(arg.substr(p.size()));
  }
  return fallback;
}

int parse_int_arg(const std::string &arg, const char *prefix, int fallback) {
  const std::string p(prefix);
  if (arg.rfind(p, 0) == 0) {
    return std::stoi(arg.substr(p.size()));
  }
  return fallback;
}

const char *kernel_name(KifmmKernel kernel) {
  return kernel == KifmmKernel::Laplace ? "laplace" : "yukawa";
}

void print_result(const std::string &case_id,
                  const KifmmHierarchyScalarConfig &config,
                  const titan::tabi::fmm_kifmm::KifmmHierarchyScalarResult
                      &result) {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << "E250_RESULT"
            << " case=" << case_id << " status=" << result.status
            << " order=" << config.hierarchy.order
            << " theta=" << config.hierarchy.theta
            << " kernel=" << kernel_name(config.kernel)
            << " kappa=" << config.kappa
            << " upward_check_scale="
            << config.hierarchy.upward_check_scale
            << " local_equiv_scale=" << config.hierarchy.local_equiv_scale
            << " local_check_scale=" << config.hierarchy.local_check_scale
            << " condition_policy=" << config.hierarchy.condition_policy
            << " source_box_policy=" << config.hierarchy.source_box_policy
            << " aggregate_rms=" << result.aggregate_relative_rms
            << " max_abs_error=" << result.max_abs_error
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
            << " mac_accepted_pairs=" << result.stats.mac_accepted_pairs
            << " check_guard_rejected_pairs="
            << result.stats.check_guard_rejected_pairs
            << " leaf_direct_pairs=" << result.stats.leaf_direct_pairs
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
            << " p2m_rhs_seconds=" << result.p2m_rhs_seconds
            << " p2m_solve_seconds=" << result.p2m_solve_seconds
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
}

int run_case(const std::string &case_id) {
  int order = 6;
  double theta = 0.6;
  KifmmKernel kernel = KifmmKernel::Laplace;
  double kappa = 0.0;
  std::vector<KifmmScalarSource> sources = separated_sources();
  std::vector<Vec3> targets = separated_targets();

  if (case_id == "scalar_laplace_all_near_q6") {
    theta = 1.0e-12;
    sources = compact_sources();
    targets = compact_targets();
  } else if (case_id == "scalar_laplace_separated_q6") {
    order = 6;
  } else if (case_id == "scalar_laplace_separated_q8") {
    order = 8;
  } else if (case_id == "scalar_yukawa_mixed_q8_kappa0p12746") {
    order = 8;
    kernel = KifmmKernel::Yukawa;
    kappa = 0.12746;
  } else {
    std::cerr << "unknown E25.0 case: " << case_id << '\n';
    return 2;
  }

  const auto result = titan::tabi::fmm_kifmm::run_hierarchical_scalar_kifmm(
      sources, targets, scalar_config(order, theta, kernel, kappa));
  print_result(case_id, scalar_config(order, theta, kernel, kappa), result);
  return result.status == "ok" ? 0 : 1;
}

int run_case_with_args(int argc, char **argv) {
  std::string case_id = "scalar_laplace_separated_q6";
  int order = 6;
  double theta = 0.6;
  double upward_check_scale = 2.40;
  double equiv_scale = 1.05;
  double local_equiv_scale = 2.40;
  double local_check_scale = 1.05;
  KifmmKernel kernel = KifmmKernel::Laplace;
  double kappa = 0.0;
  bool full_hierarchy = false;
  bool order_provided = false;
  bool theta_provided = false;
  bool kernel_provided = false;
  bool kappa_provided = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.rfind("--case=", 0) == 0) {
      case_id = arg.substr(std::string("--case=").size());
    } else if (arg.rfind("--order=", 0) == 0) {
      order = parse_int_arg(arg, "--order=", order);
      order_provided = true;
    } else if (arg.rfind("--theta=", 0) == 0) {
      theta = parse_double_arg(arg, "--theta=", theta);
      theta_provided = true;
    } else if (arg.rfind("--upward-check-scale=", 0) == 0) {
      upward_check_scale =
          parse_double_arg(arg, "--upward-check-scale=", upward_check_scale);
    } else if (arg.rfind("--equiv-scale=", 0) == 0) {
      equiv_scale = parse_double_arg(arg, "--equiv-scale=", equiv_scale);
    } else if (arg.rfind("--local-equiv-scale=", 0) == 0) {
      local_equiv_scale =
          parse_double_arg(arg, "--local-equiv-scale=", local_equiv_scale);
    } else if (arg.rfind("--local-check-scale=", 0) == 0) {
      local_check_scale =
          parse_double_arg(arg, "--local-check-scale=", local_check_scale);
    } else if (arg.rfind("--kernel=", 0) == 0) {
      const std::string value = arg.substr(std::string("--kernel=").size());
      kernel = value == "yukawa" ? KifmmKernel::Yukawa : KifmmKernel::Laplace;
      kernel_provided = true;
    } else if (arg.rfind("--kappa=", 0) == 0) {
      kappa = parse_double_arg(arg, "--kappa=", kappa);
      kappa_provided = true;
    } else if (arg == "--full-hierarchy") {
      full_hierarchy = true;
    }
  }

  if (!order_provided) {
    if (case_id.find("_q8") != std::string::npos) {
      order = 8;
    } else if (case_id.find("_q6") != std::string::npos) {
      order = 6;
    }
  }
  if (!kernel_provided && case_id.find("yukawa") != std::string::npos) {
    kernel = KifmmKernel::Yukawa;
  }
  if (!kappa_provided && case_id.find("kappa0p12746") != std::string::npos) {
    kappa = 0.12746;
  }

  std::vector<KifmmScalarSource> sources = separated_sources();
  std::vector<Vec3> targets = separated_targets();
  bool stress_case = false;
  if (case_id.find("all_near") != std::string::npos) {
    sources = compact_sources();
    targets = compact_targets();
    if (!theta_provided) {
      theta = 1.0e-12;
    }
  } else if (case_id.find("stress") != std::string::npos) {
    sources = stress_sources();
    targets = stress_targets();
    stress_case = true;
  }

  KifmmHierarchyScalarConfig config =
      scalar_config(order, theta, kernel, kappa);
  if (stress_case) {
    config.hierarchy.max_depth = 6;
    config.hierarchy.leaf_size = 8;
  }
  config.hierarchy.equiv_scale = equiv_scale;
  config.hierarchy.upward_check_scale = upward_check_scale;
  config.hierarchy.local_equiv_scale = local_equiv_scale;
  config.hierarchy.local_check_scale = local_check_scale;
  const auto result =
      full_hierarchy
          ? titan::tabi::fmm_kifmm::run_full_hierarchical_scalar_kifmm(
                sources, targets, config)
          : titan::tabi::fmm_kifmm::run_hierarchical_scalar_kifmm(
                sources, targets, config);
  print_result(case_id, config, result);
  return result.status == "ok" ? 0 : 1;
}

} // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.rfind("--case=", 0) == 0 || arg == "--full-hierarchy") {
      return run_case_with_args(argc, argv);
    }
  }

  KifmmHierarchyPolicy policy;
  policy.max_depth = 4;
  policy.leaf_size = 4;

  std::vector<Vec3> sources{{-0.9, -0.8, -0.7},
                            {-0.5, -0.2, 0.1},
                            {0.2, 0.3, -0.1},
                            {0.7, 0.8, 0.6},
                            {1.1, -0.6, 0.4},
                            {-1.2, 0.5, -0.3}};
  std::vector<Vec3> targets{
      {-0.7, 0.9, 0.2}, {0.6, -0.9, -0.4}, {1.3, 0.2, 0.7}};

  const auto tree =
      titan::tabi::fmm_kifmm::build_strict_octree(sources, targets, policy);

  require(tree.status == "ok", "strict octree did not build");
  require(tree.stats.unique_half_widths <= policy.max_depth + 1,
          "strict octree produced too many half-widths");
  require(tree.stats.max_observed_depth <= policy.max_depth,
          "strict octree exceeded max depth");

  for (const auto &node : tree.nodes) {
    for (const int child_id : node.child) {
      if (child_id < 0) {
        continue;
      }
      const auto &child = tree.nodes[static_cast<std::size_t>(child_id)];
      require(child.parent == node.id, "child parent link mismatch");
      require(std::abs(child.half_width - 0.5 * node.half_width) <=
                  1.0e-12 * std::max(1.0, node.half_width),
              "child half-width is not exactly parent/2");
    }
  }

  const auto repeat =
      titan::tabi::fmm_kifmm::build_strict_octree(sources, targets, policy);
  require(repeat.tree_hash == tree.tree_hash, "strict octree hash unstable");

  const auto worklist =
      titan::tabi::fmm_kifmm::build_hierarchy_worklist(tree, policy);
  require(worklist.worklist_hash != 0, "worklist hash missing");
  require(worklist.m2l_pairs.size() + worklist.direct_leaf_pairs.size() > 0,
          "worklist empty");
  require(worklist.policy_hash != 0, "policy hash missing");
  require(tree.stats.exact_pair_coverage ==
              static_cast<std::uint64_t>(sources.size() * targets.size()),
          "tree pair coverage is not exact");

  KifmmHierarchyScalarConfig scalar_config;
  scalar_config.hierarchy.max_depth = 4;
  scalar_config.hierarchy.leaf_size = 2;
  scalar_config.hierarchy.order = 6;
  scalar_config.hierarchy.theta = 1.0e-12;
  scalar_config.kernel = KifmmKernel::Laplace;

  const auto scalar = titan::tabi::fmm_kifmm::run_hierarchical_scalar_kifmm(
      compact_sources(), compact_targets(), scalar_config);
  require(scalar.status == "ok", "scalar hierarchy failed");
  require(scalar.aggregate_relative_rms <= 1.0e-6,
          "scalar hierarchy parity gap");
  require(scalar.output_hash != 0, "scalar output hash missing");

  KifmmHierarchyScalarConfig separated_config = scalar_config;
  separated_config.hierarchy.max_depth = 5;
  separated_config.hierarchy.leaf_size = 2;
  separated_config.hierarchy.theta = 0.6;
  separated_config.hierarchy.order = 8;
  const auto separated = titan::tabi::fmm_kifmm::run_hierarchical_scalar_kifmm(
      separated_sources(), separated_targets(), separated_config);
  require(separated.status == "ok", "separated scalar hierarchy failed");
  require(separated.stats.far_fraction >= 0.10,
          "hierarchical separated case did not use far work");
  require(separated.aggregate_relative_rms <= 1.0e-6,
          "hierarchical separated scalar parity gap");

  KifmmHierarchyScalarConfig yukawa_config = separated_config;
  yukawa_config.kernel = KifmmKernel::Yukawa;
  yukawa_config.kappa = 0.12746;
  const auto yukawa = titan::tabi::fmm_kifmm::run_hierarchical_scalar_kifmm(
      separated_sources(), separated_targets(), yukawa_config);
  require(yukawa.status == "ok", "Yukawa scalar hierarchy failed");
  require(yukawa.aggregate_relative_rms <= 1.0e-6,
          "Yukawa scalar hierarchy parity gap");

  std::cout << "FMM-E25.0 hierarchy smoke passed\n";
  return 0;
}
