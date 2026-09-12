#include "tabi/tabi_fmm_kifmm_hierarchy_real_panel.h"
#include "tabi/tabi_fmm_kifmm_sampled_panel_matvec.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Options;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::Result;
using titan::tabi::fmm_kifmm::KifmmHierarchyJcpConfig;
using titan::tabi::fmm_kifmm::KifmmHierarchyRealPanelResult;
using titan::tabi::fmm_kifmm::KifmmRealPanelCloudConfig;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::run_kifmm_hierarchy_real_panel_full_apply_sampled_audit;
using titan::tabi::fmm_kifmm::run_kifmm_hierarchy_real_panel_sampled;
using titan::tabi::fmm_kifmm::run_kifmm_lite_sampled_panel_matvec;

// E27/E28/E29/E30 smoke markers: hybrid-local-bridge,
// hybrid-source-local-bridge, hybrid-local-split-bridge,
// channel-hybrid-k123-hierarchy-k4-cartesian,
// channel-hybrid-shared-worklist-k123-hierarchy-k4-cartesian,
// dual-worklist-channel-hybrid-k123-kifmm-k4-cartesian,
// cartesian_p5_geometry0p4, full-panel-profile, build_direct_reference.
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

Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

double norm(Vec3 a) { return std::sqrt(std::max(0.0, dot(a, a))); }

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

Vec3 position(const Panel &p) { return {p.x, p.y, p.z}; }

Vec3 normal(const Panel &p) { return {p.nx, p.ny, p.nz}; }

double edge_length(Vec3 a, Vec3 b) { return norm(sub(a, b)); }

double aspect_score(const Panel &p) {
  const Vec3 a{p.v0x, p.v0y, p.v0z};
  const Vec3 b{p.v1x, p.v1y, p.v1z};
  const Vec3 c{p.v2x, p.v2y, p.v2z};
  const double e0 = edge_length(a, b);
  const double e1 = edge_length(b, c);
  const double e2 = edge_length(c, a);
  const double lo = std::max(1.0e-300, std::min({e0, e1, e2}));
  return std::max({e0, e1, e2}) / lo;
}

std::uint64_t mix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}

std::vector<Atom> read_pqr(const std::string &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to open PQR fixture: " + path);
  }
  std::vector<Atom> atoms;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) {
      continue;
    }
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) {
      fields.push_back(field);
    }
    if (fields.size() < 10) {
      throw std::runtime_error("malformed PQR row at line " +
                               std::to_string(line_number));
    }
    const std::size_t n = fields.size();
    Atom atom;
    atom.x = std::stod(fields[n - 5]);
    atom.y = std::stod(fields[n - 4]);
    atom.z = std::stod(fields[n - 3]);
    atom.q = std::stod(fields[n - 2]);
    atom.radius = std::stod(fields[n - 1]);
    if (!(atom.radius > 0.0)) {
      atom.radius = 0.6;
    }
    atoms.push_back(atom);
  }
  if (atoms.empty()) {
    throw std::runtime_error("PQR fixture contains no atoms: " + path);
  }
  return atoms;
}

std::vector<Panel> build_real_surface(const std::string &pqr_path,
                                      int surface_points) {
  Options options;
  options.surface_mode = "vdw";
  options.surface_geometry = "triangle";
  options.diagonal_mode = "solid-angle";
  options.surface_point_count = surface_points;
  Result surface_result;
  TitanExecutor executor(1);
  std::vector<Panel> panels =
      TitanPBE::Tabi::build_surface(read_pqr(pqr_path), options,
                                    &surface_result, &executor);
  if (panels.size() < 256) {
    throw std::runtime_error("surface builder produced too few panels");
  }
  return panels;
}

Vec3 centroid(const std::vector<Panel> &panels) {
  Vec3 c{};
  for (const Panel &p : panels) {
    c.x += p.x;
    c.y += p.y;
    c.z += p.z;
  }
  const double inv = 1.0 / static_cast<double>(panels.size());
  return {c.x * inv, c.y * inv, c.z * inv};
}

std::vector<double> unknowns_for(const std::vector<Panel> &panels) {
  std::vector<double> x(2 * panels.size(), 0.0);
  const Vec3 c = centroid(panels);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const Vec3 r = sub(position(panels[i]), c);
    const Vec3 n = unit(normal(panels[i]));
    x[i] = 0.12 * std::sin(0.17 * static_cast<double>(i)) + 0.004 * r.x -
           0.003 * r.y + 0.015 * dot(unit(r), n);
    x[panels.size() + i] =
        -0.10 * std::cos(0.11 * static_cast<double>(i)) + 0.003 * r.z -
        0.008 * dot(unit(r), n);
  }
  return x;
}

template <typename KeyFn>
std::vector<std::size_t> select_by_key(const std::vector<Panel> &panels,
                                       std::size_t count, bool largest,
                                       KeyFn key_fn) {
  std::vector<std::size_t> idx(panels.size());
  std::iota(idx.begin(), idx.end(), std::size_t{0});
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
    const double ka = key_fn(panels[a], a);
    const double kb = key_fn(panels[b], b);
    if (ka == kb) {
      return a < b;
    }
    return largest ? ka > kb : ka < kb;
  });
  if (idx.size() > count) {
    idx.resize(count);
  }
  std::sort(idx.begin(), idx.end());
  return idx;
}

std::vector<std::size_t> uniform_hash_sample(const std::vector<Panel> &panels,
                                             std::size_t count) {
  std::vector<std::size_t> idx(panels.size());
  std::iota(idx.begin(), idx.end(), std::size_t{0});
  std::sort(idx.begin(), idx.end(), [](std::size_t a, std::size_t b) {
    const std::uint64_t ha = mix64(static_cast<std::uint64_t>(a));
    const std::uint64_t hb = mix64(static_cast<std::uint64_t>(b));
    if (ha == hb) {
      return a < b;
    }
    return ha < hb;
  });
  if (idx.size() > count) {
    idx.resize(count);
  }
  std::sort(idx.begin(), idx.end());
  return idx;
}

std::vector<std::size_t> high_aspect_sample(const std::vector<Panel> &panels,
                                            std::size_t count) {
  return select_by_key(panels, count, true, [](const Panel &p, std::size_t) {
    return aspect_score(p);
  });
}

std::vector<std::size_t>
inward_normal_groove_sample(const std::vector<Panel> &panels,
                            std::size_t count) {
  const Vec3 c = centroid(panels);
  return select_by_key(panels, count, false, [&](const Panel &p, std::size_t) {
    return dot(unit(sub(position(p), c)), unit(normal(p)));
  });
}

std::vector<std::size_t>
closest_opposing_surface_sample(const std::vector<Panel> &panels,
                                std::size_t count) {
  const std::size_t shortlist_count =
      std::min<std::size_t>(panels.size(),
                            std::max<std::size_t>(1024, 8 * count));
  std::vector<std::size_t> targets =
      inward_normal_groove_sample(panels, shortlist_count / 2);
  std::vector<std::size_t> aspect =
      high_aspect_sample(panels, shortlist_count / 4);
  std::vector<std::size_t> hashed =
      uniform_hash_sample(panels, shortlist_count / 4);
  targets.insert(targets.end(), aspect.begin(), aspect.end());
  targets.insert(targets.end(), hashed.begin(), hashed.end());
  std::sort(targets.begin(), targets.end());
  targets.erase(std::unique(targets.begin(), targets.end()), targets.end());

  std::vector<std::size_t> search = uniform_hash_sample(panels, shortlist_count);
  search.insert(search.end(), targets.begin(), targets.end());
  std::sort(search.begin(), search.end());
  search.erase(std::unique(search.begin(), search.end()), search.end());

  std::vector<std::pair<double, std::size_t>> scored;
  scored.reserve(targets.size());
  for (const std::size_t ti : targets) {
    const Vec3 tp = position(panels[ti]);
    const Vec3 tn = unit(normal(panels[ti]));
    double best = std::numeric_limits<double>::infinity();
    for (const std::size_t sj : search) {
      if (sj == ti) {
        continue;
      }
      const Vec3 delta = sub(position(panels[sj]), tp);
      const double d = norm(delta);
      if (!(d > 0.0)) {
        continue;
      }
      const double opposing = -dot(tn, unit(normal(panels[sj])));
      const double facing = dot(tn, unit(delta));
      const double score =
          d / std::max(0.05, 0.5 + 0.5 * opposing + 0.25 * facing);
      best = std::min(best, score);
    }
    if (std::isfinite(best)) {
      scored.push_back({best, ti});
    }
  }
  std::sort(scored.begin(), scored.end(), [](const auto &a, const auto &b) {
    if (a.first == b.first) {
      return a.second < b.second;
    }
    return a.first < b.first;
  });
  std::vector<std::size_t> out;
  for (const auto &entry : scored) {
    if (out.size() == count) {
      break;
    }
    out.push_back(entry.second);
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<std::size_t> sample_for_name(const std::vector<Panel> &panels,
                                         const std::string &sample_name,
                                         std::size_t sample_count) {
  if (sample_name == "uniform_hash") {
    return uniform_hash_sample(panels, sample_count);
  }
  if (sample_name == "high_aspect") {
    return high_aspect_sample(panels, sample_count);
  }
  if (sample_name == "closest_opposing_surface") {
    return closest_opposing_surface_sample(panels, sample_count);
  }
  if (sample_name == "inward_normal_groove") {
    return inward_normal_groove_sample(panels, sample_count);
  }
  throw std::runtime_error("unknown E25.5 sample set: " + sample_name);
}

KifmmHierarchyJcpConfig hierarchy_config(int order) {
  KifmmHierarchyJcpConfig config;
  config.hierarchy.order = order;
  config.hierarchy.leaf_size = 32;
  config.hierarchy.max_depth = 6;
  config.hierarchy.theta = 0.6;
  config.hierarchy.equiv_scale = 1.05;
  config.hierarchy.upward_check_scale = 3.0;
  config.hierarchy.local_equiv_scale = 3.0;
  config.hierarchy.local_check_scale = 1.05;
  config.hierarchy.condition_policy = "qr-rank-only";
  config.hierarchy.source_box_policy = "strict-octree";
  config.eps_ratio = 80.0;
  config.kappa = 0.12746;
  config.aggregate_relative_rms_gate = 1.0e-6;
  config.max_channel_relative_rms_gate = 5.0e-6;
  return config;
}

void apply_hierarchy_overrides(
    KifmmHierarchyJcpConfig &config,
    const std::unordered_map<std::string, std::string> &args) {
  const auto get = [&](const std::string &name) {
    const auto it = args.find(name);
    return it == args.end() ? std::string{} : it->second;
  };
  if (!get("theta").empty()) {
    config.hierarchy.theta = std::stod(get("theta"));
  }
  if (!get("equiv-scale").empty()) {
    config.hierarchy.equiv_scale = std::stod(get("equiv-scale"));
  }
  if (!get("upward-check-scale").empty()) {
    config.hierarchy.upward_check_scale = std::stod(get("upward-check-scale"));
  }
  if (!get("local-equiv-scale").empty()) {
    config.hierarchy.local_equiv_scale = std::stod(get("local-equiv-scale"));
  }
  if (!get("local-check-scale").empty()) {
    config.hierarchy.local_check_scale = std::stod(get("local-check-scale"));
  }
  if (!get("condition-policy").empty()) {
    config.hierarchy.condition_policy = get("condition-policy");
  }
  if (!get("m2l-eval-policy").empty()) {
    config.hierarchy.m2l_eval_policy = get("m2l-eval-policy");
  }
  if (!get("m2l-matrix-cache-policy").empty()) {
    // Accepted values include static-translation-table-cache for E51
    // translation-invariance diagnostics and target-node-streaming-cache for
    // E49 retained memory reduction diagnostics.
    config.hierarchy.m2l_matrix_cache_policy =
        get("m2l-matrix-cache-policy");
  }
  if (!get("m2l-matrix-cache-limit-mb").empty()) {
    config.hierarchy.m2l_matrix_cache_limit_mb =
        std::stod(get("m2l-matrix-cache-limit-mb"));
  }
  if (!get("m2l-gemm-batch-policy").empty()) {
    config.hierarchy.m2l_gemm_batch_policy =
        get("m2l-gemm-batch-policy");
  }
  if (!get("m2l-parallel-policy").empty()) {
    config.hierarchy.m2l_parallel_policy = get("m2l-parallel-policy");
  }
  if (!get("m2l-parallel-workers").empty()) {
    config.hierarchy.m2l_parallel_workers =
        std::stoi(get("m2l-parallel-workers"));
  }
  if (!get("jcp-channel-policy").empty()) {
    config.hierarchy.jcp_channel_policy = get("jcp-channel-policy");
  }
  if (!get("tree-stop-policy").empty()) {
    config.hierarchy.tree_stop_policy = get("tree-stop-policy");
  }
  if (!get("tree-panel-radius-factor").empty()) {
    config.hierarchy.tree_panel_radius_factor =
        std::stod(get("tree-panel-radius-factor"));
  }
  if (!get("tree-occupancy-stop").empty()) {
    config.hierarchy.tree_occupancy_stop =
        std::stoi(get("tree-occupancy-stop"));
  }
  if (!get("real-panel-guard").empty()) {
    // E25.7 passes "strict-local-check" through to the real-panel hierarchy.
    config.real_panel_guard = get("real-panel-guard");
  }
  if (!get("guard-margin-factor").empty()) {
    config.real_panel_guard_margin_factor =
        std::stod(get("guard-margin-factor"));
  }
  if (!get("hier-debug-mode").empty()) {
    // E25.8 passes modes such as "m2l-direct-to-targets" through for
    // hierarchy propagation attribution.
    config.hier_debug_mode = get("hier-debug-mode");
  }
}

KifmmRealPanelCloudConfig lite_config() {
  KifmmRealPanelCloudConfig config;
  config.lite.order = 6;
  config.lite.leaf_size = 32;
  config.lite.opening_theta = 0.6;
  config.lite.equiv_scale = 1.05;
  config.lite.check_scale = 2.40;
  config.lite.kappa = 0.12746;
  config.lite.eps_ratio = 80.0;
  config.lite.max_condition_estimate = 1.0e12;
  config.lite.max_equivalent_relative_residual = 1.0e-10;
  config.lite.max_aggregate_relative_rms = 1.0e-6;
  config.lite.max_channel_relative_rms = 5.0e-6;
  config.lite.min_far_panel_fraction = 0.10;
  config.lite.require_far_work = true;
  config.lite.require_unit_normals = true;
  config.lite.solver = "qr";
  config.lite.condition_policy = "qr-rank-only";
  config.lite.source_box_policy = "level-cube";
  config.near_integration = "analytic_singular";
  config.quad_buffer = "off";
  return config;
}

std::unordered_map<std::string, std::string> parse_args(int argc, char **argv) {
  std::unordered_map<std::string, std::string> args;
  for (int i = 1; i < argc; ++i) {
    std::string token(argv[i]);
    if (token.rfind("--", 0) != 0) {
      continue;
    }
    token = token.substr(2);
    const std::size_t eq = token.find('=');
    if (eq == std::string::npos) {
      args[token] = "1";
    } else {
      args[token.substr(0, eq)] = token.substr(eq + 1);
    }
  }
  return args;
}

std::string get_arg(const std::unordered_map<std::string, std::string> &args,
                    const std::string &name, const std::string &fallback) {
  const auto it = args.find(name);
  return it == args.end() ? fallback : it->second;
}

void print_invalid(const std::string &name,
                   const KifmmHierarchyRealPanelResult &result) {
  std::cout << "E255_INVALID name=" << name << " status=" << result.status
            << " execution_status=" << result.execution_status << '\n';
}

double hierarchy_stage_timing_sum(const KifmmHierarchyRealPanelResult &result) {
  const auto &h = result.hierarchy;
  return h.p2m_seconds + h.m2m_seconds + h.m2l_seconds +
         h.local_solve_seconds + h.l2l_seconds + h.l2p_seconds +
         h.direct_fallback_seconds;
}

std::pair<std::string, double>
hierarchy_bottleneck(const KifmmHierarchyRealPanelResult &result) {
  const auto &h = result.hierarchy;
  const std::pair<std::string, double> stages[] = {
      {"p2m", h.p2m_seconds},
      {"m2m", h.m2m_seconds},
      {"m2l", h.m2l_seconds},
      {"local_solve", h.local_solve_seconds},
      {"l2l", h.l2l_seconds},
      {"l2p", h.l2p_seconds},
      {"hierarchy_direct_fallback", h.direct_fallback_seconds}};
  std::pair<std::string, double> best{"none", 0.0};
  for (const auto &stage : stages) {
    if (std::isfinite(stage.second) && stage.second > best.second) {
      best = stage;
    }
  }
  const double denominator = hierarchy_stage_timing_sum(result);
  const double fraction = denominator > 0.0 ? best.second / denominator : 0.0;
  return {best.first, fraction};
}

void print_profile(const std::string &fixture, const std::string &sample_name,
                   const KifmmHierarchyRealPanelResult &result,
                   bool emit_e35_attribution,
                   bool emit_e36_stage_attribution) {
  const auto &h = result.hierarchy;
  const double stage_sum = hierarchy_stage_timing_sum(result);
  const double unattributed =
      result.hierarchy_seconds > stage_sum ? result.hierarchy_seconds - stage_sum
                                           : 0.0;
  const auto bottleneck = hierarchy_bottleneck(result);
  std::cout << std::setprecision(12)
            << "E255_RESULT fixture=" << fixture << " sample=" << sample_name
            << " status=" << result.status
            << " execution_status=" << result.execution_status
            << " accuracy_status=" << result.accuracy_status
            << " accuracy_failure=" << result.accuracy_failure
            << " source_panels=" << result.source_panels
            << " target_panels=" << result.target_panels
            << " sampled_rows=" << result.sampled_rows
            << " timing_scope=" << result.timing_scope
            << " compression_status=" << result.compression_status
            << " full_apply_stage_status=" << result.full_apply_stage_status
            << " kifmm_compression_ratio="
            << result.kifmm_compression_ratio
            << " cartesian_k4_compression_ratio="
            << result.cartesian_k4_compression_ratio
            << " order=" << result.order << " theta=" << result.theta
            << " eps_ratio=" << result.eps_ratio
            << " kappa=" << result.kappa
            << " near_integration=" << result.near_integration
            << " quad_buffer=" << result.quad_buffer
            << " condition_policy=" << result.condition_policy
            << " source_box_policy=" << result.source_box_policy
            << " m2l_eval_policy=" << result.m2l_eval_policy
            << " m2l_matrix_cache_policy="
            << result.m2l_matrix_cache_policy
            << " m2l_gemm_batch_policy="
            << result.m2l_gemm_batch_policy
            << " m2l_parallel_policy=" << result.m2l_parallel_policy
            << " m2l_parallel_workers=" << result.m2l_parallel_workers
            << " tree_stop_policy=" << result.tree_stop_policy
            << " jcp_channel_policy=" << result.jcp_channel_policy
            << " hierarchy_k4_policy=" << result.hierarchy_k4_policy
            << " hierarchy_k4_evaluation_count="
            << result.hierarchy_k4_evaluation_count
            << " tree_panel_radius_factor="
            << result.tree_panel_radius_factor
            << " tree_occupancy_stop=" << result.tree_occupancy_stop
            << " real_panel_guard=" << result.real_panel_guard
            << " hier_debug_mode=" << result.hier_debug_mode
            << " guard_margin_factor="
            << result.real_panel_guard_margin_factor
            << " guard_margin_max=" << result.real_panel_guard_margin_max
            << " aggregate_gate=" << result.aggregate_relative_rms_gate
            << " channel_gate=" << result.max_channel_relative_rms_gate
            << " aggregate_rms=" << result.aggregate_relative_rms
            << " k1_rms=" << result.k1_relative_rms
            << " k2_rms=" << result.k2_relative_rms
            << " k3_rms=" << result.k3_relative_rms
            << " k4_rms=" << result.k4_relative_rms
            << " max_channel_rms=" << result.max_channel_relative_rms
            << " max_channel_name=" << result.max_channel_name
            << " k4_abs_rms=" << result.k4_abs_rms
            << " k4_reference_rms=" << result.k4_reference_rms
            << " far_aggregate_rms=" << result.far_aggregate_relative_rms
            << " far_k1_rms=" << result.far_k1_relative_rms
            << " far_k2_rms=" << result.far_k2_relative_rms
            << " far_k3_rms=" << result.far_k3_relative_rms
            << " far_k4_rms=" << result.far_k4_relative_rms
            << " fallback_k4_rms=" << result.fallback_k4_relative_rms
            << " far_direct_valid_pairs=" << result.far_direct_valid_pairs
            << " source_box_boundary_margin_min="
            << result.source_box_boundary_margin_min
            << " target_check_excess_min=" << result.target_check_excess_min
            << " target_panel_envelope_excess_min="
            << result.target_panel_envelope_excess_min
            << " target_local_check_excess_min="
            << result.target_local_check_excess_min
            << " real_guard_rejected_pairs="
            << result.real_guard_rejected_pairs
            << " real_guard_split_pairs=" << result.real_guard_split_pairs
            << " real_guard_direct_pairs=" << result.real_guard_direct_pairs
            << " negative_target_check_excess_pairs="
            << result.negative_target_check_excess_pairs
            << " negative_target_check_excess_panel_pairs="
            << result.negative_target_check_excess_panel_pairs
            << " negative_target_local_check_excess_pairs="
            << result.negative_target_local_check_excess_pairs
            << " negative_target_local_check_excess_panel_pairs="
            << result.negative_target_local_check_excess_panel_pairs
            << " child_local_check_inside_parent_local_check_min_margin="
            << result.child_local_check_inside_parent_local_check_min_margin
            << " negative_l2l_local_containment_pairs="
            << result.negative_l2l_local_containment_pairs
            << " min_leaf_half_width_over_panel_radius="
            << result.min_leaf_half_width_over_panel_radius
            << " panels_crossing_node_boundary="
            << result.panels_crossing_node_boundary
            << " tree_node_count=" << result.tree_node_count
            << " tree_leaf_count=" << result.tree_leaf_count
            << " tree_max_depth_observed="
            << result.tree_max_depth_observed
            << " tree_radius_stop_leaf_count="
            << result.tree_radius_stop_leaf_count
            << " tree_occupancy_stop_leaf_count="
            << result.tree_occupancy_stop_leaf_count
            << " tree_leaf_size_stop_count="
            << result.tree_leaf_size_stop_count
            << " tree_max_depth_stop_count="
            << result.tree_max_depth_stop_count
            << " tree_unsplittable_stop_count="
            << result.tree_unsplittable_stop_count
            << " tree_min_leaf_half_width_over_panel_radius="
            << result.tree_min_leaf_half_width_over_panel_radius
            << " tree_max_leaf_panel_radius="
            << result.tree_max_leaf_panel_radius
            << " tree_mean_leaf_source_count="
            << result.tree_mean_leaf_source_count
            << " tree_mean_leaf_target_count="
            << result.tree_mean_leaf_target_count
            << " sampled_dense_valid_pairs="
            << result.sampled_dense_reference_valid_pairs
            << " production_fallback_pairs="
            << result.production_fallback_pairs
            << " production_invalid_pairs=" << result.production_invalid_pairs
            << " nonproduction_fallback_pairs="
            << result.nonproduction_fallback_pairs
            << " far_panel_pairs=" << result.far_panel_pairs
            << " direct_panel_pairs=" << result.direct_panel_pairs
            << " coverage=" << result.exact_pair_coverage
            << " expected_coverage="
            << (result.source_panels * result.sampled_rows)
            << " full_coverage_expected=" << result.full_coverage_expected
            << " sampled_coverage_expected="
            << result.sampled_coverage_expected
            << " far_fraction=" << result.far_fraction
            << " m2m_translations=" << h.m2m_translations
            << " m2l_translations=" << h.m2l_translations
            << " l2l_translations=" << h.l2l_translations
            << " l2p_target_count=" << h.l2p_target_count
            << " p2m_leaf_solves=" << h.p2m_leaf_solves
            << " local_node_solves=" << h.local_node_solves
            << " cache_entries=" << h.equivalent_factor_cache_entries
            << " cache_hits=" << h.equivalent_factor_cache_hits
            << " cache_misses=" << h.equivalent_factor_cache_misses
            << " p2m_seconds=" << h.p2m_seconds
            << " p2m_rhs_seconds=" << h.p2m_rhs_seconds
            << " p2m_solve_seconds=" << h.p2m_solve_seconds
            << " m2m_seconds=" << h.m2m_seconds
            << " m2l_seconds=" << h.m2l_seconds
            << " m2l_target_nodes=" << result.m2l_target_nodes
            << " m2l_source_nodes=" << result.m2l_source_nodes
            << " m2l_check_points_total="
            << result.m2l_check_points_total
            << " m2l_equivalent_points=" << result.m2l_equivalent_points
            << " m2l_kernel_evals_estimated="
            << result.m2l_kernel_evals_estimated
            << " m2l_seconds_per_translation="
            << result.m2l_seconds_per_translation
            << " m2l_seconds_per_check_point="
            << result.m2l_seconds_per_check_point
            << " m2l_matrix_cache_entries="
            << result.m2l_matrix_cache_entries
            << " m2l_matrix_cache_peak_entries="
            << result.m2l_matrix_cache_peak_entries
            << " m2l_matrix_cache_hits=" << result.m2l_matrix_cache_hits
            << " m2l_matrix_cache_misses="
            << result.m2l_matrix_cache_misses
            << " m2l_matrix_cache_hit_rate="
            << result.m2l_matrix_cache_hit_rate
            << " m2l_matrix_cache_limit_bytes="
            << result.m2l_matrix_cache_limit_bytes
            << " m2l_matrix_cache_evictions="
            << result.m2l_matrix_cache_evictions
            << " m2l_matrix_cache_estimated_bytes="
            << result.m2l_matrix_cache_estimated_bytes
            << " m2l_matrix_cache_peak_estimated_bytes="
            << result.m2l_matrix_cache_peak_estimated_bytes
            << " m2l_matrix_cache_total_built_bytes="
            << result.m2l_matrix_cache_total_built_bytes
            << " m2l_matrix_cache_build_seconds="
            << result.m2l_matrix_cache_build_seconds
            << " m2l_matrix_gemm_seconds=" << result.m2l_matrix_gemm_seconds
            << " m2l_matrix_cache_key_status="
            << result.m2l_matrix_cache_key_status
            << " m2l_gemm_group_count="
            << result.m2l_gemm_group_count
            << " m2l_gemm_grouped_pair_count="
            << result.m2l_gemm_grouped_pair_count
            << " m2l_gemm_grouped_rhs_columns="
            << result.m2l_gemm_grouped_rhs_columns
            << " m2l_gemm_max_group_size="
            << result.m2l_gemm_max_group_size
            << " m2l_parallel_task_count="
            << result.m2l_parallel_task_count
            << " m2l_parallel_group_count="
            << result.m2l_parallel_group_count
            << " m2l_parallel_max_task_pairs="
            << result.m2l_parallel_max_task_pairs
            << " m2l_parallel_prepare_seconds="
            << result.m2l_parallel_prepare_seconds
            << " m2l_parallel_eval_seconds="
            << result.m2l_parallel_eval_seconds
            << " m2l_parallel_merge_seconds="
            << result.m2l_parallel_merge_seconds
            << " local_solve_seconds=" << h.local_solve_seconds
            << " l2l_seconds=" << h.l2l_seconds
            << " l2p_seconds=" << h.l2p_seconds
            << " hierarchy_direct_fallback_seconds="
            << h.direct_fallback_seconds
            << " stage_timing_sum=" << stage_sum
            << " hierarchy_unattributed_seconds=" << unattributed
            << " hierarchy_bottleneck_stage=" << bottleneck.first
            << " hierarchy_bottleneck_fraction=" << bottleneck.second
            << " hierarchy_seconds=" << result.hierarchy_seconds
            << " flat_vs_dense_aggregate="
            << result.flat_vs_dense_aggregate_relative_rms
            << " flat_vs_dense_k4="
            << result.flat_vs_dense_k4_relative_rms
            << " flat_vs_dense_max_channel="
            << result.flat_vs_dense_max_channel_relative_rms
            << " hier_vs_dense_aggregate="
            << result.hierarchy_vs_dense_aggregate_relative_rms
            << " hier_vs_dense_k4="
            << result.hierarchy_vs_dense_k4_relative_rms
            << " hier_vs_flat_aggregate="
            << result.hierarchy_vs_flat_aggregate_relative_rms
            << " hier_vs_flat_k4="
            << result.hierarchy_vs_flat_k4_relative_rms
            << " local_check_rhs_relative_rms="
            << result.local_check_rhs_relative_rms
            << " local_check_rhs_max_abs_error="
            << result.local_check_rhs_max_abs_error
            << " l2p_oracle_aggregate="
            << result.l2p_oracle_aggregate_relative_rms
            << " l2p_oracle_k4=" << result.l2p_oracle_k4_relative_rms
            << " k4_laplace_grad_relative_rms="
            << result.k4_laplace_grad_relative_rms
            << " k4_yukawa_grad_relative_rms="
            << result.k4_yukawa_grad_relative_rms
            << " k4_combined_relative_rms="
            << result.k4_combined_relative_rms
            << " split_k4_laplace="
            << result.split_k4_laplace_relative_rms
            << " split_k4_yukawa="
            << result.split_k4_yukawa_relative_rms
            << " split_k4_combined="
            << result.split_k4_combined_relative_rms
            << " split_vs_dense_k4="
            << result.split_vs_dense_k4_relative_rms
            << " split_vs_flat_k4="
            << result.split_vs_flat_k4_relative_rms
            << " hybrid_mode=" << result.hybrid_mode
            << " hierarchy_inner_status=" << result.hierarchy_inner_status
            << " hierarchy_inner_accuracy_failure="
            << result.hierarchy_inner_accuracy_failure
            << " m2e_residual=" << result.m2e_max_relative_residual
            << " local_equiv_residual="
            << result.local_equiv_max_relative_residual
            << " e2l_coeff_rms=" << result.e2l_coefficient_relative_rms
            << " e2l_coeff_status=" << result.e2l_coefficient_status
            << " k4_decomposition_status="
            << result.k4_decomposition_status
            << " l2p_gradient_method=" << result.l2p_gradient_method
            << " k1_source=" << result.k1_source
            << " k2_source=" << result.k2_source
            << " k3_source=" << result.k3_source
            << " k4_source=" << result.k4_source
            << " hierarchy_k123_max_channel_rms="
            << result.hierarchy_k123_max_channel_relative_rms
            << " hierarchy_k123_max_channel_name="
            << result.hierarchy_k123_max_channel_name
            << " hybrid_recombined_aggregate_rms="
            << result.hybrid_recombined_aggregate_relative_rms
            << " hybrid_recombined_max_channel_rms="
            << result.hybrid_recombined_max_channel_relative_rms
            << " hybrid_recombined_max_channel_name="
            << result.hybrid_recombined_max_channel_name
            << " channel_hybrid_double_count_guard="
            << result.channel_hybrid_double_count_guard
            << " channel_hybrid_recombined_hash="
            << result.channel_hybrid_recombined_hash
            << " cartesian_k4_status=" << result.cartesian_k4_status
            << " cartesian_k4_execution_status="
            << result.cartesian_k4_execution_status
            << " cartesian_k4_accuracy_status="
            << result.cartesian_k4_accuracy_status
            << " cartesian_k4_accuracy_failure="
            << result.cartesian_k4_accuracy_failure
            << " cartesian_k4_rms=" << result.cartesian_k4_relative_rms
            << " cartesian_k4_abs_rms=" << result.cartesian_k4_abs_rms
            << " cartesian_k4_reference_rms="
            << result.cartesian_k4_reference_rms
            << " cartesian_k4_far_pairs="
            << result.cartesian_k4_far_panel_pairs
            << " cartesian_k4_fallback_pairs="
            << result.cartesian_k4_fallback_panel_pairs
            << " cartesian_k4_total_seconds="
            << result.cartesian_k4_total_seconds
            << " cartesian_k4_policy_hash="
            << result.cartesian_k4_policy_hash
            << " shared_worklist_hash=" << result.shared_worklist_hash
            << " shared_theta=" << result.shared_theta
            << " shared_guard=" << result.shared_guard
            << " shared_far_pairs=" << result.shared_far_pairs
            << " shared_direct_pairs=" << result.shared_direct_pairs
            << " shared_far_fraction=" << result.shared_far_fraction
            << " cartesian_k4_shared_far_pairs="
            << result.cartesian_k4_shared_far_pairs
            << " cartesian_k4_shared_direct_pairs="
            << result.cartesian_k4_shared_direct_pairs
            << " direct_avx2_all_channel_pairs="
            << result.direct_avx2_all_channel_pairs
            << " worklist_consistency_status="
            << result.worklist_consistency_status
            << " kifmm_worklist_hash=" << result.kifmm_worklist_hash
            << " kifmm_far_pairs=" << result.kifmm_far_pairs
            << " kifmm_direct_pairs=" << result.kifmm_direct_pairs
            << " kifmm_far_fraction=" << result.kifmm_far_fraction
            << " dual_worklist_status=" << result.dual_worklist_status
            << " channel_recombination_status="
            << result.channel_recombination_status
            << " worklist_independence_status="
            << result.worklist_independence_status
            << " hybrid_vs_dense_k4="
            << result.hierarchy_vs_dense_k4_relative_rms
            << " hybrid_vs_flat_k4="
            << result.hierarchy_vs_flat_k4_relative_rms
            << " pure_hierarchy_vs_dense_k4="
            << result.hierarchy_vs_dense_k4_relative_rms
            << " flat_output_hash=" << result.flat_output_hash
            << " production_fallback_seconds="
            << result.production_fallback_seconds
            << " dense_oracle_seconds=" << result.dense_oracle_seconds
            << " total_seconds=" << result.total_seconds
            << " tree_hash=" << result.tree_hash
            << " worklist_hash=" << result.worklist_hash
            << " policy_hash=" << result.policy_hash
            << " panel_geometry_hash=" << result.panel_geometry_hash
            << " unknown_vector_hash=" << result.unknown_vector_hash
            << " sampled_target_hash=" << result.sampled_target_hash
            << " output_hash=" << result.output_hash
            << " lite_status=" << result.lite_status
            << " lite_aggregate_rms=" << result.lite_aggregate_relative_rms
            << " lite_max_channel_rms=" << result.lite_max_channel_relative_rms
            << " lite_total_seconds=" << result.lite_total_seconds << '\n';
  if (emit_e35_attribution) {
    std::cout << std::setprecision(12)
              << "E35_ATTR fixture=" << fixture << " sample=" << sample_name
              << " source_panels=" << result.source_panels
              << " target_panels=" << result.target_panels
              << " sampled_rows=" << result.sampled_rows
              << " full_coverage_expected=" << result.full_coverage_expected
              << " sampled_coverage_expected="
              << result.sampled_coverage_expected
              << " far_panel_pairs=" << result.far_panel_pairs
              << " direct_panel_pairs=" << result.direct_panel_pairs
              << " production_fallback_pairs="
              << result.production_fallback_pairs
              << " cartesian_k4_far_pairs="
              << result.cartesian_k4_far_panel_pairs
              << " cartesian_k4_fallback_pairs="
              << result.cartesian_k4_fallback_panel_pairs
              << " kifmm_compression_ratio="
              << result.kifmm_compression_ratio
              << " cartesian_k4_compression_ratio="
              << result.cartesian_k4_compression_ratio
              << " hierarchy_seconds=" << result.hierarchy_seconds
              << " cartesian_k4_total_seconds="
              << result.cartesian_k4_total_seconds
              << " production_fallback_seconds="
              << result.production_fallback_seconds
              << " dense_oracle_seconds=" << result.dense_oracle_seconds
              << " total_seconds=" << result.total_seconds << '\n';
  }
  if (emit_e36_stage_attribution) {
    std::cout << std::setprecision(12)
              << "E36_STAGE fixture=" << fixture << " sample=" << sample_name
              << " source_panels=" << result.source_panels
              << " target_panels=" << result.target_panels
              << " sampled_rows=" << result.sampled_rows
              << " far_panel_pairs=" << result.far_panel_pairs
              << " direct_panel_pairs=" << result.direct_panel_pairs
              << " far_fraction=" << result.far_fraction
              << " kifmm_compression_ratio="
              << result.kifmm_compression_ratio
              << " cartesian_k4_compression_ratio="
              << result.cartesian_k4_compression_ratio
              << " p2m_seconds=" << h.p2m_seconds
              << " p2m_rhs_seconds=" << h.p2m_rhs_seconds
              << " p2m_solve_seconds=" << h.p2m_solve_seconds
              << " m2m_seconds=" << h.m2m_seconds
              << " m2l_seconds=" << h.m2l_seconds
              << " local_solve_seconds=" << h.local_solve_seconds
              << " l2l_seconds=" << h.l2l_seconds
              << " l2p_seconds=" << h.l2p_seconds
              << " hierarchy_direct_fallback_seconds="
              << h.direct_fallback_seconds
              << " hierarchy_seconds=" << result.hierarchy_seconds
              << " stage_timing_sum=" << stage_sum
              << " hierarchy_unattributed_seconds=" << unattributed
              << " hierarchy_bottleneck_stage=" << bottleneck.first
              << " hierarchy_bottleneck_fraction=" << bottleneck.second
              << " p2m_leaf_solves=" << h.p2m_leaf_solves
              << " m2m_translations=" << h.m2m_translations
              << " m2m_parent_solves=" << h.m2m_parent_solves
              << " m2l_translations=" << h.m2l_translations
              << " local_node_solves=" << h.local_node_solves
              << " l2l_translations=" << h.l2l_translations
              << " l2p_target_count=" << h.l2p_target_count
              << " cache_entries=" << h.equivalent_factor_cache_entries
              << " cache_hits=" << h.equivalent_factor_cache_hits
              << " cache_misses=" << h.equivalent_factor_cache_misses
              << " hierarchy_status=" << h.status
              << " hierarchy_accuracy_status=" << h.accuracy_status
              << " total_seconds=" << result.total_seconds << '\n';
  }
  for (int bin = 0; bin < 5; ++bin) {
    std::cout << std::setprecision(12)
              << "E255B_ATTR row_sample=" << sample_name << " bin=" << bin
              << " mac_cell_pairs=" << result.mac_bin_cell_pairs[bin]
              << " mac_panel_pairs=" << result.mac_bin_panel_pairs[bin]
              << " min_target_check_excess="
              << result.mac_bin_min_target_check_excess[bin] << '\n';
  }
  for (const auto &field : result.field_attribution) {
    std::cout << std::setprecision(12)
              << "E258_FIELD sample=" << sample_name
              << " mode=" << result.hier_debug_mode
              << " scope=" << field.scope
              << " field_name=" << field.field_name
              << " kernel=" << field.kernel
              << " node_id=" << field.node_id
              << " node_level=" << field.node_level
              << " relative_rms=" << field.relative_rms
              << " abs_rms=" << field.abs_rms
              << " reference_rms=" << field.reference_rms
              << " max_abs_error=" << field.max_abs_error << '\n';
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    const auto args = parse_args(argc, argv);
    const std::string mode = get_arg(args, "mode", "profile");
    const std::string fixture = get_arg(args, "fixture", "1SH1");
    const std::string sample_name = get_arg(args, "sample", "uniform_hash");
    const std::size_t sample_count =
        static_cast<std::size_t>(std::stoul(get_arg(args, "sample-count", "128")));
    const int order = std::stoi(get_arg(args, "order", "8"));
    const int surface_points =
        std::stoi(get_arg(args, "surface-points", "8"));
    const std::size_t panel_limit =
        static_cast<std::size_t>(std::stoul(get_arg(args, "panel-limit", "0")));
    const std::string pqr_1sh1 =
        get_arg(args, "pqr-1sh1", "test/pbe_inputs/1sh1/1sh1.pqr");
    const std::string pqr_1a63 =
        get_arg(args, "pqr-1a63", "test/pbe_inputs/1a63/1a63.pqr");
    const std::string pqr_path = fixture == "1A63" ? pqr_1a63 : pqr_1sh1;
    const bool emit_e35_attribution =
        get_arg(args, "emit-e35-attribution", "0") == "1";
    const bool emit_e36_stage_attribution =
        get_arg(args, "emit-e36-stage-attribution", "0") == "1";

    std::vector<Panel> panels = build_real_surface(pqr_path, surface_points);
    if (panel_limit > 0 && panel_limit < panels.size()) {
      panels.resize(panel_limit);
    }
    std::vector<double> unknowns = unknowns_for(panels);
    std::vector<std::size_t> sample =
        sample_for_name(panels, sample_name, sample_count);
    KifmmHierarchyJcpConfig config = hierarchy_config(order);
    apply_hierarchy_overrides(config, args);

    if (mode == "invalid-duplicates") {
      if (sample.size() >= 2) {
        sample[1] = sample[0];
      }
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "duplicate sample rejected");
      return EXIT_SUCCESS;
    }
    if (mode == "invalid-range") {
      sample.push_back(panels.size());
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "out-of-range sample rejected");
      return EXIT_SUCCESS;
    }
    if (mode == "invalid-nonfinite-unknown") {
      unknowns[3] = std::numeric_limits<double>::quiet_NaN();
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "nonfinite unknown rejected");
      return EXIT_SUCCESS;
    }
    if (mode == "invalid-area") {
      panels[sample.front()].area = 0.0;
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "zero area rejected");
      return EXIT_SUCCESS;
    }
    if (mode == "invalid-normal") {
      panels[sample.front()].nx *= 2.0;
      panels[sample.front()].ny *= 2.0;
      panels[sample.front()].nz *= 2.0;
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "non-unit normal rejected");
      return EXIT_SUCCESS;
    }
    if (mode == "invalid-triangle") {
      panels[sample.front()].has_triangle_geometry = false;
      const auto result =
          run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample, config);
      print_invalid(mode, result);
      require(result.execution_status == "fail", "missing triangle rejected");
      return EXIT_SUCCESS;
    }

    const bool full_panel_profile = mode == "full-panel-profile";
    const auto result =
        full_panel_profile
            ? run_kifmm_hierarchy_real_panel_full_apply_sampled_audit(
                  panels, unknowns, sample, config)
            : run_kifmm_hierarchy_real_panel_sampled(panels, unknowns, sample,
                                                     config);
    print_profile(fixture, sample_name + "_" + std::to_string(sample.size()),
                  result, emit_e35_attribution, emit_e36_stage_attribution);
    require(result.execution_status == "pass", "profile execution passes");
    require(result.sampled_dense_reference_valid_pairs ==
                result.sampled_coverage_expected,
            "dense reference coverage exact");
    require(result.production_invalid_pairs == 0, "production invalid zero");
    require(result.nonproduction_fallback_pairs == 0,
            "nonproduction fallback zero");
    require(result.exact_pair_coverage == result.full_coverage_expected,
            "hierarchy full coverage exact");
    if (full_panel_profile) {
      require(result.target_panels == result.source_panels,
              "full-panel target count matches source panels");
      require(result.timing_scope == "full_panel_apply_sampled_audit",
              "full-panel timing scope reported");
    }
    require(result.hierarchy.m2l_translations > 0, "M2L active");
    require(result.output_hash != 0, "output hash nonzero");
    if (!full_panel_profile) {
      require(std::isfinite(result.far_k4_relative_rms),
              "far K4 attribution finite");
    }
    require(std::isfinite(result.source_box_boundary_margin_min),
            "source-box margin attribution finite");
    require(std::isfinite(result.target_check_excess_min),
            "target check excess attribution finite");

    KifmmRealPanelCloudConfig lite = lite_config();
    const auto lite_result =
        run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample, lite);
    std::cout << std::setprecision(12)
              << "E255_LITE fixture=" << fixture << " sample=" << sample_name
              << "_" << sample.size() << " status=" << lite_result.status
              << " aggregate_rms=" << lite_result.lite.aggregate_relative_rms
              << " max_channel_rms="
              << std::max({lite_result.lite.k1_relative_rms,
                           lite_result.lite.k2_relative_rms,
                           lite_result.lite.k3_relative_rms,
                           lite_result.lite.k4_relative_rms})
              << " total_seconds=" << lite_result.total_seconds << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
