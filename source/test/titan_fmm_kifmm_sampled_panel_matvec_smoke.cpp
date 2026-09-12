#include "tabi/tabi_fmm_kifmm_sampled_panel_matvec.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#if defined(TITAN_ENABLE_E16_CARTESIAN_COMPARATOR) &&                       \
    !defined(TITAN_ENABLE_E15_CARTESIAN_COMPARATOR)
#define TITAN_ENABLE_E15_CARTESIAN_COMPARATOR 1
#endif

#ifdef TITAN_ENABLE_E15_CARTESIAN_COMPARATOR
#include "tabi/tabi_fmm_cartesian_sampled_panel_matvec.h"
#endif

#include <algorithm>
#include <cctype>
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
#include <vector>

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Options;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::Result;
using titan::tabi::fmm_kifmm::KifmmRealPanelCloudConfig;
using titan::tabi::fmm_kifmm::KifmmSampledPanelMatvecResult;
using titan::tabi::fmm_kifmm::Vec3;
#ifdef TITAN_ENABLE_E15_CARTESIAN_COMPARATOR
using titan::tabi::fmm_kifmm::CartesianSampledPanelMatvecConfig;
using titan::tabi::fmm_kifmm::CartesianSampledPanelMatvecResult;
using titan::tabi::fmm_kifmm::run_guarded_cartesian_p4_sampled_panel_matvec;
#endif
using titan::tabi::fmm_kifmm::run_kifmm_lite_sampled_panel_matvec;

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

std::vector<Panel> build_real_surface(const std::string &pqr_path) {
  Options options;
  options.surface_mode = "vdw";
  options.surface_geometry = "triangle";
  options.diagonal_mode = "solid-angle";
  options.surface_point_count = 8;
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

KifmmRealPanelCloudConfig base_config(int order, double theta = 0.6,
                                      double min_far_panel_fraction = 0.0) {
  KifmmRealPanelCloudConfig config;
  config.lite.order = order;
  config.lite.leaf_size = 32;
  config.lite.opening_theta = theta;
  config.lite.equiv_scale = 1.05;
  config.lite.check_scale = 2.40;
  config.lite.kappa = 0.12746;
  config.lite.eps_ratio = 80.0;
  config.lite.max_condition_estimate = 1.0e12;
  config.lite.max_equivalent_relative_residual = 1.0e-10;
  config.lite.max_aggregate_relative_rms = 1.0e-6;
  config.lite.max_channel_relative_rms = 5.0e-6;
  config.lite.min_far_panel_fraction = min_far_panel_fraction;
  config.lite.require_far_work = true;
  config.lite.require_unit_normals = true;
  config.lite.solver = "qr";
  config.lite.condition_policy = "qr-rank-only";
  config.near_integration = "analytic_singular";
  config.quad_buffer = "off";
  return config;
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

std::vector<std::size_t>
inward_normal_groove_sample(const std::vector<Panel> &panels,
                            std::size_t count) {
  const Vec3 c = centroid(panels);
  return select_by_key(panels, count, false, [&](const Panel &p, std::size_t) {
    return dot(unit(sub(position(p), c)), unit(normal(p)));
  });
}

std::vector<std::size_t> high_aspect_sample(const std::vector<Panel> &panels,
                                            std::size_t count) {
  return select_by_key(panels, count, true, [](const Panel &p, std::size_t) {
    return aspect_score(p);
  });
}

std::vector<std::size_t>
closest_opposing_surface_sample(const std::vector<Panel> &panels,
                                std::size_t count) {
  const std::size_t shortlist_count =
      std::min<std::size_t>(panels.size(), std::max<std::size_t>(1024, 8 * count));
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
  for (const std::size_t i : targets) {
    const Vec3 pi = position(panels[i]);
    const Vec3 ni = unit(normal(panels[i]));
    double best_sq = std::numeric_limits<double>::infinity();
    for (const std::size_t j : search) {
      if (i == j) {
        continue;
      }
      if (dot(ni, unit(normal(panels[j]))) > -0.65) {
        continue;
      }
      const Vec3 d = sub(pi, position(panels[j]));
      best_sq = std::min(best_sq, dot(d, d));
    }
    if (std::isfinite(best_sq)) {
      scored.push_back({best_sq, i});
    }
  }
  std::sort(scored.begin(), scored.end(),
            [](const auto &a, const auto &b) {
              if (a.first == b.first) {
                return a.second < b.second;
              }
              return a.first < b.first;
            });
  std::vector<std::size_t> idx;
  idx.reserve(std::min(count, scored.size()));
  for (const auto &entry : scored) {
    if (idx.size() >= count) {
      break;
    }
    idx.push_back(entry.second);
  }
  std::sort(idx.begin(), idx.end());
  return idx;
}

double max_channel_rms(const KifmmSampledPanelMatvecResult &result) {
  return std::max({result.lite.k1_relative_rms, result.lite.k2_relative_rms,
                   result.lite.k3_relative_rms, result.lite.k4_relative_rms});
}

void print_result(const std::string &fixture, const std::string &sample,
                  const KifmmSampledPanelMatvecResult &result,
                  double theta = 0.6,
                  const std::string &condition_policy = "qr-rank-only",
                  double eps_ratio = std::numeric_limits<double>::quiet_NaN(),
                  double kappa = std::numeric_limits<double>::quiet_NaN(),
                  const std::string &near_integration = "",
                  const std::string &quad_buffer = "",
                  const std::string &fallback_kernel = "") {
  std::cout << std::setprecision(12)
            << "KIFMM_SAMPLED_PANEL_CASE fixture=" << fixture
            << " sample=" << sample << " status=" << result.status
            << " source_panels=" << result.lite.source_panels
            << " sampled_rows=" << result.sampled_target_indices.size()
            << " order=" << result.order
            << " theta=" << theta
            << " condition_policy=" << condition_policy
            << " source_box_policy=" << result.source_box_policy
            << " min_far_panel_fraction="
            << result.lite.min_far_panel_fraction;
  if (std::isfinite(eps_ratio)) {
    std::cout << " eps_ratio=" << eps_ratio;
  }
  if (std::isfinite(kappa)) {
    std::cout << " kappa=" << kappa;
  }
  if (!near_integration.empty()) {
    std::cout << " near_integration=" << near_integration;
  }
  if (!quad_buffer.empty()) {
    std::cout << " quad_buffer=" << quad_buffer;
  }
  if (!fallback_kernel.empty()) {
    std::cout << " fallback_kernel=" << fallback_kernel;
  }
  std::cout
            << " far_cell_pairs=" << result.lite.accepted_far_cell_pairs
            << " far_panel_pairs=" << result.lite.accepted_far_panel_pairs
            << " fallback_panel_pairs="
            << result.lite.direct_fallback_panel_pairs
            << " production_fallback_valid_pairs="
            << result.production_fallback_valid_pairs
            << " production_direct_reference_valid_pairs="
            << result.production_direct_reference_valid_pairs
            << " invalid_pairs=" << result.production_invalid_pairs
            << " nonproduction_fallback_pairs="
            << result.nonproduction_fallback_pairs
            << " aggregate_rms=" << result.lite.aggregate_relative_rms
            << " k1=" << result.lite.k1_relative_rms
            << " k2=" << result.lite.k2_relative_rms
            << " k3=" << result.lite.k3_relative_rms
            << " k4=" << result.lite.k4_relative_rms
            << " aggregate_abs_rms=" << result.aggregate_abs_rms
            << " k1_abs_rms=" << result.k1_abs_rms
            << " k2_abs_rms=" << result.k2_abs_rms
            << " k3_abs_rms=" << result.k3_abs_rms
            << " k4_abs_rms=" << result.k4_abs_rms
            << " k1_reference_rms=" << result.k1_reference_rms
            << " k2_reference_rms=" << result.k2_reference_rms
            << " k3_reference_rms=" << result.k3_reference_rms
            << " k4_reference_rms=" << result.k4_reference_rms
            << " max_equiv_residual="
            << result.lite.max_equivalent_relative_residual
            << " max_condition=" << result.lite.max_condition_estimate
            << " condition_computed="
            << (result.condition_computed ? "true" : "false")
            << " condition_gate_enabled="
            << (result.condition_gate_enabled ? "true" : "false")
            << " direct_reference_seconds="
            << result.direct_reference_seconds
            << " tree_build_seconds=" << result.tree_build_seconds
            << " worklist_build_seconds=" << result.worklist_build_seconds
            << " equivalent_build_solve_seconds="
            << result.equivalent_build_solve_seconds
            << " equivalent_matrix_build_seconds="
            << result.equivalent_matrix_build_seconds
            << " equivalent_condition_seconds="
            << result.equivalent_condition_seconds
            << " equivalent_factor_seconds=" << result.equivalent_factor_seconds
            << " equivalent_rhs_build_seconds="
            << result.equivalent_rhs_build_seconds
            << " equivalent_solve_seconds=" << result.equivalent_solve_seconds
            << " equivalent_factor_cache_hits="
            << result.equivalent_factor_cache_hits
            << " equivalent_factor_cache_misses="
            << result.equivalent_factor_cache_misses
            << " cache_hit_rate=" << result.equivalent_cache_hit_rate
            << " cache_entries=" << result.equivalent_system_cache_entries
            << " unique_half_widths="
            << result.unique_equivalent_half_widths
            << " unique_kappa_h_values="
            << result.unique_equivalent_kappa_h_values
            << " source_tree_levels=" << result.source_tree_levels
            << " unique_source_node_half_widths="
            << result.unique_source_node_half_widths
            << " source_root_half_width=" << result.source_root_half_width
            << " level_cube_geometry="
            << (result.level_cube_geometry ? "true" : "false")
            << " half_width_min=" << result.equivalent_half_width_min
            << " half_width_max=" << result.equivalent_half_width_max
            << " laplace_systems=" << result.equivalent_laplace_systems
            << " yukawa_systems=" << result.equivalent_yukawa_systems
            << " systems_per_source_node="
            << result.equivalent_systems_per_source_node
            << " far_eval_seconds=" << result.far_eval_seconds
            << " fallback_eval_seconds=" << result.fallback_eval_seconds
            << " total_seconds=" << result.total_seconds
            << " unique_far_source_nodes=" << result.unique_far_source_nodes
            << " cached_source_nodes=" << result.cached_source_nodes
            << " source_tree_bytes=" << result.source_tree_bytes
            << " target_tree_bytes=" << result.target_tree_bytes
            << " worklist_bytes=" << result.worklist_bytes
            << " equivalent_charge_storage_bytes="
            << result.equivalent_charge_storage_bytes
            << " cached_channel_object_estimated_bytes="
            << result.cached_channel_object_estimated_bytes
            << " max_equivalent_rank=" << result.max_equivalent_rank
            << " min_equivalent_sigma_min="
            << result.min_equivalent_sigma_min
            << " max_equivalent_sigma_max="
            << result.max_equivalent_sigma_max
            << " max_equivalent_rhs_norm=" << result.max_equivalent_rhs_norm
            << " max_equivalent_charge_norm="
            << result.max_equivalent_charge_norm
            << " max_equivalent_charge_to_rhs_norm="
            << result.max_equivalent_charge_to_rhs_norm
            << " panel_geometry_hash=" << result.panel_geometry_hash
            << " unknown_vector_hash=" << result.unknown_vector_hash
            << " sampled_target_hash=" << result.sampled_target_hash
            << " worklist_hash=" << result.lite.worklist_hash
            << " policy_hash=" << result.policy_hash << '\n';
}

void require_common_clean(const std::string &name,
                          const KifmmSampledPanelMatvecResult &result,
                          std::size_t source_count) {
  require(result.production_invalid_pairs == 0,
          name + " has no invalid production near kernel calls");
  require(result.nonproduction_fallback_pairs == 0,
          name + " has no nonproduction fallback pairs");
  require(!result.used_nonproduction_fallback,
          name + " does not flag nonproduction fallback");
  require(result.production_fallback_valid_pairs ==
              result.lite.direct_fallback_panel_pairs,
          name + " production fallback count matches fallback work");
  require(result.production_direct_reference_valid_pairs ==
              result.sampled_target_indices.size() * source_count,
          name + " dense direct reference covers sampled rows");
  require(result.lite.worklist_hash != 0, name + " worklist hash recorded");
  require(result.policy_hash != 0, name + " policy hash recorded");
}

void require_candidate(const std::string &fixture, const std::string &sample,
                       const KifmmSampledPanelMatvecResult &result,
                       std::size_t source_count) {
  print_result(fixture, sample, result);
  const std::string name = fixture + "/" + sample;
  require(result.finite, name + " finite");
  require(result.status == "pass", name + " status pass");
  require_common_clean(name, result, source_count);
  require(result.lite.accepted_far_panel_pairs > 0,
          name + " has accepted far work");
  require(result.lite.direct_fallback_panel_pairs > 0,
          name + " has production direct fallback work");
  require(result.lite.aggregate_relative_rms <= 1.0e-6,
          name + " aggregate RMS gate");
  require(max_channel_rms(result) <= 5.0e-6, name + " channel RMS gate");
  require(result.lite.max_equivalent_relative_residual <= 1.0e-10,
          name + " equivalent residual gate");
  require(result.lite.max_condition_estimate <= 1.0e12,
          name + " condition gate");
}

bool passes_candidate_gate(const KifmmSampledPanelMatvecResult &result) {
  return result.finite && result.status == "pass" &&
         result.lite.accepted_far_panel_pairs > 0 &&
         result.lite.direct_fallback_panel_pairs > 0 &&
         result.production_invalid_pairs == 0 &&
         result.nonproduction_fallback_pairs == 0 &&
         result.lite.aggregate_relative_rms <= 1.0e-6 &&
         max_channel_rms(result) <= 5.0e-6 &&
         result.lite.max_equivalent_relative_residual <= 1.0e-10 &&
         result.lite.max_condition_estimate <= 1.0e12;
}

void require_diagnostic_clean(const std::string &fixture,
                              const std::string &sample,
                              const KifmmSampledPanelMatvecResult &result,
                              std::size_t source_count) {
  print_result(fixture, sample, result);
  require_common_clean(fixture + "/" + sample, result, source_count);
  require(result.lite.accepted_far_panel_pairs > 0,
          fixture + "/" + sample + " has accepted far work");
  require(result.lite.direct_fallback_panel_pairs > 0,
          fixture + "/" + sample + " has production fallback work");
  require(result.status == "pass" ||
              result.status == "aggregate_rms_gate_failed" ||
              result.status == "channel_rms_gate_failed",
          fixture + "/" + sample + " has interpretable diagnostic status");
  require(result.lite.max_equivalent_relative_residual <= 1.0e-10,
          fixture + "/" + sample + " equivalent residual gate");
  require(result.lite.max_condition_estimate <= 1.0e12,
          fixture + "/" + sample + " condition gate");
}

KifmmSampledPanelMatvecResult
run_case(const std::vector<Panel> &panels, const std::vector<double> &unknowns,
         const std::vector<std::size_t> &sample, int order,
         double theta = 0.6,
         const std::string &condition_policy = "qr-rank-only",
         const std::string &source_box_policy = "tight-bbox",
         double min_far_panel_fraction = 0.0) {
  KifmmRealPanelCloudConfig config =
      base_config(order, theta, min_far_panel_fraction);
  config.lite.condition_policy = condition_policy;
  config.lite.source_box_policy = source_box_policy;
  return run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample, config);
}

void require_repeat_hash(const std::string &name,
                         const KifmmSampledPanelMatvecResult &a,
                         const KifmmSampledPanelMatvecResult &b) {
  require(a.lite.worklist_hash == b.lite.worklist_hash,
          name + " worklist hash is deterministic");
  require(a.policy_hash == b.policy_hash,
          name + " policy hash is deterministic");
  require(a.panel_geometry_hash == b.panel_geometry_hash,
          name + " panel geometry hash is deterministic");
  require(a.sampled_target_hash == b.sampled_target_hash,
          name + " sampled target hash is deterministic");
  require(a.unique_source_node_half_widths == b.unique_source_node_half_widths,
          name + " source half-width count is deterministic");
  require(a.lite.accepted_far_panel_pairs == b.lite.accepted_far_panel_pairs,
          name + " far pair count is deterministic");
  require(a.lite.direct_fallback_panel_pairs ==
              b.lite.direct_fallback_panel_pairs,
          name + " fallback pair count is deterministic");
}

void run_e14_far_fraction_floor_smoke(const std::string &pqr_1sh1) {
  const std::vector<Panel> panels = build_real_surface(pqr_1sh1);
  const std::vector<double> unknowns = unknowns_for(panels);
  const auto sample = uniform_hash_sample(panels, 128);
  KifmmRealPanelCloudConfig baseline_config = base_config(6);
  baseline_config.lite.source_box_policy = "level-cube";
  baseline_config.lite.min_far_panel_fraction = 0.0;
  const KifmmSampledPanelMatvecResult baseline =
      run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample,
                                          baseline_config);
  print_result("1SH1", "e14_far_floor_off", baseline, 0.6,
               baseline_config.lite.condition_policy);
  require(baseline.status == "pass", "far floor off preserves pass");

  KifmmRealPanelCloudConfig q4_config = baseline_config;
  q4_config.lite.order = 4;
  const KifmmSampledPanelMatvecResult q4 =
      run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample, q4_config);
  print_result("1SH1", "e14_q4_failure_timing", q4, 0.6,
               q4_config.lite.condition_policy);
  require(q4.status == "aggregate_rms_gate_failed",
          "q4 keeps expected aggregate RMS failure");
  require(q4.total_seconds > 0.0, "q4 failure row records elapsed time");

  KifmmRealPanelCloudConfig guarded_config = baseline_config;
  guarded_config.lite.min_far_panel_fraction = 0.99;
  const KifmmSampledPanelMatvecResult guarded =
      run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample,
                                          guarded_config);
  print_result("1SH1", "e14_far_floor_forced", guarded, 0.6,
               guarded_config.lite.condition_policy);
  require(!guarded.finite, "fallback-dominated row is not finite/pass");
  require(guarded.status == "fallback_dominated",
          "high far-fraction floor reports fallback_dominated");
  require(guarded.lite.accepted_far_panel_pairs ==
              baseline.lite.accepted_far_panel_pairs,
          "far floor does not rebuild worklist geometry");
  require(guarded.policy_hash != baseline.policy_hash,
          "far floor changes policy hash");
  require(guarded.total_seconds > 0.0,
          "fallback-dominated failure row records elapsed time");
  std::cout << "PASS titan_fmm_kifmm_sampled_panel_matvec_smoke E14\n";
}

#ifdef TITAN_ENABLE_E15_CARTESIAN_COMPARATOR
void run_e15_cartesian_comparator_case(const std::string &fixture_label,
                                       const std::string &sample_label,
                                       const std::vector<Panel> &panels,
                                       const std::vector<double> &unknowns,
                                       const std::vector<std::size_t> &sample,
                                       bool enforce_smoke_accuracy,
                                       const CartesianSampledPanelMatvecConfig
                                           *cartesian_override = nullptr,
                                       int kifmm_order = 6) {
  KifmmRealPanelCloudConfig kifmm_config =
      base_config(kifmm_order, 0.6, 0.10);
  kifmm_config.lite.source_box_policy = "level-cube";
  const KifmmSampledPanelMatvecResult kifmm =
      run_kifmm_lite_sampled_panel_matvec(panels, unknowns, sample,
                                          kifmm_config);
  print_result(fixture_label, "e15_kifmm_" + sample_label, kifmm, 0.6,
               kifmm_config.lite.condition_policy,
               kifmm_config.lite.eps_ratio, kifmm_config.lite.kappa,
               kifmm_config.near_integration, kifmm_config.quad_buffer,
               "production_near_kernel");
  require(kifmm.status == "pass", "E15 KIFMM baseline passes");

  CartesianSampledPanelMatvecConfig cart_config;
  if (cartesian_override != nullptr)
    cart_config = *cartesian_override;
  cart_config.eps_ratio = kifmm_config.lite.eps_ratio;
  cart_config.kappa = kifmm_config.lite.kappa;
  require(std::abs(cart_config.eps_ratio - kifmm_config.lite.eps_ratio) <
              1.0e-14,
          "Cartesian and KIFMM eps_ratio match");
  require(std::abs(cart_config.kappa - kifmm_config.lite.kappa) < 1.0e-14,
          "Cartesian and KIFMM kappa match");
  require(cart_config.near_integration == kifmm_config.near_integration,
          "Cartesian and KIFMM near integration match");
  require(cart_config.quad_buffer == kifmm_config.quad_buffer,
          "Cartesian and KIFMM quad buffer match");
  const CartesianSampledPanelMatvecResult cart =
      run_guarded_cartesian_p4_sampled_panel_matvec(panels, unknowns, sample,
                                                    cart_config);
  require(cart.timing_scope == "sampled_rows_only",
          "Cartesian p4 comparator reports true sampled-row timing scope");
  require(cart.pair_count_scope == "sampled_rows_only",
          "Cartesian p4 comparator reports true sampled-row pair counts");
  std::cout << "CARTESIAN_SAMPLE_RESULT"
            << " fixture=" << fixture_label
            << " sample=" << sample_label
            << " status=" << cart.status
            << " execution_status=" << cart.execution_status
            << " accuracy_status=" << cart.accuracy_status
            << " accuracy_failure=" << cart.accuracy_failure
            << " order=" << cart.order
            << " theta=" << cart_config.theta
            << " tail_tol=" << cart_config.tail_tol
            << " geometry_theta=" << cart_config.strict_geometry_theta
            << " cache_mode=" << cart_config.cache_mode
            << " eps_ratio=" << cart_config.eps_ratio
            << " kappa=" << cart_config.kappa
            << " near_integration=" << cart_config.near_integration
            << " quad_buffer=" << cart_config.quad_buffer
            << " fallback_kernel=" << cart_config.fallback_kernel
            << " fallback_kernel_effective="
            << cart.fallback_kernel_effective
            << " timing_scope=" << cart.timing_scope
            << " pair_count_scope=" << cart.pair_count_scope
            << " source_panels=" << cart.source_panels
            << " sampled_rows=" << cart.sampled_rows
            << " aggregate_rms=" << cart.aggregate_relative_rms
            << " max_channel_rms=" << cart.max_channel_relative_rms
            << " max_channel_name=" << cart.max_channel_name
            << " k4_abs_rms=" << cart.k4_abs_rms
            << " k4_reference_rms=" << cart.k4_reference_rms
            << " k1=" << cart.k1_relative_rms
            << " k2=" << cart.k2_relative_rms
            << " k3=" << cart.k3_relative_rms
            << " k4=" << cart.k4_relative_rms
            << " linf_max_abs_tail=" << cart.linf_max_abs_tail
            << " linf_max_relative_tail=" << cart.linf_max_relative_tail
            << " geometry_guard_rejected_cells="
            << cart.geometry_guard_rejected_cells
            << " geometry_guard_direct_fallback_pairs="
            << cart.geometry_guard_direct_fallback_pairs
            << " geometry_ratio_max=" << cart.geometry_ratio_max
            << " fallback_refined_pairs=" << cart.fallback_refined_pairs
            << " fallback_base_near_pairs=" << cart.fallback_base_near_pairs
            << " sampled_far_panel_pairs=" << cart.far_panel_pairs
            << " sampled_fallback_panel_pairs="
            << cart.fallback_panel_pairs
            << " valid_pairs=" << cart.production_direct_reference_valid_pairs
            << " invalid_pairs=" << cart.invalid_pairs
            << " nonproduction_fallback_pairs="
            << cart.nonproduction_fallback_pairs
            << " sampled_target_hash=" << cart.sampled_target_hash
            << " policy_hash=" << cart.policy_hash
            << " sampled_apply_seconds=" << cart.apply_seconds
            << " sampled_setup_seconds=" << cart.setup_seconds
            << " direct_reference_seconds=" << cart.direct_reference_seconds
            << " total_seconds=" << cart.total_seconds << '\n';
  require(cart.execution_status == "pass",
          "Cartesian p4 sampled comparator executes finitely");
  require(cart.status != "pass" || cart.accuracy_status == "pass",
          "Cartesian pass status implies accuracy pass");
  require(cart.accuracy_status == "pass" || cart.accuracy_status == "fail",
          "Cartesian comparator reports evaluated accuracy status");
  require(std::isfinite(cart.k4_abs_rms) &&
              std::isfinite(cart.k4_reference_rms),
          "Cartesian comparator reports finite K4 absolute/reference RMS");
  require(cart.sampled_target_hash == kifmm.sampled_target_hash,
          "Cartesian and KIFMM use identical sampled targets");
  require(cart.production_direct_reference_valid_pairs ==
              panels.size() * sample.size(),
          "Cartesian comparator uses complete dense sampled reference");
  require(cart.nonproduction_fallback_pairs == 0,
          "Cartesian comparator does not use nonproduction fallback");
  if (enforce_smoke_accuracy) {
    require(cart.aggregate_relative_rms <= 5.0e-6,
            "Cartesian p4 sampled rows remain within smoke tolerance");
  }
}

void run_e15_cartesian_comparator_smoke(const std::string &pqr_1sh1) {
  const std::vector<Panel> panels = build_real_surface(pqr_1sh1);
  const std::vector<double> unknowns = unknowns_for(panels);
  const auto sample = uniform_hash_sample(panels, 128);
  run_e15_cartesian_comparator_case("1SH1", "uniform_hash_128", panels,
                                    unknowns, sample, true);
  std::cout << "PASS titan_fmm_kifmm_sampled_panel_matvec_smoke E15\n";
}
#endif

} // namespace

int main(int argc, char **argv) {
  try {
    const std::string pqr_1sh1 =
        argc > 1 ? argv[1] : "test/pbe_inputs_public/1sh1_full_amber_ph7.pqr";
    const std::string pqr_1a63 =
        argc > 2 ? argv[2] : "test/pbe_inputs_public/1a63_full_amber_ph7.pqr";
    bool full = false;
    bool e9_profile = false;
    bool e14_far_floor_smoke = false;
    bool e15_cartesian_comparator_smoke = false;
    bool e15_profile = false;
    bool e17_profile = false;
    bool e18_profile = false;
    std::string e9_fixture = "1SH1";
    std::string e9_sample = "uniform_hash";
    int e9_order = 6;
    int kifmm_order = 6;
    double e9_theta = 0.6;
    int cartesian_order = 4;
    double cartesian_theta = 0.6;
    double cartesian_tail_tol = 0.002;
    double cartesian_geometry_theta = 0.5;
    std::string cartesian_cache_mode = "green-jet-compressed";
    std::string condition_policy = "qr-rank-only";
    std::string source_box_policy = "tight-bbox";
    double min_far_panel_fraction = 0.0;
    std::size_t requested_sample_count = 0;
    for (int i = 3; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--full") {
        full = true;
      } else if (arg == "--e14-far-floor-smoke") {
        e14_far_floor_smoke = true;
      } else if (arg == "--e15-cartesian-comparator-smoke" ||
                 arg == "--e16-cartesian-comparator-smoke") {
        e15_cartesian_comparator_smoke = true;
      } else if (arg == "--e15-profile" || arg == "--e16-profile") {
        e15_profile = true;
      } else if (arg == "--e17-profile") {
        e17_profile = true;
      } else if (arg == "--e18-profile") {
        e18_profile = true;
      } else if (arg == "--e9-profile" || arg == "--e10-factor-profile") {
        e9_profile = true;
      } else if (arg.rfind("--sample-count=", 0) == 0) {
        requested_sample_count =
            static_cast<std::size_t>(std::stoull(arg.substr(15)));
      } else if (arg.rfind("--fixture=", 0) == 0) {
        e9_fixture = arg.substr(10);
      } else if (arg.rfind("--sample=", 0) == 0) {
        e9_sample = arg.substr(9);
      } else if (arg.rfind("--order=", 0) == 0) {
        e9_order = std::stoi(arg.substr(8));
      } else if (arg.rfind("--kifmm-order=", 0) == 0) {
        kifmm_order = std::stoi(arg.substr(14));
      } else if (arg.rfind("--orders=", 0) == 0) {
        const std::string orders = arg.substr(9);
        const std::size_t comma = orders.find(',');
        e9_order = std::stoi(comma == std::string::npos
                                 ? orders
                                 : orders.substr(0, comma));
      } else if (arg.rfind("--theta=", 0) == 0) {
        e9_theta = std::stod(arg.substr(8));
      } else if (arg.rfind("--cartesian-order=", 0) == 0) {
        cartesian_order = std::stoi(arg.substr(18));
      } else if (arg.rfind("--cartesian-theta=", 0) == 0) {
        cartesian_theta = std::stod(arg.substr(18));
      } else if (arg.rfind("--cartesian-tail-tol=", 0) == 0) {
        cartesian_tail_tol = std::stod(arg.substr(21));
      } else if (arg.rfind("--cartesian-geometry-theta=", 0) == 0) {
        cartesian_geometry_theta = std::stod(arg.substr(27));
      } else if (arg.rfind("--cartesian-cache-mode=", 0) == 0) {
        cartesian_cache_mode = arg.substr(23);
      } else if (arg.rfind("--min-far-panel-fraction=", 0) == 0) {
        min_far_panel_fraction = std::stod(arg.substr(25));
      } else if (arg.rfind("--condition-policy=", 0) == 0) {
        condition_policy = arg.substr(19);
      } else if (arg.rfind("--source-box-policy=", 0) == 0) {
        source_box_policy = arg.substr(20);
      }
    }
    if (source_box_policy != "tight-bbox" &&
        source_box_policy != "level-cube") {
      throw std::runtime_error("invalid_source_box_policy");
    }
    const std::size_t sample_count =
        requested_sample_count > 0 ? requested_sample_count
                                   : (full ? 512 : 32);

    if (e14_far_floor_smoke) {
      run_e14_far_fraction_floor_smoke(pqr_1sh1);
      return EXIT_SUCCESS;
    }
    if (e15_cartesian_comparator_smoke) {
#ifdef TITAN_ENABLE_E15_CARTESIAN_COMPARATOR
      run_e15_cartesian_comparator_smoke(pqr_1sh1);
      return EXIT_SUCCESS;
#else
      throw std::runtime_error(
          "E15 Cartesian comparator smoke was built without "
          "TITAN_ENABLE_E15_CARTESIAN_COMPARATOR");
#endif
    }
    if (e15_profile || e17_profile || e18_profile) {
#ifdef TITAN_ENABLE_E15_CARTESIAN_COMPARATOR
      if (kifmm_order != 6 && kifmm_order != 8)
        throw std::runtime_error("E18 KIFMM order must be 6 or 8");
      std::string fixture_label = e9_fixture;
      std::transform(fixture_label.begin(), fixture_label.end(),
                     fixture_label.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::toupper(c));
                     });
      const std::string pqr_path =
          fixture_label == "1A63" ? pqr_1a63 : pqr_1sh1;
      const std::vector<Panel> panels = build_real_surface(pqr_path);
      const std::vector<double> unknowns = unknowns_for(panels);
      std::vector<std::size_t> sample;
      if (e9_sample == "uniform_hash") {
        sample = uniform_hash_sample(panels, sample_count);
      } else if (e9_sample == "inward_normal_groove") {
        sample = inward_normal_groove_sample(panels, sample_count);
      } else if (e9_sample == "closest_opposing_surface") {
        sample = closest_opposing_surface_sample(panels, sample_count);
      } else if (e9_sample == "high_aspect") {
        sample = high_aspect_sample(panels, sample_count);
      } else {
        throw std::runtime_error("unknown E15 sample set: " + e9_sample);
      }
      require(!sample.empty(), "E15 sample is non-empty");
      require(sample.size() <= sample_count,
              "E15 sample does not exceed request");
      CartesianSampledPanelMatvecConfig cartesian_override;
      cartesian_override.order = cartesian_order;
      cartesian_override.theta = cartesian_theta;
      cartesian_override.tail_tol = cartesian_tail_tol;
      cartesian_override.strict_geometry_theta = cartesian_geometry_theta;
      cartesian_override.cache_mode = cartesian_cache_mode;
      run_e15_cartesian_comparator_case(
          fixture_label, e9_sample + "_" + std::to_string(sample.size()),
          panels, unknowns, sample, false,
          (e17_profile || e18_profile) ? &cartesian_override : nullptr,
          kifmm_order);
      std::cout << "PASS titan_fmm_kifmm_sampled_panel_matvec_smoke "
                << (e18_profile   ? "E18_PROFILE"
                    : e17_profile ? "E17_PROFILE"
                                  : "E15_PROFILE")
                << '\n';
      return EXIT_SUCCESS;
#else
      throw std::runtime_error(
          "E15 Cartesian comparator profile was built without "
          "TITAN_ENABLE_E15_CARTESIAN_COMPARATOR");
#endif
    }

    if (e9_profile) {
      std::string fixture_label = e9_fixture;
      std::transform(fixture_label.begin(), fixture_label.end(),
                     fixture_label.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::toupper(c));
                     });
      const std::string pqr_path =
          fixture_label == "1A63" ? pqr_1a63 : pqr_1sh1;
      const std::vector<Panel> panels = build_real_surface(pqr_path);
      const std::vector<double> unknowns = unknowns_for(panels);
      std::vector<std::size_t> sample;
      if (e9_sample == "uniform_hash") {
        sample = uniform_hash_sample(panels, sample_count);
      } else if (e9_sample == "inward_normal_groove") {
        sample = inward_normal_groove_sample(panels, sample_count);
      } else if (e9_sample == "closest_opposing_surface") {
        sample = closest_opposing_surface_sample(panels, sample_count);
      } else if (e9_sample == "high_aspect") {
        sample = high_aspect_sample(panels, sample_count);
      } else {
        throw std::runtime_error("unknown E9 sample set: " + e9_sample);
      }
      require(!sample.empty(), "E9 sample is non-empty");
      require(sample.size() <= sample_count, "E9 sample does not exceed request");
      const KifmmSampledPanelMatvecResult result =
          run_case(panels, unknowns, sample, e9_order, e9_theta,
                   condition_policy, source_box_policy,
                   min_far_panel_fraction);
      print_result(fixture_label, e9_sample, result, e9_theta,
                   condition_policy);
      require_common_clean(fixture_label + "/" + e9_sample, result,
                           panels.size());
      require(result.lite.accepted_far_panel_pairs > 0,
              "E9 row has accepted far work");
      require(result.lite.direct_fallback_panel_pairs > 0,
              "E9 row has production fallback work");
      require(result.status == "pass" ||
                  result.status == "aggregate_rms_gate_failed" ||
                  result.status == "channel_rms_gate_failed",
              "E9 row has interpretable status");
      require(result.lite.max_equivalent_relative_residual <= 1.0e-10,
              "E9 equivalent residual gate");
      if (condition_policy == "svd-full-debug" ||
          condition_policy == "svd-values") {
        require(result.lite.max_condition_estimate <= 1.0e12,
                "E9 condition gate");
      }
      std::cout << "PASS titan_fmm_kifmm_sampled_panel_matvec_smoke E9\n";
      return EXIT_SUCCESS;
    }

    const std::vector<Panel> panels_1sh1 = build_real_surface(pqr_1sh1);
    const std::vector<double> x_1sh1 = unknowns_for(panels_1sh1);

    const auto uniform_1sh1 = uniform_hash_sample(panels_1sh1, sample_count);
    require(uniform_1sh1.size() == sample_count, "1SH1 uniform sample size");

    const KifmmSampledPanelMatvecResult u4 =
        run_case(panels_1sh1, x_1sh1, uniform_1sh1, 4, 0.6,
                 "qr-rank-only", "tight-bbox", min_far_panel_fraction);
    print_result("1SH1", "uniform_hash_q4_diagnostic", u4);
    require_common_clean("1SH1/uniform q4", u4, panels_1sh1.size());

    const KifmmSampledPanelMatvecResult u6 =
        run_case(panels_1sh1, x_1sh1, uniform_1sh1, 6, 0.6,
                 "qr-rank-only", "tight-bbox", min_far_panel_fraction);
    require_diagnostic_clean("1SH1", "uniform_hash_q6", u6,
                             panels_1sh1.size());
    const KifmmSampledPanelMatvecResult u6_repeat =
        run_case(panels_1sh1, x_1sh1, uniform_1sh1, 6, 0.6,
                 "qr-rank-only", "tight-bbox", min_far_panel_fraction);
    require_repeat_hash("1SH1 uniform q6", u6, u6_repeat);

    const KifmmSampledPanelMatvecResult u8 =
        run_case(panels_1sh1, x_1sh1, uniform_1sh1, 8, 0.6,
                 "qr-rank-only", "tight-bbox", min_far_panel_fraction);
    require_candidate("1SH1", "uniform_hash_q8_confirm", u8,
                      panels_1sh1.size());
    require(passes_candidate_gate(u6) ||
                u8.lite.aggregate_relative_rms < u6.lite.aggregate_relative_rms,
            "1SH1 q8 improves any q6 diagnostic gap");

    double groove_rms = -1.0;
    double opposing_rms = -1.0;
    double aspect_rms = -1.0;
    double aspect_q8_rms = -1.0;
    if (full) {
      const auto groove_1sh1 =
          inward_normal_groove_sample(panels_1sh1, sample_count);
      const auto opposing_1sh1 =
          closest_opposing_surface_sample(panels_1sh1, sample_count);
      require(groove_1sh1.size() == sample_count, "1SH1 groove sample size");
      require(opposing_1sh1.size() == sample_count,
              "1SH1 opposing sample size");
      const KifmmSampledPanelMatvecResult g6 =
          run_case(panels_1sh1, x_1sh1, groove_1sh1, 6, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_diagnostic_clean("1SH1", "inward_normal_groove_q6", g6,
                               panels_1sh1.size());
      const KifmmSampledPanelMatvecResult o6 =
          run_case(panels_1sh1, x_1sh1, opposing_1sh1, 6, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_diagnostic_clean("1SH1", "closest_opposing_surface_q6", o6,
                               panels_1sh1.size());
      const KifmmSampledPanelMatvecResult g8 =
          run_case(panels_1sh1, x_1sh1, groove_1sh1, 8, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_candidate("1SH1", "inward_normal_groove_q8_confirm", g8,
                        panels_1sh1.size());
      const KifmmSampledPanelMatvecResult o8 =
          run_case(panels_1sh1, x_1sh1, opposing_1sh1, 8, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_candidate("1SH1", "closest_opposing_surface_q8_confirm", o8,
                        panels_1sh1.size());
      groove_rms = g6.lite.aggregate_relative_rms;
      opposing_rms = o6.lite.aggregate_relative_rms;

      const std::vector<Panel> panels_1a63 = build_real_surface(pqr_1a63);
      const std::vector<double> x_1a63 = unknowns_for(panels_1a63);
      const auto aspect_1a63 = high_aspect_sample(panels_1a63, 128);
      require(!aspect_1a63.empty(), "1A63 high-aspect sample is non-empty");
      const KifmmSampledPanelMatvecResult a6 =
          run_case(panels_1a63, x_1a63, aspect_1a63, 6, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_diagnostic_clean("1A63", "high_aspect_q6", a6,
                               panels_1a63.size());
      const KifmmSampledPanelMatvecResult a8 =
          run_case(panels_1a63, x_1a63, aspect_1a63, 8, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require_candidate("1A63", "high_aspect_q8_confirm", a8,
                        panels_1a63.size());
      aspect_rms = a6.lite.aggregate_relative_rms;
      aspect_q8_rms = a8.lite.aggregate_relative_rms;
    }

    std::vector<std::size_t> duplicate = uniform_1sh1;
    if (!duplicate.empty()) {
      duplicate.back() = duplicate.front();
      const KifmmSampledPanelMatvecResult dup =
          run_case(panels_1sh1, x_1sh1, duplicate, 6, 0.6,
                   "qr-rank-only", "tight-bbox", min_far_panel_fraction);
      require(!dup.finite, "duplicate sample indices are rejected");
      require(dup.status == "invalid_target_indices",
              "duplicate sample status");
    }

    std::vector<double> bad_unknowns = x_1sh1;
    bad_unknowns[0] = std::numeric_limits<double>::quiet_NaN();
    const KifmmSampledPanelMatvecResult bad =
        run_kifmm_lite_sampled_panel_matvec(
            panels_1sh1, bad_unknowns, uniform_1sh1,
            base_config(6, 0.6, min_far_panel_fraction));
    require(!bad.finite, "non-finite unknowns are rejected");
    require(bad.status == "invalid_unknowns", "invalid unknown status");

    std::cout << std::setprecision(12)
              << "KIFMM_SAMPLED_PANEL_SUMMARY mode="
              << (full ? "full" : "quick")
              << " rows=" << sample_count
              << " 1sh1_uniform_rms=" << u6.lite.aggregate_relative_rms
              << " 1sh1_groove_rms=" << groove_rms
              << " 1sh1_opposing_rms=" << opposing_rms
              << " 1a63_aspect_rms=" << aspect_rms
              << " 1sh1_uniform_q8_rms=" << u8.lite.aggregate_relative_rms
              << " 1a63_aspect_q8_rms=" << aspect_q8_rms
              << " 1sh1_hash=" << u6.lite.worklist_hash << '\n';
    std::cout << "PASS titan_fmm_kifmm_sampled_panel_matvec_smoke\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
