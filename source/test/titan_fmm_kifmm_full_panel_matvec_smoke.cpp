#include "tabi/tabi_fmm_kifmm_full_panel_matvec.h"
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
#include <vector>

using TitanPBE::Tabi::Atom;
using TitanPBE::Tabi::Options;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::Result;
using titan::tabi::fmm_kifmm::apply_prepared_kifmm_lite_full_panel_operator;
using titan::tabi::fmm_kifmm::assemble_kifmm_lite_tabi_system_matvec;
using titan::tabi::fmm_kifmm::KifmmFullPanelPreparedOperator;
using titan::tabi::fmm_kifmm::KifmmFullPanelMatvecResult;
using titan::tabi::fmm_kifmm::KifmmJcpRow;
using titan::tabi::fmm_kifmm::KifmmRealPanelCloudConfig;
using titan::tabi::fmm_kifmm::prepare_kifmm_lite_full_panel_operator;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::run_kifmm_lite_full_panel_matvec;

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

std::vector<double> e20_vector_zero(const std::vector<Panel> &panels) {
  return std::vector<double>(2 * panels.size(), 0.0);
}

std::vector<double> e20_vector_sparse_unit(const std::vector<Panel> &panels) {
  std::vector<double> out(2 * panels.size(), 0.0);
  if (!panels.empty()) {
    out[panels.size() / 3] = 1.0;
    out[panels.size() + (2 * panels.size()) / 3] = -0.5;
  }
  return out;
}

std::vector<double> e20_vector_rademacher(const std::vector<Panel> &panels,
                                          std::uint64_t seed) {
  std::vector<double> out(2 * panels.size(), 0.0);
  std::uint64_t state = seed == 0 ? 1 : seed;
  for (double &value : out) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    value = (state & 1ull) ? 0.25 : -0.25;
  }
  return out;
}

std::vector<double> e20_linear_combination(const std::vector<double> &x,
                                           const std::vector<double> &y,
                                           double a, double b) {
  require(x.size() == y.size(), "linear combination size mismatch");
  std::vector<double> out(x.size(), 0.0);
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = a * x[i] + b * y[i];
  }
  return out;
}

bool finite_rows(const std::vector<KifmmJcpRow> &rows) {
  return std::all_of(rows.begin(), rows.end(), [](const KifmmJcpRow &row) {
    return std::isfinite(row.phi_row) && std::isfinite(row.phin_row);
  });
}

double e20_linearity_relative_error(const std::vector<KifmmJcpRow> &lhs,
                                    const std::vector<KifmmJcpRow> &ax,
                                    const std::vector<KifmmJcpRow> &ay,
                                    double a, double b) {
  require(lhs.size() == ax.size() && lhs.size() == ay.size(),
          "linearity row size mismatch");
  double num = 0.0;
  double den = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double ref_phi = a * ax[i].phi_row + b * ay[i].phi_row;
    const double ref_phin = a * ax[i].phin_row + b * ay[i].phin_row;
    const double d_phi = lhs[i].phi_row - ref_phi;
    const double d_phin = lhs[i].phin_row - ref_phin;
    num += d_phi * d_phi + d_phin * d_phin;
    den += ref_phi * ref_phi + ref_phin * ref_phin;
  }
  return std::sqrt(num / std::max(den, 1.0e-300));
}

KifmmRealPanelCloudConfig base_config(int order) {
  KifmmRealPanelCloudConfig config;
  config.lite.order = order;
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
      const double score = d / std::max(0.05, 0.5 + 0.5 * opposing + 0.25 * facing);
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
  out.reserve(std::min(count, scored.size()));
  for (const auto &entry : scored) {
    if (out.size() == count) {
      break;
    }
    out.push_back(entry.second);
  }
  std::sort(out.begin(), out.end());
  return out;
}

double max_channel_rms(const KifmmFullPanelMatvecResult &result) {
  return std::max({result.k1_relative_rms, result.k2_relative_rms,
                   result.k3_relative_rms, result.k4_relative_rms});
}

std::string max_channel_name(const KifmmFullPanelMatvecResult &result) {
  const double values[4] = {result.k1_relative_rms, result.k2_relative_rms,
                            result.k3_relative_rms, result.k4_relative_rms};
  const char *names[4] = {"K1", "K2", "K3", "K4"};
  int best = 0;
  for (int i = 1; i < 4; ++i) {
    if (values[i] > values[best]) {
      best = i;
    }
  }
  return names[best];
}

void print_result(const std::string &fixture, const std::string &sample,
                  const KifmmFullPanelMatvecResult &result,
                  const KifmmFullPanelMatvecResult &repeat,
                  const KifmmRealPanelCloudConfig &config) {
  const std::size_t total_pairs =
      result.lite.accepted_far_panel_pairs +
      result.lite.direct_fallback_panel_pairs;
  const double far_fraction =
      total_pairs == 0
          ? 0.0
          : static_cast<double>(result.lite.accepted_far_panel_pairs) /
                static_cast<double>(total_pairs);
  std::cout << std::setprecision(12)
            << "KIFMM_FULL_PANEL_CASE fixture=" << fixture
            << " sample=" << sample << " status=" << result.status
            << " source_panels=" << result.lite.source_panels
            << " target_panels=" << result.lite.target_panels
            << " sampled_rows=" << result.sampled_target_indices.size()
            << " order=" << result.order
            << " theta=" << config.lite.opening_theta
            << " full_apply_policy=" << result.kifmm_full_apply_policy
            << " condition_policy=" << config.lite.condition_policy
            << " source_box_policy=" << config.lite.source_box_policy
            << " eps_ratio=" << config.lite.eps_ratio
            << " kappa=" << config.lite.kappa
            << " near_integration=" << config.near_integration
            << " quad_buffer=" << config.quad_buffer
            << " full_pair_coverage=" << result.full_pair_coverage
            << " production_fallback_valid_pairs="
            << result.production_fallback_valid_pairs
            << " sampled_direct_valid_pairs="
            << result.sampled_direct_reference_valid_pairs
            << " sampled_invalid_pairs=" << result.sampled_invalid_pairs
            << " nonproduction_fallback_pairs="
            << result.nonproduction_fallback_pairs
            << " far_panel_pairs=" << result.lite.accepted_far_panel_pairs
            << " fallback_panel_pairs="
            << result.lite.direct_fallback_panel_pairs
            << " far_fraction=" << far_fraction
            << " aggregate_rms=" << result.aggregate_relative_rms
            << " max_channel_rms=" << max_channel_rms(result)
            << " max_channel_name=" << max_channel_name(result)
            << " k1=" << result.k1_relative_rms
            << " k2=" << result.k2_relative_rms
            << " k3=" << result.k3_relative_rms
            << " k4=" << result.k4_relative_rms
            << " max_equiv_residual="
            << result.lite.max_equivalent_relative_residual
            << " max_condition=" << result.lite.max_condition_estimate
            << " cache_entries=" << result.equivalent_system_cache_entries
            << " cache_hit_rate=" << result.equivalent_cache_hit_rate
            << " unique_far_source_nodes=" << result.unique_far_source_nodes
            << " unique_half_widths=" << result.unique_equivalent_half_widths
            << " unique_kappa_h_values="
            << result.unique_equivalent_kappa_h_values
            << " tree_build_seconds=" << result.tree_build_seconds
            << " worklist_build_seconds=" << result.worklist_build_seconds
            << " equivalent_build_solve_seconds="
            << result.equivalent_build_solve_seconds
            << " kifmm_channel_build_seconds="
            << result.kifmm_channel_build_seconds
            << " full_far_eval_seconds=" << result.full_far_eval_seconds
            << " full_fallback_eval_seconds="
            << result.full_fallback_eval_seconds
            << " kifmm_far_eval_seconds=" << result.kifmm_far_eval_seconds
            << " kifmm_fallback_eval_seconds="
            << result.kifmm_fallback_eval_seconds
            << " kifmm_target_tasks=" << result.kifmm_target_tasks
            << " kifmm_parallel_workers=" << result.kifmm_parallel_workers
            << " kifmm_parallel_speedup_vs_serial="
            << result.kifmm_parallel_speedup_vs_serial
            << " sampled_dense_audit_seconds="
            << result.sampled_dense_audit_seconds
            << " setup_seconds=" << result.setup_seconds
            << " apply_seconds=" << result.apply_seconds
            << " first_apply_seconds=" << result.first_apply_seconds
            << " repeat_apply_seconds=" << repeat.apply_seconds
            << " repeat_apply_min_seconds="
            << std::min(result.apply_seconds, repeat.apply_seconds)
            << " repeat_apply_max_seconds="
            << std::max(result.apply_seconds, repeat.apply_seconds)
            << " total_seconds=" << result.total_seconds
            << " panel_geometry_hash=" << result.panel_geometry_hash
            << " unknown_vector_hash=" << result.unknown_vector_hash
            << " sampled_target_hash=" << result.sampled_target_hash
            << " full_output_hash_phi=" << result.full_output_hash_phi
            << " full_output_hash_phin=" << result.full_output_hash_phin
            << " full_channel_output_hash=" << result.full_channel_output_hash
            << " worklist_hash=" << result.lite.worklist_hash
            << " policy_hash=" << result.policy_hash << '\n';
}

void require_candidate(const std::string &name,
                       const KifmmFullPanelMatvecResult &result,
                       std::size_t source_count) {
  require(result.finite, name + " finite");
  require(result.status == "pass", name + " status pass");
  require(result.lite.source_panels == source_count, name + " source count");
  require(result.lite.target_panels == source_count, name + " target count");
  require(result.full_pair_coverage == source_count * source_count,
          name + " full pair coverage exact");
  require(result.sampled_direct_reference_valid_pairs ==
              result.sampled_target_indices.size() * source_count,
          name + " sampled dense direct coverage exact");
  require(result.sampled_invalid_pairs == 0, name + " sampled invalid zero");
  require(result.nonproduction_fallback_pairs == 0,
          name + " nonproduction fallback zero");
  require(result.lite.accepted_far_panel_pairs > 0,
          name + " accepted far work exists");
  require(result.lite.direct_fallback_panel_pairs > 0,
          name + " production fallback work exists");
  require(result.production_fallback_valid_pairs ==
              result.lite.direct_fallback_panel_pairs,
          name + " fallback coverage exact");
  require(result.aggregate_relative_rms <= 1.0e-6,
          name + " aggregate gate");
  require(max_channel_rms(result) <= 5.0e-6, name + " channel gate");
  require(result.lite.max_equivalent_relative_residual <= 1.0e-10,
          name + " equivalent residual gate");
  require(result.panel_geometry_hash != 0, name + " panel hash");
  require(result.unknown_vector_hash != 0, name + " unknown hash");
  require(result.sampled_target_hash != 0, name + " sample hash");
  require(result.full_output_hash_phi != 0, name + " phi output hash");
  require(result.full_output_hash_phin != 0, name + " phin output hash");
  require(result.full_channel_output_hash != 0, name + " channel output hash");
  require(result.lite.worklist_hash != 0, name + " worklist hash");
  require(result.policy_hash != 0, name + " policy hash");
  require(result.setup_seconds >= 0.0, name + " setup seconds finite");
  require(result.apply_seconds > 0.0, name + " apply seconds positive");
  require(result.first_apply_seconds > 0.0, name + " first apply seconds positive");
}

int run_profile(const std::string &pqr_1sh1, const std::string &pqr_1a63,
                const std::string &fixture, const std::string &sample_name,
                std::size_t sample_count, int order) {
  const std::string &pqr_path = fixture == "1A63" ? pqr_1a63 : pqr_1sh1;
  const std::vector<Panel> panels = build_real_surface(pqr_path);
  const std::vector<double> unknowns = unknowns_for(panels);
  std::vector<std::size_t> sample;
  if (sample_name == "uniform_hash") {
    sample = uniform_hash_sample(panels, sample_count);
  } else if (sample_name == "high_aspect") {
    sample = high_aspect_sample(panels, sample_count);
  } else if (sample_name == "closest_opposing_surface") {
    sample = closest_opposing_surface_sample(panels, sample_count);
  } else if (sample_name == "inward_normal_groove") {
    sample = inward_normal_groove_sample(panels, sample_count);
  } else {
    throw std::runtime_error("unknown E19 sample set: " + sample_name);
  }
  KifmmRealPanelCloudConfig config = base_config(order);
  const KifmmFullPanelMatvecResult result =
      run_kifmm_lite_full_panel_matvec(panels, unknowns, sample, config);
  require_candidate(fixture + "/" + sample_name, result, panels.size());

  const KifmmFullPanelMatvecResult repeat =
      run_kifmm_lite_full_panel_matvec(panels, unknowns, sample, config);
  print_result(fixture, sample_name + "_" + std::to_string(sample.size()),
               result, repeat, config);
  require(result.lite.worklist_hash == repeat.lite.worklist_hash,
          "repeat worklist hash stable");
  require(result.policy_hash == repeat.policy_hash, "repeat policy hash stable");
  require(result.panel_geometry_hash == repeat.panel_geometry_hash,
          "repeat panel hash stable");
  require(result.unknown_vector_hash == repeat.unknown_vector_hash,
          "repeat unknown hash stable");
  require(result.sampled_target_hash == repeat.sampled_target_hash,
          "repeat sampled hash stable");
  require(result.full_output_hash_phi == repeat.full_output_hash_phi,
          "repeat phi output hash stable");
  require(result.full_output_hash_phin == repeat.full_output_hash_phin,
          "repeat phin output hash stable");
  require(result.full_channel_output_hash == repeat.full_channel_output_hash,
          "repeat channel output hash stable");

  std::vector<std::size_t> duplicate = sample;
  if (!duplicate.empty()) {
    duplicate.back() = duplicate.front();
    const KifmmFullPanelMatvecResult dup =
        run_kifmm_lite_full_panel_matvec(panels, unknowns, duplicate, config);
    require(!dup.finite, "duplicate sampled targets rejected");
    require(dup.status == "invalid_target_indices",
            "duplicate sample status");
  }
  std::vector<std::size_t> out_of_range = sample;
  out_of_range.push_back(panels.size());
  const KifmmFullPanelMatvecResult oor =
      run_kifmm_lite_full_panel_matvec(panels, unknowns, out_of_range, config);
  require(!oor.finite, "out-of-range sampled targets rejected");
  require(oor.status == "invalid_target_indices",
          "out-of-range sample status");
  std::cout << "PASS titan_fmm_kifmm_full_panel_matvec_smoke E19_PROFILE\n";
  return EXIT_SUCCESS;
}

double far_fraction_for(const KifmmFullPanelMatvecResult &result) {
  const std::size_t total_pairs =
      result.lite.accepted_far_panel_pairs +
      result.lite.direct_fallback_panel_pairs;
  return total_pairs == 0
             ? 0.0
             : static_cast<double>(result.lite.accepted_far_panel_pairs) /
                   static_cast<double>(total_pairs);
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
  throw std::runtime_error("unknown E20 sample set: " + sample_name);
}

int run_e20_profile(const std::string &pqr_1sh1,
                    const std::string &pqr_1a63,
                    const std::string &fixture,
                    const std::string &sample_name,
                    std::size_t sample_count, int order) {
  const std::string &pqr_path = fixture == "1A63" ? pqr_1a63 : pqr_1sh1;
  const std::vector<Panel> panels = build_real_surface(pqr_path);
  const std::vector<std::size_t> sample =
      sample_for_name(panels, sample_name, sample_count);
  KifmmRealPanelCloudConfig config = base_config(order);
  // E20 probes linearity of the accelerated operator itself. Dense-direct audit
  // metrics are still computed, but random vectors are allowed to exceed the
  // E19 solution-like audit gate without aborting the diagnostic.
  config.lite.max_aggregate_relative_rms = 1.0e100;
  config.lite.max_channel_relative_rms = 1.0e100;

  KifmmFullPanelPreparedOperator prepared =
      prepare_kifmm_lite_full_panel_operator(panels, sample, config);
  const std::vector<double> x = unknowns_for(panels);
  const std::vector<double> y = e20_vector_rademacher(panels, 0xE20A63ull);
  const std::vector<double> z = e20_linear_combination(x, y, 0.37, -1.21);
  const std::vector<double> sparse = e20_vector_sparse_unit(panels);
  const std::vector<double> zero = e20_vector_zero(panels);

  const KifmmFullPanelMatvecResult ax =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, x);
  const KifmmFullPanelMatvecResult ay =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, y);
  const KifmmFullPanelMatvecResult az =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, z);
  const KifmmFullPanelMatvecResult ax_repeat =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, x);
  const KifmmFullPanelMatvecResult asparse =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, sparse);
  const KifmmFullPanelMatvecResult azero =
      apply_prepared_kifmm_lite_full_panel_operator(prepared, zero);

  require(finite_rows(ax.lite.rows), "x finite rows");
  require(finite_rows(ay.lite.rows), "y finite rows");
  require(finite_rows(az.lite.rows), "z finite rows");
  require(finite_rows(ax_repeat.lite.rows), "x repeat finite rows");
  require(finite_rows(asparse.lite.rows), "sparse finite rows");
  require(finite_rows(azero.lite.rows), "zero finite rows");

  const double linearity_error =
      e20_linearity_relative_error(az.lite.rows, ax.lite.rows, ay.lite.rows,
                                   0.37, -1.21);
  const bool repeat_output_hash_stable =
      ax.full_output_hash_phi == ax_repeat.full_output_hash_phi &&
      ax.full_output_hash_phin == ax_repeat.full_output_hash_phin &&
      ax.full_channel_output_hash == ax_repeat.full_channel_output_hash;
  const bool repeat_worklist_hash_stable =
      ax.lite.worklist_hash == ay.lite.worklist_hash &&
      ax.lite.worklist_hash == az.lite.worklist_hash &&
      ax.lite.worklist_hash == ax_repeat.lite.worklist_hash;
  const bool repeat_policy_hash_stable =
      ax.policy_hash == ay.policy_hash && ax.policy_hash == az.policy_hash &&
      ax.policy_hash == ax_repeat.policy_hash;
  const bool coverage_ok =
      az.full_pair_coverage == panels.size() * panels.size() &&
      az.sampled_direct_reference_valid_pairs == sample.size() * panels.size();
  const bool prepared_context_reused =
      prepared.frozen_worklist && prepared.frozen_equivalent_system_cache &&
      prepared.prepared_equivalent_system_cache_entries > 0 &&
      prepared.prepared_far_source_nodes > 0 && az.tree_build_seconds == 0.0 &&
      az.worklist_build_seconds == 0.0 &&
      az.equivalent_factor_cache_misses == 0;
  const bool e20_pass = linearity_error <= 1.0e-10 &&
                        repeat_output_hash_stable &&
                        repeat_worklist_hash_stable &&
                        repeat_policy_hash_stable && coverage_ok &&
                        az.sampled_invalid_pairs == 0 &&
                        az.nonproduction_fallback_pairs == 0 &&
                        far_fraction_for(az) >= 0.10 &&
                        prepared_context_reused;

  std::cout << std::setprecision(12)
            << "KIFMM_E20_CASE fixture=" << fixture
            << " sample=" << sample_name
            << " sampled_rows=" << sample.size()
            << " order=" << order
            << " status=" << (e20_pass ? "pass" : "fail")
            << " operator_status_x=" << ax.status
            << " operator_status_y=" << ay.status
            << " operator_status_z=" << az.status
            << " linearity_relative_error=" << linearity_error
            << " prepare_build_count=" << prepared.prepare_build_count
            << " apply_count=" << prepared.apply_count
            << " prepared_context_reused="
            << (prepared_context_reused ? 1 : 0)
            << " frozen_worklist=" << (prepared.frozen_worklist ? 1 : 0)
            << " frozen_equivalent_system_cache="
            << (prepared.frozen_equivalent_system_cache ? 1 : 0)
            << " prepared_worklist_pairs=" << prepared.prepared_worklist_pairs
            << " prepared_far_source_nodes="
            << prepared.prepared_far_source_nodes
            << " prepared_equivalent_system_cache_entries="
            << prepared.prepared_equivalent_system_cache_entries
            << " repeat_output_hash_stable="
            << (repeat_output_hash_stable ? 1 : 0)
            << " repeat_worklist_hash_stable="
            << (repeat_worklist_hash_stable ? 1 : 0)
            << " repeat_policy_hash_stable="
            << (repeat_policy_hash_stable ? 1 : 0)
            << " full_pair_coverage=" << az.full_pair_coverage
            << " sampled_direct_valid_pairs="
            << az.sampled_direct_reference_valid_pairs
            << " sampled_invalid_pairs=" << az.sampled_invalid_pairs
            << " nonproduction_fallback_pairs="
            << az.nonproduction_fallback_pairs
            << " far_fraction=" << far_fraction_for(az)
            << " tree_build_seconds=" << az.tree_build_seconds
            << " worklist_build_seconds=" << az.worklist_build_seconds
            << " equivalent_factor_cache_misses="
            << az.equivalent_factor_cache_misses
            << " equivalent_factor_cache_hits="
            << az.equivalent_factor_cache_hits
            << " setup_seconds=" << az.setup_seconds
            << " apply_seconds=" << az.apply_seconds
            << " total_seconds=" << az.total_seconds
            << " worklist_hash=" << ax.lite.worklist_hash
            << " policy_hash=" << ax.policy_hash
            << " panel_geometry_hash=" << ax.panel_geometry_hash
            << " sampled_target_hash=" << ax.sampled_target_hash
            << " x_output_hash_phi=" << ax.full_output_hash_phi
            << " y_output_hash_phi=" << ay.full_output_hash_phi
            << " z_output_hash_phi=" << az.full_output_hash_phi
            << " x_output_hash_phin=" << ax.full_output_hash_phin
            << " y_output_hash_phin=" << ay.full_output_hash_phin
            << " z_output_hash_phin=" << az.full_output_hash_phin
            << " x_channel_hash=" << ax.full_channel_output_hash
            << " y_channel_hash=" << ay.full_channel_output_hash
            << " z_channel_hash=" << az.full_channel_output_hash
            << " sparse_channel_hash=" << asparse.full_channel_output_hash
            << " zero_channel_hash=" << azero.full_channel_output_hash
            << '\n';
  return e20_pass ? EXIT_SUCCESS : EXIT_FAILURE;
}

double full_row_relative_difference(const KifmmFullPanelMatvecResult &a,
                                    const KifmmFullPanelMatvecResult &b) {
  require(a.lite.rows.size() == b.lite.rows.size(),
          "full row relative difference size mismatch");
  double num = 0.0;
  double den = 0.0;
  for (std::size_t i = 0; i < a.lite.rows.size(); ++i) {
    const double d_phi = a.lite.rows[i].phi_row - b.lite.rows[i].phi_row;
    const double d_phin = a.lite.rows[i].phin_row - b.lite.rows[i].phin_row;
    num += d_phi * d_phi + d_phin * d_phin;
    den += a.lite.rows[i].phi_row * a.lite.rows[i].phi_row +
           a.lite.rows[i].phin_row * a.lite.rows[i].phin_row;
  }
  return std::sqrt(num / std::max(den, 1.0e-300));
}

double full_channel_relative_difference(const KifmmFullPanelMatvecResult &a,
                                        const KifmmFullPanelMatvecResult &b) {
  require(a.lite.channel_rows.size() == b.lite.channel_rows.size(),
          "full channel relative difference size mismatch");
  double num = 0.0;
  double den = 0.0;
  for (std::size_t i = 0; i < a.lite.channel_rows.size(); ++i) {
    const double dk1 = a.lite.channel_rows[i].k1 - b.lite.channel_rows[i].k1;
    const double dk2 = a.lite.channel_rows[i].k2 - b.lite.channel_rows[i].k2;
    const double dk3 = a.lite.channel_rows[i].k3 - b.lite.channel_rows[i].k3;
    const double dk4 = a.lite.channel_rows[i].k4 - b.lite.channel_rows[i].k4;
    num += dk1 * dk1 + dk2 * dk2 + dk3 * dk3 + dk4 * dk4;
    den += a.lite.channel_rows[i].k1 * a.lite.channel_rows[i].k1 +
           a.lite.channel_rows[i].k2 * a.lite.channel_rows[i].k2 +
           a.lite.channel_rows[i].k3 * a.lite.channel_rows[i].k3 +
           a.lite.channel_rows[i].k4 * a.lite.channel_rows[i].k4;
  }
  return std::sqrt(num / std::max(den, 1.0e-300));
}

int run_e22_profile(const std::string &pqr_1sh1,
                    const std::string &pqr_1a63,
                    const std::string &fixture,
                    const std::string &sample_name,
                    std::size_t sample_count, int order,
                    std::size_t workers) {
  const std::string &pqr_path = fixture == "1A63" ? pqr_1a63 : pqr_1sh1;
  const std::vector<Panel> panels = build_real_surface(pqr_path);
  const std::vector<double> unknowns = unknowns_for(panels);
  const std::vector<std::size_t> sample =
      sample_for_name(panels, sample_name, sample_count);

  KifmmRealPanelCloudConfig serial_config = base_config(order);
  serial_config.full_apply_policy = "serial";
  KifmmFullPanelPreparedOperator serial_prepared =
      prepare_kifmm_lite_full_panel_operator(panels, sample, serial_config);
  KifmmFullPanelMatvecResult serial =
      apply_prepared_kifmm_lite_full_panel_operator(serial_prepared, unknowns,
                                                    true);
  require_candidate(fixture + "/" + sample_name + "/serial", serial,
                    panels.size());

  KifmmRealPanelCloudConfig parallel_config = base_config(order);
  parallel_config.full_apply_policy = "target-node-parallel";
  KifmmFullPanelPreparedOperator parallel_prepared =
      prepare_kifmm_lite_full_panel_operator(panels, sample, parallel_config);
  TitanExecutor executor(std::max<std::size_t>(1, workers));
  KifmmFullPanelMatvecResult parallel =
      apply_prepared_kifmm_lite_full_panel_operator(
          parallel_prepared, unknowns, &executor, true);
  require_candidate(fixture + "/" + sample_name + "/parallel", parallel,
                    panels.size());

  const double row_diff = full_row_relative_difference(serial, parallel);
  const double channel_diff =
      full_channel_relative_difference(serial, parallel);
  const bool hash_stable =
      serial.panel_geometry_hash == parallel.panel_geometry_hash &&
      serial.sampled_target_hash == parallel.sampled_target_hash &&
      serial.lite.worklist_hash == parallel.lite.worklist_hash &&
      serial.full_output_hash_phi == parallel.full_output_hash_phi &&
      serial.full_output_hash_phin == parallel.full_output_hash_phin &&
      serial.full_channel_output_hash == parallel.full_channel_output_hash;
  const double speedup =
      parallel.apply_seconds > 0.0 ? serial.apply_seconds / parallel.apply_seconds
                                   : 0.0;
  parallel.kifmm_parallel_speedup_vs_serial = speedup;

  print_result(fixture, sample_name + "_" + std::to_string(sample.size()) +
                            "_serial",
               serial, serial, serial_config);
  print_result(fixture, sample_name + "_" + std::to_string(sample.size()) +
                            "_target_node_parallel",
               parallel, parallel, parallel_config);

  std::cout << std::setprecision(12)
            << "KIFMM_E22_PARITY fixture=" << fixture
            << " sample=" << sample_name
            << " sampled_rows=" << sample.size()
            << " order=" << order
            << " workers=" << workers
            << " serial_status=" << serial.status
            << " parallel_status=" << parallel.status
            << " serial_apply_seconds=" << serial.apply_seconds
            << " parallel_apply_seconds=" << parallel.apply_seconds
            << " speedup=" << speedup
            << " row_relative_diff=" << row_diff
            << " channel_relative_diff=" << channel_diff
            << " hash_stable=" << (hash_stable ? 1 : 0)
            << " serial_worklist_hash=" << serial.lite.worklist_hash
            << " parallel_worklist_hash=" << parallel.lite.worklist_hash
            << " serial_policy_hash=" << serial.policy_hash
            << " parallel_policy_hash=" << parallel.policy_hash
            << " serial_far_panel_pairs="
            << serial.lite.accepted_far_panel_pairs
            << " parallel_far_panel_pairs="
            << parallel.lite.accepted_far_panel_pairs
            << " serial_fallback_panel_pairs="
            << serial.lite.direct_fallback_panel_pairs
            << " parallel_fallback_panel_pairs="
            << parallel.lite.direct_fallback_panel_pairs
            << '\n';
  require(row_diff <= 1.0e-12, "serial/parallel row parity");
  require(channel_diff <= 1.0e-12, "serial/parallel channel parity");
  require(hash_stable, "serial/parallel output hashes stable");
  std::cout << "PASS titan_fmm_kifmm_full_panel_matvec_smoke E22_PROFILE\n";
  return EXIT_SUCCESS;
}

double full_vector_relative_difference(const std::vector<double> &a,
                                       const std::vector<double> &b) {
  require(a.size() == b.size(), "full vector relative difference size mismatch");
  double num = 0.0;
  double den = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double d = a[i] - b[i];
    num += d * d;
    den += a[i] * a[i];
  }
  return std::sqrt(num / std::max(den, 1.0e-300));
}

int run_e23_profile(const std::string &pqr_1sh1,
                    const std::string &pqr_1a63,
                    const std::string &fixture,
                    const std::string &sample_name,
                    std::size_t sample_count, int order,
                    std::size_t workers) {
  const std::string &pqr_path = fixture == "1A63" ? pqr_1a63 : pqr_1sh1;
  const std::vector<Panel> panels = build_real_surface(pqr_path);
  const std::vector<double> unknowns = unknowns_for(panels);
  const std::vector<std::size_t> sample =
      sample_for_name(panels, sample_name, sample_count);

  KifmmRealPanelCloudConfig serial_config = base_config(order);
  serial_config.full_apply_policy = "serial";
  KifmmFullPanelPreparedOperator serial_prepared =
      prepare_kifmm_lite_full_panel_operator(panels, sample, serial_config);
  const KifmmFullPanelMatvecResult serial =
      apply_prepared_kifmm_lite_full_panel_operator(serial_prepared, unknowns,
                                                    true);
  require_candidate(fixture + "/" + sample_name + "/e23_serial", serial,
                    panels.size());

  KifmmRealPanelCloudConfig parallel_config = base_config(order);
  parallel_config.full_apply_policy = "target-node-parallel";
  KifmmFullPanelPreparedOperator parallel_prepared =
      prepare_kifmm_lite_full_panel_operator(panels, sample, parallel_config);
  TitanExecutor executor(std::max<std::size_t>(1, workers));
  const KifmmFullPanelMatvecResult parallel =
      apply_prepared_kifmm_lite_full_panel_operator(
          parallel_prepared, unknowns, &executor, true);
  require_candidate(fixture + "/" + sample_name + "/e23_parallel", parallel,
                    panels.size());

  std::vector<double> serial_matvec;
  std::vector<double> parallel_matvec;
  assemble_kifmm_lite_tabi_system_matvec(serial, unknowns, 0.5, 0.5,
                                         serial_matvec);
  assemble_kifmm_lite_tabi_system_matvec(parallel, unknowns, 0.5, 0.5,
                                         parallel_matvec);
  require(serial_matvec.size() == 2 * panels.size(),
          "E23 serial bridge output size");
  require(parallel_matvec.size() == 2 * panels.size(),
          "E23 parallel bridge output size");
  require(std::all_of(serial_matvec.begin(), serial_matvec.end(),
                      [](double v) { return std::isfinite(v); }),
          "E23 serial bridge finite");
  require(std::all_of(parallel_matvec.begin(), parallel_matvec.end(),
                      [](double v) { return std::isfinite(v); }),
          "E23 parallel bridge finite");
  const double bridge_diff =
      full_vector_relative_difference(serial_matvec, parallel_matvec);
  const double row_diff = full_row_relative_difference(serial, parallel);
  const double channel_diff =
      full_channel_relative_difference(serial, parallel);
  const bool hash_stable =
      serial.panel_geometry_hash == parallel.panel_geometry_hash &&
      serial.sampled_target_hash == parallel.sampled_target_hash &&
      serial.lite.worklist_hash == parallel.lite.worklist_hash;

  std::cout << std::setprecision(12)
            << "KIFMM_E23_BRIDGE fixture=" << fixture
            << " sample=" << sample_name
            << " sampled_rows=" << sample.size()
            << " order=" << order
            << " workers=" << workers
            << " serial_status=" << serial.status
            << " parallel_status=" << parallel.status
            << " bridge_relative_diff=" << bridge_diff
            << " row_relative_diff=" << row_diff
            << " channel_relative_diff=" << channel_diff
            << " worklist_hash_stable=" << (hash_stable ? 1 : 0)
            << " serial_worklist_hash=" << serial.lite.worklist_hash
            << " parallel_worklist_hash=" << parallel.lite.worklist_hash
            << " serial_policy_hash=" << serial.policy_hash
            << " parallel_policy_hash=" << parallel.policy_hash
            << " serial_cache_entries="
            << serial.equivalent_system_cache_entries
            << " parallel_cache_entries="
            << parallel.equivalent_system_cache_entries
            << " serial_far_panel_pairs="
            << serial.lite.accepted_far_panel_pairs
            << " parallel_far_panel_pairs="
            << parallel.lite.accepted_far_panel_pairs
            << " serial_fallback_panel_pairs="
            << serial.lite.direct_fallback_panel_pairs
            << " parallel_fallback_panel_pairs="
            << parallel.lite.direct_fallback_panel_pairs
            << '\n';
  require(bridge_diff <= 1.0e-12, "E23 serial/parallel bridge parity");
  require(row_diff <= 1.0e-12, "E23 serial/parallel row parity");
  require(channel_diff <= 1.0e-12, "E23 serial/parallel channel parity");
  require(hash_stable, "E23 serial/parallel hashes stable");
  std::cout << "PASS titan_fmm_kifmm_full_panel_matvec_smoke E23_PROFILE\n";
  return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 3) {
      throw std::runtime_error(
          "usage: titan_fmm_kifmm_full_panel_matvec_smoke <1sh1.pqr> <1a63.pqr> [--e19-profile|--e20-profile|--e22-profile|--e23-profile]");
    }
    bool e19_profile = false;
    bool e20_profile = false;
    bool e22_profile = false;
    bool e23_profile = false;
    std::string fixture = "1SH1";
    std::string sample = "uniform_hash";
    std::size_t sample_count = 128;
    int order = 6;
    std::size_t workers = 4;
    for (int i = 3; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--e19-profile") {
        e19_profile = true;
      } else if (arg == "--e20-profile") {
        e20_profile = true;
      } else if (arg == "--e22-profile") {
        e22_profile = true;
      } else if (arg == "--e23-profile") {
        e23_profile = true;
      } else if (arg.rfind("--fixture=", 0) == 0) {
        fixture = arg.substr(10);
      } else if (arg.rfind("--sample=", 0) == 0) {
        sample = arg.substr(9);
      } else if (arg.rfind("--sample-count=", 0) == 0) {
        sample_count = static_cast<std::size_t>(std::stoull(arg.substr(15)));
      } else if (arg.rfind("--order=", 0) == 0) {
        order = std::stoi(arg.substr(8));
      } else if (arg.rfind("--workers=", 0) == 0) {
        workers = static_cast<std::size_t>(std::stoull(arg.substr(10)));
      } else {
        throw std::runtime_error("unknown argument: " + arg);
      }
    }
    const int profile_count = (e19_profile ? 1 : 0) + (e20_profile ? 1 : 0) +
                              (e22_profile ? 1 : 0) + (e23_profile ? 1 : 0);
    if (profile_count != 1) {
      throw std::runtime_error(
          "smoke requires exactly one of --e19-profile, --e20-profile, --e22-profile, or --e23-profile");
    }
    if (!(fixture == "1SH1" || fixture == "1A63")) {
      throw std::runtime_error("fixture must be 1SH1 or 1A63");
    }
    if (!(order == 6 || order == 8)) {
      throw std::runtime_error("order must be 6 or 8");
    }
    if (e20_profile) {
      return run_e20_profile(argv[1], argv[2], fixture, sample, sample_count,
                             order);
    }
    if (e22_profile) {
      return run_e22_profile(argv[1], argv[2], fixture, sample, sample_count,
                             order, workers);
    }
    if (e23_profile) {
      return run_e23_profile(argv[1], argv[2], fixture, sample, sample_count,
                             order, workers);
    }
    return run_profile(argv[1], argv[2], fixture, sample, sample_count, order);
  } catch (const std::exception &ex) {
    std::cerr << "FAIL: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
