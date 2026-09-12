#include "tabi/tabi_fmm_kifmm_real_panel_cloud.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
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
using titan::tabi::fmm_kifmm::KifmmRealPanelCloudResult;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::run_kifmm_lite_real_panel_matvec;

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
  const std::vector<Atom> atoms = read_pqr(pqr_path);
  std::vector<Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, options, &surface_result, &executor);
  if (panels.size() < 128) {
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

double edge_length(Vec3 a, Vec3 b) { return norm(sub(a, b)); }

double aspect_score(const Panel &p) {
  const Vec3 a{p.v0x, p.v0y, p.v0z};
  const Vec3 b{p.v1x, p.v1y, p.v1z};
  const Vec3 c{p.v2x, p.v2y, p.v2z};
  const double e0 = edge_length(a, b);
  const double e1 = edge_length(b, c);
  const double e2 = edge_length(c, a);
  const double lo = std::max(1.0e-300, std::min({e0, e1, e2}));
  const double hi = std::max({e0, e1, e2});
  return hi / lo;
}

template <typename KeyFn>
std::vector<std::size_t> select_by_key(const std::vector<Panel> &panels,
                                       std::size_t count, bool largest,
                                       KeyFn key_fn) {
  std::vector<std::size_t> idx(panels.size());
  std::iota(idx.begin(), idx.end(), 0);
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

std::vector<std::size_t>
select_nearest(const std::vector<Panel> &panels, Vec3 anchor,
               std::size_t count) {
  return select_by_key(panels, count, false, [&](const Panel &p, std::size_t) {
    const Vec3 d = sub(position(p), anchor);
    return dot(d, d);
  });
}

std::vector<Panel> gather(const std::vector<Panel> &panels,
                          const std::vector<std::vector<std::size_t>> &groups) {
  std::vector<Panel> out;
  std::vector<std::size_t> used;
  for (const auto &group : groups) {
    for (const std::size_t id : group) {
      if (std::find(used.begin(), used.end(), id) == used.end()) {
        used.push_back(id);
        out.push_back(panels[id]);
      }
    }
  }
  return out;
}

std::vector<double> unknowns_for(const std::vector<Panel> &panels) {
  std::vector<double> x(2 * panels.size(), 0.0);
  const Vec3 c = centroid(panels);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const Vec3 r = sub(position(panels[i]), c);
    const Vec3 n = normal(panels[i]);
    x[i] = 0.15 * std::sin(0.37 * static_cast<double>(i)) + 0.01 * r.x -
           0.006 * r.y + 0.02 * dot(unit(r), unit(n));
    x[panels.size() + i] =
        -0.12 * std::cos(0.23 * static_cast<double>(i)) + 0.004 * r.z -
        0.01 * dot(unit(r), unit(n));
  }
  return x;
}

KifmmRealPanelCloudConfig base_config() {
  KifmmRealPanelCloudConfig config;
  config.lite.order = 6;
  config.lite.leaf_size = 8;
  config.lite.opening_theta = 0.6;
  config.lite.equiv_scale = 1.05;
  config.lite.check_scale = 2.40;
  config.lite.kappa = 0.12746;
  config.lite.eps_ratio = 80.0;
  config.lite.max_condition_estimate = 1.0e12;
  config.lite.max_equivalent_relative_residual = 1.0e-10;
  config.lite.max_aggregate_relative_rms = 1.0e-6;
  config.lite.max_channel_relative_rms = 5.0e-6;
  config.lite.require_unit_normals = true;
  config.lite.solver = "qr";
  config.near_integration = "analytic_singular";
  config.quad_buffer = "off";
  return config;
}

double max_channel_rms(const KifmmRealPanelCloudResult &r) {
  return std::max({r.lite.k1_relative_rms, r.lite.k2_relative_rms,
                   r.lite.k3_relative_rms, r.lite.k4_relative_rms});
}

void print_result(const std::string &name,
                  const KifmmRealPanelCloudResult &result) {
  std::cout << std::setprecision(12)
            << "KIFMM_PROTEIN_PANEL_CASE name=" << name
            << " status=" << result.status
            << " panels=" << result.lite.source_panels
            << " leaves=" << result.lite.source_leaves
            << " far_cell_pairs=" << result.lite.accepted_far_cell_pairs
            << " far_panel_pairs=" << result.lite.accepted_far_panel_pairs
            << " fallback_panel_pairs="
            << result.lite.direct_fallback_panel_pairs
            << " production_fallback_valid_pairs="
            << result.production_fallback_valid_pairs
            << " production_direct_reference_valid_pairs="
            << result.production_direct_reference_valid_pairs
            << " invalid_near_pairs=" << result.production_invalid_pairs
            << " nonproduction_fallback_pairs="
            << result.nonproduction_fallback_pairs
            << " aggregate_rms=" << result.lite.aggregate_relative_rms
            << " max_channel_rms=" << max_channel_rms(result)
            << " max_equiv_residual="
            << result.lite.max_equivalent_relative_residual
            << " max_condition=" << result.lite.max_condition_estimate
            << " hash=" << result.lite.worklist_hash
            << " policy_hash=" << result.policy_hash << '\n';
}

void require_common(const std::string &name,
                    const KifmmRealPanelCloudResult &result) {
  require(result.nonproduction_fallback_pairs == 0,
          name + " uses no nonproduction fallback");
  require(!result.used_nonproduction_fallback,
          name + " did not flag nonproduction fallback");
  require(result.production_invalid_pairs == 0,
          name + " has no invalid production near kernel calls");
  require(result.production_fallback_valid_pairs ==
              result.lite.direct_fallback_panel_pairs,
          name + " production near count matches fallback work");
  require(result.production_direct_reference_valid_pairs > 0,
          name + " direct reference pairs are recorded");
  require(result.policy_hash != 0, name + " policy hash is recorded");
}

void require_direct_reference(const std::string &name,
                              const KifmmRealPanelCloudResult &result,
                              std::size_t n) {
  print_result(name, result);
  require(result.finite, name + " finite");
  require(result.status == "direct_only_reference",
          name + " is a direct-only diagnostic");
  require_common(name, result);
  require(result.lite.accepted_far_cell_pairs == 0,
          name + " has zero far cell pairs");
  require(result.production_fallback_valid_pairs == n * n,
          name + " dense production fallback covers all pairs");
  require(result.production_direct_reference_valid_pairs == n * n,
          name + " dense direct reference covers all pairs");
  require(result.lite.aggregate_relative_rms == 0.0,
          name + " direct-only RMS is zero");
}

void require_parity_case(const std::string &name,
                         const KifmmRealPanelCloudResult &result) {
  print_result(name, result);
  require(result.finite, name + " finite");
  require(result.status == "pass", name + " status pass");
  require_common(name, result);
  require(result.lite.accepted_far_cell_pairs > 0,
          name + " has accepted far work");
  require(result.lite.direct_fallback_panel_pairs > 0,
          name + " keeps direct fallback work");
  require(result.lite.aggregate_relative_rms <= 1.0e-6,
          name + " aggregate RMS gate");
  require(max_channel_rms(result) <= 5.0e-6, name + " channel RMS gate");
  require(result.lite.max_equivalent_relative_residual <= 1.0e-10,
          name + " equivalent residual gate");
  require(result.lite.max_condition_estimate <= 1.0e12,
          name + " condition gate");
}

KifmmRealPanelCloudResult run_case(const std::vector<Panel> &panels,
                                   KifmmRealPanelCloudConfig config) {
  return run_kifmm_lite_real_panel_matvec(panels, unknowns_for(panels),
                                          config);
}

} // namespace

int main(int argc, char **argv) {
  try {
    const std::string pqr_path =
        argc > 1 ? argv[1] : "test/pbe_inputs_public/1sh1_full_amber_ph7.pqr";
    const std::vector<Panel> full_panels = build_real_surface(pqr_path);
    require(full_panels.size() >= 256,
            "protein-derived surface has enough panels");
    const Vec3 c = centroid(full_panels);

    const auto min_x =
        select_by_key(full_panels, 14, false,
                      [](const Panel &p, std::size_t) { return p.x; });
    const auto max_x =
        select_by_key(full_panels, 14, true,
                      [](const Panel &p, std::size_t) { return p.x; });
    const auto min_y =
        select_by_key(full_panels, 10, false,
                      [](const Panel &p, std::size_t) { return p.y; });
    const auto radial_inward =
        select_by_key(full_panels, 16, false,
                      [&](const Panel &p, std::size_t) {
                        return dot(unit(sub(position(p), c)), unit(normal(p)));
                      });
    const auto high_aspect =
        select_by_key(full_panels, 24, true,
                      [](const Panel &p, std::size_t) {
                        return aspect_score(p);
                      });

    const std::vector<Panel> all_near_cloud =
        gather(full_panels, {select_nearest(full_panels, position(full_panels[min_x[0]]), 8)});
    KifmmRealPanelCloudConfig direct_config = base_config();
    direct_config.lite.require_far_work = false;
    direct_config.lite.leaf_size = static_cast<int>(all_near_cloud.size());
    const KifmmRealPanelCloudResult all_near =
        run_case(all_near_cloud, direct_config);
    require_direct_reference("all_near_protein_patch", all_near,
                             all_near_cloud.size());

    const std::vector<Panel> separated_cloud =
        gather(full_panels, {min_x, max_x});
    const KifmmRealPanelCloudResult separated =
        run_case(separated_cloud, base_config());
    require_parity_case("separated_surface_patches", separated);

    const std::vector<Panel> mixed_cloud = gather(
        full_panels,
        {select_nearest(full_panels, position(full_panels[min_y[0]]), 16),
         max_x});
    const KifmmRealPanelCloudResult mixed = run_case(mixed_cloud, base_config());
    require_parity_case("mixed_near_far_patch", mixed);

    const std::vector<Panel> groove_cloud =
        gather(full_panels, {radial_inward, max_x});
    const KifmmRealPanelCloudResult groove =
        run_case(groove_cloud, base_config());
    require_parity_case("inward_normal_groove_patch", groove);

    const std::vector<Panel> aspect_cloud =
        gather(full_panels, {high_aspect, min_x});
    const KifmmRealPanelCloudResult aspect =
        run_case(aspect_cloud, base_config());
    require_parity_case("high_aspect_triangle_stress", aspect);

    const KifmmRealPanelCloudResult mixed_repeat =
        run_case(mixed_cloud, base_config());
    require(mixed_repeat.lite.worklist_hash == mixed.lite.worklist_hash,
            "mixed protein worklist hash is deterministic");
    require(mixed_repeat.lite.accepted_far_panel_pairs ==
                mixed.lite.accepted_far_panel_pairs,
            "mixed protein far panel count is deterministic");
    require(mixed_repeat.production_fallback_valid_pairs ==
                mixed.production_fallback_valid_pairs,
            "mixed protein fallback count is deterministic");
    require(mixed_repeat.policy_hash == mixed.policy_hash,
            "mixed protein policy hash is deterministic");

    KifmmRealPanelCloudConfig order4 = base_config();
    order4.lite.order = 4;
    const KifmmRealPanelCloudResult order4_mixed = run_case(mixed_cloud, order4);
    print_result("mixed_near_far_patch_q4_diagnostic", order4_mixed);
    require(order4_mixed.lite.accepted_far_cell_pairs ==
                mixed.lite.accepted_far_cell_pairs,
            "q4 uses the same mixed work partition");
    require(order4_mixed.production_invalid_pairs == 0,
            "q4 diagnostic has valid production fallback");

    KifmmRealPanelCloudConfig order8 = base_config();
    order8.lite.order = 8;
    const KifmmRealPanelCloudResult order8_mixed = run_case(mixed_cloud, order8);
    require_parity_case("mixed_near_far_patch_q8_confirm", order8_mixed);
    require(order8_mixed.lite.aggregate_relative_rms <=
                std::max(1.0e-12, 2.0 * mixed.lite.aggregate_relative_rms),
            "q8 mixed accuracy is consistent with q6");

    std::cout << std::setprecision(12)
              << "KIFMM_PROTEIN_PANEL_SUMMARY panels_total="
              << full_panels.size()
              << " separated_rms=" << separated.lite.aggregate_relative_rms
              << " mixed_rms=" << mixed.lite.aggregate_relative_rms
              << " groove_rms=" << groove.lite.aggregate_relative_rms
              << " aspect_rms=" << aspect.lite.aggregate_relative_rms
              << " mixed_hash=" << mixed.lite.worklist_hash << '\n';
    std::cout << "PASS titan_fmm_kifmm_protein_panel_cloud_smoke\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
