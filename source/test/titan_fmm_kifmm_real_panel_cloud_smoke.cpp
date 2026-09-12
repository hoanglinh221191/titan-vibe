#include "tabi/tabi_fmm_kifmm_real_panel_cloud.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using TitanPBE::Tabi::Panel;
using titan::tabi::fmm_kifmm::KifmmLiteConfig;
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

double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

Panel make_panel(Vec3 center, Vec3 normal, double area, std::size_t id) {
  constexpr double pi = 3.141592653589793238462643383279502884;
  normal = unit(normal);
  Vec3 tangent = cross(normal, {0.0, 0.0, 1.0});
  if (norm(tangent) < 1.0e-8) {
    tangent = cross(normal, {0.0, 1.0, 0.0});
  }
  tangent = unit(tangent);
  const Vec3 bitangent = unit(cross(normal, tangent));
  const double side = std::sqrt(4.0 * area / std::sqrt(3.0));
  const double r = side / std::sqrt(3.0);

  Panel panel;
  panel.x = center.x;
  panel.y = center.y;
  panel.z = center.z;
  panel.nx = normal.x;
  panel.ny = normal.y;
  panel.nz = normal.z;
  panel.area = area;
  panel.solid_angle = 2.0 * pi;
  panel.diagonal = 0.5;
  panel.atom_index = id;
  const Vec3 v0{center.x + r * tangent.x,
                center.y + r * tangent.y,
                center.z + r * tangent.z};
  const Vec3 v1{center.x + r * (-0.5 * tangent.x + 0.8660254037844386 * bitangent.x),
                center.y + r * (-0.5 * tangent.y + 0.8660254037844386 * bitangent.y),
                center.z + r * (-0.5 * tangent.z + 0.8660254037844386 * bitangent.z)};
  const Vec3 v2{center.x + r * (-0.5 * tangent.x - 0.8660254037844386 * bitangent.x),
                center.y + r * (-0.5 * tangent.y - 0.8660254037844386 * bitangent.y),
                center.z + r * (-0.5 * tangent.z - 0.8660254037844386 * bitangent.z)};
  panel.v0x = v0.x;
  panel.v0y = v0.y;
  panel.v0z = v0.z;
  panel.v1x = v1.x;
  panel.v1y = v1.y;
  panel.v1z = v1.z;
  panel.v2x = v2.x;
  panel.v2y = v2.y;
  panel.v2z = v2.z;
  panel.v0_index = 3 * id;
  panel.v1_index = 3 * id + 1;
  panel.v2_index = 3 * id + 2;
  panel.has_triangle_geometry = true;
  return panel;
}

std::vector<Panel> cluster(Vec3 center, double spacing, int count,
                           std::size_t offset) {
  std::vector<Panel> panels;
  panels.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const int ix = i % 3;
    const int iy = (i / 3) % 3;
    const int iz = i / 9;
    const Vec3 p{center.x + spacing * (ix - 1),
                 center.y + spacing * (iy - 1),
                 center.z + spacing * (iz - 0.5)};
    const Vec3 n = unit({0.31 + 0.07 * i, -0.22 + 0.05 * i,
                         0.63 - 0.02 * i});
    panels.push_back(make_panel(p, n, 0.08 + 0.004 * (i % 5),
                                offset + static_cast<std::size_t>(i)));
  }
  return panels;
}

std::vector<double> unknowns_for(std::size_t n) {
  std::vector<double> x(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = ((i % 2) == 0 ? 0.31 : -0.24) + 0.009 * static_cast<double>(i);
    x[n + i] = ((i % 3) == 0 ? -0.41 : 0.29) - 0.007 * static_cast<double>(i);
  }
  return x;
}

KifmmRealPanelCloudConfig base_config() {
  KifmmRealPanelCloudConfig config;
  config.lite.order = 8;
  config.lite.leaf_size = 8;
  config.lite.opening_theta = 0.6;
  config.lite.equiv_scale = 1.05;
  config.lite.check_scale = 2.40;
  config.lite.kappa = 0.12746;
  config.lite.eps_ratio = 80.0;
  config.lite.max_aggregate_relative_rms = 1.0e-6;
  config.lite.max_channel_relative_rms = 5.0e-6;
  config.lite.require_unit_normals = true;
  config.lite.solver = "qr";
  config.near_integration = "analytic_singular";
  config.quad_buffer = "off";
  return config;
}

void print_result(const KifmmRealPanelCloudResult &result,
                  const std::string &name) {
  std::cout << std::setprecision(10) << "KIFMM_REAL_PANEL_CASE name=" << name
            << " aggregate=" << result.lite.aggregate_relative_rms
            << " k1=" << result.lite.k1_relative_rms
            << " k2=" << result.lite.k2_relative_rms
            << " k3=" << result.lite.k3_relative_rms
            << " k4=" << result.lite.k4_relative_rms
            << " far_pairs=" << result.lite.accepted_far_cell_pairs
            << " fallback_pairs=" << result.lite.direct_fallback_panel_pairs
            << " production_fallback_valid_pairs="
            << result.production_fallback_valid_pairs
            << " production_direct_reference_valid_pairs="
            << result.production_direct_reference_valid_pairs
            << " production_invalid_pairs="
            << result.production_invalid_pairs
            << " nonproduction_fallback_pairs="
            << result.nonproduction_fallback_pairs
            << " policy_hash=" << result.policy_hash
            << " status=" << result.status << '\n';
}

void require_common_clean(const KifmmRealPanelCloudResult &result,
                          const std::string &name) {
  require(result.nonproduction_fallback_pairs == 0,
          name + " has no nonproduction fallback pairs");
  require(!result.used_nonproduction_fallback, name + " no nonproduction fallback");
  require(result.production_invalid_pairs == 0,
          name + " has no invalid production near kernels");
  require(result.policy_hash != 0, name + " has a policy hash");
}

void require_direct_reference(const KifmmRealPanelCloudResult &result,
                              const std::string &name) {
  print_result(result, name);
  require(result.finite, name + " finite direct reference");
  require(result.status == "direct_only_reference",
          name + " direct-only reference status");
  require(result.lite.status == result.status,
          name + " lite status mirrors top-level status");
  require_common_clean(result, name);
  require(result.lite.aggregate_relative_rms == 0.0,
          name + " aggregate RMS is zero");
  require(result.production_direct_reference_valid_pairs > 0,
          name + " direct reference pair count is recorded");
}

void require_parity(const KifmmRealPanelCloudResult &result, double gate,
                    const std::string &name) {
  print_result(result, name);
  require(result.finite, name + " finite");
  require(result.status == "pass", name + " status pass");
  require_common_clean(result, name);
  require(result.lite.aggregate_relative_rms <= gate,
          name + " aggregate RMS gate");
  require(result.lite.k1_relative_rms <= gate, name + " K1 RMS gate");
  require(result.lite.k2_relative_rms <= gate, name + " K2 RMS gate");
  require(result.lite.k3_relative_rms <= gate, name + " K3 RMS gate");
  require(result.lite.k4_relative_rms <= gate, name + " K4 RMS gate");
}

} // namespace

int main() {
  const KifmmLiteConfig default_lite_config;
  require(default_lite_config.max_channel_relative_rms == 5.0e-6,
          "default channel RMS gate is E6 report gate");

  const KifmmRealPanelCloudConfig config = base_config();

  const auto near_panels = cluster({0.0, 0.0, 0.0}, 0.025, 8, 0);
  KifmmRealPanelCloudConfig all_near_config = config;
  all_near_config.lite.require_far_work = false;
  const KifmmRealPanelCloudResult all_near =
      run_kifmm_lite_real_panel_matvec(near_panels,
                                       unknowns_for(near_panels.size()),
                                       all_near_config);
  require_direct_reference(all_near, "all-near");
  require(all_near.lite.accepted_far_cell_pairs == 0,
          "all-near has no far cell pairs");
  require(all_near.production_fallback_valid_pairs ==
              near_panels.size() * near_panels.size(),
          "all-near production fallback covers every pair");
  require(all_near.production_direct_reference_valid_pairs ==
              near_panels.size() * near_panels.size(),
          "all-near direct reference covers every pair");

  auto separated_panels = cluster({-8.0, 0.0, 0.0}, 0.045, 10, 100);
  const auto separated_tail = cluster({8.0, 0.2, -0.1}, 0.045, 10, 200);
  separated_panels.insert(separated_panels.end(), separated_tail.begin(),
                          separated_tail.end());
  const KifmmRealPanelCloudResult separated =
      run_kifmm_lite_real_panel_matvec(separated_panels,
                                       unknowns_for(separated_panels.size()),
                                       config);
  require_parity(separated, 1.0e-4, "separated");
  require(separated.lite.accepted_far_panel_pairs > 0,
          "separated has accepted far work");
  require(separated.production_fallback_valid_pairs > 0,
          "separated still uses production near fallback");

  auto mixed_panels = cluster({-8.0, 0.0, 0.0}, 0.045, 12, 300);
  const auto mixed_tail = cluster({3.0, 0.1, 0.0}, 0.035, 8, 400);
  mixed_panels.insert(mixed_panels.end(), mixed_tail.begin(), mixed_tail.end());
  const KifmmRealPanelCloudResult mixed =
      run_kifmm_lite_real_panel_matvec(mixed_panels,
                                       unknowns_for(mixed_panels.size()),
                                       config);
  require_parity(mixed, 1.0e-4, "mixed");
  require(mixed.lite.accepted_far_panel_pairs > 0, "mixed far pairs");
  require(mixed.lite.direct_fallback_panel_pairs > 0, "mixed fallback pairs");
  require(mixed.production_fallback_valid_pairs ==
              mixed.lite.direct_fallback_panel_pairs,
          "mixed production near count matches fallback pairs");
  const KifmmRealPanelCloudResult mixed_repeat =
      run_kifmm_lite_real_panel_matvec(mixed_panels,
                                       unknowns_for(mixed_panels.size()),
                                       config);
  require(mixed_repeat.lite.worklist_hash == mixed.lite.worklist_hash,
          "mixed worklist hash is deterministic");
  require(mixed_repeat.policy_hash == mixed.policy_hash,
          "mixed policy hash is deterministic");

  KifmmRealPanelCloudConfig looser_channel_gate = config;
  looser_channel_gate.lite.max_channel_relative_rms = 1.0e-4;
  const KifmmRealPanelCloudResult looser =
      run_kifmm_lite_real_panel_matvec(mixed_panels,
                                       unknowns_for(mixed_panels.size()),
                                       looser_channel_gate);
  require(looser.lite.worklist_hash == mixed.lite.worklist_hash,
          "accuracy gate change preserves worklist hash");
  require(looser.policy_hash != mixed.policy_hash,
          "accuracy gate change alters policy hash");

  std::vector<Panel> piercing = separated_panels;
  piercing[10].v0x = separated_panels[0].x;
  piercing[10].v0y = separated_panels[0].y;
  piercing[10].v0z = separated_panels[0].z;
  const KifmmRealPanelCloudResult vertex_guard =
      run_kifmm_lite_real_panel_matvec(piercing,
                                       unknowns_for(piercing.size()), config);
  require_parity(vertex_guard, 1.0e-4, "vertex-check-fallback");
  require(vertex_guard.lite.direct_fallback_panel_pairs >
              separated.lite.direct_fallback_panel_pairs,
          "target vertex inside check surface forces fallback");

  std::vector<Panel> shared_vertex = separated_panels;
  shared_vertex[10].v0_index = shared_vertex[0].v0_index;
  const KifmmRealPanelCloudResult adjacency_guard =
      run_kifmm_lite_real_panel_matvec(shared_vertex,
                                       unknowns_for(shared_vertex.size()),
                                       config);
  require_parity(adjacency_guard, 1.0e-4, "shared-vertex-fallback");
  require(adjacency_guard.lite.direct_fallback_panel_pairs >
              separated.lite.direct_fallback_panel_pairs,
          "shared vertex forces production fallback");

  KifmmRealPanelCloudConfig too_tight = config;
  too_tight.lite.max_channel_relative_rms = 1.0e-12;
  const KifmmRealPanelCloudResult gate_fail =
      run_kifmm_lite_real_panel_matvec(mixed_panels,
                                       unknowns_for(mixed_panels.size()),
                                       too_tight);
  require(!gate_fail.finite, "too-tight channel RMS gate fails in library");
  require(gate_fail.status == "channel_rms_gate_failed",
          "too-tight channel RMS status");

  std::vector<Panel> bad_normal = near_panels;
  bad_normal.front().nx *= 2.0;
  bad_normal.front().ny *= 2.0;
  bad_normal.front().nz *= 2.0;
  const KifmmRealPanelCloudResult nonunit =
      run_kifmm_lite_real_panel_matvec(bad_normal,
                                       unknowns_for(bad_normal.size()),
                                       all_near_config);
  require(!nonunit.finite, "non-unit real panel normal is rejected");
  require(nonunit.status == "invalid_panel", "non-unit normal status");

  std::vector<Panel> invalid = near_panels;
  invalid.front().area = -1.0;
  const KifmmRealPanelCloudResult bad =
      run_kifmm_lite_real_panel_matvec(invalid, unknowns_for(invalid.size()),
                                       all_near_config);
  require(!bad.finite, "invalid real panel is rejected");
  require(bad.status == "invalid_panel", "invalid panel status");

  std::vector<Panel> invalid_vertex = near_panels;
  invalid_vertex.front().v0x = std::numeric_limits<double>::quiet_NaN();
  const KifmmRealPanelCloudResult bad_vertex =
      run_kifmm_lite_real_panel_matvec(invalid_vertex,
                                       unknowns_for(invalid_vertex.size()),
                                       all_near_config);
  require(!bad_vertex.finite, "non-finite vertex is rejected");
  require(bad_vertex.status == "invalid_panel", "non-finite vertex status");

  std::vector<Panel> degenerate = near_panels;
  degenerate.front().v1x = degenerate.front().v0x;
  degenerate.front().v1y = degenerate.front().v0y;
  degenerate.front().v1z = degenerate.front().v0z;
  degenerate.front().v2x = degenerate.front().v0x;
  degenerate.front().v2y = degenerate.front().v0y;
  degenerate.front().v2z = degenerate.front().v0z;
  const KifmmRealPanelCloudResult bad_triangle =
      run_kifmm_lite_real_panel_matvec(degenerate,
                                       unknowns_for(degenerate.size()),
                                       all_near_config);
  require(!bad_triangle.finite, "degenerate triangle is rejected");
  require(bad_triangle.status == "invalid_panel", "degenerate triangle status");

  std::cout << std::setprecision(10)
            << "KIFMM_REAL_PANEL all_near_rms="
            << all_near.lite.aggregate_relative_rms
            << " separated_rms=" << separated.lite.aggregate_relative_rms
            << " mixed_rms=" << mixed.lite.aggregate_relative_rms
            << " mixed_hash=" << mixed.lite.worklist_hash
            << " policy_hash=" << mixed.policy_hash
            << " production_fallback_valid_pairs="
            << mixed.production_fallback_valid_pairs
            << '\n';
  std::cout << "PASS titan_fmm_kifmm_real_panel_cloud_smoke\n";
  return EXIT_SUCCESS;
}
