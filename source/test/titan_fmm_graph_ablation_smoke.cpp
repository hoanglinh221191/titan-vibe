#include "tabi/tabi_fmm_m2l_prototype.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::FmmM2lAnalyticContext;
using titan::tabi::fmm_m2l::KernelKind;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_graph_ablation_smoke: " << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<PrototypePanel> make_panels(std::size_t count) {
  std::vector<PrototypePanel> panels(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double t = static_cast<double>(i);
    panels[i].x = 0.17 * t + 0.31 * std::sin(0.73 * t);
    panels[i].y = 0.27 * std::cos(0.41 * t) + 0.03 * (i % 5);
    panels[i].z = 0.19 * std::sin(0.29 * t) + 0.07 * (i % 7);
    panels[i].nx = 0.0;
    panels[i].ny = 0.0;
    panels[i].nz = 1.0;
    panels[i].area = 0.4;
  }
  return panels;
}

std::vector<std::size_t> brute_knn(const std::vector<PrototypePanel> &panels,
                                   std::size_t target, std::size_t cap) {
  std::vector<std::pair<double, std::size_t>> ranked;
  ranked.reserve(panels.size());
  for (std::size_t source = 0; source < panels.size(); ++source) {
    const double dx = panels[source].x - panels[target].x;
    const double dy = panels[source].y - panels[target].y;
    const double dz = panels[source].z - panels[target].z;
    ranked.push_back({dx * dx + dy * dy + dz * dz, source});
  }
  std::sort(ranked.begin(), ranked.end());
  std::vector<std::size_t> row;
  for (std::size_t i = 0; i < std::min(cap, ranked.size()); ++i)
    row.push_back(ranked[i].second);
  std::sort(row.begin(), row.end());
  return row;
}

} // namespace

int main() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.12;
  config.max_leaf_size = 6;
  config.expansion_order = 4;
  config.separation_theta = 0.6;
  const auto panels = make_panels(47);
  FmmM2lAnalyticContext context(config, panels, 1.0 / 80.0,
                                config.kappa);

  constexpr std::size_t cap = 8;
  const auto knn = context.distance_knn_panel_graph(cap);
  expect(knn.size() == panels.size(), "kNN graph row count mismatch");
  for (std::size_t target = 0; target < panels.size(); ++target) {
    expect(knn[target] == brute_knn(panels, target, cap),
           "tree kNN differs from brute force at target " +
               std::to_string(target));
    expect(std::binary_search(knn[target].begin(), knn[target].end(), target),
           "kNN row omits self edge");
  }

  const auto near_graph = context.near_panel_graph(cap);
  expect(near_graph.size() == panels.size(), "FMM-near graph row mismatch");
  for (std::size_t target = 0; target < panels.size(); ++target) {
    expect(!near_graph[target].empty() && near_graph[target].size() <= cap,
           "FMM-near graph degree is invalid");
    expect(std::binary_search(near_graph[target].begin(),
                              near_graph[target].end(), target),
           "FMM-near row omits self edge");
  }

  const auto groups = context.distance_knn_overlap_panel_groups(6, 8);
  std::vector<int> owned_count(panels.size(), 0);
  for (const auto &group : groups) {
    expect(group.factor_panels.size() ==
               group.owned_panels.size() + group.halo_count,
           "factor/halo accounting mismatch");
    for (const std::size_t panel : group.owned_panels) {
      expect(panel < owned_count.size(), "owned panel out of range");
      ++owned_count[panel];
      expect(std::find(group.factor_panels.begin(), group.factor_panels.end(),
                       panel) != group.factor_panels.end(),
             "factor group omits owned panel");
    }
  }
  expect(std::all_of(owned_count.begin(), owned_count.end(),
                     [](int count) { return count == 1; }),
         "kNN overlap ownership is not an exact partition");

  std::cout << "PASS titan_fmm_graph_ablation_smoke panels=" << panels.size()
            << " cap=" << cap << " groups=" << groups.size() << '\n';
  return 0;
}
