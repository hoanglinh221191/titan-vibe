#include "tabi/tabi_fmm_kifmm_surface.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

namespace {

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmSurfacePoint;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::instantiate_surface;
using titan::tabi::fmm_kifmm::kifmm_cube_surface_unique_count;
using titan::tabi::fmm_kifmm::make_cube_surface_rule;

void fail(const std::string &message) {
  std::cerr << message << '\n';
  std::exit(1);
}

void require_true(bool ok, const std::string &message) {
  if (!ok) {
    fail(message);
  }
}

bool nearly_equal(double a, double b) {
  return std::abs(a - b) <= 1.0e-14;
}

std::string key_for_unit_point(const KifmmSurfacePoint &point) {
  auto key_component = [](double value) -> int {
    if (nearly_equal(value, -1.0)) {
      return -1;
    }
    if (nearly_equal(value, 0.0)) {
      return 0;
    }
    if (nearly_equal(value, 1.0)) {
      return 1;
    }
    return 1000 + static_cast<int>(std::llround(value * 1000000.0));
  };
  std::ostringstream out;
  out << key_component(point.x.x) << ',' << key_component(point.x.y) << ','
      << key_component(point.x.z);
  return out.str();
}

double dot(Vec3 a, Vec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

void check_rule(int order, std::size_t expected_count) {
  const auto rule = make_cube_surface_rule(order, 1.05, 1.60);
  const auto repeat = make_cube_surface_rule(order, 1.05, 1.60);
  require_true(kifmm_cube_surface_unique_count(order) == expected_count,
               "wrong closed-form surface count");
  require_true(rule.unit_equiv.size() == expected_count,
               "unexpected equivalent surface count");
  require_true(rule.unit_check.size() == expected_count,
               "unexpected check surface count");
  require_true(repeat.unit_equiv.size() == rule.unit_equiv.size(),
               "repeat surface count changed");

  std::set<std::string> seen;
  for (std::size_t i = 0; i < rule.unit_equiv.size(); ++i) {
    const auto &equiv = rule.unit_equiv[i];
    const auto &check = rule.unit_check[i];
    const auto &again = repeat.unit_equiv[i];
    require_true(nearly_equal(equiv.x.x, again.x.x) &&
                     nearly_equal(equiv.x.y, again.x.y) &&
                     nearly_equal(equiv.x.z, again.x.z),
                 "surface order is not deterministic");
    require_true(std::isfinite(equiv.x.x) && std::isfinite(equiv.x.y) &&
                     std::isfinite(equiv.x.z),
                 "non-finite surface coordinate");
    require_true(nearly_equal(equiv.weight, 1.0), "unexpected unit weight");

    const int normal_axes =
        (std::abs(equiv.outward_normal.x) > 0.5 ? 1 : 0) +
        (std::abs(equiv.outward_normal.y) > 0.5 ? 1 : 0) +
        (std::abs(equiv.outward_normal.z) > 0.5 ? 1 : 0);
    require_true(normal_axes == 1, "normal is not axis-aligned");
    require_true(dot(equiv.x, equiv.outward_normal) > 0.0,
                 "normal is not outward");

    const std::string key = key_for_unit_point(equiv);
    require_true(seen.insert(key).second, "duplicate cube-surface key");
    require_true(nearly_equal(check.outward_normal.x, equiv.outward_normal.x) &&
                     nearly_equal(check.outward_normal.y, equiv.outward_normal.y) &&
                     nearly_equal(check.outward_normal.z, equiv.outward_normal.z),
                 "check/equiv normal ordering mismatch");
    require_true(dot(check.x, check.x) > dot(equiv.x, equiv.x),
                 "check surface is not outside equivalent surface");
  }

  const KifmmBox box{{2.0, -3.0, 4.5}, 2.0, 7, 42};
  const auto placed_equiv = instantiate_surface(rule.unit_equiv, box);
  const auto placed_check = instantiate_surface(rule.unit_check, box);
  require_true(placed_equiv.size() == rule.unit_equiv.size(),
               "placed equivalent surface size mismatch");
  for (std::size_t i = 0; i < placed_equiv.size(); ++i) {
    const Vec3 rel_equiv{placed_equiv[i].x.x - box.center.x,
                         placed_equiv[i].x.y - box.center.y,
                         placed_equiv[i].x.z - box.center.z};
    const Vec3 rel_check{placed_check[i].x.x - box.center.x,
                         placed_check[i].x.y - box.center.y,
                         placed_check[i].x.z - box.center.z};
    require_true(dot(rel_check, rel_check) > dot(rel_equiv, rel_equiv),
                 "placed check surface is not outside placed equivalent surface");
  }
}

} // namespace

int main() {
  check_rule(4, 56);
  check_rule(6, 152);
  check_rule(8, 296);
  std::cout << "PASS titan_fmm_kifmm_surface_smoke\n";
  return 0;
}
