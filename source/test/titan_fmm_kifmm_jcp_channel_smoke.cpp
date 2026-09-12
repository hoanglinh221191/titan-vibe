#include "tabi/tabi_fmm_kifmm_fields.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmJcpChannelRow;
using titan::tabi::fmm_kifmm::KifmmJcpChannels;
using titan::tabi::fmm_kifmm::KifmmJcpRow;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::Vec3;
using titan::tabi::fmm_kifmm::assemble_direct_jcp_channel_row;
using titan::tabi::fmm_kifmm::assemble_direct_jcp_row;
using titan::tabi::fmm_kifmm::assemble_kifmm_jcp_channel_row;
using titan::tabi::fmm_kifmm::assemble_kifmm_jcp_channels;
using titan::tabi::fmm_kifmm::assemble_kifmm_jcp_row;
using titan::tabi::fmm_kifmm::make_cube_surface_rule;

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

double norm(Vec3 a) {
  return std::sqrt(dot(a, a));
}

double distance(Vec3 a, Vec3 b) {
  return norm({a.x - b.x, a.y - b.y, a.z - b.z});
}

Vec3 unit(Vec3 a) {
  const double n = std::max(norm(a), 1.0e-300);
  return {a.x / n, a.y / n, a.z / n};
}

std::vector<KifmmPanelSource> deterministic_sources(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {{box.center.x - 0.18 * h, box.center.y + 0.07 * h,
        box.center.z - 0.10 * h},
       unit({0.62, -0.37, 0.44}),
       0.79,
       0.31,
       -0.47},
      {{box.center.x + 0.15 * h, box.center.y - 0.13 * h,
        box.center.z + 0.09 * h},
       unit({-0.35, 0.81, -0.27}),
       0.96,
       -0.52,
       0.58},
      {{box.center.x + 0.03 * h, box.center.y + 0.17 * h,
        box.center.z - 0.16 * h},
       unit({0.21, 0.44, 0.87}),
       0.68,
       0.49,
       -0.26},
      {{box.center.x - 0.11 * h, box.center.y - 0.05 * h,
        box.center.z + 0.14 * h},
       unit({-0.67, -0.28, 0.49}),
       0.88,
       -0.34,
       0.39},
  };
}

std::vector<Vec3> external_targets(const KifmmBox &box) {
  const double h = box.half_width;
  return {
      {box.center.x + 7.80 * h, box.center.y + 0.22 * h,
       box.center.z - 0.31 * h},
      {box.center.x - 7.90 * h, box.center.y - 0.12 * h,
       box.center.z + 0.37 * h},
      {box.center.x + 0.19 * h, box.center.y + 7.70 * h,
       box.center.z + 0.27 * h},
      {box.center.x - 0.29 * h, box.center.y - 7.85 * h,
       box.center.z - 0.16 * h},
      {box.center.x + 0.24 * h, box.center.y + 0.36 * h,
       box.center.z + 8.05 * h},
      {box.center.x - 0.34 * h, box.center.y - 0.18 * h,
       box.center.z - 7.70 * h},
  };
}

std::vector<Vec3> target_normals() {
  return {unit({0.96, 0.15, -0.12}),  unit({-0.92, -0.08, 0.22}),
          unit({0.12, 0.98, 0.09}),   unit({-0.14, -0.95, -0.18}),
          unit({0.08, 0.15, 0.99}),   unit({-0.11, -0.19, -0.97})};
}

double rel_rms_rows(const std::vector<KifmmJcpRow> &approx,
                    const std::vector<KifmmJcpRow> &direct) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const double d_phi = approx[i].phi_row - direct[i].phi_row;
    const double d_phin = approx[i].phin_row - direct[i].phin_row;
    err_sq += d_phi * d_phi + d_phin * d_phin;
    ref_sq += direct[i].phi_row * direct[i].phi_row +
              direct[i].phin_row * direct[i].phin_row;
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

double rel_rms_channel(const std::vector<KifmmJcpChannelRow> &approx,
                       const std::vector<KifmmJcpChannelRow> &direct,
                       char channel) {
  double err_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    double a = 0.0;
    double d = 0.0;
    if (channel == '1') {
      a = approx[i].k1;
      d = direct[i].k1;
    } else if (channel == '2') {
      a = approx[i].k2;
      d = direct[i].k2;
    } else if (channel == '3') {
      a = approx[i].k3;
      d = direct[i].k3;
    } else {
      a = approx[i].k4;
      d = direct[i].k4;
    }
    const double diff = a - d;
    err_sq += diff * diff;
    ref_sq += d * d;
  }
  return std::sqrt(err_sq / std::max(ref_sq, 1.0e-300));
}

} // namespace

int main() {
  const KifmmBox box{{0.15, -0.25, 0.40}, 1.0, 0, 7};
  const auto rule = make_cube_surface_rule(6, 1.05, 2.40);
  const auto sources = deterministic_sources(box);
  const auto targets = external_targets(box);
  const auto normals = target_normals();
  double target_min_distance = std::numeric_limits<double>::infinity();
  double target_max_distance = 0.0;
  for (const Vec3 target : targets) {
    const double d = distance(target, box.center);
    target_min_distance = std::min(target_min_distance, d);
    target_max_distance = std::max(target_max_distance, d);
  }
  constexpr double kappa = 0.12746;
  constexpr double eps_ratio = 80.0;
  constexpr double inv_eps_ratio = 1.0 / eps_ratio;

  const KifmmJcpChannels channels = assemble_kifmm_jcp_channels(
      rule, box, sources, kappa, eps_ratio, inv_eps_ratio, "qr");
  require(channels.finite, "KIFMM JCP channel solves finite");
  require(channels.status == "pass", "KIFMM JCP channel status pass");

  std::vector<KifmmJcpRow> approx;
  std::vector<KifmmJcpRow> direct;
  std::vector<KifmmJcpChannelRow> approx_channels;
  std::vector<KifmmJcpChannelRow> direct_channels;
  for (std::size_t i = 0; i < targets.size(); ++i) {
    approx.push_back(assemble_kifmm_jcp_row(rule, box, channels, kappa,
                                            eps_ratio, inv_eps_ratio,
                                            targets[i], normals[i]));
    direct.push_back(assemble_direct_jcp_row(sources, kappa, eps_ratio,
                                             inv_eps_ratio, targets[i],
                                             normals[i]));
    approx_channels.push_back(assemble_kifmm_jcp_channel_row(
        rule, box, channels, kappa, eps_ratio, inv_eps_ratio, targets[i],
        normals[i]));
    direct_channels.push_back(assemble_direct_jcp_channel_row(
        sources, kappa, eps_ratio, inv_eps_ratio, targets[i], normals[i]));
  }

  const double rms = rel_rms_rows(approx, direct);
  const double k1_rms = rel_rms_channel(approx_channels, direct_channels, '1');
  const double k2_rms = rel_rms_channel(approx_channels, direct_channels, '2');
  const double k3_rms = rel_rms_channel(approx_channels, direct_channels, '3');
  const double k4_rms = rel_rms_channel(approx_channels, direct_channels, '4');
  std::cout << std::setprecision(10)
            << "KIFMM_JCP order=" << rule.order
            << " equiv_scale=" << rule.equiv_scale
            << " check_scale=" << rule.check_scale
            << " external_target_min_distance=" << target_min_distance
            << " external_target_max_distance=" << target_max_distance
            << " kernel=laplace+yukawa kappa=" << kappa
            << " solver=qr"
            << " k1_rel_rms=" << k1_rms << " k2_rel_rms=" << k2_rms
            << " k3_rel_rms=" << k3_rms << " k4_rel_rms=" << k4_rms
            << " aggregate_row_rel_rms=" << rms
            << " max_equiv_residual=" << channels.max_relative_residual
            << " max_condition=" << channels.max_condition_estimate
            << " equivalent_system_builds="
            << channels.equivalent_system_builds << '\n';
  require(std::isfinite(rms) && rms <= 1.0e-4,
          "synthetic K1-K4 row RMS gate");
  require(std::isfinite(k1_rms) && std::isfinite(k2_rms) &&
              std::isfinite(k3_rms) && std::isfinite(k4_rms),
          "synthetic K1-K4 channel RMS finite");
  std::cout << "PASS titan_fmm_kifmm_jcp_channel_smoke\n";
  return EXIT_SUCCESS;
}
