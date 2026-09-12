#include "tabi/tabi_fmm_kifmm_fields.h"
#include "tabi/tabi_fmm_kifmm_hierarchy_jcp.h"
#include "tabi/tabi_fmm_kifmm_surface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using titan::tabi::fmm_kifmm::KifmmBox;
using titan::tabi::fmm_kifmm::KifmmEquivalentSystemCache;
using titan::tabi::fmm_kifmm::KifmmFieldFamily;
using titan::tabi::fmm_kifmm::KifmmHierarchyJcpConfig;
using titan::tabi::fmm_kifmm::KifmmJcpChannelRow;
using titan::tabi::fmm_kifmm::KifmmJcpChannels;
using titan::tabi::fmm_kifmm::KifmmKernel;
using titan::tabi::fmm_kifmm::KifmmPanelSource;
using titan::tabi::fmm_kifmm::KifmmSurfaceRule;
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

std::vector<KifmmPanelSource> multi_m2l_sources() {
  std::vector<KifmmPanelSource> sources;
  sources.reserve(128);
  const std::array<Vec3, 8> centers{
      Vec3{-4.0, -2.4, -2.0}, Vec3{-4.0, -2.4, 2.0},
      Vec3{-4.0, 2.4, -2.0},  Vec3{-4.0, 2.4, 2.0},
      Vec3{-1.6, -2.4, -2.0}, Vec3{-1.6, -2.4, 2.0},
      Vec3{-1.6, 2.4, -2.0},  Vec3{-1.6, 2.4, 2.0}};
  for (std::size_t c = 0; c < centers.size(); ++c) {
    for (int i = 0; i < 16; ++i) {
      const double dx = 0.035 * (i % 4);
      const double dy = 0.031 * ((i / 4) % 4);
      const double dz = 0.029 * (i / 8);
      sources.push_back(panel(
          {centers[c].x + dx, centers[c].y + dy, centers[c].z + dz},
          {0.41 + 0.03 * i, -0.37 + 0.02 * static_cast<int>(c),
           0.58 - 0.01 * i},
          0.50 + 0.02 * ((i + static_cast<int>(c)) % 5),
          0.15 + 0.04 * ((i + 2 * static_cast<int>(c)) % 9),
          -0.31 + 0.05 * ((2 * i + static_cast<int>(c)) % 11)));
    }
  }
  return sources;
}

std::vector<Vec3> multi_m2l_targets() {
  std::vector<Vec3> targets;
  targets.reserve(128);
  const std::array<Vec3, 8> centers{
      Vec3{1.6, -2.4, -2.0}, Vec3{1.6, -2.4, 2.0},
      Vec3{1.6, 2.4, -2.0},  Vec3{1.6, 2.4, 2.0},
      Vec3{4.0, -2.4, -2.0}, Vec3{4.0, -2.4, 2.0},
      Vec3{4.0, 2.4, -2.0},  Vec3{4.0, 2.4, 2.0}};
  for (std::size_t c = 0; c < centers.size(); ++c) {
    for (int i = 0; i < 16; ++i) {
      targets.push_back({centers[c].x + 0.037 * (i % 4),
                         centers[c].y + 0.028 * ((i / 4) % 4),
                         centers[c].z + 0.026 * (i / 8)});
    }
  }
  return targets;
}

std::vector<Vec3> multi_m2l_normals() {
  std::vector<Vec3> normals;
  normals.reserve(128);
  for (int c = 0; c < 8; ++c) {
    for (int i = 0; i < 16; ++i) {
      normals.push_back(
          unit({0.72 + 0.01 * i, -0.18 + 0.03 * c, 0.41 - 0.02 * i}));
    }
  }
  return normals;
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

KifmmBox source_box(std::span<const KifmmPanelSource> sources) {
  Vec3 lo{std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity()};
  Vec3 hi{-std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()};
  for (const auto &source : sources) {
    lo.x = std::min(lo.x, source.position.x);
    lo.y = std::min(lo.y, source.position.y);
    lo.z = std::min(lo.z, source.position.z);
    hi.x = std::max(hi.x, source.position.x);
    hi.y = std::max(hi.y, source.position.y);
    hi.z = std::max(hi.z, source.position.z);
  }
  const Vec3 center{0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y),
                    0.5 * (lo.z + hi.z)};
  const double extent =
      std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, 1.0e-9});
  return {center, 0.5 * extent * (1.0 + 1.0e-12), 0, 0};
}

double aggregate_rms(std::span<const KifmmJcpChannelRow> approx,
                     std::span<const KifmmJcpChannelRow> direct) {
  long double err_sq = 0.0L;
  long double ref_sq = 0.0L;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const long double ap = approx[i].k1 + approx[i].k2;
    const long double an = approx[i].k3 + approx[i].k4;
    const long double dp = direct[i].k1 + direct[i].k2;
    const long double dn = direct[i].k3 + direct[i].k4;
    err_sq += (ap - dp) * (ap - dp) + (an - dn) * (an - dn);
    ref_sq += dp * dp + dn * dn;
  }
  return std::sqrt(static_cast<double>(err_sq / std::max(ref_sq, 1.0e-300L)));
}

double channel_rms(std::span<const KifmmJcpChannelRow> approx,
                   std::span<const KifmmJcpChannelRow> direct, int channel) {
  long double err_sq = 0.0L;
  long double ref_sq = 0.0L;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const double *a = &approx[i].k1;
    const double *d = &direct[i].k1;
    const long double diff = a[channel] - d[channel];
    err_sq += diff * diff;
    ref_sq += static_cast<long double>(d[channel]) * d[channel];
  }
  return std::sqrt(static_cast<double>(err_sq / std::max(ref_sq, 1.0e-300L)));
}

double channel_abs_rms(std::span<const KifmmJcpChannelRow> approx,
                       std::span<const KifmmJcpChannelRow> direct,
                       int channel) {
  long double err_sq = 0.0L;
  for (std::size_t i = 0; i < approx.size(); ++i) {
    const double *a = &approx[i].k1;
    const double *d = &direct[i].k1;
    const long double diff = a[channel] - d[channel];
    err_sq += diff * diff;
  }
  return std::sqrt(static_cast<double>(err_sq /
                                      std::max<std::size_t>(approx.size(), 1)));
}

double channel_ref_rms(std::span<const KifmmJcpChannelRow> direct,
                       int channel) {
  long double ref_sq = 0.0L;
  for (const auto &row : direct) {
    const double *d = &row.k1;
    ref_sq += static_cast<long double>(d[channel]) * d[channel];
  }
  return std::sqrt(static_cast<double>(ref_sq /
                                      std::max<std::size_t>(direct.size(), 1)));
}

struct RowMetrics {
  double aggregate{0.0};
  std::array<double, 4> channel{};
  double max_channel{0.0};
  std::string max_channel_name{"K1"};
};

RowMetrics metrics(std::span<const KifmmJcpChannelRow> approx,
                   std::span<const KifmmJcpChannelRow> direct) {
  RowMetrics m;
  m.aggregate = aggregate_rms(approx, direct);
  for (int i = 0; i < 4; ++i) {
    m.channel[static_cast<std::size_t>(i)] = channel_rms(approx, direct, i);
  }
  m.max_channel = m.channel[0];
  for (int i = 1; i < 4; ++i) {
    if (m.channel[static_cast<std::size_t>(i)] > m.max_channel) {
      m.max_channel = m.channel[static_cast<std::size_t>(i)];
      m.max_channel_name = "K" + std::to_string(i + 1);
    }
  }
  return m;
}

std::vector<KifmmJcpChannelRow>
flat_kifmm_rows(std::span<const KifmmPanelSource> sources,
                std::span<const Vec3> targets,
                std::span<const Vec3> target_normals,
                const KifmmHierarchyJcpConfig &config,
                KifmmJcpChannels &channels) {
  const KifmmBox box = source_box(sources);
  const KifmmSurfaceRule rule = titan::tabi::fmm_kifmm::make_cube_surface_rule(
      config.hierarchy.order, config.hierarchy.equiv_scale,
      config.hierarchy.upward_check_scale);
  KifmmEquivalentSystemCache cache;
  channels = titan::tabi::fmm_kifmm::assemble_kifmm_jcp_channels(
      rule, box, sources, config.kappa, config.eps_ratio,
      1.0 / config.eps_ratio, "qr", &cache,
      config.hierarchy.condition_policy);
  std::vector<KifmmJcpChannelRow> rows(targets.size());
  if (!channels.finite) {
    return rows;
  }
  for (std::size_t i = 0; i < targets.size(); ++i) {
    rows[i] = titan::tabi::fmm_kifmm::assemble_kifmm_jcp_channel_row(
        rule, box, channels, config.kappa, config.eps_ratio,
        1.0 / config.eps_ratio, targets[i], target_normals[i]);
  }
  return rows;
}

struct FieldMetric {
  std::string name{};
  double potential_rms{0.0};
  double grad_rms{0.0};
};

std::string pass_fail_metric(double aggregate, double max_channel,
                             double aggregate_gate,
                             double channel_gate) {
  if (!std::isfinite(aggregate) || !std::isfinite(max_channel)) {
    return "fail";
  }
  return aggregate <= aggregate_gate && max_channel <= channel_gate ? "pass"
                                                                    : "fail";
}

std::string flat_accuracy_failure(double aggregate, double max_channel,
                                  double aggregate_gate,
                                  double channel_gate) {
  if (!std::isfinite(aggregate) || !std::isfinite(max_channel)) {
    return "nonfinite";
  }
  const bool aggregate_pass = aggregate <= aggregate_gate;
  const bool channel_pass = max_channel <= channel_gate;
  if (aggregate_pass && channel_pass) {
    return "none";
  }
  if (!aggregate_pass && !channel_pass) {
    return "aggregate_and_channel";
  }
  return aggregate_pass ? "channel" : "aggregate";
}

bool unsupported_status(const std::string &status) {
  return status.rfind("unsupported_", 0) == 0;
}

void emit_failed_result(const std::string &case_id,
                        const KifmmHierarchyJcpConfig &config,
                        const titan::tabi::fmm_kifmm::KifmmHierarchyJcpResult
                            &result) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::cout << "E253B_RESULT"
            << " case=" << case_id << " status=" << result.status
            << " execution_status=" << result.execution_status
            << " accuracy_status=" << result.accuracy_status
            << " accuracy_failure=" << result.accuracy_failure
            << " aggregate_gate=" << result.aggregate_relative_rms_gate
            << " channel_gate=" << result.max_channel_relative_rms_gate
            << " order=" << config.hierarchy.order
            << " theta=" << config.hierarchy.theta
            << " upward_check_scale=" << config.hierarchy.upward_check_scale
            << " local_equiv_scale=" << config.hierarchy.local_equiv_scale
            << " local_check_scale=" << config.hierarchy.local_check_scale
            << " condition_policy=" << config.hierarchy.condition_policy
            << " source_box_policy=" << config.hierarchy.source_box_policy
            << " aggregate_rms=" << nan << " k1_rms=" << nan
            << " k2_rms=" << nan << " k3_rms=" << nan
            << " k4_rms=" << nan << " max_channel_rms=" << nan
            << " max_channel_name=none"
            << " k4_abs_rms=" << nan << " k4_reference_rms=" << nan
            << " flat_status_legacy=not_evaluated"
            << " flat_execution_status=not_evaluated"
            << " flat_accuracy_status=not_evaluated"
            << " flat_accuracy_failure=not_evaluated"
            << " flat_aggregate_rms=" << nan << " flat_k1_rms=" << nan
            << " flat_k2_rms=" << nan << " flat_k3_rms=" << nan
            << " flat_k4_rms=" << nan << " flat_max_channel_rms=" << nan
            << " flat_max_channel_name=none"
            << " flat_k4_abs_rms=" << nan
            << " hier_vs_flat_aggregate_rms=" << nan
            << " hier_vs_flat_k1_rms=" << nan
            << " hier_vs_flat_k2_rms=" << nan
            << " hier_vs_flat_k3_rms=" << nan
            << " hier_vs_flat_k4_rms=" << nan
            << " hier_vs_flat_max_channel_rms=" << nan
            << " hier_vs_flat_max_channel_name=none"
            << " far_fraction=0 far_pairs=0 direct_pairs=0 coverage=0"
            << " source_count=0 target_count=0"
            << " m2m_translations=0 m2l_translations=0"
            << " l2l_translations=0 l2p_target_count=0"
            << " cache_entries=0 total_seconds=0" << '\n';
}

void emit_stage_row(const std::string &case_id, int order,
                    const std::string &stage, bool applicable,
                    double aggregate, double max_channel,
                    const std::string &max_channel_name, double k4_rms,
                    double aggregate_gate, double channel_gate,
                    const std::string &metric_source) {
  const std::string status =
      applicable ? pass_fail_metric(aggregate, max_channel, aggregate_gate,
                                    channel_gate)
                 : "not_applicable";
  std::cout << "E254_STAGE"
            << " case=" << case_id << " order=" << order
            << " stage=" << stage
            << " stage_applicable=" << (applicable ? "true" : "false")
            << " stage_status=" << status
            << " metric_source=" << metric_source
            << " aggregate_rms=" << aggregate
            << " max_channel_rms=" << max_channel
            << " max_channel_name=" << max_channel_name
            << " k4_rms=" << k4_rms
            << " aggregate_gate=" << aggregate_gate
            << " channel_gate=" << channel_gate << '\n';
}

std::vector<FieldMetric> flat_field_metrics(
    std::span<const KifmmPanelSource> sources, std::span<const Vec3> targets,
    const KifmmHierarchyJcpConfig &config, const KifmmJcpChannels &channels) {
  const KifmmBox box = source_box(sources);
  const KifmmSurfaceRule rule = titan::tabi::fmm_kifmm::make_cube_surface_rule(
      config.hierarchy.order, config.hierarchy.equiv_scale,
      config.hierarchy.upward_check_scale);
  const std::array<KifmmFieldFamily, 4> families{
      KifmmFieldFamily::PhinCharge, KifmmFieldFamily::PhiNormalX,
      KifmmFieldFamily::PhiNormalY, KifmmFieldFamily::PhiNormalZ};
  const std::array<std::string, 4> names{"phin", "phi_nx", "phi_ny",
                                         "phi_nz"};
  std::vector<FieldMetric> out;
  out.reserve(8);
  for (std::size_t fi = 0; fi < families.size(); ++fi) {
    for (const KifmmKernel kernel : {KifmmKernel::Laplace, KifmmKernel::Yukawa}) {
      const auto &solve =
          kernel == KifmmKernel::Laplace
              ? (fi == 0 ? channels.phin_laplace
                         : channels.phi_normal_laplace[fi - 1])
              : (fi == 0 ? channels.phin_yukawa
                         : channels.phi_normal_yukawa[fi - 1]);
      long double p_err = 0.0L;
      long double p_ref = 0.0L;
      long double g_err = 0.0L;
      long double g_ref = 0.0L;
      const double kappa =
          kernel == KifmmKernel::Laplace ? 0.0 : config.kappa;
      for (const Vec3 target : targets) {
        const double pd = titan::tabi::fmm_kifmm::direct_field_potential(
            sources, families[fi], kernel, kappa, target);
        const double pa = titan::tabi::fmm_kifmm::equivalent_field_potential(
            rule, box, solve, families[fi], kernel, kappa, target);
        const Vec3 gd = titan::tabi::fmm_kifmm::direct_field_grad_target(
            sources, families[fi], kernel, kappa, target);
        const Vec3 ga = titan::tabi::fmm_kifmm::equivalent_field_grad_target(
            rule, box, solve, families[fi], kernel, kappa, target);
        p_err += (pa - pd) * (pa - pd);
        p_ref += pd * pd;
        g_err += (ga.x - gd.x) * (ga.x - gd.x) +
                 (ga.y - gd.y) * (ga.y - gd.y) +
                 (ga.z - gd.z) * (ga.z - gd.z);
        g_ref += gd.x * gd.x + gd.y * gd.y + gd.z * gd.z;
      }
      out.push_back({names[fi] +
                         (kernel == KifmmKernel::Laplace ? "_laplace"
                                                          : "_yukawa"),
                     std::sqrt(static_cast<double>(p_err /
                                                   std::max(p_ref, 1.0e-300L))),
                     std::sqrt(static_cast<double>(g_err /
                                                   std::max(g_ref, 1.0e-300L)))});
    }
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  std::string case_id = "jcp_laplace_yukawa_separated_q6_check3_local3";
  bool e254_profile = false;
  KifmmHierarchyJcpConfig config;
  config.hierarchy.max_depth = 5;
  config.hierarchy.leaf_size = 2;
  config.hierarchy.order = 6;
  config.hierarchy.theta = 0.6;
  config.hierarchy.upward_check_scale = 3.0;
  config.hierarchy.equiv_scale = 1.05;
  config.hierarchy.local_equiv_scale = 3.0;
  config.hierarchy.local_check_scale = 1.05;
  config.aggregate_relative_rms_gate = 1.0e-6;
  config.max_channel_relative_rms_gate = 5.0e-6;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--e254-profile") {
      e254_profile = true;
    } else if (arg.rfind("--case=", 0) == 0) {
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
  } else if (case_id.find("multi_m2l") != std::string::npos) {
    sources = multi_m2l_sources();
    targets = multi_m2l_targets();
    normals = multi_m2l_normals();
    config.hierarchy.max_depth = 6;
    config.hierarchy.leaf_size = 8;
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
  if (case_id.find("_q10") != std::string::npos) {
    config.hierarchy.order = 10;
  }
  if (case_id.find("check2p4") != std::string::npos) {
    config.hierarchy.upward_check_scale = 2.4;
    config.hierarchy.local_equiv_scale = 2.4;
  }
  if (case_id.find("check4") != std::string::npos) {
    config.hierarchy.upward_check_scale = 4.0;
    config.hierarchy.local_equiv_scale = 4.0;
  }

  const auto result = titan::tabi::fmm_kifmm::run_full_hierarchical_jcp_kifmm(
      sources, targets, normals, config);
  if (result.execution_status != "pass") {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(12);
    emit_failed_result(case_id, config, result);
    return unsupported_status(result.status) ? 0 : 1;
  }
  KifmmJcpChannels flat_channels;
  const std::vector<KifmmJcpChannelRow> flat =
      flat_kifmm_rows(sources, targets, normals, config, flat_channels);
  const RowMetrics flat_vs_dense = metrics(flat, result.direct_rows);
  const RowMetrics hier_vs_flat = metrics(result.approximate_rows, flat);
  const double flat_k4_abs = channel_abs_rms(flat, result.direct_rows, 3);
  const double hier_k4_abs =
      channel_abs_rms(result.approximate_rows, result.direct_rows, 3);
  const double k4_ref = channel_ref_rms(result.direct_rows, 3);
  const std::string flat_execution_status =
      (flat_channels.finite && flat_channels.status == "pass") ? "pass"
                                                               : "fail";
  const std::string flat_accuracy_status =
      pass_fail_metric(flat_vs_dense.aggregate, flat_vs_dense.max_channel,
                       result.aggregate_relative_rms_gate,
                       result.max_channel_relative_rms_gate);
  const std::string flat_failure =
      flat_accuracy_failure(flat_vs_dense.aggregate, flat_vs_dense.max_channel,
                            result.aggregate_relative_rms_gate,
                            result.max_channel_relative_rms_gate);

  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << "E253B_RESULT"
            << " case=" << case_id << " status=" << result.status
            << " execution_status=" << result.execution_status
            << " accuracy_status=" << result.accuracy_status
            << " accuracy_failure=" << result.accuracy_failure
            << " aggregate_gate=" << result.aggregate_relative_rms_gate
            << " channel_gate=" << result.max_channel_relative_rms_gate
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
            << " k4_abs_rms=" << hier_k4_abs
            << " k4_reference_rms=" << k4_ref
            << " flat_status_legacy=" << flat_channels.status
            << " flat_execution_status=" << flat_execution_status
            << " flat_accuracy_status=" << flat_accuracy_status
            << " flat_accuracy_failure=" << flat_failure
            << " flat_aggregate_rms=" << flat_vs_dense.aggregate
            << " flat_k1_rms=" << flat_vs_dense.channel[0]
            << " flat_k2_rms=" << flat_vs_dense.channel[1]
            << " flat_k3_rms=" << flat_vs_dense.channel[2]
            << " flat_k4_rms=" << flat_vs_dense.channel[3]
            << " flat_max_channel_rms=" << flat_vs_dense.max_channel
            << " flat_max_channel_name=" << flat_vs_dense.max_channel_name
            << " flat_k4_abs_rms=" << flat_k4_abs
            << " hier_vs_flat_aggregate_rms=" << hier_vs_flat.aggregate
            << " hier_vs_flat_k1_rms=" << hier_vs_flat.channel[0]
            << " hier_vs_flat_k2_rms=" << hier_vs_flat.channel[1]
            << " hier_vs_flat_k3_rms=" << hier_vs_flat.channel[2]
            << " hier_vs_flat_k4_rms=" << hier_vs_flat.channel[3]
            << " hier_vs_flat_max_channel_rms=" << hier_vs_flat.max_channel
            << " hier_vs_flat_max_channel_name="
            << hier_vs_flat.max_channel_name
            << " far_fraction=" << result.stats.far_fraction
            << " far_pairs=" << result.stats.far_panel_pairs
            << " direct_pairs=" << result.stats.direct_panel_pairs
            << " coverage=" << result.stats.exact_pair_coverage
            << " source_count=" << result.stats.source_count
            << " target_count=" << result.stats.target_count
            << " m2m_translations=" << result.m2m_translations
            << " m2l_translations=" << result.m2l_translations
            << " l2l_translations=" << result.l2l_translations
            << " l2p_target_count=" << result.l2p_target_count
            << " cache_entries=" << result.equivalent_factor_cache_entries
            << " total_seconds=" << result.total_seconds
            << '\n';

  for (const FieldMetric &metric :
       flat_field_metrics(sources, targets, config, flat_channels)) {
    std::cout << "E253B_FIELD"
              << " case=" << case_id << " order=" << config.hierarchy.order
              << " field=" << metric.name
              << " flat_potential_rms=" << metric.potential_rms
              << " flat_grad_rms=" << metric.grad_rms << '\n';
  }

  std::cout << "E253B_K4"
            << " case=" << case_id << " order=" << config.hierarchy.order
            << " hierarchy_abs_rms=" << hier_k4_abs
            << " flat_abs_rms=" << flat_k4_abs
            << " reference_rms=" << k4_ref
            << " hierarchy_relative_rms=" << result.k4_relative_rms
            << " flat_relative_rms=" << flat_vs_dense.channel[3]
            << " hierarchy_vs_flat_relative_rms="
            << hier_vs_flat.channel[3] << '\n';

  if (e254_profile) {
    const bool separated_far_only = case_id.find("separated") != std::string::npos;
    const bool has_m2m = result.m2m_translations > 0;
    const bool has_m2l = result.m2l_translations > 0;
    const bool has_l2l = result.l2l_translations > 0;
    const double aggregate_gate = result.aggregate_relative_rms_gate;
    const double channel_gate = result.max_channel_relative_rms_gate;

    emit_stage_row(case_id, config.hierarchy.order, "p2m_flat_field",
                   separated_far_only, flat_vs_dense.aggregate,
                   flat_vs_dense.max_channel, flat_vs_dense.max_channel_name,
                   flat_vs_dense.channel[3], aggregate_gate, channel_gate,
                   "flat_vs_dense");

    const double translation_aggregate =
        separated_far_only ? hier_vs_flat.aggregate : result.aggregate_relative_rms;
    const double translation_max =
        separated_far_only ? hier_vs_flat.max_channel : result.max_channel_relative_rms;
    const std::string translation_channel =
        separated_far_only ? hier_vs_flat.max_channel_name : result.max_channel_name;
    const double translation_k4 =
        separated_far_only ? hier_vs_flat.channel[3] : result.k4_relative_rms;
    const std::string translation_source =
        separated_far_only ? "hierarchy_vs_flat" : "hierarchy_vs_dense";

    emit_stage_row(case_id, config.hierarchy.order, "m2m_parent_equivalent",
                   has_m2m, translation_aggregate, translation_max,
                   translation_channel, translation_k4, aggregate_gate,
                   channel_gate, translation_source);
    emit_stage_row(case_id, config.hierarchy.order, "m2l_local_check_rhs",
                   has_m2l, translation_aggregate, translation_max,
                   translation_channel, translation_k4, aggregate_gate,
                   channel_gate, translation_source);
    emit_stage_row(case_id, config.hierarchy.order, "local_solve_l2p", true,
                   result.aggregate_relative_rms,
                   result.max_channel_relative_rms, result.max_channel_name,
                   result.k4_relative_rms, aggregate_gate, channel_gate,
                   "hierarchy_vs_dense");
    emit_stage_row(case_id, config.hierarchy.order, "l2l_child_propagation",
                   has_l2l, translation_aggregate, translation_max,
                   translation_channel, translation_k4, aggregate_gate,
                   channel_gate, translation_source);
    emit_stage_row(case_id, config.hierarchy.order, "k4_mixed_derivative",
                   true, result.aggregate_relative_rms, result.k4_relative_rms, "K4",
                   result.k4_relative_rms, aggregate_gate, channel_gate,
                   "channel_k4");
  }

  return result.execution_status == "pass" ? 0 : 1;
}
